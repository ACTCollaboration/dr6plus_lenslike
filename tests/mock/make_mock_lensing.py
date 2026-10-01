"""Build the mock lensing bandpowers of the verification tests.

The mock sky is the DR6 simulation / normalization fiducial cosmology
cosmo2017_10K_acc3 (parameters in cosmo2017/inputs/cosmo2017_10K_acc3_params.ini).
The mock bandpowers are the binned fiducial C_L^kk of that cosmology,

    b_i = sum_L B_iL C_L^kk,   C_L^kk = [L(L+1)]^2 C_L^phiphi / 4,

with B the 18-bin DR6 binning matrix (column j <-> L = j) and C_L^phiphi read
from the 2017 CAMB output cosmo2017/inputs/cosmo2017_10K_acc3_lenspotentialCls.dat
(verbatim copy of the file packaged with the likelihood corrections). This
reproduces the packaged data/v1.0/clkk_bandpowers_fiducial.txt; the bin sums
use math.fsum (exactly rounded, independent of BLAS and thread count), which
agrees with the packaged file bit for bit in most bins and to one unit in the
last place in the others (the packaged file came from a different summation).

Usage (repo root, cluster env loaded):  python tests/mock/make_mock_lensing.py
"""
import datetime
import hashlib
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
DDIR = os.path.join(REPO, "src", "dr6plus_lenslike", "data", "v1.0")
OUT = os.path.join(HERE, "cosmo2017")

INPUTS = {
    "lenspotential_cls": os.path.join(OUT, "inputs", "cosmo2017_10K_acc3_lenspotentialCls.dat"),
    "binning_matrix": os.path.join(DDIR, "binning_matrix_act.txt"),
    "params_ini": os.path.join(OUT, "inputs", "cosmo2017_10K_acc3_params.ini"),
}
# the likelihood's own copy of the spectra, when present (a symlinked data dir)
PACKAGED_CLS = os.path.join(DDIR, "like_corrs", "cosmo2017_10K_acc3_lenspotentialCls.dat")
PACKAGED = os.path.join(DDIR, "clkk_bandpowers_fiducial.txt")
MOCK = os.path.join(OUT, "clkk_bandpowers_mock.txt")


def md5(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def fiducial_clkk():
    """C_L^kk on L = 0..Lmax from the CAMB lenspotentialCls file (PP column)."""
    L, pp = np.loadtxt(INPUTS["lenspotential_cls"], unpack=True, usecols=[0, 5])
    kk = np.zeros(int(L.max()) + 1)
    kk[L.astype(int)] = pp * 2. * np.pi / 4.
    return kk


def binned(binmat, kk):
    x = kk[:binmat.shape[1]]
    return np.array([math.fsum(binmat[i] * x) for i in range(binmat.shape[0])])


def main():
    B = np.loadtxt(INPUTS["binning_matrix"])
    mock = binned(B, fiducial_clkk())
    np.savetxt(MOCK, mock)

    ref = np.loadtxt(PACKAGED)
    rel = np.abs(mock / ref - 1.)
    check = {
        "packaged_file": os.path.relpath(PACKAGED, REPO),
        "n_bins": int(mock.size),
        "n_bitwise_equal": int(np.sum(mock == ref)),
        "max_rel_diff": float(rel.max()),
        "passed": bool(rel.max() <= 2 * np.finfo(float).eps),
    }
    if os.path.exists(PACKAGED_CLS):
        check["input_cls_identical_to_likelihood_copy"] = md5(PACKAGED_CLS) == md5(INPUTS["lenspotential_cls"])
    prov = {
        "product": "mock lensing bandpowers, 18 DR6 bins (likelihood band 40<L<1100 = bins 2..13)",
        "cosmology": "cosmo2017_10K_acc3 (DR6 simulation and normalization fiducial)",
        "written": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "script": os.path.relpath(__file__, REPO),
        "script_md5": md5(__file__),
        "inputs": {k: {"path": os.path.relpath(v, REPO), "md5": md5(v)} for k, v in INPUTS.items()},
        "output": {"path": os.path.relpath(MOCK, REPO), "md5": md5(MOCK)},
        "check_vs_packaged": check,
    }
    with open(os.path.join(OUT, "PROVENANCE_lensing.json"), "w") as f:
        json.dump(prov, f, indent=1)
    print(json.dumps(check, indent=1))
    if not check["passed"]:
        raise SystemExit("mock does not reproduce the packaged fiducial bandpowers")


if __name__ == "__main__":
    main()
