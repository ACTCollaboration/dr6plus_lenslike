"""Tests of the per-variant CMB-correction path (per_variant_corrs) for
dr6plus_variant0, dr6plus_optimal, actbase (9h) and act_planck (9a): products in
data/v1.0/like_corrs_<variant>/ from src/export_per_variant_corrs.py and
src/build_variant_corrs_from_stages.py (profile-hardened, ph, norm matrices;
N1-deconvolved form).

 1. Blinding: measured bandpowers refused without mock=True.
 2. lens_only analytic_marg: covariance = packaged base + Eq. 34 addition on
    the 40 < L < 1100 band, theory = B clkk + Delta_b, Hartlap on 599 sims;
    lnlike reproduced by hand and independent of the input CMB spectra. With
    mock=True the sky is fiducial (CMB = normalization fiducial), so Delta_b = 0
    and lnlike at the fiducial clkk is zero.
 3. Full mode: base covariance (not CMB marginalized); the CMB term vanishes
    at the fiducial CMB and equals Mb @ (C - C_fid) for a perturbed CMB.

Run with the cluster module environment loaded (see repo README).
"""
import os
import sys

import numpy as np
import pytest

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, REPO_ROOT)

from dr6plus_lenslike.dr6plus_lenslike import load_data, generic_lnlike  # noqa: E402

DDIR = os.path.join(REPO_ROOT, "dr6plus_lenslike", "data", "v1.0")
BASE = {"dr6plus_variant0": "covmat_clkk_dr6plus_variant0.txt",
        "dr6plus_optimal": "covmat_clkk_hilcTP_nightday_glsloo.txt",
        "actbase": "covmat_clkk_actbase.txt",
        "act_planck": "covmat_clkk_act_planck.txt"}
VARIANTS = [v for v in BASE if os.path.exists(os.path.join(DDIR, f"like_corrs_{v}", "response_binned.npz"))]
BAND = np.s_[2:-4]
TRIM = 2998

pytestmark = pytest.mark.skipif(not VARIANTS, reason="no per-variant products packaged")


def _fid_cls():
    ls, tt, ee, bb, te = np.loadtxt(os.path.join(DDIR, "like_corrs", "cosmo2017_10K_acc3_lensedCls.dat"),
                                    unpack=True)
    fac = 2 * np.pi / (ls * (ls + 1))
    ell = np.arange(ls.size + 2)
    out = {}
    for k, c in (("tt", tt), ("ee", ee), ("bb", bb), ("te", te)):
        a = np.zeros(ell.size)
        a[ls.astype(int)] = c * fac
        out[k] = a
    Lp, pp = np.loadtxt(os.path.join(DDIR, "like_corrs", "cosmo2017_10K_acc3_lenspotentialCls.dat"),
                        unpack=True, usecols=[0, 5])
    kk = np.zeros(int(Lp.max()) + 1)
    kk[Lp.astype(int)] = pp * 2 * np.pi / 4.
    return ell, out, np.arange(kk.size), kk


@pytest.mark.parametrize("V", VARIANTS)
def test_blinded(V):
    with pytest.raises(ValueError, match="BLINDED"):
        load_data(V, lens_only=True, analytic_marg=True, like_corrections=False)


@pytest.mark.parametrize("V", VARIANTS)
def test_lens_only_analytic_marg(V):
    LC = os.path.join(DDIR, f"like_corrs_{V}")
    d = load_data(V, lens_only=True, analytic_marg=True, like_corrections=False, mock=True)
    cov = np.loadtxt(os.path.join(LC, "covmat_cmbmarg_analytic.txt"))[BAND, BAND]
    assert np.allclose(d["cov"], cov, rtol=0, atol=0)
    hart = (599 - 12 - 2.) / (599 - 1.)
    assert np.allclose(d["cinv"], np.linalg.inv(cov) * hart, rtol=1e-12)
    # mock: fiducial sky, no Eq. 35 recentering (the file Delta_b encodes the
    # measured CMB spectra and is nonzero)
    delta = np.loadtxt(os.path.join(LC, "lens_only_delta.txt"))[BAND]
    assert np.any(delta != 0)
    assert np.array_equal(d["lens_only_delta"], np.zeros(delta.size))
    ell, cl, Lk, kk = _fid_cls()
    lnl, th = generic_lnlike(d, Lk, kk, ell, cl["tt"], cl["ee"], cl["te"], cl["bb"], return_theory=True)
    bclkk = d["binmat_act"] @ np.pad(kk, (0, max(0, TRIM + 2 - kk.size)))[:TRIM + 2]
    assert np.allclose(th, bclkk, rtol=1e-12)
    r = d["data_binned_clkk"] - bclkk
    assert np.isclose(lnl, -0.5 * r @ d["cinv"] @ r, rtol=1e-12)
    assert abs(lnl) < 1e-10
    lnl2 = generic_lnlike(d, Lk, kk, ell, 1.05 * cl["tt"], cl["ee"], cl["te"], cl["bb"])
    assert lnl2 == lnl


@pytest.mark.parametrize("V", VARIANTS)
def test_full_mode(V):
    LC = os.path.join(DDIR, f"like_corrs_{V}")
    d = load_data(V, lens_only=False, like_corrections=True, mock=True)
    base = np.loadtxt(os.path.join(DDIR, BASE[V]))[BAND, BAND]
    assert np.allclose(d["cov"], base, rtol=0, atol=0)
    assert d.get("binned_corr") and d["Mb_binned"].shape == (3, 12, TRIM + 2)
    ell, cl, Lk, kk = _fid_cls()
    _, th0 = generic_lnlike(d, Lk, kk, ell, cl["tt"], cl["ee"], cl["te"], cl["bb"], return_theory=True)
    b0 = d["binmat_act"] @ np.pad(kk, (0, max(0, TRIM + 2 - kk.size)))[:TRIM + 2]
    assert np.allclose(th0, b0, rtol=1e-10)
    _, th1 = generic_lnlike(d, Lk, kk, ell, 1.02 * cl["tt"], cl["ee"], cl["te"], cl["bb"], return_theory=True)
    Mb = np.load(os.path.join(LC, "response_binned.npz"))["Mb"][0][BAND]
    exp = Mb[:, :TRIM + 1] @ (0.02 * cl["tt"][:TRIM + 1])
    assert np.allclose(th1 - th0, exp, rtol=1e-8)
