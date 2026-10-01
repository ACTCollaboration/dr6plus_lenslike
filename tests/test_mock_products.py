"""Tests of the verification mock products in tests/mock/ and of the mock_file
option of load_data.

 1. The lensing mock (tests/mock/cosmo2017/clkk_bandpowers_mock.txt) is the
    binned cosmo2017 C_L^kk and reproduces the packaged
    clkk_bandpowers_fiducial.txt to float64 rounding.
 2. mock_file routes the blinded variants to that file, and lnlike at the
    fiducial C_L^kk is zero in lens_only analytic_marg mode.
 3. mock_file without mock=True is refused.
 4. The primary-CMB and BAO mocks (make_mock_cmb_bao.py) match their recorded
    md5s and keep the real covariances.
 5. Tier 0: every likelihood of the verification runs, on the mock files and at
    the mock point in the chains' parametrization (theta_fid.yaml), through
    cobaya: CMB and BAO chi2 ~ 0; lensing chi2 = the recorded 2017-file vs
    current-CAMB theory offset (LENS_CHI2); cobaya passes mock_file through
    (data vector bit-equal to the mock file).

Run with the cluster module environment loaded (see repo README); test 5 makes
one CAMB call at P-ACT accuracy (about 10 s).
"""
import contextlib
import importlib.util
import io
import os
import sys

import numpy as np
import pytest

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from dr6plus_lenslike.dr6plus_lenslike import load_data, generic_lnlike  # noqa: E402

DDIR = os.path.join(REPO_ROOT, "src", "dr6plus_lenslike", "data", "v1.0")
MOCK = os.path.join(REPO_ROOT, "tests", "mock", "cosmo2017", "clkk_bandpowers_mock.txt")
VARIANTS = ("dr6plus_variant0", "dr6plus_optimal", "actbase", "act_planck")
BAND = np.s_[2:-4]
TRIM = 2998

pytestmark = pytest.mark.skipif(not os.path.exists(MOCK), reason="run tests/mock/make_mock_lensing.py first")


def _generator():
    spec = importlib.util.spec_from_file_location(
        "make_mock_lensing", os.path.join(REPO_ROOT, "tests", "mock", "make_mock_lensing.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_lensing_mock_reproduces_packaged_fiducial():
    gen = _generator()
    B = np.loadtxt(gen.INPUTS["binning_matrix"])
    rebuilt = gen.binned(B, gen.fiducial_clkk())
    mock = np.loadtxt(MOCK)
    assert np.array_equal(rebuilt, mock)
    ref = np.loadtxt(os.path.join(DDIR, "clkk_bandpowers_fiducial.txt"))
    assert np.max(np.abs(mock / ref - 1)) <= 2 * np.finfo(float).eps


@pytest.mark.parametrize("V", VARIANTS)
def test_mock_file_lens_only(V):
    with contextlib.redirect_stdout(io.StringIO()):
        d = load_data(V, lens_only=True, analytic_marg=True, like_corrections=False,
                      mock=True, mock_file=MOCK)
    assert np.array_equal(d["data_binned_clkk"], np.loadtxt(MOCK)[BAND])
    assert not np.any(d["lens_only_delta"])
    gen = _generator()
    kk = gen.fiducial_clkk()
    Lk = np.arange(kk.size)
    ell = np.arange(TRIM + 2)
    zero = np.zeros(ell.size)
    lnl = generic_lnlike(d, Lk, kk, ell, zero, zero, zero, zero)
    assert abs(lnl) < 1e-10


@pytest.mark.parametrize("name,key", [("ACTbase", "actbase"), ("ACT_Planck", "act_planck")])
def test_mixed_case_variant_names(name, key):
    # names are lowercased; 'act_planck' must not switch on the Planck lensing bandpowers
    with contextlib.redirect_stdout(io.StringIO()):
        d = load_data(name, lens_only=True, analytic_marg=True, like_corrections=False, mock=True)
        e = load_data(key, lens_only=True, analytic_marg=True, like_corrections=False, mock=True)
    assert not d["include_planck"]
    assert d["data_binned_clkk"].size == 12
    assert np.array_equal(d["cov"], e["cov"])


def test_mock_file_requires_mock():
    with pytest.raises(ValueError, match="mock_file requires mock=True"):
        load_data("dr6plus_variant0", lens_only=True, analytic_marg=True,
                  like_corrections=False, mock_file=MOCK)


def _cmb_bao_generator():
    spec = importlib.util.spec_from_file_location(
        "make_mock_cmb_bao", os.path.join(REPO_ROOT, "tests", "mock", "make_mock_cmb_bao.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PROV_CMB_BAO = os.path.join(REPO_ROOT, "tests", "mock", "cosmo2017", "PROVENANCE_cmb_bao.json")
# Tier-0 lensing chi2 at the mock point (ACT reference theory block, CAMB 1.6.2), 2026-09-29
LENS_CHI2 = {"lens_joint_dr6plus_variant0": 0.02864, "lens_only_dr6plus_variant0": 0.01853,
             "lens_joint_dr6plus_optimal": 0.04787, "lens_only_dr6plus_optimal": 0.03045,
             "lens_joint_actbase": 0.03485, "lens_only_actbase": 0.02250,
             "lens_joint_act_planck": 0.04774, "lens_only_act_planck": 0.03046}


@pytest.mark.skipif(not os.path.exists(PROV_CMB_BAO), reason="run tests/mock/make_mock_cmb_bao.py first")
def test_cmb_bao_mocks_recorded_and_real_covariances():
    import hashlib
    import json
    import sacc
    prov = json.load(open(PROV_CMB_BAO))
    for rel, h in prov["outputs"].items():
        assert hashlib.md5(open(os.path.join(REPO_ROOT, rel), "rb").read()).hexdigest() == h, rel
    assert prov["swap_test_passed"]
    g = _cmb_bao_generator()
    ins = prov["inputs"]
    mock, real = sacc.Sacc.load_fits(g.ACT_MOCK), sacc.Sacc.load_fits(ins["act_sacc"]["path"])
    assert np.array_equal(mock.covariance.covmat, real.covariance.covmat)
    assert len(mock.mean) == len(real.mean) and not np.array_equal(mock.mean, real.mean)
    mean_md5 = hashlib.md5(np.ascontiguousarray(mock.mean, dtype="<f8").tobytes()).hexdigest()
    assert mean_md5 == prov["act_mock_mean_md5"]
    same = {"plik/c_matrix_plik_v22.dat": os.path.join(g.PLIK_MOCK_DIR, "c_matrix_plik_v22.dat"),
            "bao_cov": os.path.join(g.BAO_MOCK_DIR, g.BAO_COV)}
    for k, p in same.items():
        assert hashlib.md5(open(p, "rb").read()).hexdigest() == ins[k]["md5"], k


@pytest.mark.skipif(not os.path.exists(PROV_CMB_BAO), reason="run tests/mock/make_mock_cmb_bao.py first")
def test_tier0_all_likelihoods_at_mock_point():
    os.environ.setdefault("COBAYA_NOMPI", "1")
    import yaml
    from cobaya.model import get_model
    g = _cmb_bao_generator()
    fid = yaml.safe_load(open(os.path.join(g.OUT, "theta_fid.yaml")))
    params = {k: v for k, v in fid.items() if k != "logA"}
    params["As"] = float(1e-10 * np.exp(fid["logA"]))
    lik = g.likelihood_block(mock=True)
    for V in VARIANTS:
        common = dict(variant=V, mock=True, mock_file=MOCK, trim_lmax=TRIM, apply_hartlap=True, lmax=4000)
        lik[f"lens_joint_{V}"] = {"class": "dr6plus_lenslike.ACTDR6LensLike", **common, "lens_only": False}
        lik[f"lens_only_{V}"] = {"class": "dr6plus_lenslike.ACTDR6LensLike", **common, "lens_only": True,
                                 "analytic_marg": True, "no_like_corrections": True, "cov_file": None}
    with contextlib.redirect_stdout(io.StringIO()):
        m = get_model({"theory": g.theory_block(), "likelihood": lik, "params": params})
        ll = m.loglikes({}, as_dict=True, return_derived=False)
    chi2 = {k: -2 * v for k, v in ll.items()}
    for k in ("act_dr6_cmbonly.PlanckActCut", "act_dr6_cmbonly.ACTDR6CMBonly", "bao.desi_dr2.desi_bao_all"):
        assert chi2[k] < 1e-3, (k, chi2[k])
    for V in VARIANTS:
        for mode in ("joint", "only"):
            name = f"lens_{mode}_{V}"
            assert abs(chi2[name] - LENS_CHI2[name]) < 1e-3, (name, chi2[name])
            assert np.array_equal(m.likelihood[name].data["data_binned_clkk"], np.loadtxt(MOCK)[BAND])
