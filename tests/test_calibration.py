"""Tests for dr6plus_lenslike/calibration.py."""

import os
import numpy as np
import pytest

from dr6plus_lenslike.calibration import load_calibration_weights, compute_dCl_from_calibration

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "src", "dr6plus_lenslike", "data", "v1.0")
DATA_DIR = os.path.normpath(DATA_DIR)


# ---------------------------------------------------------------------------
# Test A: shape and non-negativity of loaded weights
# ---------------------------------------------------------------------------
def test_load_calibration_weights_shape_and_nonneg():
    """Weight files should load to (4, 5001) arrays with all non-negative values."""
    w_T, w_E = load_calibration_weights(DATA_DIR)

    assert w_T.shape == (4, 5001), f"w_T shape {w_T.shape} != (4, 5001)"
    assert w_E.shape == (4, 5001), f"w_E shape {w_E.shape} != (4, 5001)"

    assert np.all(w_T >= 0), "w_T contains negative values"
    assert np.all(w_E >= 0), "w_E contains negative values"


# ---------------------------------------------------------------------------
# Test B: zero calibration → zero ΔC
# ---------------------------------------------------------------------------
def test_zero_calibration_gives_zero_dC():
    """With delta_c = delta_p = 0, all four dC outputs must be identically zero."""
    N = 200
    w_T = np.ones((4, N))
    w_E = np.ones((4, N))
    C_TT = np.ones(N) * 5.0
    C_EE = np.ones(N) * 3.0
    C_TE = np.ones(N) * 2.0

    dC_TT, dC_EE, dC_TE, dC_BB = compute_dCl_from_calibration(
        np.zeros(4), np.zeros(4), w_T, w_E, C_TT, C_EE, C_TE
    )

    np.testing.assert_array_equal(dC_TT, 0.0)
    np.testing.assert_array_equal(dC_EE, 0.0)
    np.testing.assert_array_equal(dC_TE, 0.0)
    np.testing.assert_array_equal(dC_BB, 0.0)


# ---------------------------------------------------------------------------
# Test C: dC_BB is always zero
# ---------------------------------------------------------------------------
def test_dC_BB_always_zero():
    """dC_BB must be zero even when delta_c is nonzero."""
    N = 100
    w_T = np.random.rand(4, N) + 0.1
    w_E = np.random.rand(4, N) + 0.1
    delta_c = np.array([0.01, -0.02, 0.005, 0.015])
    delta_p = np.zeros(4)
    C_TT = np.random.rand(N) + 1.0
    C_EE = np.random.rand(N) + 0.5
    C_TE = np.random.rand(N) * 0.5

    _, _, _, dC_BB = compute_dCl_from_calibration(
        delta_c, delta_p, w_T, w_E, C_TT, C_EE, C_TE
    )

    np.testing.assert_array_equal(dC_BB, 0.0)


# ---------------------------------------------------------------------------
# Test D: uniform gain + uniform weights → analytic result
# ---------------------------------------------------------------------------
def test_uniform_gain_uniform_weights():
    """With uniform weights and identical gain for all arrays, eff = delta."""
    N = 100
    delta = 0.02
    delta_c = np.full(4, delta)
    delta_p = np.zeros(4)

    # Uniform weights: every array contributes equally at every ℓ.
    w_T = np.ones((4, N))
    w_E = np.ones((4, N))

    C_TT = np.linspace(10.0, 1.0, N)
    C_EE = np.linspace(5.0, 0.5, N)
    C_TE = np.linspace(3.0, 0.1, N)

    dC_TT, dC_EE, dC_TE, dC_BB = compute_dCl_from_calibration(
        delta_c, delta_p, w_T, w_E, C_TT, C_EE, C_TE
    )

    # eff_T = eff_E = delta everywhere, so:
    #   dC_TT = 2 * delta * C_TT_fid
    #   dC_TE = (delta + delta) * C_TE_fid = 2 * delta * C_TE_fid
    np.testing.assert_allclose(dC_TT, 2 * delta * C_TT, atol=1e-12)
    np.testing.assert_allclose(dC_TE, 2 * delta * C_TE, atol=1e-12)


# ---------------------------------------------------------------------------
# Test E: pol efficiency enters E but not T
# ---------------------------------------------------------------------------
def test_poleff_enters_E_not_T():
    """Pol efficiency affects EE and TE but leaves TT unchanged."""
    N = 100
    delta = 0.03
    delta_c = np.zeros(4)
    delta_p = np.full(4, delta)

    w_T = np.ones((4, N))
    w_E = np.ones((4, N))

    C_TT = np.linspace(8.0, 1.0, N)
    C_EE = np.linspace(4.0, 0.4, N)
    C_TE = np.linspace(2.0, 0.1, N)

    dC_TT, dC_EE, dC_TE, dC_BB = compute_dCl_from_calibration(
        delta_c, delta_p, w_T, w_E, C_TT, C_EE, C_TE
    )

    # eff_T = 0 (no gain deviation), eff_E = delta (pol efficiency only).
    np.testing.assert_array_equal(dC_TT, 0.0)   # T unaffected by pol efficiency
    np.testing.assert_allclose(dC_EE, 2 * delta * C_EE, atol=1e-12)
    np.testing.assert_allclose(dC_TE, delta * C_TE, atol=1e-12)   # one E factor


# ---------------------------------------------------------------------------
# Test F: live sanity check with real weights and fiducial spectra
# ---------------------------------------------------------------------------
def test_live_sanity_real_data():
    """1% uniform gain change should produce a nonzero dC_TT in a physically
    reasonable range when using the real weight files and fiducial spectra."""
    w_T, w_E = load_calibration_weights(DATA_DIR)

    # Load fiducial CMB spectra: columns are ℓ, TT, EE, BB, TE in Dℓ (µK²).
    fid_file = os.path.join(DATA_DIR, "like_corrs", "cosmo2017_10K_acc3_lensedCls.dat")
    fid = np.loadtxt(fid_file)
    ells = fid[:, 0].astype(int)   # ℓ values starting at 2

    # Build Cℓ arrays of length 5001 (ℓ = 0..5000), zero-padded at low ℓ.
    N_ell = 5001
    C_TT = np.zeros(N_ell)
    C_EE = np.zeros(N_ell)
    C_TE = np.zeros(N_ell)

    for i, ell in enumerate(ells):
        if ell >= N_ell:
            break
        # Convert Dℓ → Cℓ: Cℓ = Dℓ × 2π / (ℓ(ℓ+1))
        factor = 2 * np.pi / (ell * (ell + 1))
        C_TT[ell] = fid[i, 1] * factor   # TT
        C_EE[ell] = fid[i, 2] * factor   # EE
        C_TE[ell] = fid[i, 4] * factor   # TE (column 4; column 3 is BB)

    delta_c = np.full(4, 0.01)   # 1% uniform gain
    delta_p = np.zeros(4)

    dC_TT, dC_EE, dC_TE, dC_BB = compute_dCl_from_calibration(
        delta_c, delta_p, w_T, w_E, C_TT, C_EE, C_TE
    )

    # dC_TT should be nonzero in the main survey range ℓ=600..3000.
    assert np.any(dC_TT[600:3001] != 0.0), "dC_TT is zero over ℓ=600..3000"

    # Sanity range for 1% gain in Cℓ µK² units: roughly 2% of C_TT.
    # At ℓ~600, C_TT ≈ D_TT × 2π/(ℓ(ℓ+1)) ≈ 3000 × 1.7e-5 ≈ 0.05 µK²,
    # so dC_TT ≈ 2% × 0.05 ≈ 1e-3 µK².  Upper bound set to 1e-2 for headroom.
    max_abs = np.max(np.abs(dC_TT[600:3001]))
    assert 1e-10 < max_abs < 1e-2, (
        f"max|dC_TT[600:3000]| = {max_abs:.3e} outside sanity range [1e-10, 1e-2]"
    )

    # BB must be identically zero.
    np.testing.assert_array_equal(dC_BB, 0.0)
