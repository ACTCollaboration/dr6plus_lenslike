"""
Calibration-induced CMB spectrum corrections for the DR6+ lensing likelihood.

Overview
--------
The ACT DR6+ maps are made from four detector arrays:
    pa5a = pa5_f090,  pa5b = pa5_f150,
    pa6a = pa6_f090,  pa6b = pa6_f150.

Two types of calibration nuisance affect the maps:

1. **Gain** (c_a, map-level): a multiplicative factor applied to each array's
   raw map before coadding.  Because the same factor multiplies both T and E
   pixels, it affects all three CMB power spectra: TT, TE, and EE.

2. **Polarisation efficiency** (p_a, E-only): an additional multiplicative
   factor that affects only the polarisation response of the detector.  A
   miscalibrated pol-efficiency rescales the E map by an extra (1+p_a) factor
   (to first order), so it enters EE twice, TE once, but not TT at all.

Coadding weights
----------------
The final T or E map is a noise-weighted coadd over the four arrays.  The
files `noise_{arr}_{T,E}_weights.txt` give, for each multipole ℓ from 0 to
5000, the noise-weighted contribution of array `arr` to the coadded map.
Concretely, if the coadded T map is

    T_coadd(ℓ) = Σ_a  w_{T,a}(ℓ) / [Σ_a w_{T,a}(ℓ)]  ×  T_a(ℓ)

then the effective fractional calibration change in the coadded map is

    eff_T(ℓ) = Σ_a  w_{T,a}(ℓ) × δc_a  /  Σ_a w_{T,a}(ℓ)

and similarly for E (replacing δc with δc+δp).  Stacking the four arrays
into a (4, N_ell) matrix makes this a single matrix–vector operation.

Units
-----
All ΔC values returned by `compute_dCl_from_calibration` are in Cℓ (not Dℓ)
units, because `get_corrected_clkk` in the main likelihood expects Cℓ.

MCMC interface
--------------
`delta_c` and `delta_p` are MCMC nuisance parameters supplied at each step by
the Cobaya sampler.  This module does not load or store calibration samples;
it only evaluates the linearised CMB-spectrum response given a particular
(delta_c, delta_p) pair.

T² term
-------
The separate term that scales theory clkk by R(delta_c) uses
`response_cal_matrix.txt` (18×9).  That matrix is loaded by the likelihood's
`load_data()` method; this module handles only the CMB spectrum side.
"""

import os
import numpy as np


# ---------------------------------------------------------------------------
# Array ordering used throughout this module.
# ---------------------------------------------------------------------------
_ARRAYS = ["pa5a", "pa5b", "pa6a", "pa6b"]


def load_calibration_weights(data_dir):
    """Load noise-coadding weight files and return stacked arrays.

    Parameters
    ----------
    data_dir : str
        Directory containing the eight weight files
        `noise_{arr}_{T,E}_weights.txt` (e.g. `dr6plus_lenslike/data/v1.0/`).

    Returns
    -------
    w_T : ndarray, shape (4, 5001)
        T coadding weights.  Row order: [pa5a, pa5b, pa6a, pa6b].
        Each row is the per-ℓ noise weight of that array in the coadded T map
        for ℓ = 0, 1, …, 5000.
    w_E : ndarray, shape (4, 5001)
        E coadding weights.  Same row order and ℓ indexing.
    """
    rows_T = []
    rows_E = []

    for arr in _ARRAYS:
        # Each file contains 5001 lines, one weight per ℓ from ℓ=0 to ℓ=5000.
        # The weight gives the ℓ-by-ℓ contribution of this array to the coadd.
        rows_T.append(
            np.loadtxt(os.path.join(data_dir, f"noise_{arr}_T_weights.txt"))
        )
        rows_E.append(
            np.loadtxt(os.path.join(data_dir, f"noise_{arr}_E_weights.txt"))
        )

    # Stack into (4, 5001): rows are arrays, columns are ℓ values.
    # The shape (4, N_ell) lets us compute weighted sums over arrays with a
    # single matrix operation: (w * delta[:, None]).sum(axis=0).
    w_T = np.vstack(rows_T)   # shape (4, 5001)
    w_E = np.vstack(rows_E)   # shape (4, 5001)

    return w_T, w_E


def compute_dCl_from_calibration(
    delta_c, delta_p, w_T, w_E, C_TT_fid, C_EE_fid, C_TE_fid
):
    """Compute calibration-induced changes to the CMB power spectra.

    This mirrors `_calibrate_spectra` from `act_dr6_mflike` exactly.

    Parameters
    ----------
    delta_c : array_like, shape (4,)
        Fractional gain deviations for [pa5a, pa5b, pa6a, pa6b].
        ``delta_c[a] = sampled_gain_a − fiducial_gain_a`` (dimensionless).
        Gain c_a is map-level: the T and E maps of array a are both multiplied
        by (1 + c_a), so gain enters TT, TE, and EE.
    delta_p : array_like, shape (4,)
        Pol-efficiency deviations for [pa5a, pa5b, pa6a, pa6b].
        ``delta_p[a] = sampled_poleff_a − fiducial_poleff_a`` (dimensionless).
        Pass ``np.zeros(4)`` to vary gain only and leave pol efficiency fixed.
        Pol efficiency p_a is E-only: the E map of array a is multiplied by
        an additional (1 + p_a), so pol efficiency enters EE and TE but not TT.
    w_T : ndarray, shape (4, N_ell)
        T coadding weights from `load_calibration_weights`.
    w_E : ndarray, shape (4, N_ell)
        E coadding weights from `load_calibration_weights`.
    C_TT_fid : ndarray, shape (N_ell,)
        Fiducial TT power spectrum in Cℓ units (µK²).
    C_EE_fid : ndarray, shape (N_ell,)
        Fiducial EE power spectrum in Cℓ units (µK²).
    C_TE_fid : ndarray, shape (N_ell,)
        Fiducial TE power spectrum in Cℓ units (µK²).

    Returns
    -------
    dC_TT : ndarray, shape (N_ell,)
        Calibration-induced change ΔC_TT in Cℓ units.
    dC_EE : ndarray, shape (N_ell,)
        Calibration-induced change ΔC_EE in Cℓ units.
    dC_TE : ndarray, shape (N_ell,)
        Calibration-induced change ΔC_TE in Cℓ units.
    dC_BB : ndarray, shape (N_ell,)
        Calibration-induced change ΔC_BB in Cℓ units.  Always zero: lensing
        B-modes are small and BB is not used in the norm correction.
    """
    delta_c = np.asarray(delta_c, dtype=float)   # (4,)
    delta_p = np.asarray(delta_p, dtype=float)   # (4,)

    # ------------------------------------------------------------------
    # Sum weights over arrays at each ℓ.
    # Used as the denominator of the noise-weighted average.
    # ------------------------------------------------------------------
    w_T_sum = w_T.sum(axis=0)   # (N_ell,)
    w_E_sum = w_E.sum(axis=0)   # (N_ell,)

    # Replace zeros with 1.0 to avoid division-by-zero at low ℓ where no
    # array contributes (e.g. ℓ < survey minimum).  The resulting eff values
    # at those ℓ are zeroed out explicitly below.
    w_T_sum_safe = np.where(w_T_sum > 0, w_T_sum, 1.0)
    w_E_sum_safe = np.where(w_E_sum > 0, w_E_sum, 1.0)

    # ------------------------------------------------------------------
    # Effective fractional T gain at each ℓ.
    # The coadded T map is T_coadd = Σ_a w_a T_a / Σ_a w_a.
    # Gain rescales T_a → (1+c_a) T_a, so to first order
    #     δT_coadd / T_coadd = Σ_a w_a δc_a / Σ_a w_a = eff_T.
    # ------------------------------------------------------------------
    eff_T = (w_T * delta_c[:, None]).sum(axis=0) / w_T_sum_safe   # (N_ell,)

    # ------------------------------------------------------------------
    # Effective fractional E gain at each ℓ.
    # Gain (c_a) AND pol efficiency (p_a) both rescale E:
    #     E_a → (1+c_a)(1+p_a) E_a ≈ (1 + c_a + p_a) E_a  (first order).
    # So the combined per-array deviation entering E is delta_c + delta_p.
    # ------------------------------------------------------------------
    eff_E = (w_E * (delta_c + delta_p)[:, None]).sum(axis=0) / w_E_sum_safe   # (N_ell,)

    # ------------------------------------------------------------------
    # Zero out effective gains where no array contributes at this ℓ.
    # Without this, the safe-denominator trick above would give spurious
    # nonzero eff values at masked multipoles.
    # ------------------------------------------------------------------
    eff_T[w_T_sum == 0] = 0.0
    eff_E[w_E_sum == 0] = 0.0

    # ------------------------------------------------------------------
    # TT = T × T:  both factors pick up (1+eff_T), so the fractional
    # change in the power spectrum is (1+eff_T)² − 1 ≈ 2 eff_T.
    # Multiply by C_TT_fid to get the absolute change in Cℓ.
    # ------------------------------------------------------------------
    dC_TT = 2 * eff_T * C_TT_fid

    # ------------------------------------------------------------------
    # EE = E × E:  same argument with eff_E.
    # ------------------------------------------------------------------
    dC_EE = 2 * eff_E * C_EE_fid

    # ------------------------------------------------------------------
    # TE = T × E:  one T factor and one E factor, so the first-order
    # fractional change is eff_T + eff_E.
    # ------------------------------------------------------------------
    dC_TE = (eff_T + eff_E) * C_TE_fid

    # ------------------------------------------------------------------
    # BB:  lensing B-modes are much smaller than T/E, and BB is not used
    # in get_corrected_clkk's normalization correction.  Return zeros.
    # ------------------------------------------------------------------
    dC_BB = np.zeros_like(C_TT_fid)

    return dC_TT, dC_EE, dC_TE, dC_BB
