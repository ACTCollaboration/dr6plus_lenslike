# ACT DR4+Planck CMB Covariance Verification - Findings

## Objective
Verify whether using ACT DR4+Planck 2015 CMB covariance (instead of ACT DR6) in the analytic Eq. 34 formula recovers the larger error bar increases (~5.5%) from the chain-based approach.

## Approach Taken
1. Loaded ACT DR4 CMB bandpowers and covariance (260 bins from 2 patches)
2. Extracted diagonal variance for TT, TE, EE spectra
3. Interpolated variance to full ell resolution
4. Applied Eq. 34 modified covariance formula using M matrices

## Critical Finding: Unit Mismatch

The implementation **failed due to a severe unit mismatch** between ACT DR4 data and the M matrix expectations:

### Scale of the Problem
```
DR6 approach:   cov_add/cov_original ~ 0.024  (2.4% increase) ✓
DR4 approach:   cov_add/cov_original ~ 10^12  (astronomical!) ✗
```

### Root Cause Analysis

1. **Lensing covariance** (covmat_act.txt):
   - Diagonal values: ~10^-18 to 10^-15
   - Units: Clκκ (lensing convergence power spectrum)

2. **ACT DR4 CMB variance**:
   - Diagonal values: ~10^-6
   - File format: likely Cl in dimensionless units or (µK)^2
   - **10^12 times larger than expected!**

3. **M matrices (dAL/dC)**:
   - Calibrated for Dl in specific units (µK^2)
   - Expect CMB variance in matching units
   - **Unit conversion from DR4 format was not applied**

### Data Format Issues

**ACT DR4** (`cl_cmb_ap.dat`):
```
ell=600.5:  Cl=0.0436,  sigma=0.0031
ell=650.5:  Cl=0.0290,  sigma=0.0021
```
- Values around 0.01-0.05 suggest Cl units (dimensionless or normalized)
- 260 bins from 2 independent patches (deep + wide)
- Fortran binary covariance matrix

**Expected by M matrices**:
- Dl in µK^2 (typically 1000-5000 µK^2 for TT at ell~600)
- Variance in (µK^2)^2 units
- Proper conversion: Dl = Cl × ell(ell+1)/(2π)

## Comparison with Reference

**DR6 analytic approach** (new_covmat.ipynb):
- Uses ACT DR6 SACC file with foreground-marginalized bandpowers
- Properly binned using Planck plik_lite covariance and weights
- Applies full 9-term cross-spectrum covariance
- **Result: 2.35% max diagonal increase**

**DR4 approach** (this notebook):
- Used diagonal-only interpolated variance
- Did not apply unit conversion
- **Result: 10^14% increase (unphysical)**

## Why DR6 Gives Smaller Increase than Chain-Based (~5.5%)

The key insight from project-context.md Section 11:

1. **Full covariance structure matters**:
   - Off-diagonal correlations within spectra
   - Cross-spectrum correlations (TT-TE, TT-EE, TE-EE)
   - Proper binning using CMB bandpower windows

2. **DR6 uses recent data with**:
   - Different noise levels
   - Different sky coverage
   - Different foreground treatment

3. **Chain-based approach** implicitly includes:
   - Full parameter space exploration
   - Non-linear effects
   - Different effective covariance from MCMC sampling

## Correct Path Forward

To properly verify if earlier data recovers the 5.5% increase:

### Option 1: Use Planck plik_lite Directly (Recommended)
```python
# Load Planck plik_lite covariance (already used in project)
planck_dir = '/scratch/jiaqu/likelihood_data/data/planck_2018_pliklite_native/'
cov_planck = np.loadtxt(planck_dir + 'c_matrix_plik_v22.dat')
weights = np.loadtxt(planck_dir + 'bweight.dat')

# Bin M matrices using Planck weights (see project-context.md)
# Apply full covariance with ell cuts (TT<1000, TE/EE<600)
# This matches the pattern in new_covmat.ipynb Section 11
```

**Advantages**:
- Known data format (used in project)
- Proper units and binning
- Full cross-spectrum covariance
- Well-documented in project-context.md

### Option 2: Fix ACT DR4 Unit Conversion
```python
# Convert ACT DR4 from Cl to Dl variance
ell = ell_act[idx]
factor = ell * (ell + 1) / (2 * np.pi)
var_Dl = var_Cl * factor**2  # Variance transforms as square

# But still need to verify absolute units (µK^2 vs dimensionless)
```

**Challenges**:
- Unclear absolute units in DR4 data
- Need to verify calibration factors
- May still have normalization issues

## Output Files

Generated (with caveats):
- `verify_dr4_planck_act_dr4_cov_structure.png` - ACT DR4 covariance structure (valid)
- `verify_dr4_planck_variance_interpolation.png` - Interpolated variance (valid plot)
- `verify_dr4_planck_covmat_comparison.png` - Comparison (shows unit issue)
- `verify_dr4_planck_diagonal_increase.npy` - Results (unphysical due to units)

## Conclusion

**The ACT DR4 approach as implemented does NOT provide a valid verification** due to critical unit mismatches.

**Key learnings**:
1. CMB data format and units vary between releases (DR4 vs DR6)
2. The M matrices (dAL/dC) are calibrated for specific units
3. Proper binning and unit conversion are essential
4. The 2.35% increase from DR6 analytic is likely correct
5. To understand the 5.5% chain-based result, need to investigate:
   - Which CMB dataset was used in chains
   - What covariance was actually applied
   - Possible non-linear effects in MCMC

**Recommended next step**: Implement Option 1 using Planck plik_lite covariance, following the exact pattern from new_covmat.ipynb Section 11, to properly test if Planck-only (without DR6) recovers the larger increase.
