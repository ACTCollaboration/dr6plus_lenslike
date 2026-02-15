# DR6+ Lensing Likelihood - Project Context

## Overview
Python package for ACT DR6+ CMB lensing likelihood computations. Used for cosmological parameter estimation with Cobaya.

## Key Files
- `dr6plus_lenslike/dr6plus_lenslike.py` — Main likelihood implementation
- `tests/test_dr6plus_lenslike.py` — Unit tests (uses mocked data)
- `base_config/` — YAML configuration files for Cobaya runs
- `runs/` — Run configurations

## Variants
Supported variants: `act_baseline`, `act_extended`, `actplanck_baseline`, `actplanck_extended`, `act_polonly`, `act_cibdeproj`, `act_cinpaint`, `spt3g`, `actspt3g_baseline`, `actspt3g_extended`, `actplanckspt3g_baseline`, `actplanckspt3g_extended`, `dr6plus_fiducial_baseline`, `dr6plus_fiducial_extended`, `day_baseline`, `day_extended`

## Two Operating Modes

### 1. Full likelihood corrections (`lens_only=False`)
- Uses regular covariance matrices
- Requires `tt, te, ee, bb, pp` spectra from theory
- Applies likelihood corrections via `get_corrected_clkk()`:
  - Normalization corrections (`dAL_dC`)
  - N1 bias corrections for both kk and CMB spectra

### 2. Lens-only mode (`lens_only=True`)
- Uses CMB-marginalized covariance matrices (`*_cmbmarg.txt`)
- Only requires `pp` spectrum from theory
- No likelihood corrections applied
- Alternative marginalization over 2pt CMB power spectra

## Current Work
- Implementing and testing the `lens_only=True` pathway
- Ensuring backward compatibility with existing functionality

## Data Files (v1.0)
- Bandpowers: `clkk_bandpowers_*.txt`
- Binning matrices: `binning_matrix_*.txt`
- Covariance matrices: `covmat_*.txt` (regular), `covmat_*_cmbmarg.txt` (CMB-marginalized)
- Likelihood corrections: `like_corrs/` subdirectory

## External Data

### ACT DR6 CMB-only
- SACC file: `/home/jiaqu/DR6-ACT-lite/act_dr6_cmbonly/data/act_dr6_cmb_sacc.fits`
- Contains foreground-marginalized TT, TE, EE bandpowers (ℓ≥600)
- Reference: `act_dr6_cmbonly.ACTDR6CMBonly` likelihood

### Planck plik_lite
- Location: `/scratch/jiaqu/likelihood_data/data/planck_2018_pliklite_native/`
- Files: `cl_cmb_plik_v22.dat`, `c_matrix_plik_v22.dat`, `blmin.dat`, `blmax.dat`
- Contains TT, TE, EE bandpowers (ℓ=30-2508)
- Reference: `act_dr6_cmbonly.PlanckActCut` likelihood

### Planck+ACT Combination (PlanckActCut)
ell cuts for combining Planck (low ℓ) with ACT (high ℓ):
- TT: Planck ℓ<1000, ACT ℓ≥600
- TE: Planck ℓ<600, ACT ℓ≥600
- EE: Planck ℓ<600, ACT ℓ≥600

Covariance cutting pattern (from `PlanckActCut.py`):
```python
cov[to_cut, :] = 0.0
cov[:, to_cut] = 0.0
cov[to_cut, to_cut] = 1e10
```

### Planck plik_lite Binning (from Cobaya `planck_pliklite.py`)
**Key insight**: All spectra (TT, TE, EE) use the **same** bin indices 0:nbin, not sequential slices:
- TT: blmin/blmax[0:215], ℓ=32-2492 (215 bins)
- TE: blmin/blmax[0:199], ℓ=32-1988 (199 bins)
- EE: blmin/blmax[0:199], ℓ=32-1988 (199 bins)

Binning operation (theory Dℓ → binned Cℓ):
```python
# weights include 2π/(ℓ(ℓ+1)) factor, so output is Cℓ
weights_raw = np.loadtxt('bweight.dat')
ls = np.arange(len(weights_raw)) + 30  # bin_lmin_offset=30
weights = np.hstack((np.zeros(30), weights_raw * 2*np.pi / (ls*(ls+1))))
# For bin i:
binned_Cl = np.dot(Dl[blmin[i]:blmax[i]+1], weights[blmin[i]:blmax[i]+1])
```

Data file `cl_cmb_plik_v22.dat` stores Cℓ (not Dℓ). Convert to Dℓ: `Dl = Cl * ell*(ell+1)/(2π)`

## Notebooks
- `notebooks/explore_dr6_cmbonly.ipynb` — Planck+ACT combined CMB power spectra exploration
  - Loads both Planck plik_lite and ACT DR6 data
  - Applies ell cuts following PlanckActCut pattern
  - Visualizes bandpowers and covariance matrices with cuts

## Equation (35): CMB-Marginalized Lensing Bandpowers

The lensing bandpower with CMB marginalization:
```
Ĉᵢᵠᵠ = Mᵢᵠ·ᴸ Cᴸᵠᵠ|fid + Mᵢˣ·ˡ (Ĉₗˣ - Cₗˣ|fid)
```

### Components

| Symbol | Description | Shape | File |
|--------|-------------|-------|------|
| Mᵢᵠ·ᴸ | Lensing binning matrix | (nbins, Lmax) | `binning_matrix_act.txt` |
| Mᵢˣ·ˡ | dAL/dC response matrix | (4, 4001, 3001) | `norm_correction_matrix_Lmin0_Lmax4000.npy` |
| Cᴸᵠᵠ\|fid | Fiducial clkk | (4001,) | `cosmo2017_10K_acc3_lenspotentialCls.dat` |
| Cₗˣ\|fid | Fiducial CMB Dℓ | (3001,) | `cosmo2017_10K_acc3_lensedCls.dat` |
| Ĉₗˣ | Measured Planck+ACT Dℓ | interpolated | From notebook |

### Likelihood Correction Files
Location: `/home/jiaqu/act_dr6_lenslike/data/v1.1/like_corrs/`
- `norm_correction_matrix_Lmin0_Lmax4000.npy` — dAL/dC for ACT (shape: 4, 4001, 3001)
- `P18_norm_correction_matrix_Lmin0_Lmax3000.npy` — dAL/dC for Planck
- `cosmo2017_10K_acc3_lensedCls.dat` — Fiducial CMB Dℓ (TT, EE, BB, TE)
- `cosmo2017_10K_acc3_lenspotentialCls.dat` — Fiducial φφ spectrum
- `n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax4000.txt` — Fiducial AL (normalization)
- `N1der_*.txt` — N1 bias derivatives

### Key Insights (verified 2026-01-07)

1. **dAL/dC expects Dℓ (not Cℓ)** in µK² units
2. **Spectrum order in dAL/dC**: [TT, EE, BB, TE] (index 0,1,2,3)
3. **ell range**: Input 0-3000, output L 0-4000
4. **Correction is significant only at low L** (L < 200)

### Correction Formula (from `get_corrected_clkk`)
```python
# For each spectrum X in [TT, EE, BB, TE]:
cldiff = Dl_measured - Dl_fiducial  # Dℓ units!
raw = dAL_dC[X] @ cldiff
norm_corr += -2 * raw / fAL

# Final corrected clkk:
nclkk = clkk_theory + norm_corr * clkk_fid + N1_corrections
```

### Test Results with Planck+ACT Data
Using combined Planck+ACT CMB vs fiducial:
- L=100: -0.32% correction (TT: -0.12%, EE: -0.07%, TE: -0.13%)
- L=500: ~0% (high L insensitive to CMB variations)
- L=1000: ~0%

The Planck+ACT measured CMB is ~8% lower than fiducial at low ℓ, causing the negative correction.

## Equation 35: Constant Correction for lens_only Mode (Section 10)

For `lens_only=True`, the lensing bandpowers need a constant correction computed from measured CMB:

$$\mathrm{const} = 2[d\ln R/dC] \Delta C + [dN_1/dC^j] \Delta C^j - [dN_1/dC^{\kappa\kappa}] C^{\kappa\kappa,\mathrm{fid}}$$

**Difference from Eq. A3 (Section 9):**
| Term | Eq. A3 | Eq. 35 |
|------|--------|--------|
| Normalization | `norm_corr * clkk_fid` | `norm_corr` (raw) |
| N1 κκ | Not computed | `-dN1_kk @ clkk_fid` |

**Computed values (L=100):**
- Normalization: -0.32% of clkk_fid (dominant)
- N1 CMB: ~-0.002% of clkk_fid
- N1 κκ: ~-0.01% of clkk_fid (NEW in Eq. 35)

**Output files (in notebooks/):**
- `eq35_constant_correction.txt`: L, norm_term, n1_cmb_term, n1_kk_term, total
- `eq35_constant_correction.npy`: NumPy dict with same data + clkk_fid

## Equation 34: Modified Covariance for lens_only Mode (Section 11)

For `lens_only=True`, the lensing covariance needs modification to account for CMB uncertainties:

$$\bar{\Sigma}_{ij} = \Sigma_{ij} + M_i^{X,\ell} \, \mathrm{cov}_{\mathrm{CMB}}^{X\ell;Y\ell'} \, M_j^{Y,\ell'}$$

where $M_i^{X,\ell} = -2 \frac{dA_L/dC^X}{f_{A_L}} + \frac{dN_1}{dC^X}$

### Implementation (notebook Section 11)

**Correct approach**: Transform M matrices to CMB bin space, then use binned Planck covariance directly.

1. **M matrix binning** (ℓ dimension):
   - M_binned has shape (18, 3001) after lensing L binning
   - Transform to (18, n_cmb_bins) using Planck weights:
   ```python
   M_doubly_binned[:, b] = sum(M_binned[:, ell] * weight[ell]) / sum(weight)
   ```

2. **Planck covariance blocks** (after ell cuts):
   - TT: 114 bins (ℓ ≤ 1000)
   - TE: 69 bins (ℓ ≤ 600)
   - EE: 69 bins (ℓ ≤ 600)
   - Full cross-covariances included (TT-TE, TT-EE, TE-EE)

3. **Covariance addition** (9 terms):
   ```python
   cov_add = M_TT @ cov_TT_TT @ M_TT.T + M_TT @ cov_TT_TE @ M_TE.T + ...
   ```

**Results:**
- Max diagonal increase: ~4.86% at lowest L bin
- Effect concentrated at low L (as expected)
- Modified covariance remains positive definite

**Output files (in notebooks/):**
- `eq34_modified_covariance.npy`: Modified covariance with metadata
- `eq34_modified_covariance.png`: Diagnostic plots

**Production files (in dr6plus_lenslike/data/v1.0/):**
- `covmat_act_cmbmarg_analytic.txt`: Full 18x18 modified covariance (same format as covmat_act.txt)
- `covmat_act_cmbmarg_analytic_info.npy`: Metadata including original/added covariance, diagonal increase %

**IMPORTANT**: The 2:-6 bin cuts are applied in dr6plus_lenslike.py during load_data(), NOT in the covariance file.

**IMPORTANT**: The previous interpolation-based approach was incorrect. Binned covariance should NOT be interpolated to full ℓ resolution.

### Section 12: Verification of Double-binning vs Interpolated Approaches

Notebook Section 12 compares two approaches for Eq. 34 modified covariance:

1. **Approach 1 (Double-binning)**: M matrices binned along ℓ using ACT bandpower windows, with full binned CMB covariance including cross-spectrum terms.

2. **Approach 2 (Interpolated)**: Unbinned M matrices with diagonal-only interpolated variance from ACT bins.

**Key finding**: Approach 1 is more rigorous as it properly accounts for full CMB covariance structure including cross-spectrum correlations. Approach 2 provides a simpler approximation but ignores off-diagonal covariance.

**Output files:**
- `eq34_approach_comparison.png`: Side-by-side comparison plot

## Window Function Format Differences: DR4 vs DR6

**Critical finding from 2026-02-12 investigation**:

ACT DR4 and DR6 use **fundamentally different window function representations** that require different mathematical treatments in Eq. 34 modified covariance calculations.

### DR6 (SACC format)
- **File**: `/home/jiaqu/DR6-ACT-lite/act_dr6_cmbonly/data/act_dr6_cmb_sacc.fits`
- **Format**: Compact, sparse windows
- **Bin support**: ~50 ells per bin (e.g., Bin 0: ℓ ∈ [626, 675])
- **Normalization**: ∑w_ℓ = 1.0 (exact)
- **Window shape**: (44 TT bins, 6501 theory ells)
- **CMB smoothness**: VALID assumption (Δℓ = 50)
- **Correct windowing**: Simple sum over bin support
  ```python
  M_doubly[:, b] = sum(M[:, ells_in_bin])  # ~50 ells
  ```
- **Result**: 2.35% max diagonal increase (correct)

### DR4 (pyactlike format)
- **File**: `/home/jiaqu/pyactlike/pyactlike/data/coadd_bpwf_*.npz`
- **Format**: Extended, dense windows
- **Bin support**: ~3000 ells per bin (e.g., Bin 0: ℓ ∈ [0, 7923], 3162 non-zero weights)
- **Normalization**: ∑w_ℓ ≈ 1.16 (approximately normalized)
- **Window shape**: (52 bins, 7924 theory ells)
- **CMB smoothness**: INVALID assumption (Δℓ = 3000)
- **Correct windowing**: Matrix multiplication
  ```python
  M_doubly = M @ window.T  # Weighted average over all ells
  ```
- **Result**: 0.025% max diagonal increase (correct for this format)

### Why Different Methods are Required

The "simple sum" method assumes:
1. All C_ℓ within a bin are approximately equal (smoothness)
2. When bin-averaged C changes by δC, each C_ℓ in the bin changes by δC
3. Total lensing response: ∑_{ℓ ∈ bin} M_ℓ

This is valid for DR6's compact 50-ell bins but **catastrophically fails** for DR4's extended 3000-ell bins.

**Tested scenarios**:
- DR6 + simple sum: 2.35% ✓ (correct)
- DR6 + M @ window.T: 0.000% ✗ (wrong - underestimates)
- DR4 + simple sum: 2.6 million % ✗ (wrong - massive overestimate)
- DR4 + M @ window.T: 0.025% ✓ (correct)

**Implication**: Cannot directly compare DR4 and DR6 Eq. 34 results. The window formats encode bin-averaging differently, leading to different numerical values even for identical underlying physics.

## 2026-02-15 Update: DR4 Binning Methodology — RESOLVED

### Source Code Analysis

Examined `pyactlike/like.py` (lines 195-213) and `notebooks/plot_spectra.ipynb`:

```python
# DR4 binning operation (from like.py)
# Step 1: Convert theory Dℓ → Cℓ
cltt[1:tt_lmax] = dell_tt / l_list / (l_list + 1.0) * 2.0 * np.pi

# Step 2: Apply window via matrix multiplication
cth_tt = win_func_w[2*bmax : 3*bmax, 1:lmax_win] @ cltt[1:lmax_win]
```

### DR4 vs DR6 Window Comparison

| Property | DR4 | DR6 |
|----------|-----|-----|
| Array shape | (n_bins, n_ell) | (n_ell, n_bins) |
| Normalization | sum(W) ≈ 1.00 | sum(W) = 1.00 |
| Effective bin width | ~50 ells (>0.01 threshold) | 50 ells exact |
| Weight distribution | Varying with small tails | Uniform within support |

**Key finding**: DR4 bins are ~50 ells wide (same as DR6!) when using threshold >0.01. The "170 ells" estimate used threshold >0.001 which includes negligible tails.

### Physical Derivation: Simple Sum is Correct

For Eq. 34 modified covariance `Σ̄ = Σ + M·Cov_CMB·M^T`:

1. **Binning**: `Ĉ_b = Σ_ℓ W_{b,ℓ} × C_ℓ` with `Σ W = 1`
2. **Smoothness**: If C_ℓ ≈ constant in bin, then `Ĉ_b = c`
3. **Fluctuation**: When `Ĉ_b → Ĉ_b + δĈ`, each C_ℓ shifts by δĈ
4. **Response**: `δĈ_κκ = Σ_ℓ M_ℓ × δC_ℓ = δĈ × Σ_ℓ M_ℓ`

**Conclusion**: `M_binned = Σ_ℓ M_ℓ` (simple sum over bin support)

The `M @ W.T` method computes a weighted average ≈ M at single ell, which is WRONG for Eq. 34.

### Bug Fix for verify_dr4_covmat_v2.ipynb

**Before (WRONG)**:
```python
M_doubly_dr4[key] = M_trunc @ window_trunc.T
```

**After (CORRECT)**:
```python
for b in range(n_cmb_bins):
    sig_ells = np.where(window_trunc[b, :] > 0.01)[0]
    M_doubly_dr4[key][:, b] = np.sum(M_trunc[:, sig_ells], axis=1)
```

### Test Results

| Notebook | Method | Result |
|----------|--------|--------|
| DR6 (new_covmat.ipynb) | Simple sum | 2.35% ✓ |
| DR4 (verify_dr4_covmat_v2.ipynb) | Simple sum | **394%** ✗ (way too large) |
| DR4 (verify_dr4_covmat_v2.ipynb) | M @ W.T | 0.04% (too small) |
| Chain-based reference | MCMC | ~5.5% |

### Two-Patch Issue

DR4 has two patches (deep + wide), which changes the Eq. 34 summation structure:

**Eq. 34**: $\bar{\Sigma}_{ij} = \Sigma_{ij} + \sum_{X,Y} \sum_{\beta,\beta'} \tilde{M}_i^{X,\beta} \text{Cov}(C_\beta^X, C_{\beta'}^Y) \tilde{M}_j^{Y,\beta'}$

| Dataset | Spectra | Covariance blocks |
|---------|---------|-------------------|
| DR6 | 3 (TT, TE, EE) | 9 |
| DR4 | 6 (TT/TE/EE × 2 patches) | 36 |

The 36 DR4 blocks include:
- 6 auto-blocks (TT_p1-TT_p1, etc.)
- 6 cross-patch same-spectrum (TT_p1-TT_p2, etc.)
- 24 cross-spectrum blocks

**Open question**: Should cross-patch covariances contribute? If patches observe different sky, their cosmic variance is independent.

## Next Steps
1. ~~Verify which windowing method is physically correct~~ — **DONE** (simple sum for DR6)
2. **[BLOCKED]** DR4 simple sum gives 394% — need to understand two-patch structure
3. ~~Planck plik_lite gives 305-3043%~~ — **RESOLVED**: Use M at bin center, Cℓ² cov, gives 1.2%
4. Test full `analytic_marg=True` pipeline end-to-end

## Planck Eq. 34 Implementation — RESOLVED (2026-02-15)

### Paper Quote (Planck lensing)
> "We evaluate the CMB power correction using the plik_lite band powers... To relate plik_lite bins to the M_i^{X,ℓ'} bins, we assume that the underlying CMB power spectra are represented only by modes that are smooth over Δℓ = 50. The plik_lite bandpower covariance cov_CMB is similarly used to calculate Eq. (34). The increase in the diagonal of the covariance is about 6% at its largest."

### Root Cause Analysis (debug_planck_eq34.py)

The original bug (305-3043%) came from **TWO incorrect approaches**:

1. **Unit conversion bug**: Converting Planck covariance from Cℓ² to Dℓ² multiplies by (ℓ(ℓ+1)/(2π))² ≈ 10^10

2. **M summation bug**: Summing M over ~50 ells (per "smoothness" assumption) multiplies result by N² where N is bin width

### Correct Approach: M at bin center, NO unit conversion

| Method | Result | Status |
|--------|--------|--------|
| M at bin center, Cℓ² cov (H1) | **1.2%** | ✓ CORRECT |
| M summed, Cℓ² cov (H7) | 305% | ✗ Too high |
| M summed, Dℓ² cov (H4) | 10^14% | ✗ WAY too high |
| Super-bins, Dℓ² cov (H5) | 10^13% | ✗ WAY too high |

### Why H1 Works

1. **M at bin center**: For narrow bins (~5-10 ells), using M at the bin center ell gives correct response
2. **Cℓ² covariance directly**: plik_lite covariance is stored in Cℓ² units, and using it directly (without Dℓ conversion) gives reasonable results
3. **No double-counting**: Using M at center avoids multiplying by bin width

### Gap Between 1.2% and Paper's 6%

The 5× gap is explained by:
- Our M matrices use ACT lmin=600, missing ℓ < 600 contribution
- Planck paper likely uses Planck-specific M matrices with lmin=100
- Low ℓ contributes more to lensing normalization correction (dAL/dC peaks at low ℓ)

### Recommended Implementation

```python
# For each plik_lite bin b with center ell_b >= 600:
M_doubly[:, b] = M_binned[:, ell_center_b]  # M at bin center, NOT summed

# Use Cℓ² covariance directly (NO unit conversion!)
cov_add = M_doubly @ cov_planck_Cl2 @ M_doubly.T
```

### Debug Files
- `notebooks/debug_planck_eq34.py` - Systematic hypothesis testing (10 approaches)
- `notebooks/verify_planck_covmat.ipynb` - Original (broken) attempts
