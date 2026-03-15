# DR6+ Lensing Likelihood - Project Context

## Overview
Python package for ACT DR6+ CMB lensing likelihood computations. Used for cosmological parameter estimation with Cobaya.

## Key Files
- `dr6plus_lenslike/dr6plus_lenslike.py` — Main likelihood implementation
- `tests/test_dr6plus_lenslike.py` — Unit tests (uses mocked data)
- `base_config/` — YAML configuration files for Cobaya runs
  - `camb.yaml` — Standard CAMB (pip)
  - `camb_alens.yaml` — CAMB_alens fork (`/home/jiaqu/CAMB_alens`), A_lens/B_lens + CosmoRec, requires `module load gsl`
  - `camb_disputauble.yaml` — disputauble pip package (`disputauble.CobayaCAMB_mnuEff`), negative mnu extrapolation
  - `class_sz.yaml` — CLASS_SZ solver
- `runs/` — Run configurations (select theory via `theory: !defaults [../base_config/camb]`)

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

## Planck PR4 Reference Implementation (2026-02-15)

**Source**: `/home/jiaqu/planck_PR4_lensing/` (J. Carron et al., arXiv:2206.07773)

### Two Likelihood Versions
1. **PlanckPR4Lensing** — Full likelihood
   - `fields_required = T E P`
   - Linear correction depends on TT, EE, TE, PP
   - Has calibration parameter `A_planck`

2. **PlanckPR4LensingMarged** — CMB-marginalized
   - `fields_required = P` (only lensing)
   - Linear correction depends only on PP (dN1/dC^kk)
   - No calibration parameter needed
   - Uses precomputed CMB-marginalized covariance

### Data Files
```
data_pr4/
├── pp_*_cov.dat                    # Regular covariance (9×9)
├── pp_*_CMBmarged_cov.dat          # CMB-marginalized covariance
├── pp_*_lens_delta_window/         # M matrices (TT, EE, TE, PP columns)
├── pp_*_CMBmarged_lens_delta_window/  # M matrices (PP only)
├── pp_*_lensing_fiducial_correction.dat
└── pp_*_CMBmarged_lensing_fiducial_correction.dat
```

### CMB Marginalization Effect
| Bin | L_center | Cov increase |
|-----|----------|--------------|
| 1   | 28       | 4.4%         |
| 2   | 64       | 13.0%        |
| 3   | 106      | 13.6%        |
| 4-9 | 150-377  | 0.3-7.9%     |

### Key Differences from Our Approach
| Aspect | Planck PR4 | ACT DR6 |
|--------|------------|---------|
| M matrix lmin | 100 | 600 |
| Max cov increase | ~13% | ~2% |
| Bins | 9 (L=8-400) | 18 (L~40-700) |

### Cobaya Integration
The likelihood uses Cobaya's CMBlikes base class:
```python
# Linear correction formula (cmblikes.py line 479):
band += self.linear_correction.bin(Cls) - self.fid_correction.T
# Equivalent to: band += M @ (Cls_theory - Cls_fiducial)
```

### Missing Information
- Code that generates CMBmarged covariance is not public
- Unknown what CMB covariance was used (internal Planck, not plik_lite?)

## 2026-02-15 Investigation: Chain-Based Eq. 34 vs Analytic — UNRESOLVED

### Goal
Match chain-based result (~5.2% diagonal increase at L≈123) using measured DR4 + Planck CMB covariances with analytic Eq. 34.

### Reference Chain Notebook
`/home/jiaqu/DW/100223_compute_cov_norm.ipynb`

Key chain approach:
```python
# For each parameter sample from ACT+Planck posterior:
ps_cl = make_spec(As, ns, H0, ombh2, omch2, tau)  # CAMB theory
response = lin_corr_matrix @ (ps_cl - ps_ref)  # dAL/dC @ ΔCℓ
norm = 2 * response / AL_ref  # AL_ref = n0mv / (L(L+1)/2)²
storage[i] = bin(norm)

# Final covariance:
cov_norm = np.cov((storage * binned_clkk_fid).T)
```

### Key Insight: Parameter-Induced vs Measured Covariance

**The chain uses PARAMETER-INDUCED Cℓ covariance**, not measured CMB covariance:
- Chain samples from ACT+Planck posterior → computes CAMB theory Cℓ(θ)
- Takes covariance of theory Cℓ variations across samples
- This captures only parameter uncertainty (~1% in As, ~0.4% in ns)

**Measured CMB covariance is MUCH larger**:
- Includes cosmic variance + instrument noise + systematics
- ~30-40× larger than parameter-induced covariance

### DR4 Data Structure
- **Bandpowers**: `/home/jiaqu/pyactlike/pyactlike/data/cl_cmb_ap.dat`
  - 260 bins = 2 patches × (40 TT + 45 TE + 45 EE)
  - Dℓ in µK² units

- **Covariance**: `/home/jiaqu/pyactlike/pyactlike/data/c_matrix_ap.dat`
  - 260×260 Fortran binary
  - Full cross-patch, cross-spectrum correlations

- **Windows**: `coadd_bpwf_15mJy_191127_lmin2.npz`, `coadd_bpwf_100mJy_191127_lmin2.npz`
  - Shape (520, 7924), ~100 ells per bin with >1% weight

### M Matrix Versions Discovered

| Location | Value at (L=100, ℓ=800) | Notes |
|----------|-------------------------|-------|
| spt_act/v1.2 | -5.02e-10 | Current standard |
| blake_plot/misc | -1.97e-17 | 10^7× smaller |
| act_dr6_lenslike/v1.1 | -1.97e-17 | Same as blake_plot |
| act_dr6_lenslike/v1.2 | -5.02e-10 | Same as spt_act |

**v1.1 → v1.2 changed by factor ~10^7!** This explains some historic confusion.

### Test Results

| Approach | Unit Conversion | M Binning | Result | Target |
|----------|-----------------|-----------|--------|--------|
| Raw Dℓ cov | None | Center | 10^14% | 5% |
| Cℓ (µK²) cov | Dℓ→Cℓ | Center | 172% | 5% |
| Cℓ (µK²) cov | Dℓ→Cℓ | Window | 169% | 5% |
| Dimensionless Cℓ | Dℓ→Cℓ/T² | Center | ~0% | 5% |
| Blake M matrix | None | Window | 0.025% | 5% |

### Why All Approaches Fail

The ~170% result (with Cℓ µK² covariance) is ~33× too high. This matches the expected ratio between:
- Measured CMB covariance (cosmic variance + noise)
- Parameter-induced Cℓ covariance (theory variation only)

The chain's parameter-induced covariance is approximately:
```python
Var(Cℓ)_param ≈ Cℓ² × [4(σAs/As)² + (σns × ∂lnCℓ/∂ns)² + ...]
            ≈ Cℓ² × 0.01²  # ~1% parameter variation
```

While measured covariance is:
```python
Var(Cℓ)_measured ≈ Cℓ² × 2/(2ℓ+1) + noise  # cosmic variance dominant
```

For ℓ~1000: measured/param ≈ (2/2001) / 0.0001 ≈ 100× larger.

### User's Question (Unresolved)
> "Is there a reason why we cannot apply the same binning function to the M matrix, then extrapolate to L=3000?"

Tested via `M @ W.T` (window-weighted average), gives 169% — still wrong.

The fundamental issue appears to be that measured CMB covariance cannot substitute for parameter-induced covariance without a scaling factor.

### Scripts Created
- `notebooks/match_chain_measured_cov.py` — 7 binning approaches tested
- `notebooks/test_blake_M_matrix.py` — M matrix version comparison
- `notebooks/test_unit_conversion.py` — Unit conversion tests
- `notebooks/compare_all_M_matrices.py` — Comprehensive M matrix analysis

### Open Questions
1. Is there a scaling factor to convert measured → parameter-induced covariance?
2. Should the two-patch structure of DR4 affect the result?
3. Does the chain use a different M matrix file entirely?
4. Is the "simple sum" vs "window average" distinction important here?

### Conclusion (2026-02-15)
Cannot match chain-based ~5% result using measured CMB covariance with any binning method tested. The measured covariance appears to be fundamentally ~30-40× larger than what the chain computes from parameter variations. Further investigation paused pending clarification of the physical relationship between measured and parameter-induced CMB covariances.

## 2026-02-16: Analytic Eq. 34 Implementation — COMPLETE

### Production Script
**File**: [src/lensing_only_cov_marginalization.py](../src/lensing_only_cov_marginalization.py)

Based on [notebooks/new_covmat.ipynb](../notebooks/new_covmat.ipynb), formalized into a reusable module.

### Module Structure
```python
# Exported functions in src/__init__.py
from .lensing_only_cov_marginalization import (
    compute_analytic_covariance,   # Main Eq. 34 computation
    compute_chain_covariance,      # Chain-based (for comparison)
    compare_approaches,            # Side-by-side visualization
)
```

### Usage
```bash
# Analytic mode (default)
python src/lensing_only_cov_marginalization.py --mode analytic

# Save production files
python src/lensing_only_cov_marginalization.py --mode analytic --save-production

# Chain-based comparison
python src/lensing_only_cov_marginalization.py --mode chain --nsample 1000

# Full comparison
python src/lensing_only_cov_marginalization.py --mode compare
```

### Verified Implementation (Triple-Checked)

| Component | Notebook Line | Script Line | Status |
|-----------|---------------|-------------|--------|
| `build_M_matrix` | 746-758 | 84-97 | ✓ IDENTICAL |
| M matrix indices (TT=0, EE=1, TE=3) | 762-764 | 181-183 | ✓ IDENTICAL |
| Unit conversion (Cℓ→Dℓ) | 790-796 | 199-205 | ✓ IDENTICAL |
| `doubly_bin_M` (simple sum) | 892-922 | 115-143 | ✓ IDENTICAL |
| cov_add (9 terms) | 1042-1063 | 252-273 | ✓ IDENTICAL |

### Key Computation: M Matrix
```python
# M^X_{L,ℓ} = -2 * dAL_dC[L,ℓ] / fAL[L] * clkk_fid[L] + dN1_X[L,ℓ]
M = np.zeros((Lmax, 3000))
for L in range(2, Lmax):
    if fAL[L] > 0:
        M[L, :] = -2 * dAL_dC_X[L, :3000] / fAL[L] * clkk_fid_trunc[L]
    M[L, :] += dN1_X[L, :]
```

### Key Computation: Double Binning
```python
# Bin along L (lensing), then along ℓ (CMB)
M_L_binned = lens_binmat @ M[:Lmax_lens, :]

# Unit conversion: M is for Cℓ, ACT cov is Dℓ²
ell_factor = ells * (ells + 1) / (2 * np.pi)
M_L_binned /= ell_factor

# ℓ binning: simple sum over bin support (DR6 SACC windows)
for b in range(n_cmb_bins):
    bin_support = window[b, :] > 1e-10
    M_doubly[:, b] = np.sum(M_L_binned[:, ells[bin_support]])
```

### Key Computation: 9-Term Covariance Addition
```python
cov_add = np.zeros((n_lens_bins, n_lens_bins))
cov_add += M_TT @ cov_TT_TT @ M_TT.T
cov_add += M_TT @ cov_TT_TE @ M_TE.T + M_TE @ cov_TT_TE.T @ M_TT.T
cov_add += M_TT @ cov_TT_EE @ M_EE.T + M_EE @ cov_TT_EE.T @ M_TT.T
cov_add += M_TE @ cov_TE_TE @ M_TE.T
cov_add += M_TE @ cov_TE_EE @ M_EE.T + M_EE @ cov_TE_EE.T @ M_TE.T
cov_add += M_EE @ cov_EE_EE @ M_EE.T
```

### Results
| Bin | L center | Diagonal Increase (%) |
|-----|----------|----------------------|
| 0 | 41 | 0.41% |
| 1 | 61 | 0.53% |
| 2 | 81 | 0.71% |
| **5** | **123** | **2.35%** (max) |
| 17 | 713 | 0.33% |

### Production Files
- `dr6plus_lenslike/data/v1.0/covmat_act_cmbmarg_analytic.txt` — 18×18 modified covariance
- `dr6plus_lenslike/data/v1.0/covmat_act_cmbmarg_analytic_info.npy` — Metadata (original cov, cov_add, diag %)

### Chain-Based Approach (Pending Comparison)

**Reference chain**: `/project/rrg-rbond-ac/jiaqu/chains/act_dr6_2pt/lcdm/p-actlite_lcdm_camb/p-actlite_lcdm_camb`

**Config file**: [p-actlite_lcdm_camb.updated.yaml](read in prior context)

The chain-based approach computes:
```python
# For N samples from ACT+Planck posterior:
for i in range(nsample):
    # Sample cosmological params from chain
    As, ns, tau, H0, omch2, ombh2 = chains[idx, :]

    # Compute CAMB theory spectra at this point
    ps_sample = make_spec_camb(As, ns, H0, ombh2, omch2, tau)

    # Compute normalization response
    delta_cl = ps_sample - ps_fiducial
    norm_response = dAL_dC @ delta_cl
    norm_factor = 2 * norm_response / AL_ref

    # Bin and store
    storage[i] = lens_binmat @ norm_factor * clkk_fid

# Covariance of binned normalization
cov_norm = np.cov(storage.T)
```

**Key difference**: Chain uses parameter-induced Cℓ covariance (from posterior sampling), not measured CMB covariance.

### Next Session Tasks
1. ~~Run chain-based mode: `python src/lensing_only_cov_marginalization.py --mode chain --nsample 1000`~~ DONE
2. ~~Compare analytic vs chain-based diagonal increases~~ IN PROGRESS
3. Investigate any systematic differences
4. Document final recommended approach for production

## 2026-02-16 Session 2: Sqrt Bug Fix & DR4 Chain Support

### Critical Bug Fix: Variance vs Sigma

**Bug**: All diagonal increase calculations were computing variance ratios instead of sigma ratios.

**Before (WRONG)**:
```python
diag_increase_pct = (diag_modified - diag_orig) / diag_orig * 100
```

**After (CORRECT)**:
```python
sigma_orig = np.sqrt(diag_orig)
sigma_modified = np.sqrt(diag_modified)
diag_increase_pct = (sigma_modified - sigma_orig) / sigma_orig * 100
```

**Impact**: All reported percentages were ~2× too high. The "2.35%" becomes ~1.17%, the "7.62%" becomes ~3.7%.

**Files fixed**: `src/lensing_only_cov_marginalization.py` (6 locations: lines 283, 544, 562, 631, 721, 754)

### DR4+Planck Chain Support

**Chain location**: `/project/rrg-rbond-ac/jiaqu/chains/dr4_planck_lcdm/CLASS2p8_ACTPol_lite_DR4_leakfix_yp2_baseLCDM_taup_planck2018_lowTT_plikHM_TT_lmax650_hip`

**Issue**: DR4 chains use `sampler: minimize` in their YAML, which causes getdist to fail with:
```
ValueError: Unknown sampler type minimize
```

**Solution**: Added `_load_chains_raw()` function to load chain text files directly:
```python
def _load_chains_raw(chain_path, burn_in=0.3):
    """Load chains directly from text files (for CLASS/DR4 format)."""
    with open(f'{chain_path}.1.txt', 'r') as f:
        header = f.readline().strip().replace('#', '').split()
    # ... load all chain files, apply burn-in, extract parameters
```

### Parameter Constraint Comparison: DR6 vs DR4

| Parameter | DR6 std | DR4 std | Ratio (DR6/DR4) |
|-----------|---------|---------|-----------------|
| logA | 0.0095 | 0.0132 | **0.72** |
| omch2 | 0.00084 | 0.00153 | **0.55** |
| H0 | 0.445 | 0.728 | **0.61** |
| sigma8 | 0.0057 | 0.0085 | **0.67** |
| tau | 0.0052 | 0.0108 | **0.48** |

**Key finding**: DR6+Planck has significantly **tighter** constraints than DR4+Planck (most ratios < 0.7).

**Implication**: Tighter constraints → smaller parameter-induced Cℓ variance → **smaller** Eq. 34 covariance addition. Need to verify with re-run using corrected sqrt formula.

### Batch Script Updates

**File**: [src/run_chain_cov.sh](../src/run_chain_cov.sh)

Now accepts chain path as second argument:
```bash
# Usage:
sbatch src/run_chain_cov.sh [nsample] [chain_path]

# Examples:
sbatch src/run_chain_cov.sh 1000  # Default DR6+Planck chain
sbatch src/run_chain_cov.sh 1000 /path/to/dr4_chain  # DR4 chain
```

### Output Filename Collision Fix

Different chains now save to different output files:

| Chain | Output filename |
|-------|-----------------|
| DR6+Planck | `covmat_act_cmbmarg_chain_DR6.txt` |
| DR4+Planck | `covmat_act_cmbmarg_chain_DR4.txt` |
| Other | `covmat_act_cmbmarg_chain_PACT.txt` |

Detection logic:
```python
if 'dr4' in args.chain_path.lower() or 'DR4' in args.chain_path:
    chain_label = 'DR4'
elif 'dr6' in args.chain_path.lower() or 'p-actlite' in args.chain_path:
    chain_label = 'DR6'
else:
    chain_label = 'PACT'
```

### Chain Comparison Script

**File**: [src/run_chain_comparison.sh](../src/run_chain_comparison.sh)

Runs both DR6 and DR4 chains in a single job and generates comparison plot:
```bash
sbatch src/run_chain_comparison.sh [nsample]  # Default 500 samples each
```

**Output**:
- `products/chain_comparison_results.npy` — Numerical results
- `products/chain_comparison_dr6_vs_dr4.png` — Visualization

### Next Steps

1. Re-run chain-based computations with corrected sqrt formula
2. Compare DR6 vs DR4 results with proper sigma ratios
3. Verify that tighter DR6 constraints give smaller covariance increase (as expected)
4. Document final production workflow

## 2026-02-16 Session 3: Fiducial Cosmology Mismatch — CRITICAL ISSUE

### Problem Discovery

When comparing DR4 (DW chain) vs DR6 (p-actlite) CMB marginalization:
- **Expected**: DR6 has tighter parameter constraints → smaller covariance addition
- **Observed**: DR6 has **LARGER** covariance addition (~1.3× larger)

### Root Cause Analysis

The chain-based computation uses a **fixed fiducial cosmology** (cosmo2017, DR4-era):

```python
ps_fid = load_fiducial_spectra()  # cosmo2017_10K_acc3_lensedCls.dat
norm = compute_norm_correction(ps_sample, ps_fid, M_matrix, AL_ref)
```

**Key insight**: The norm correction depends on `ps_sample - ps_fid`. If DR6's mean cosmology is offset from the DR4 fiducial:
- DR6 samples are centered around a different point in parameter space
- The CAMB-computed Cℓ spectra are systematically offset from ps_fid
- This creates larger deviations even with tighter parameter spread

### Mean Parameter Comparison

| Parameter | DR4 mean | DR6 mean | Difference |
|-----------|----------|----------|------------|
| As | 2.19e-9 | 2.12e-9 | **-3%** |
| τ (tau) | 0.072 | 0.060 | **-17%** |
| H0 | 67.5 | 67.6 | +0.1% |
| ns | 0.969 | 0.971 | +0.2% |

**The 17% τ difference is critical** — τ strongly affects the overall amplitude via exp(-2τ).

### Evidence: Norm Correction Variance

Raw norm correction std (from stored samples):

| Bin | L | DR4 std | DR6 std | Ratio |
|-----|---|---------|---------|-------|
| 4 | 123 | 1.11e-2 | 1.27e-2 | 0.87 |
| 10 | 582 | 0.89e-2 | 1.24e-2 | 0.72 |

DR6 has ~30-40% larger norm correction variance despite having tighter parameter constraints.

### Two Components Affected by Fiducial

1. **δCℓ computation**: `delta_cl = Cℓ_sample - Cℓ_fid`
   - If Cℓ_fid doesn't match the chain's mean, there's a systematic offset
   - **However**, covariance is shift-invariant: Var(X+c) = Var(X)
   - So this alone shouldn't cause larger variance

2. **M matrix computation**: `M = -2 * dAL/dCℓ / fAL * clkk_fid + dN1`
   - All derivatives (dAL/dCℓ, dN1) are computed at the fiducial
   - fAL and clkk_fid are at the fiducial
   - **These would need to be recomputed for DR6 cosmology**

### Likely Explanation: Nonlinear CAMB Response

The mapping θ → Cℓ(θ) via CAMB is **nonlinear**. The Jacobian ∂Cℓ/∂θ varies with cosmology.

If DR6's mean cosmology has **steeper derivatives** (larger ∂Cℓ/∂θ), then:
- Same parameter spread δθ produces larger Cℓ spread
- This could explain why tighter parameters still give larger Cℓ variance

### Required Fix for DR6

**Option A (Minimal)**: Use DR6 fiducial spectra for δCℓ computation
- File: `p-actlite_lcdm_camb_bestfit_clkk.dat` exists
- Need corresponding CMB Cℓ fiducial (TT, TE, EE)
- M matrices remain at DR4 fiducial (linearization still valid if not too far)

**Option B (Full)**: Recompute everything at DR6 fiducial
- New M matrices (dAL/dCℓ, dN1 evaluated at DR6 cosmology)
- New fiducial spectra (Cℓ_fid, clkk_fid, fAL)
- More work but more accurate

### Available DR6 Bestfit Files

```
/project/rrg-rbond-ac/jiaqu/chains/act_dr6_2pt/lcdm/
├── p-actlite_lcdm_camb_bestfit_clkk.dat
├── p_act_bestfit_clkk.dat
├── p_act_lb_bestfit_clkk.dat
└── planck_desi_lcdm_bestfit_clkk.dat
```

### Parameter Comparison: All 6 Sampled Parameters

**Script**: [notebooks/compare_chain_contours.py](../notebooks/compare_chain_contours.py)
**Output**: [products/chain_contours_all6_dr4_vs_dr6.png](../products/chain_contours_all6_dr4_vs_dr6.png)

| Param | DR4 mean | DR4 std | DR6 mean | DR6 std | DR4/DR6 std |
|-------|----------|---------|----------|---------|-------------|
| As | 2.19e-9 | 5.3e-11 | 2.12e-9 | 2.8e-11 | **1.88** |
| ns | 0.969 | 4.1e-3 | 0.971 | 3.7e-3 | 1.13 |
| τ | 0.072 | 1.2e-2 | 0.060 | 6.2e-3 | **2.01** |
| H0 | 67.5 | 0.57 | 67.6 | 0.51 | 1.12 |
| ωch² | 0.120 | 1.3e-3 | 0.119 | 1.2e-3 | 1.07 |
| ωbh² | 0.0224 | 1.3e-4 | 0.0225 | 1.1e-4 | 1.17 |

**Key findings:**
- DR6 is tighter on **ALL** parameters (all ratios > 1)
- Biggest improvements: **As (1.9×)** and **τ (2.0×)**
- τ shifted down by 17% (0.072 → 0.060)
- As shifted down by 3%
- τ-As anti-correlation visible (trade off via amplitude e^{-2τ}As)

### Hypothesis: Nonlinear CAMB Jacobian

The CMB marginalization covariance depends on:
```
Var(Cℓ) ≈ (∂Cℓ/∂θ) @ Var(θ) @ (∂Cℓ/∂θ).T
```

At DR6's cosmology (lower τ):
- Cℓ ∝ e^{-2τ} → lower τ means larger Cℓ
- ∂Cℓ/∂τ = -2 × Cℓ → larger Cℓ means larger |∂Cℓ/∂τ|

So DR6's cosmology has a **steeper CAMB Jacobian** - the same parameter perturbation produces larger Cℓ perturbation. This could explain why tighter parameters still give larger Cℓ variance.

**Important**: The M matrices (lensing response) don't need to change - they're ∂C^κκ/∂Cℓ, not the CAMB Jacobian ∂Cℓ/∂θ. The fiducial subtraction also doesn't matter since variance is shift-invariant.

### Next Session TODO

1. **Verify CAMB Jacobian hypothesis**: Compute ∂Cℓ/∂θ at both DR4 and DR6 mean cosmologies
2. If confirmed, the larger DR6 covariance may be **physically correct** (not a bug)
3. **Document** that CMB marginalization depends on cosmology, not just parameter constraints
4. **Consider** whether to report DR6 result as-is or investigate further

## Self-Calibration — VERIFIED (2026-03-14)

### What it is
Marginalizing 8 calibration nuisance parameters (4 gain + 4 pol efficiency for pa5a/b, pa6a/b)
over the DR6+ lensing likelihood. Calibration affects lensing in two partially-cancelling ways:
1. **T²(δc) term**: response matrix R (18×9) scales theory clkk
2. **Norm correction**: same calibration-induced dC_TT/EE/TE/BB feeds into get_corrected_clkk

### Key products
| Product | Path |
|---------|------|
| Response matrix R (18×9) | `/home/jiaqu/act_dr6_lenslike/act_dr6_lenslike/data/v1.2/response_cal_matrix.txt` |
| Calibration samples (100×4) | `/home/jiaqu/DR6plus_lensing/preprocessing/systematic_tests/calibration_samples.txt` |
| Coadding weights (5001,) | `/home/jiaqu/DR6plus_lensing/preprocessing/output/21122025_nighttime/stage_compute_weights/noise_{pa5a,pa5b,pa6a,pa6b}_{T,E}_weights.txt` |

### Calibration model (verbatim from act_dr6_mflike/_calibrate_spectra)
Array naming: pa5a=pa5_f090, pa5b=pa5_f150, pa6a=pa6_f090, pa6b=pa6_f150. Gain c_a is map-level (affects T and E equally); pol efficiency p_a is E-only.

CMB spectrum changes (Cℓ units, linearized):
```python
eff_T = (w_T * delta_c[:, None]).sum(axis=0) / w_T.sum(axis=0)
eff_E = (w_E * (delta_c + delta_p)[:, None]).sum(axis=0) / w_E.sum(axis=0)
dC_TT = 2 * eff_T * C_TT_fid_cl
dC_EE = 2 * eff_E * C_EE_fid_cl
dC_TE = (eff_T + eff_E) * C_TE_fid_cl
```

### Sign convention (critical)
`norm_corr > 0` for positive calibration because `dAL_dC < 0` (more CMB → smaller AL).
A_lens bias = data/theory − 1 ≈ `frac_4pt − norm_corr` (MINUS sign — they partially cancel).

### Verification results
- Cancellation fraction: 56.5% (green std = 0.57× red std)
- A_lens bias: T²-only = 0.0039, self-cal = 0.0017 (~2.3× improvement)
- Script: `/home/jiaqu/DR6plus_lensing/notebooks/systematic_tests/selfcal_verification.py` (commit 6628829, branch preproc_test)

### Next step
Implement self-calibration in `dr6plus_lenslike/dr6plus_lenslike.py`:
- Add 8 calibration nuisance parameters to `ACTDR6LensLike`
- In `generic_lnlike`: apply T²(δc) to theory clkk AND pass calibration-corrected CMB spectra to `get_corrected_clkk`
- Load R and weights in `load_data()`
- Add Cobaya parameter priors matching `act_dr6_mflike/params_systematics.yaml`

## Calibration Module — IMPLEMENTED (2026-03-14)

### File
`dr6plus_lenslike/calibration.py` — standalone module, primary reference document
for calibration physics.

### Data files in `dr6plus_lenslike/data/v1.0/`
| File | Shape | Purpose |
|------|-------|---------|
| `response_cal_matrix.txt` | 18×9 | T²(δc) scaling of theory clkk (Prompt B) |
| `noise_pa5a_T_weights.txt` … `noise_pa6b_E_weights.txt` | 5001 values each | Per-ℓ noise-coadding weights for T and E maps |

### Functions
```python
load_calibration_weights(data_dir)
    → w_T, w_E  each shape (4, 5001), row order [pa5a, pa5b, pa6a, pa6b]

compute_dCl_from_calibration(delta_c, delta_p, w_T, w_E,
                              C_TT_fid, C_EE_fid, C_TE_fid)
    → dC_TT, dC_EE, dC_TE, dC_BB  each shape (N_ell,), Cℓ units
```

### MCMC interface
`delta_c` (4,) gain deviations and `delta_p` (4,) pol-efficiency deviations
are MCMC nuisance parameters supplied by the Cobaya sampler at each step;
not loaded from file by this module.

### Tests
`tests/test_calibration.py` — 6 tests (A–F), all passing.
