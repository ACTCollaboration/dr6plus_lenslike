# Agent Manager Status

## Current Focus
Foreground-marginalized CMB bandpowers integration for `lens_only` mode.

## Session History

### 2026-01-06: Initial setup
- Initialized agent management files
- Reviewed codebase structure and existing `lens_only` implementation
- Identified key code paths and existing test coverage

### 2026-01-06: Foreground marginalization focus
- User clarified: focus on foreground marginalized step
- Source data: `/home/jiaqu/DR6-ACT-lite/` contains DR6 foreground-marginalized CMB bandpowers
- Task: Read README, extract bandpowers and fiducial data as needed

### 2026-01-07: Planck+ACT cut version notebook
- Updated `notebooks/explore_dr6_cmbonly.ipynb` to use Planck+ACT combined data
- Installed Planck plik_lite data via `cobaya-install act_dr6_cmbonly.PlanckActCut`
- **Data sources:**
  - ACT DR6: `/home/jiaqu/DR6-ACT-lite/act_dr6_cmbonly/data/act_dr6_cmb_sacc.fits`
  - Planck plik_lite: `/scratch/jiaqu/likelihood_data/data/planck_2018_pliklite_native/`
- **ell cuts (following PlanckActCut.yaml):**
  - TT: Planck ℓ<1000, ACT ℓ≥600
  - TE: Planck ℓ<600, ACT ℓ≥600
  - EE: Planck ℓ<600, ACT ℓ≥600
- Implemented covariance cutting following `PlanckActCut.py` pattern:
  - Zero out rows/columns for cut bins
  - Set diagonal to 1e10 (infinite variance)
- Refactored to use shared `ELL_CUTS` dict for consistency between bandpowers and covariance
- Generated combined spectra plot: `planck_act_combined_spectra.png`

### 2026-01-07: Fixed Planck binning implementation
- **Bug identified**: Notebook used sequential blmin/blmax slices for TT/TE/EE (wrong)
- **Fix**: All spectra use same bin indices 0:nbin (following Cobaya `planck_pliklite.py`)
- **Key insight from Cobaya**:
  - TT bins 0-214 use blmin/blmax[0:215], ℓ=32-2492
  - TE bins 0-198 use blmin/blmax[0:199], ℓ=32-1988
  - EE bins 0-198 use blmin/blmax[0:199], ℓ=32-1988
- **Binning operation**: `weights[ell]` indexed by ell, includes 2π/(ℓ(ℓ+1)) factor
- Commit 843430c

### 2026-01-07: Aligned bandpower cut with covariance cut
- **Issue**: Bandpower used `ell_avg < lmax`, covariance used `blmax > lmax` (1-bin mismatch)
- **Fix**: Changed bandpower cut to use `blmax <= lmax` (matches PlanckActCut.py)
- **Result**: Perfect consistency - TE/EE now 69 bins in both (was 70 bandpower, 69 covariance)
- Commit ab6d42a

### 2026-01-07: Analyzed equation (35) likelihood corrections
- **Goal**: Understand how to construct CMB-marginalized lensing bandpowers
- **Key finding**: dAL/dC matrix expects **Dℓ** (not Cℓ) in µK² units
- **Spectrum order**: [TT, EE, BB, TE] in dAL/dC array
- **Test result**: Planck+ACT vs fiducial gives -0.32% correction at L=100
- **Correction mostly affects low L** (L < 200), negligible at higher L
- **Documentation**: Full details added to project-context.md

### 2026-01-07: Implemented Equation A3 corrections in notebook
- **Commit**: 8051734
- Added Section 9 to `notebooks/explore_dr6_cmbonly.ipynb`
- **Equation A3**: `C_L^{κκ,th} ≈ C_L^{κκ} + 2[d ln R/dC] ΔC · C_L^{κκ}|fid + [dN1/dC] ΔC`
- **Key implementation details**:
  - Normalization correction: `-2 * (dAL_dC @ ΔDℓ) / fAL` (sign verified against paper)
  - N1 CMB correction: `dN1_xx @ ΔDℓ` for TT, EE, TE
  - BB excluded (no measurement available)
- **Loaded matrices from** `/home/jiaqu/act_dr6_lenslike/data/v1.1/like_corrs/`:
  - `n0mv_fiducial_*.txt` — fiducial A_L normalization
  - `norm_correction_matrix_*.npy` — dA_L/dC (shape 4×4001×3001)
  - `N1der_*.txt` — N1 derivatives (shape 3000×3000)
- **Plots**: 4-panel figure with correction breakdown by spectrum
- **Numerical summary**: Corrections at L=100, 500, 1000
- **Verified**: Notebook implementation matches `get_corrected_clkk()` in dr6plus_lenslike.py
  - Same formula: `-2 * (dAL_dC @ cldiff) / fAL` for normalization
  - Same formula: `dN1_X @ cldiff` for N1 CMB correction
  - Only difference: notebook uses measured Planck+ACT CMB instead of theory CMB

### 2026-01-07: Implemented Equation 35 constant correction in notebook
- **Commit**: 4e855fb
- Added Section 10 to `notebooks/explore_dr6_cmbonly.ipynb`
- **Equation 35 constant**: `const = 2[d ln R/dC] ΔC + [dN₁/dCʲ] ΔCʲ - [dN₁/dC^{κκ}] C^{κκ,fid}`
- **Key differences from Eq. A3**:
  - Normalization: raw `norm_corr` (NOT multiplied by clkk_fid)
  - N1 κκ: new term `-dN1_kk @ clkk_fid`
- **Output files**:
  - `notebooks/eq35_constant_correction.txt`
  - `notebooks/eq35_constant_correction.npy`
  - `notebooks/eq35_constant_correction.png`

### 2026-01-07: Attempted Equation 34 modified covariance (BLOCKED)
- **Goal**: Compute CMB-marginalized lensing covariance: Σ̄ᵢⱼ = Σᵢⱼ + Mᵢˣ'ℓ cov_CMBˣℓ;Yℓ' Mⱼʸ'ℓ'
- **Commits**: d948837, a4b6dac (attempted fixes)
- **FUNDAMENTAL PROBLEM**: Dimension mismatch between matrices
  - **M matrices** (dAL_dC, dN1): Operate on **unbinned** CMB spectra at full ℓ resolution (ℓ=0-3000)
  - **Planck/ACT covariances**: **Binned** covariance (613 bins for Planck, 127 for ACT)
  - These are **incompatible** - cannot simply interpolate binned covariance to full ℓ
- **Why interpolation fails**: Binned covariance represents variance of bin-averaged Cℓ, not per-ℓ variance. The variance of ∑wᵢCℓᵢ ≠ interpolated version of var(bin-avg)
- **Status**: Section 11 added to notebook but produces unphysical results (% increase ~10^15)

### 2026-01-07: Equation 34 modified covariance (PARTIAL)
- **Commit**: d56813f
- **Reference**: ACT DR6 lensing paper states they use plik_lite bandpower covariance directly
- **Key assumption**: "CMB power spectra are represented only by modes that are smooth over Δℓ = 50"
- **Correct approach**: Transform M matrices to CMB bin space (NOT interpolate covariance)
- **Formula**:
  ```
  M̃ᵢᵇ = [Σ(ℓ∈b) Mᵢℓ] / [Σ(ℓ∈b) wℓ]
  ```
- **Results (Planck only)**:
  | Bin | L_center | % increase |
  |-----|----------|------------|
  | 0   | ~14      | 4.86%      |
  | 1   | ~30      | 0.02%      |
  | 2+  | >50      | <0.001%    |
- **ISSUE IDENTIFIED (2026-01-12)**: Only Planck covariance is used!
  - Planck covers low ℓ (TT: ℓ≤1000, TE/EE: ℓ≤600)
  - ACT covariance (`act_cov_full`, 127×127) is loaded but NOT used
  - This explains why corrections are small at high L — ACT uncertainties not propagated
- **Output files**:
  - `notebooks/eq34_modified_covariance.npy` (incomplete — Planck only)
  - `notebooks/eq34_modified_covariance.png`

### 2026-01-12: CRITICAL CORRECTION — Planck not needed
- **User clarification**: ACT lensing is reconstructed from ℓ_CMB ∈ [600, 3000] only
- **Norm correction matrices** also use `lmin=600, lmax=3000`
- **Implication**: We do NOT need Planck covariance at all — only ACT CMB covariance
- **Previous approach was wrong**: Eq. 35 in old notebook incorrectly used Planck low-ℓ data
- **New notebook**: `notebooks/new_covmat.ipynb` — user is consolidating correct approach here

### 2026-01-12: Reference files for corrections
- **Norm correction generation**: `/home/jiaqu/DW/280922_AL_finitediffs.ipynb`
  - Uses finite differences with `lmin=600, lmax=3000`
  - Generates `norm_correction_matrix_*.npy` files
- **N1 derivative generation**: `/home/jiaqu/DR6_lenslike/n1_correction.py`
  - Uses `lmin=600, lmax=3000` by default
  - Output: `N1der_{ESTIMATOR}_lmin600_lmax3000_full.txt`
  - Core functions in `/home/jiaqu/so-lenspipe/solenspipe/biastheory.py`
- **Planck variant** (NOT needed for ACT): `/home/jiaqu/DR6_lenslike/n1correction_planck.py` uses `lmin=100, lmax=2048`

### 2026-01-12: Eq. 34 Debugging — ROOT CAUSE IDENTIFIED

**Initial problem**: Diagonal increase was ~10^5% (way too large)

**Fix 1 applied**: Simple sum over bin support (not window-weighted) for M rebinning
- Per DR6 paper: "CMB power spectra are smooth over Δℓ = 50"
- Changed `doubly_bin_M` to sum M values over bin support without weighting

**Fix 2 applied**: Unit conversion M^Cℓ → M^Dℓ
- M matrices built in Cℓ units (from dAL_dC which expects Cℓ)
- ACT covariance is in Dℓ² units
- Conversion: M^Dℓ = M^Cℓ × 2π/(ℓ(ℓ+1))

**Result after fixes**: Diagonal increase now ~0.03% (way too small vs paper's ~6%)

**Diagnostic findings** (added print statements to trace magnitudes):
```
M_TT raw (Cℓ units): max ~10^-5
M_TT after L-binning: max ~10^-5
M_TT after Dℓ conversion: max ~10^-12 (divided by ell_factor ~10^5-10^6)
M_TT_doubly (after ℓ-binning): max ~10^-12 (NO accumulation!)
cov_TT_TT diagonal: 1-1000 µK⁴
cov_add diagonal: ~10^-22
lens_cov_orig diagonal: ~10^-17 to 10^-16
```

**ROOT CAUSE: Sign cancellation in M matrix**
- M values oscillate in sign within each ACT CMB bin
- When summing M over bin support (~100-150 ℓ per ACT bin), positive and negative values cancel
- Result: M_doubly ≈ M_unbinned (no accumulation from sum)

**Why paper gets 6% but we get 0.03%:**
- Paper uses **Planck plik_lite covariance** with narrow bins (Δℓ ≈ 10-50)
- We use **ACT CMB covariance** with wide bins (Δℓ ≈ 100-150)
- Wider bins → more ℓ values per bin → more cancellation → smaller effect

**Key insight from paper (DR6 lensing):**
> "We evaluate the CMB power correction using the plik_lite band powers, which are calculated from the full plik high-ℓ likelihood by marginalizing over the foreground model"

The paper explicitly uses Planck (not ACT) CMB covariance for the marginalization!

### 2026-01-12: Eq. 34 RESOLVED — Correct path to correction matrices

**Root cause**: Wrong path to likelihood correction matrices (was using old v1.1, needed v1.2)

**Fixed path**: `/home/jiaqu/spt_act_likelihood/act_dr6_spt_lenslike/data/v1.2/like_corrs/`

**Results now reasonable**: Diagonal increase ~1-2% at low L bins (see user's plot)

**Verification**: Compared double-binning vs interpolated approaches — both give similar results

### 2026-01-12: Production covariance file saved

**Output file**: `dr6plus_lenslike/data/v1.0/covmat_act_cmbmarg_analytic.txt`
- Full 18×18 covariance (no bin cuts applied — cuts done in likelihood)
- Uses double-binning approach with ACT CMB covariance
- Metadata saved to `covmat_act_cmbmarg_analytic_info.npy`

**Key implementation details**:
- M matrices binned along ℓ using ACT bandpower windows
- Full cross-spectrum CMB covariance included (TT-TT, TT-TE, TT-EE, TE-TE, TE-EE, EE-EE)
- Bin cuts (2:-6 for baseline) applied in `dr6plus_lenslike.py`, NOT in covariance file

### 2026-01-12: Added `analytic_marg` option to likelihood

**Commit**: 6291965

**Changes to `dr6plus_lenslike/dr6plus_lenslike.py`**:
1. Added `analytic_marg = False` class attribute (line 617)
2. Added validation: raises `ValueError` if `analytic_marg=True` but `lens_only=False`
3. Passed `analytic_marg` through to `load_data()`
4. Modified covariance loading (lines 447-450):
   - `analytic_marg=True` → loads `covmat_act_cmbmarg_analytic.txt`
   - `analytic_marg=False` → loads `covmat_act_cmbmarg.txt` (existing behavior)

**Status**: Covariance selection working.

### 2026-01-13: Added `get_lens_only_corrected_clkk` function

**Commit**: 3368dbc

**Changes to `dr6plus_lenslike/dr6plus_lenslike.py`**:
1. Added `get_lens_only_corrected_clkk(data_dict, clkk)` function (lines 146-158)
   - Formula: `nclkk = clkk + dN1_kk @ clkk + lens_only_const`
2. Data loading when `analytic_marg=True` (lines 470-478):
   - `dN1_kk` matrix from `like_corrs/N1der_KK_lmin600_lmax3000_full.txt`
   - `lens_only_const` from `lens_only_const.npy` (key: `eq35_const`)
   - Stored `d['analytic_marg'] = analytic_marg` (line 336)
3. Modified `generic_lnlike()` (lines 613-622):
   - Three-way logic: `likelihood_corrections` → `analytic_marg` → raw `cl_kk`

**Required data file**: `dr6plus_lenslike/data/v1.0/lens_only_const.npy` (copy from `notebooks/eq35_constant_correction.npy`)

### 2026-02-12: Verification of CMB covariance precision effect

**Goal**: Verify that the difference between analytic (blue, ~2.5% max error bar increase) and chain-based (orange, ~5.5% max) approaches is due to CMB precision.

**Hypothesis**: Using ACT DR4+Planck CMB covariance (instead of ACT DR6) in the analytic Eq. 34 formula should recover the larger orange error bars from the chain-based approach.

**Data sources**:
- ACT DR4: `/home/jiaqu/pyactlike/pyactlike/data/` (260 bins TT+TE+EE, covariance)
- Planck 2015: `/home/jiaqu/pyactlike/pyactlike/data/planck2015.dat` (8000 lines)
- Current analytic (blue): Uses ACT DR6 from `/home/jiaqu/DR6-ACT-lite/`
- Chain-based (orange): Sampled from Planck+DR4 MCMC chains

**Approach**: Create new notebook `notebooks/verify_dr4_planck_covmat.ipynb` that:
1. Loads ACT DR4 + Planck 2015 CMB bandpowers and covariance
2. Applies same Eq. 34 modified covariance formula as `new_covmat.ipynb`
3. Compares diagonal increase to both blue (DR6 analytic) and orange (chain-based) bars
4. If successful, should recover ~5.5% max increase (matching orange)

### 2026-02-12: ACT DR4+Planck CMB Covariance Verification — COMPLETED

**Task**: Verify if using ACT DR4 (instead of DR6) CMB covariance in the analytic Eq. 34 formula recovers the larger error bar increases (~5.5%) seen in chain-based approaches.

**Result**: ❌ **DR4 gives SMALLER increases than DR6** (opposite of expected)

**Key Findings**:
- DR6 analytic (reference): max 2.35% diagonal increase
- **DR4 analytic (this work): max 0.135% diagonal increase** (17× smaller!)
- Chain-based target: ~5.5% increase
- **Conclusion**: ACT DR4 does NOT explain the chain-based increases. It makes things worse.

**Work Completed**:
- Created `verify_dr4_covmat_v2.ipynb` following exact structure of `new_covmat.ipynb`
- Loaded ACT DR4 data (260 bins from 2 patches: deep + wide)
- Debugged and fixed THREE major unit conversion issues:
  1. **DR4 uses Cℓ (dimensionless), DR6 uses Dℓ (µK²)**
     - Conversion: Dℓ = Cℓ × ℓ(ℓ+1)/(2π)
     - For variance: var(Dℓ) = var(Cℓ) × [ℓ(ℓ+1)/(2π)]²
  2. **M matrices need clkk_fid normalization**
     - Correct formula: M = -2 × (dAL/dC) / fAL × clkk_fid
     - Missing clkk_fid (~10⁻⁷) caused 10⁷ factor error
  3. **M matrices expect Cℓ variance, not Dℓ variance**
     - Even after converting DR4 to Dℓ for plotting, must convert BACK to Cℓ variance for Eq. 34
     - This was the final 10¹¹ factor at ℓ~2000
- Generated theory vs data plots (look good)
- Applied Eq. 34 modified covariance with properly constructed M matrices

**Output Files**:
- `verify_dr4_v2_theory_vs_data.png` — Theory vs data comparison
- `verify_dr4_v2_combined_spectra.png` — Combined patches
- `verify_dr4_v2_comparison.png` — Shows DR4 < DR6 << chain target
- `verify_dr4_v2_results.npy` — Numerical data
- `verify_dr4_v2_RESULTS.md` — Full documentation

**Implications**:
The chain-based ~5.5% increase remains **unexplained**. Possible reasons:
1. Different CMB dataset used in chains (maybe Planck-only, not ACT)
2. Non-linear parameter space effects in MCMC sampling
3. Off-diagonal covariance terms not captured by diagonal-only Eq. 34
4. Different ℓ cuts or foreground treatment in chain analysis
5. Chain-based approach may include additional systematic uncertainties

**Next Steps**: Need to investigate:
- What CMB dataset the chains actually used
- Whether chains used different ℓ cuts or foreground models
- Whether off-diagonal covariance terms are significant
- Other effects beyond CMB covariance propagation

### 2026-02-12: EXHAUSTIVE INVESTIGATION — Why DR4 gives small covariance blow-up

**Previous conclusion (INCORRECT)**: DR4 windows have "extended support over 3000 ells" requiring different treatment.

**NEW INVESTIGATION (same session)**: User questioned whether DR4 bins are actually simple averages like DR6.

#### Key Finding: DR4 windows are NOT truly extended!

**Analysis of DR4 window function (`coadd_bpwf_15mJy_191127_lmin2.npz`)**:

```
Bin 0 weight distribution:
  ell < 200:    0.11% of total weight
  ell 200-500:  99.59% of total weight  ← CONCENTRATED HERE
  ell 500-1000: 0.30% of total weight
  ell > 1000:   0.00% of total weight

Weight at ell=0:   4.80e-07
Weight at ell=352: 2.41e-02 (peak)
Ratio: 2e-05 (ell=0 is negligible)
```

**The "3000 ell support" was misleading** — it's just tiny numerical leakage (values ~10⁻⁷), not meaningful weight.

**90% of weight concentrated in ~90-200 ells** (comparable to DR6's ~50 ells)

#### However: Windows are NOT simple averages

```
std/mean of weights within significant region = 1.6
(simple average would have std/mean ≈ 0)
```

The weights vary significantly and even have **negative values at edges** (not flat/uniform).

#### The 170× discrepancy between methods

For DR4 bin 0 with mock M ~ 1/ell²:
- Simple sum (significant ells): 8.45e-04
- M @ window.T (weighted):       4.98e-06
- **Ratio: 170×**

The `M @ window.T` method effectively computes `M_at_bin_center × sum(w)` ≈ single ell value.
The simple sum method computes total response over ~170 ells.

**Status**: Investigation interrupted. Need to verify if DR4 should use simple sum like DR6.

### 2026-02-15: DR4 Binning Methodology — ROOT CAUSE IDENTIFIED

**Investigation source**: Examined `pyactlike/like.py` and `notebooks/plot_spectra.ipynb`

#### DR4 Binning Operation (from like.py lines 195-213)

```python
# Step 1: Convert theory Dℓ → Cℓ (dimensionless)
cltt[1:tt_lmax] = dell_tt / l_list / (l_list + 1.0) * 2.0 * np.pi

# Step 2: Apply window via matrix multiplication
cth_tt = win_func_w[2*bmax : 3*bmax, 1:lmax_win] @ cltt[1:lmax_win]
# Output: binned Cℓ (dimensionless)
```

**Key insight**: Binning is `C_binned = Window @ C_ell` (matrix multiplication)

#### Window Properties Verified

| Property | DR4 Value | DR6 Value |
|----------|-----------|-----------|
| Normalization | sum(W) ≈ 1.00 | sum(W) = 1.00 |
| Effective bin width | ~50 ells (threshold >0.01) | 50 ells |
| Array shape | (n_bins, n_ell) | (n_ell, n_bins) |
| Significant ells | ~170-350 at >0.001 threshold | 50 exact |

**CRITICAL FINDING**: DR4 bins are ~50 ells wide (same as DR6!) when using appropriate threshold (>0.01 of max weight). The earlier "100-150 ells" estimate used threshold >0.001 which includes tails.

#### Physical Derivation: Simple Sum is Correct

For Eq. 34 modified covariance:
1. Binning relation: `Ĉ_b = Σ_ℓ W_{b,ℓ} × C_ℓ`
2. For smooth CMB (C_ℓ = constant in bin): `Ĉ_b = c × Σ W = c × 1.0 = c`
3. When `Ĉ_b` fluctuates by `δĈ_b`, each `C_ℓ` fluctuates by `δĈ_b` (smoothness)
4. Lensing response: `δĈ_κκ = Σ_ℓ M_ℓ × δC_ℓ = δĈ_b × Σ_ℓ M_ℓ`

**Therefore: M_binned = Σ_ℓ M_ℓ (SIMPLE SUM over bin support)**

The `M @ W.T` method computes weighted average ≈ M at single representative ell, NOT total response.

#### Bug in verify_dr4_covmat_v2.ipynb

Current code (WRONG):
```python
M_doubly_dr4[key] = M_trunc @ window_trunc.T  # Weighted average
```

Correct code:
```python
for b in range(n_cmb_bins):
    sig_ells = np.where(window_trunc[b, :] > 0.01)[0]  # bin support
    M_doubly_dr4[key][:, b] = np.sum(M_trunc[:, sig_ells], axis=1)
```

**Expected result after fix**: DR4 should give ~2% diagonal increase (similar to DR6), not 0.04%.

### 2026-02-15: Simple Sum Test Result — STILL TOO LARGE

**Result**: DR4 with simple sum gives **394% max diagonal increase** (expected ~2-5%)

**Analysis**: The simple sum approach overcounts because:
1. DR4 has **two patches** (deep + wide), contributing 36 covariance blocks vs DR6's 9
2. The equation sums over all (X,Y) spectrum pairs AND all (β,β') bin pairs
3. Cross-patch terms may be non-negligible

**Eq. 34 structure**:
$$\bar{\Sigma}_{ij} = \Sigma_{ij} + \sum_{X,Y} \sum_{\beta,\beta'} \tilde{M}_i^{X,\beta} \text{Cov}(C_\beta^X, C_{\beta'}^Y) \tilde{M}_j^{Y,\beta'}$$

| Dataset | Spectra | Covariance blocks | Total CMB bins |
|---------|---------|-------------------|----------------|
| DR6 | 3 (TT, TE, EE) | 9 | ~127 |
| DR4 | 6 (TT/TE/EE × 2 patches) | 36 | 260 |

**Open questions**:
1. Are cross-patch covariances (TT_p1-TT_p2, etc.) significant or near-zero?
2. Should patches be combined before applying Eq. 34, not summed separately?
3. Is the M matrix the same for both patches, or should it differ?

**Status**: Paused for investigation. Need to check cross-patch covariance magnitudes.

### 2026-02-15: Planck plik_lite Eq. 34 — RESOLVED

**Commits**: af777b6, 1aa5ad9

**Root cause of 305-3043% bug**: TWO incorrect approaches combined:
1. **Unit conversion bug**: Converting Planck cov from Cℓ² to Dℓ² multiplies by (ℓ(ℓ+1)/(2π))² ≈ 10^10
2. **M summation bug**: Summing M over ~50 ells multiplies result by N²

**Correct approach**: M at bin center, use Cℓ² covariance directly (NO unit conversion)

| Method | Result | Status |
|--------|--------|--------|
| M at bin center, Cℓ² cov | **1.2%** | ✓ CORRECT |
| M summed, Cℓ² cov | 305% | ✗ Too high |
| M summed, Dℓ² cov | 10^14% | ✗ WAY too high |

**Gap to paper's 6%**: Explained by ACT M matrices using lmin=600 vs Planck's lmin=100

**Files created**:
- `notebooks/debug_planck_eq34.py` — Systematic hypothesis testing (10 approaches)
- Updated `notebooks/verify_planck_covmat.ipynb` with corrected implementation

### 2026-02-15: Planck PR4 Lensing CMB Marginalization Investigation

**Source**: `/home/jiaqu/planck_PR4_lensing/` (Planck PR4 NPIPE lensing likelihood by J. Carron et al.)

#### Two Likelihood Versions
1. `PlanckPR4Lensing` — Regular (requires TT, EE, TE, PP from theory)
2. `PlanckPR4LensingMarged` — CMB-marginalized (only requires PP)

#### Covariance Increase from CMB Marginalization (Eq. 34)

| Bin | L_center | % increase |
|-----|----------|------------|
| 1   | 28       | 4.4%       |
| 2   | 64       | 13.0%      |
| 3   | 106      | 13.6%      |
| 4   | 150      | 7.9%       |
| 5   | 195      | 3.5%       |
| 6   | 240      | 1.3%       |
| 7   | 286      | 0.6%       |
| 8   | 331      | 0.5%       |
| 9   | 377      | 0.3%       |

**Peak increase**: ~13% at L~100 (vs our ~2% with ACT M matrices)

#### Linear Correction Structure
- **Regular**: M matrices for TT, EE, TE, PP (~148KB per bin window file)
- **CMBmarged**: Only PP correction (~26KB per bin window file)
- CMB contribution (TT+EE+TE) is **92-99%** of total fiducial correction

#### M Matrix Properties
- Non-zero at ℓ ∈ [100, 2048] (Planck CMB reconstruction range)
- Values ~10⁻¹² (derivative dC_L^φφ/dD_ℓ^XX)
- Oscillatory (positive and negative values)

#### Key Files
- `data_pr4/pp_*_cov.dat` — Regular covariance (9×9)
- `data_pr4/pp_*_CMBmarged_cov.dat` — CMB-marginalized covariance
- `data_pr4/pp_*_lens_delta_window/` — M matrices for linear corrections
- `data_pr4/pp_*_lensing_fiducial_correction.dat` — Precomputed fiducial term

#### Comparison to Our Approach

| Aspect | Planck PR4 | Our ACT approach |
|--------|------------|------------------|
| M matrix lmin | 100 | 600 |
| Max cov increase | ~13% | ~2% |
| CMB covariance | Internal (unknown) | plik_lite / ACT DR6 |

#### Gap Explanation
The 13% vs 2% gap is likely due to:
1. **ℓ range**: Planck lmin=100 captures more CMB modes than ACT lmin=600
2. **CMB covariance**: Planck may use internal covariance (not plik_lite)

#### Missing Information
- Code that **generates** the CMBmarged files is not in the repository
- README states "band-powers likelihood code has been left out"
- Unclear what exact CMB covariance was used

## Outstanding Tasks
1. ~~Explore `/home/jiaqu/DR6-ACT-lite/`~~ — **DONE**
2. ~~Debug Eq. 34 modified covariance~~ — **DONE** (fixed path issue)
3. ~~Save production covariance file~~ — **DONE** (`covmat_act_cmbmarg_analytic.txt`)
4. ~~Add `analytic_marg` option for covariance selection~~ — **DONE** (commit 6291965)
5. ~~Add `get_lens_only_corrected_clkk` function~~ — **DONE** (commit 3368dbc)
6. ~~Identify correct windowing method for Eq. 34~~ — **DONE** (simple sum for DR6)
7. **[ACTIVE]** Match chain-based approach analytically — see Worker Prompt below
8. ~~Planck plik_lite gives 305-3043%~~ — **RESOLVED** (M at bin center gives 1.2%)
9. ~~Investigate what CMB dataset was used in chain-based approach~~ — ACT DR4+Planck chains
10. Copy `lens_only_const.npy` data file and test the full pipeline
11. Compare likelihood results between `analytic_marg=True` vs `False`
12. Add/improve tests for `lens_only=True` pathway

## 2026-02-15 Investigation Summary

**Problem:** Analytic Eq. 34 gives wrong results compared to chain-based (5.2% at L≈123)
- M at bin center: 1.2-1.8% (wrong magnitude)
- M summed over bins: 69-3043% (way too large)

**Key findings:**
- plik_lite/DR4 CMB covariance is essentially diagonal (correlations ~0.04 adjacent, ~0 otherwise)
- M matrix (dAL/dC) dominated by high ℓ (1800-3000), where dAL/dC is 17× larger than at ℓ<1000
- The 5%-95% cumulative weight of M falls in ℓ ∈ [1800, 2990]

**Reference chain notebook:** `/home/jiaqu/DW/100223_compute_cov_norm.ipynb`
- Uses ACT M matrices with lmin=600
- Samples from ACT DR4+Planck chains → CAMB Cℓ → M @ ΔCℓ → bin to lensing L
- Result: ~5.2% diagonal increase at L≈123

**Next step:** Find correct binning/extrapolation to match chain result using measured CMB covariances (DR4 + Planck 2015 bandpowers and covariance matrices)

---

## Worker Prompt: Match Chain-Based Eq. 34 Using Measured CMB Covariance

Read `.agent/worker-instructions.md` and `.agent/project-context.md`.

### Task
Find the correct binning/extrapolation method to make analytic Eq. 34 match the chain-based result using **measured** DR4 + Planck 2015 CMB bandpowers and covariances.

### Reference
Chain-based approach: `/home/jiaqu/DW/100223_compute_cov_norm.ipynb`
- Uses ACT M matrices with lmin=600 (NOT lmin=100)
- Samples from ACT DR4 + Planck chains → CAMB Cℓ → M @ ΔCℓ → bin to lensing L
- Result: ~5.2% diagonal increase at L≈123

### Data to Load

1. **DR4 CMB bandpowers and covariance:**
   - `/home/jiaqu/pyactlike/pyactlike/data/cl_cmb_ap.dat` (260 bins: 2 patches × (40 TT + 45 TE + 45 EE))
   - `/home/jiaqu/pyactlike/pyactlike/data/c_matrix_ap.dat` (260×260 Fortran binary)
   - Window functions: `coadd_bpwf_15mJy_191127_lmin2.npz`, `coadd_bpwf_100mJy_191127_lmin2.npz`

2. **Planck 2015 CMB bandpowers and covariance:**
   - `/home/jiaqu/pyactlike/pyactlike/data/planck2015.dat`
   - Note: This is different from plik_lite - check the chain notebook for exact files used

3. **ACT M matrices (lmin=600):**
   - `/home/jiaqu/spt_act_likelihood/act_dr6_spt_lenslike/data/v1.2/like_corrs/`
   - `norm_correction_matrix_Lmin0_Lmax4000.npy` — dAL/dC shape (4, 4001, 3001)
   - `n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax4000.txt` — fAL normalization
   - `N1der_*.txt` — N1 derivatives
   - `cosmo2017_10K_acc3_lensedCls.dat` — fiducial CMB Dℓ
   - `cosmo2017_10K_acc3_lenspotentialCls.dat` — fiducial φφ spectrum

### Key Question
The M matrices operate on unbinned Cℓ (ℓ=600-3000), but we only have binned CMB covariance. What is the correct transformation?

Options to test:
1. Interpolate binned covariance to unbinned ℓ space
2. Use window functions to properly propagate bin → ℓ covariance
3. Derive unbinned Cℓ covariance from chains (compute CAMB Cℓ for each sample, take covariance)

### Validation
The analytic result MUST match chain-based:
- Peak at L≈123 (not L≈1400)
- ~5% diagonal increase (not 1% or 3000%)

### Output
- Script: `notebooks/match_chain_measured_cov.py`
- Diagnostic plots comparing approaches
- Document which binning method works

Add any new source files and commit changes to all tracked files.
Report what was changed and where.

## ~~CRITICAL DEBUG: Planck plik_lite Eq. 34 Investigation~~ — RESOLVED 2026-02-15

**Solution**: Use M at bin center (not summed), Cℓ² covariance directly (no unit conversion)
**Result**: 1.2% diagonal increase (gap to paper's 6% explained by lmin=600 vs lmin=100)
**Files**: `notebooks/debug_planck_eq34.py`, `notebooks/verify_planck_covmat.ipynb`

See "2026-02-15: Planck plik_lite Eq. 34 — RESOLVED" in Session History for details.

### 2026-02-15: ROOT CAUSE IDENTIFIED — DR4 Analytic Failure

**Investigation files created:**
- `notebooks/test_chain_vs_analytic.py` — Systematic comparison of M binning strategies
- `notebooks/investigate_M_ell_range.py` — Analysis of dAL/dC ℓ-dependence
- `notebooks/test_chain_vs_analytic.png` — Strategy comparison plot
- `notebooks/M_matrix_structure.png` — M matrix structure visualization
- `notebooks/mock_parameter_cov.png` — Mock parameter covariance test

**Key finding: M matrix dominated by high ℓ (1800-3000)**

| ℓ range | dAL/dC max value |
|---------|------------------|
| 0-600 | 0 (lmin=600) |
| 600-1000 | 1.92e-07 |
| 1000-2000 | 1.53e-06 |
| 2000-3000 | **3.21e-05** (17× larger) |

**Why chain-based works:**
1. Samples cosmological parameters (As, ns, etc.)
2. Parameter variations create **correlated** ΔCℓ across ALL ℓ
3. When As changes by 1%, Cℓ changes by ~1% at EVERY ℓ
4. The high-ℓ modes (ℓ~2900 where dAL/dC peaks) are fully captured

**Why analytic with binned covariance fails:**
1. Binned CMB covariance captures bin-to-bin measurement correlations
2. Does NOT capture the cross-ℓ correlations from parameter variations
3. M binning (center/sum/weighted) all fail because they don't account for parameter-space structure

**Mock test with parameter-derived covariance:**
Using Cov(Cℓ,Cℓ') = Cℓ×Cℓ'×(σAs/As)² (1% As uncertainty):
- Mock analytic: **4.10% at L≈232** (close to reference!)
- Chain-based: **5.20% at L≈123**

**Conclusion:** The analytic Eq. 34 approach requires **parameter-derived CMB covariance**, not measurement covariance. The binned CMB covariance (DR4/DR6/Planck) represents measurement uncertainty, not the parameter-induced correlations that dominate the lensing normalization response.

**Solutions:**
1. Continue using chain-based approach (already works)
2. Construct parameter-derived CMB covariance from Fisher matrix
3. Use Planck parameter covariance → CMB covariance → Eq. 34

### 2026-02-16: Chain-based Covariance Comparison — RESOLVED

**Investigation**: Comparing chain-based covariance using `p-actlite_lcdm_camb` (DR6+Planck) vs reference (DR4+Planck)

**Confirmed**: `p-actlite_lcdm_camb` DOES include Planck data:
- `planck_2018_lowl.TT` + `planck_2018_lowl.EE_sroll2` (low-ℓ)
- `act_dr6_cmbonly.PlanckActCut` (Planck high-ℓ cut at TT<1000, TE/EE<600)
- `act_dr6_cmbonly.ACTDR6CMBonly` (ACT DR6, ℓ=600-8500)

**DR4+Planck chains downloaded**:
- Path: `/project/rrg-rbond-ac/jiaqu/chains/dr4_planck_lcdm/`
- File: `CLASS2p8_ACTPol_lite_DR4_leakfix_yp2_baseLCDM_taup_planck2018_lowTT_plikHM_TT_lmax650_hip`
- Likelihoods: ACT DR4 + Planck 2018 low-ℓ TT + Planck plik TT (ℓmax=650)

**Parameter constraint comparison (DR6 vs DR4)**:

| Parameter | DR6+Planck σ | DR4+Planck σ | DR6/DR4 |
|-----------|-------------|-------------|---------|
| logA      | 0.0185      | 0.0257      | **0.72** |
| ns        | 0.0062      | 0.0064      | 0.97    |
| ombh2     | 0.00021     | 0.00019     | 1.06    |
| omch2     | 0.0015      | 0.0027      | **0.55** |
| H0        | 0.68        | 1.12        | **0.61** |
| sigma8    | 0.0089      | 0.0132      | **0.67** |
| tau       | 0.0064      | 0.0134      | **0.48** |

**Key finding**: DR6+Planck has **tighter** constraints than DR4+Planck (most ratios < 1)!

**Paradox**: If DR6 has tighter constraints, why does it give LARGER covariance increase (7.62% vs 5.20%)?

**Possible explanations**:
1. The reference 5.20% may be from a different chain (not this DR4 chain)
2. Parameter degeneracies matter — specific combinations affect lensing differently
3. The sqrt fix changed the numbers — need to re-run with corrected formula

**Next step**: Re-run chain-based covariance with corrected sqrt formula for both DR6 and DR4 using batch job

## 2026-03-14: Self-Calibration End-to-End Test Configs — COMPLETED

### Run configs created and committed (branch dev)
- `runs/self_cal_test/evaluate_selfcal_false.yaml`
- `runs/self_cal_test/evaluate_selfcal_true.yaml`
- `runs/self_cal_test/mcmc_selfcal_false.yaml`
- `runs/self_cal_test/mcmc_selfcal_true.yaml`
- `runs/self_cal_test/cobaya.sh`
- `runs/self_cal_test/submit_selfcal_false.sh`
- `runs/self_cal_test/submit_selfcal_true.sh`

### Source fixes committed (commits 5d5ea00, 84c0f1a)
1. **`dr6plus_lenslike.py`**: `do_norm_corr=` → `fid_norm=` (latent bug fix in `get_corrected_clkk` call)
2. **`dr6plus_lenslike.py`**: added `Clens = params_values.get('Clens', 1.0); cl_kk = cl_kk * Clens`
3. **`dr6plus_lenslike.py`**: added `get_allow_agnostic()` returning `True`

### Evaluate test results (2026-03-15) — PASSED
- Output: `/scratch/jiaqu/chains/selfcal_test/`
- `chi2__dr6plus_lenslike.ACTDR6LensLike` = **19.729483** in BOTH selfcal_false and selfcal_true
- mflike chi2 differs (expected: evaluate uses random ref draws for mflike nuisances)
- **Conclusion**: selfcal=True with all params=0 gives identical lensing loglike to selfcal=False ✓

### MCMC run status (2026-03-15)
- `mcmc_selfcal_true`: **completed**, Rminus1=0.0097 (~979k steps) — converged
- `mcmc_selfcal_false`: **still running** (~861k steps, Rminus1=NaN), lock at 21:44
- `mcmc_selfcal_false.yaml` has uncommitted change: `force: false` → `resume: True`

### Uncommitted change
- `runs/self_cal_test/mcmc_selfcal_false.yaml`: `resume: True` (user set to resume an interrupted run)

## 2026-03-15: selfcal_linked run created

**Commit**: worker output (see git log)

- `runs/self_cal_test/mcmc_selfcal_linked.yaml` — selfcal=True but c_pa5a/b/c_pa6a/b derived from `cal_dr6_pa5_f090` etc., p_pa5a/b/p_pa6a/b derived from `calE_dr6_pa5_f090` etc.
- `runs/self_cal_test/submit_selfcal_linked.sh` — SLURM submission for linked run
- Output: `/scratch/jiaqu/chains/selfcal_test/mcmc_selfcal_linked`

**Design rule (permanent)**: In all future run configs with `selfcal: true`, the 8 selfcal params MUST be derived from the mflike cal params — never independently sampled. Independently sampling them double-counts calibration uncertainty and widens A_lens instead of tightening it.

## Notes
- Existing tests in `tests/test_dr6plus_lenslike.py` use mocked data
- `lens_only=True` currently loads CMB-marginalized covariance and skips likelihood corrections
- M matrices operate on Cℓ (not Dℓ)

## 2026-06-01: Lensing + DESI DR2 BAO fg-marg runs (manager-direct)

Extend the widened-prior A_fg test by adding fiducial DESI DR2 BAO, to measure the **sigma8** degradation (BAO breaks the sigma8-Omega_m degeneracy that lensing-only leaves open, so sigma8 — not just S8_lens — becomes meaningful). Done directly by manager at user request.

- BAO likelihood: `bao.desi_dr2.desi_bao_all_fd` (fiducial, full tracer set), added to the existing extended class_sz+lensing configs. class_sz already computes background distances (skip_chi/hubble=0); `cobaya-run --test` confirmed the combined model initializes and A_fg routes correctly.
- Reference for BAO setup: `/home/jiaqu/mnu_de/.agent/project-context.md` (sec "BAO: DESI DR2") and `mnu_de/base/bao_all_fd.yaml`.

| YAML                              | A_fg prior  | Job ID  | output (/scratch/jiaqu/chains/fg_marg_test/) |
|-----------------------------------|-------------|---------|----------------------------------------------|
| mcmc_control_extended_bao.yaml    | (no marg)   | 1692859 | control_extended_bao   |
| mcmc_afg_extended_bao.yaml        | U[-2, 2]    | 1692860 | afg_extended_bao       |
| mcmc_afg_extended_u5_bao.yaml     | U[-5, 5]    | 1692862 | afg_extended_u5_bao    |
| mcmc_afg_extended_u10_bao.yaml    | U[-10, 10]  | 1692863 | afg_extended_u10_bao   |

Submitted 2026-06-01 ~14:55. Follow-up: once converged, make the triangle plot (sigma8 vs A_fg) with degradation % in legend, analogous to `notebooks/fg_marg_S825_Afg_triangle.py`.

## 2026-05-21: Foreground bias marginalization (A_fg) — dispatched

### Task
Add a single nuisance parameter `A_fg` to the lensing likelihood that scales a fiducial foreground bias template added to the theory `clkk` bandpowers. Defensive marginalization following MacCrann et al. 2023 (arXiv:2304.05196).

### Design decisions
- **Model**: `C_L^model = C_L^κκ(θ) + A_fg · T(L)`, template on theory side.
- **Default**: `A_fg = 0.0` (opt-in via YAML; existing analyses unchanged).
- **Variant scope**: `act_baseline` / `act_extended` only; `ValueError` otherwise.
- **Template choice**: `total_mv_prh` (index 9, bias-hardened MV total) — matches DR6 baseline profile hardening.
- **Binning**: fine-L template (length 4501) truncated to L=0..lmax-1 and binned via `binning_matrix_act` at load time → cached bandpower-space vector of length `nbins_act`.
- **Compatible with** `lens_only=True/False` (additive on theory side) and `selfcal=True` (applied after T² rescaling — A_fg is a residual contamination, not a calibration effect).
- **Class attributes added**: `fg_marg: bool=False`, `fg_template_file: str=None`, `fg_template_index: int=9`.

### Template source
`/scratch/kaper/fg_dr6plus/agora_coadded_maps/coadded_fgsonly_nonoise_all_fgs_len_fluxcutTrue_50.0mJy_clsubTrue_clmodelcandidatesModelMap_minsnr5.0_psrcnemosubFalse_psrcdorysubTrue_psrcmodelpsrcmodel_minsnr4.0_3000.npy`
- Shape: `(10, 4501)`, dtype float64
- Row order: `total_qe, primary_qe, secondary_qe, trispectrum_qe, total_prh, primary_prh, secondary_prh, trispectrum_prh, total_mv_qe, total_mv_prh`
- Worker copies to `dr6plus_lenslike/data/v1.0/fg_template_act_baseline.npy`

### Dispatched commit scope
Prompt instructs worker to commit ONLY:
- `dr6plus_lenslike/dr6plus_lenslike.py` (A_fg additions only — the existing pending diff for forecast infra `n_drop_high`/`cov_file` is *not* part of this task)
- `dr6plus_lenslike/data/v1.0/fg_template_act_baseline.npy` (copied template)
- `tests/test_dr6plus_lenslike.py` (5 new tests F–J)

Other pending modified files left untouched by worker:
- `dr6plus_lenslike/dr6plus_lenslike.py` already has uncommitted forecast infra changes — worker layers A_fg on top and commits both blocks together (acceptable bundling within a single file).
- `runs/self_cal_test/mcmc_selfcal_{false,linked}.yaml` — explicitly excluded from A_fg commit.

### Open follow-ups (after A_fg lands)
1. Run an evaluate sanity check at `A_fg = 0` vs `fg_marg=False` to confirm identical loglike.
2. Run evaluate at `A_fg = 1` to measure the bandpower shift.
3. Draft an MCMC config with `A_fg ~ U[0, 2]` (separate worker prompt).

## 2026-05-21: A_fg pp→kk unit fix (commit adae9d6)

**Bug**: Agora .npy stores ΔC_L^{φφ}; my first commit added it directly to `bclkk` (C_L^{κκ} units) → bandpower template was ~10¹¹× too small (effectively a no-op). At `A_fg=1` lnlike would have been indistinguishable from `A_fg=0`.

**Fix**: `load_data` applies `pp_to_kk(template, L) = template * (L(L+1))² / 4` on the fine-L template *before* the binning matrix.

**How the bug was caught**: Visual comparison of `template/C_L^{φφ}_fid` vs MacCrann Fig.1: ratio of `-2.75%` at L=1000 matched the Agora-BH curve in φφ units, while ratio against `C_L^{κκ}_fid` was ~10⁻⁷ — clear unit mismatch.

**How tests still passed despite the bug**: The closed-form check in `test_fg_J` (`lnlike(A) - lnlike(0) = r·Σ⁻¹·t - ½ t·Σ⁻¹·t`) is algebraically valid for any cached `t`, so it doesn't catch scale errors. Useful follow-up test: assert `|template_bandpower / clkk_data|` is in `[1e-4, 1e-1]` (right physical scale).

## 2026-05-21: A_fg blind-mode sanity test plan

### Design
Four chains on fiducial data (no foreground bias injected), all `lens_only=True` for speed:
| Run | variant | fg_marg | Notes |
|---|---|---|---|
| 1 | `dr6plus_fiducial_baseline` | `false` | control |
| 2 | `dr6plus_fiducial_extended` | `false` | control |
| 3 | `dr6plus_fiducial_baseline` | `true`  | A_fg ~ U[-2, 2] |
| 4 | `dr6plus_fiducial_extended` | `true`  | A_fg ~ U[-2, 2] |

### Prior choice
**Sanity-test prior is U[-2, 2]** (not the spec's production U[0, 2]). Reasons:
1. Mode is interior, not on a boundary — visually unambiguous.
2. Posterior mean recovers 0 on fiducial-as-data → no "data weak, prior dragging" caveat to explain.
3. If a future injection test recovers the wrong sign, U[-2, 2] surfaces it immediately; U[0, 2] hides it at the boundary.

### Variant-guard change
`fg_marg=True` currently requires `variant in {'act_baseline', 'act_extended'}`. The blind test needs `dr6plus_fiducial_baseline`/`_extended`. The Agora template + M matrices apply equally to fiducial-data variants (same binning matrix and lensing reconstruction config), so widening is safe.

**New allowed set**: `{'act_baseline', 'act_extended', 'dr6plus_fiducial_baseline', 'dr6plus_fiducial_extended'}`.

### Expected from data (lens_only=True, baseline)
- `σ(A_fg) ≈ 7` (likelihood-only). Data weakly constraining.
- With `U[-2, 2]`: posterior mode at 0, mean ≈ 0, width ≈ prior width (~1.15).
- Cosmology shifts: `≲ 0.1σ` between control and marg runs.

### Reference paths
- Fiducial bandpowers: `dr6plus_lenslike/data/v1.0/clkk_act_fiducial.txt` (loaded automatically by `dr6plus_fiducial_*` variant)
- Base templates: `base_config/lensing_fiducial.yaml`, `base_config/lcdm_param.yaml`, `base_config/mcmc_sampler.yaml`, `base_config/class_sz.yaml`
- Existing analogue: `runs/self_cal_test/` (selfcal MCMC scaffolding)
- Output dir: `/scratch/jiaqu/chains/fg_marg_test/`

## 2026-05-21: A_fg blind sanity scaffolding committed and dispatched

### Commits
- `7e1c612` Add A_fg foreground bias marginalization (+ bundled forecast infra n_drop_high/cov_file)
- `adae9d6` Fix A_fg template units: apply pp_to_kk before binning
- `ab7555d` Add A_fg blind sanity test scaffolding (widened variant guard, runs/fg_marg_test/)

### SLURM jobs submitted (1h walltime each)
| Job ID  | Script                          | YAML                            |
|---------|---------------------------------|---------------------------------|
| 1619613 | submit_control_baseline.sh      | mcmc_control_baseline.yaml      |
| 1619614 | submit_control_extended.sh      | mcmc_control_extended.yaml      |
| 1619615 | submit_afg_baseline.sh          | mcmc_afg_baseline.yaml          |
| 1619616 | submit_afg_extended.sh          | mcmc_afg_extended.yaml          |

Chains write to `/scratch/jiaqu/chains/fg_marg_test/<run_name>`. Submission queued at 2026-05-21 12:40-ish.

### Post-completion checklist
1. **A_fg posterior**: read out 1D marginal from `afg_baseline` and `afg_extended`. Verify mode ≈ 0 (interior of U[-2,2]) and mean ≈ 0. If mode pegs at ±2, the chain is broken (or A_fg has a sign bug we missed).
2. **Cosmology comparison**: compare 1D marginals `control_baseline` vs `afg_baseline` (same for extended). Expect:
   - Means shift ≲ 0.1σ.
   - σ(σ8), σ(Ωm h²), σ(H0) slightly inflated under marg.
3. **σ(A_fg) from chain**: should be close to the analytic prediction (~0.58 with U[-2,2] dominating ~7 likelihood width → posterior ≈ prior in width).
4. **Cosmology error inflation table**: write a small comparison plot to `notebooks/fg_marg_comparison.png` (separate worker prompt after chains finish).

## 2026-05-28: Widened-prior extended A_fg runs (manager-direct, user override)

Collaborators were surprised S8_lens degrades so little under `A_fg ~ U[-2,2]`. User asked to test wider priors on the extended cosmology. Done directly by manager at user's explicit instruction ("DISPATCH AND DO THEM NOW").

New files (copies of `mcmc_afg_extended.yaml`, identical except prior/proposal/output):
| YAML                          | A_fg prior  | proposal | output dir                                  | Job ID  |
|-------------------------------|-------------|----------|---------------------------------------------|---------|
| mcmc_afg_extended_u5.yaml     | U[-5, 5]    | 1.5      | /scratch/jiaqu/chains/fg_marg_test/afg_extended_u5  | 1661090 |
| mcmc_afg_extended_u10.yaml    | U[-10, 10]  | 3.0      | /scratch/jiaqu/chains/fg_marg_test/afg_extended_u10 | 1661091 |

Plus `submit_afg_extended_u5.sh` / `submit_afg_extended_u10.sh`. Variant `dr6plus_fiducial_extended`, `lens_only`, `mock_covmat.npy` — unchanged so they compare apples-to-apples against the U[-2,2] extended run. Proposal widened from 0.5 since both priors are prior-dominated (likelihood-only σ(A_fg)≈7). Submitted 2026-05-28.

### Open follow-ups
- Compare S8_lens 1D marginal across U[-2,2] / U[-5,5] / U[-10,10] extended runs once chains finish (this answers the collaborators' question).
- Selfcal mock test fix (pre-existing test_selfcal_C/D/E failures from `(18,9)` vs `(10,9)` `response_cal_matrix` shape mismatch — unrelated to A_fg).
- Commit `.agent/manager-status.md` history updates from this session.
- The runs/self_cal_test/mcmc_selfcal_{false,linked}.yaml WIP edits are still uncommitted (intentionally — out of scope for this task).
