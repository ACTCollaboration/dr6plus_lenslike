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
- **Possible solutions** (require manager decision):
  1. **Get unbinned CMB covariance**: Would need analytic model (Knox formula) or separate calculation
  2. **Modify M matrices**: Transform to work with binned CMB (requires bandpower window convolution)
  3. **Skip this correction**: Use existing `covmat_*_cmbmarg.txt` files which may already include this effect
  4. **Different formulation**: Work entirely in binned space from the start
- **Status**: Section 11 added to notebook but produces unphysical results (% increase ~10^15)

## Session Handoff Notes (2026-01-07)
**Where we left off:**
- Completed Eq. 35 constant correction in notebook (commit 4e855fb) ✓
- **BLOCKED** on Eq. 34 modified covariance (binned vs unbinned mismatch)

**Decision needed from manager:**
- How to handle the CMB covariance for Eq. 34? (see options above)
- May need to consult with domain expert about correct formulation

**Key files:**
- Notebook: `notebooks/explore_dr6_cmbonly.ipynb` (Section 9 + 10 working, Section 11 broken)
- Correction output: `notebooks/eq35_constant_correction.npy`
- Likelihood code: `dr6plus_lenslike/dr6plus_lenslike.py`

## Outstanding Tasks
1. ~~Explore `/home/jiaqu/DR6-ACT-lite/`~~ — **DONE** (commit ebb4ecb)
2. ~~Update notebook to Planck+ACT cut version~~ — **DONE**
3. ~~Implement Equation A3 corrections in notebook~~ — **DONE** (commit 8051734)
4. ~~Implement Equation 35 constant correction in notebook~~ — **DONE** (commit 4e855fb)
5. **[NEXT]** Integrate precomputed correction into dr6plus_lenslike likelihood code
6. Add/improve tests for `lens_only=True` pathway
7. Ensure backward compatibility

## Notes
- Existing tests in `tests/test_dr6plus_lenslike.py` use mocked data
- `lens_only=True` currently loads CMB-marginalized covariance and skips likelihood corrections
- Foreground marginalization is the key focus for the CMB 2pt marginalization approach
- Planck+ACT notebook now serves as reference for data extraction and ell cuts
