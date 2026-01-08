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

## Next Steps
1. Implement equation (35) constant correction in `lens_only=True` mode
2. Load pre-computed correction from notebook output files
3. Apply correction to lensing bandpowers in likelihood
