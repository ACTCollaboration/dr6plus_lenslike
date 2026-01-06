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
- DR6 CMB-only SACC file: `/home/jiaqu/DR6-ACT-lite/act_dr6_cmbonly/data/act_dr6_cmb_sacc.fits`
  - Contains foreground-marginalized TT, TE, EE bandpowers
  - Used for lensing likelihood corrections (normalization, N1 bias)

## Notebooks
- `notebooks/explore_dr6_cmbonly.ipynb` — Explores DR6 CMB SACC data products (tracers, bandpowers, covariance, windows)
