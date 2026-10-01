# DR6+ Lensing Likelihood

Likelihood for the Atacama Cosmology Telescope DR6+ CMB lensing data.

## Installation

### From source (development)

```bash
git clone <repo-url>
cd dr6plus_lenslike
pip install -e .
```

### Dependencies

Requires [Cobaya](https://cobaya.readthedocs.io/) and at least one Boltzmann solver (CAMB or CLASS).

## Theory Codes

Three CAMB variants are supported. Each has a corresponding config in `base_config/`:

| Config file | Theory code | Description |
|-------------|-------------|-------------|
| `camb.yaml` | Standard CAMB (pip) | Default Boltzmann solver |
| `camb_alens.yaml` | CAMB_alens (`/home/jiaqu/CAMB_alens`) | Fork with A_lens, B_lens parameters and CosmoRec. Requires `module load gsl` |
| `camb_disputauble.yaml` | [disputauble](https://github.com/NoahSailer/disputauble) | CAMB wrapper allowing negative effective neutrino mass extrapolation ([arXiv:2504.16932](https://arxiv.org/abs/2504.16932)) |
| `class_sz.yaml` | CLASS_SZ | CLASS-based solver with SZ extensions |

### Installing theory codes

**Standard CAMB** (already available):
```bash
pip install camb
```

**CAMB_alens** (custom fork, already at `/home/jiaqu/CAMB_alens`):
```bash
# Requires module load gsl for the CosmoRec/GSL dependency (libgsl.so.25)
cd /home/jiaqu/CAMB_alens
pip install -e .
```

**disputauble** (already installed):
```bash
pip install -v git+https://github.com/NoahSailer/disputauble --user
```

### Selecting a theory code

Run configs in `runs/` select their theory code via `!defaults`:

```yaml
# Use standard CAMB
theory: !defaults [../base_config/camb]

# Use CAMB_alens (A_lens support)
theory: !defaults [../base_config/camb_alens]

# Use disputauble (negative mnu)
theory: !defaults [../base_config/camb_disputauble]

# Use CLASS_SZ
theory: !defaults [../base_config/class_sz]
```

## Project Structure

```
base_config/          # Composable YAML building blocks
  camb.yaml           # Standard CAMB theory settings
  camb_alens.yaml     # CAMB_alens theory settings
  camb_disputauble.yaml  # disputauble theory settings
  class_sz.yaml       # CLASS_SZ theory settings
  lensing.yaml        # Lensing likelihood settings
  lcdm_param.yaml     # LCDM parameter priors
  mcmc_sampler.yaml   # Sampler settings
runs/                 # Complete run configurations
  cobaya.sh           # SLURM submission script
  act_lcdm.yaml       # Example: lensing + CLASS_SZ
src/dr6plus_lenslike/  # Main likelihood package
  dr6plus_lenslike.py # Likelihood implementation
  data/               # Bandpowers, covariance matrices, corrections
```

## Running

### Module deck (ComputeCanada/Alliance)

```bash
module load StdEnv/2023
module load aocl-lapack/5.1
module load openblas
module load gsl          # Required for CAMB_alens (CosmoRec/GSL)
module load openmpi
module load fftw
module load cfitsio
module load python
```

### Submit a run

```bash
sbatch runs/cobaya.sh
```

Edit the last line of `runs/cobaya.sh` to point to the desired run config.
