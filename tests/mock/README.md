# Verification mocks: cosmo2017

Noiseless mock data vectors for the reference (verification) chains of the DR6+
lensing likelihood: lensing alone, lensing + BAO, and lensing + primary CMB + BAO.
Every mock is evaluated at one cosmology, `cosmo2017_10K_acc3`, the fiducial of
the DR6 lensing simulations and of the lensing normalization, so the mock sky has
the fiducial CMB spectra C_fid and the mock lensing bandpowers carry no
normalization or N1 response. This sets to zero only the lens-only Eq. 35
recentering, M_b (C_data - C_fid), since for the mock sky C_data = C_fid. The likelihood
corrections themselves stay on: in joint (primary CMB + lensing) runs the theory
keeps the normalization and N1 correction M_b (C_th - C_fid), which accounts for
the difference between the sampled cosmology and C_fid at every step of the
chain and vanishes only at the fiducial point itself. This file records how each
mock was made, from which inputs, and with which theory settings, so that the
products can be regenerated and checked independently.

## Contents

| Path | Content | Generator |
|---|---|---|
| `cosmo2017/inputs/cosmo2017_10K_acc3_params.ini` | CAMB parameter file of the fiducial (defines the cosmology) | copied |
| `cosmo2017/inputs/cosmo2017_10K_acc3_lenspotentialCls.dat` | 2017 CAMB lensing-potential spectrum of the fiducial | copied |
| `cosmo2017/inputs/act_dr6_theory_camb.yaml` | ACT DR6 reference CAMB block | copied |
| `cosmo2017/inputs/act_dr6_likelihood_p_act_lite.yaml` | ACT DR6 P-ACT lite likelihood block | copied |
| `cosmo2017/clkk_bandpowers_mock.txt` | lensing mock, 18 DR6 bins | `make_mock_lensing.py` |
| `cosmo2017/cmb/ACTDR6CMBonly/v1.0/dr6_data_cmbonly_mock.fits` | ACT DR6 CMB-only mock (sacc) | `make_mock_cmb_bao.py` |
| `cosmo2017/cmb/planck_2018_pliklite_mock/` | Planck 2018 plik-lite mock (all 613 bins) | `make_mock_cmb_bao.py` |
| `cosmo2017/bao/desi_gaussian_bao_ALL_GCcomb_mean_mock.txt` | DESI DR2 BAO mock (real DR2 covariance beside it) | `make_mock_cmb_bao.py` |
| `cosmo2017/theta_fid.yaml` | the mock point in the chains' parametrization | `make_mock_cmb_bao.py` |
| `cosmo2017/PROVENANCE_lensing.json`, `PROVENANCE_cmb_bao.json` | input, output and script md5s, settings, checks | both |

The low-ell CMB likelihoods (`planck_2018_lowl.TT`, `planck_2018_lowl.EE_sroll2`)
are not mocked: the verification runs use the real low-ell data.

**Status.** The lensing-likelihood side referred to below (the `mock_file` option,
the mock-mode Eq. 35 change, and `tests/test_mock_products.py`) is not yet in the
repository. Until it is, `mock: true` uses the packaged
`clkk_bandpowers_fiducial.txt`, which equals `clkk_bandpowers_mock.txt` to one unit
in the last place. The two generators need only the files in this folder and the
tracked `binning_matrix_act.txt` / `clkk_bandpowers_fiducial.txt`.

## Copied inputs

All four files are verbatim copies (md5 identical to their sources):

| File | md5 | Source |
|---|---|---|
| `cosmo2017_10K_acc3_params.ini` | f70b827851598d5f409703cef54e4af9 | `actsims/data/cosmo2017_10K_acc3_params.ini` |
| `cosmo2017_10K_acc3_lenspotentialCls.dat` | a8293b03926e3687a1fa9c68bd1daf88 | `dr6plus_lenslike/data/v1.0/like_corrs/` (a link to the ACT DR6 lensing likelihood data, `act_dr6_spt_lenslike/data/v1.2/like_corrs/`) |
| `act_dr6_theory_camb.yaml` | bb2bea5ed1a121349074193ae0088094 | DR6-ACT-lite `yamls/theories/camb.yaml`, commit 0e0cd2c |
| `act_dr6_likelihood_p_act_lite.yaml` | d553cc8fef2391838511ae500bc963f5 | DR6-ACT-lite `yamls/likelihoods/p_act_lite.yaml`, commit 0e0cd2c |

## Mock cosmology

The six base parameters are those of the `.ini` and are common to all mocks.
The neutrino and helium sector follows the code that produced each mock:

| Parameter | Lensing mock (2017 CAMB file) | CMB and BAO mocks (CAMB 1.6.2) |
|---|---|---|
| H0 | 67.02393 | 67.02393 |
| ombh2 | 0.02219218 | 0.02219218 |
| omch2 | 0.1203058 | 0.1203058 |
| A_s | 2.15086031154146e-9 (ln(10^10 A_s) = 3.068453) | same |
| n_s | 0.9625356 | 0.9625356 |
| tau | 0.06574325 | 0.06574325 |
| Sum m_nu | one massive state, omnuh2 = 0.00064 (0.0595 eV) | one massive state, 0.06 eV |
| N_eff | 3.046 | 3.044 |
| YHe | 0.2448964 (fixed) | 0.24578 (BBN) |

The differences in the last three rows change sigma8 by 3e-5 and the lensing
chi2 at the mock point by less than 0.002 (see the last section): the mocks
describe the same cosmology.

## Lensing bandpowers

**Source spectrum.** `cosmo2017_10K_acc3_lenspotentialCls.dat` (L = 2 to 10000;
columns L, TT, EE, BB, TE, PP, TP, EP, with PP = [L(L+1)]^2 C_L^phiphi / 2 pi) is
the output of a 2017 Fortran CAMB run with the `.ini` above: `accuracy_boost = 3`,
`l_accuracy_boost = 3`, `high_accuracy_default = T`, `do_nonlinear = 3` (no
`halofit_version` set, so the halofit default of that CAMB release),
`lensing_method = 1`, RECFAST recombination. It is the spectrum behind the DR6
lensing simulations and normalization.

**Procedure** (`make_mock_lensing.py`, no CAMB call):

    C_L^kk = PP_L x 2 pi / 4 = [L(L+1)]^2 C_L^phiphi / 4
    b_i    = sum_{L=0}^{2999} B_iL C_L^kk

with B the DR6 binning matrix `dr6plus_lenslike/data/v1.0/binning_matrix_act.txt`
(18 x 3000, column j is L = j; md5 6c8f94dd257020394e9d7f12f0c7ed0b). The sums
use `math.fsum` (exactly rounded), so the output does not depend on the BLAS
library or the thread count. The likelihood uses bins 2 to 13 (40 < L < 1100).

**Check.** The mock (md5 c3dc66b233bb602a54c7fdab4cc83c48) reproduces the
packaged `dr6plus_lenslike/data/v1.0/clkk_bandpowers_fiducial.txt`: 10 of the 18
bins bit for bit and the other 8 to one unit in the last place (2.2e-16
relative). An exactly rounded sum cannot do better, so the residual comes from
the summation order of the original, not from the content.

**Use in a run** (lensing likelihood options):

```yaml
mock: true
mock_file: <repo>/tests/mock/cosmo2017/clkk_bandpowers_mock.txt
```

In lens-only mode with `mock: true` the Eq. 35 recentering is zero (the mock sky
has the fiducial CMB) and the CMB-marginalized covariance is kept. In joint mode
(`lens_only: false`) nothing changes with `mock: true`: the theory carries the
normalization and N1 corrections M_b (C_th - C_fid) evaluated at each sampled
cosmology, with the uninflated covariance, exactly as for data.

## Primary CMB and BAO

**Theory.** The ACT DR6 reference CAMB block, `inputs/act_dr6_theory_camb.yaml`,
with no additions:

| extra_arg | value |
|---|---|
| kmax | 10 |
| k_per_logint | 130 |
| nonlinear | true |
| lens_potential_accuracy | 8 |
| lens_margin | 2050 |
| lAccuracyBoost | 1.2 |
| min_l_logl_sampling | 6000 |
| DoLateRadTruncation | false |
| recombination_model | CosmoRec |

Everything the block does not set takes the CAMB default, recorded in
`PROVENANCE_cmb_bao.json` as used: `nnu = 3.044`, `num_massive_neutrinos = 1`,
`halofit_version = mead2020`, BBN helium. The multipole range is set by the
likelihood requirements (`lmax_theory = 9000` for ACT). CAMB is version 1.6.2, as
resolved by cobaya from its packages path
(`/scratch/jiaqu/likelihood_data/code/CAMB`, built with CosmoRec).

**Likelihoods.** The high-ell entries of `inputs/act_dr6_likelihood_p_act_lite.yaml`
(`act_dr6_cmbonly.ACTDR6CMBonly` and `act_dr6_cmbonly.PlanckActCut`) and
`bao.desi_dr2.desi_bao_all`, with the nuisance parameters at A_act = A_planck = 1
and P_act = 1.

**Real data files replaced** (all from the cobaya packages path
`/scratch/jiaqu/likelihood_data/data/`; md5s of every file in `PROVENANCE_cmb_bao.json`):

| Likelihood | File | md5 |
|---|---|---|
| ACT DR6 CMB-only | `ACTDR6CMBonly/v1.0/dr6_data_cmbonly.fits` (public v1.0) | b615072fd09d5f5c614505a97eaacf33 |
| Planck plik-lite | `planck_2018_pliklite_native/cl_cmb_plik_v22.dat` | 690b6b708572bf00486db3250ea988d8 |
| DESI DR2 BAO | `bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt` | 9c2d2c029e7cd58562cb5b9bf2abf5ff |

**Procedure** (`make_mock_cmb_bao.py`, one cobaya model evaluated at the mock
point). Each mock value is the likelihood's own prediction, written into a copy
of the real file; covariances, windows and binning are unchanged.

- ACT: every point of the sacc file, b = W D_ell (A_act = P_act = 1), with W the
  point's bandpower window and D_ell = ell(ell+1) C_ell / 2 pi.
- Planck: every bin of TT (215), TE (199) and EE (199),
  b = sum over [blmin, blmax] of w_ell D_ell, with the plik-lite weights
  (A_planck = 1).
- BAO: every row (D_V/r_d, D_M/r_d, D_H/r_d at the DR2 redshifts), from the cobaya
  BAO likelihood's `theory_fun`; the real DR2 covariance is copied beside it.

**Derived at the mock point:** cosmomc_theta = 0.01040734, sigma8 = 0.82196,
r_drag = 147.2184 Mpc, YHe = 0.24578. `theta_fid.yaml` gives the point in the
chains' parametrization (cosmomc_theta, logA).

**Checks.**

- Swap test: the likelihoods read from the mock files give chi2 = 6e-27 (Planck),
  0 (ACT), 3e-21 (BAO) at the mock point; the real data give 219.8, 174.7, 36.6.
- The ACT covariance in the mock sacc file equals the real one exactly.
- The sacc FITS writer stores Python object ids and hash-ordered table columns, so
  the bytes of `dr6_data_cmbonly_mock.fits` change on every regeneration while its
  content does not. Its reproducible checksum is the md5 of the mean vector
  (float64, little endian), `act_mock_mean_md5` = b30d95344a68c835a17b8ce4192fd49c.
  The Planck and BAO files are byte-reproducible (md5 c649673b17d144d54577973bcddb3ddd
  and 7388e1ba63ddc09d0c68e83263c7de5f).

**Use in a run.** Include the same theory block (`inputs/act_dr6_theory_camb.yaml`)
and override the three likelihoods (the rest of each block as in the P-ACT lite
yaml):

```yaml
act_dr6_cmbonly.ACTDR6CMBonly:
  data_folder: <repo>/tests/mock/cosmo2017/cmb/ACTDR6CMBonly
  version: v1.0
  input_file: dr6_data_cmbonly_mock.fits
act_dr6_cmbonly.PlanckActCut:
  path: <repo>/tests/mock/cosmo2017/cmb/planck_2018_pliklite_mock
bao.desi_dr2.desi_bao_all:
  path: <repo>/tests/mock/cosmo2017/bao
  measurements_file: desi_gaussian_bao_ALL_GCcomb_mean_mock.txt
  cov_file: desi_gaussian_bao_ALL_GCcomb_cov.txt
```

`base_config/camb_pact.yaml` is not equivalent to the reference block: it sets
N_eff = 3.046 with three degenerate massive states, which changes the high-ell
CMB chi2 at the mock point by 0.26 (Planck) + 0.90 (ACT).

## Tier 0 and the known lensing offset

`tests/test_mock_products.py::test_tier0_all_likelihoods_at_mock_point` evaluates
every verification likelihood on the mock files at `theta_fid.yaml`, through
cobaya, with the reference theory block:

| Likelihood | chi2 at the mock point |
|---|---|
| Planck plik-lite / ACT CMB-only / DESI DR2 | 1.4e-5 / 4.3e-5 / 3.6e-6 |
| Lensing dr6plus_variant0, joint / lens-only | 0.0286 / 0.0185 |
| Lensing dr6plus_optimal, joint / lens-only | 0.0479 / 0.0304 |

The CMB and BAO residuals come from the cosmomc_theta to H0 round trip (67.02455
vs 67.02393). The lensing chi2 is a theory-code offset: the lensing mock is the
2017 CAMB spectrum, while the chains predict C_L^kk with CAMB 1.6.2 at the same
parameters. In joint mode the prediction also includes the correction
M_b (C_th - C_fid), which at the mock point is small but not zero, because the
CAMB 1.6.2 spectra C_th differ slightly from the 2017 lensed spectra that define
C_fid; this, together with the different covariances of the two modes, is why
the joint chi2 exceeds the lens-only one. In lens-only mode the offset splits as
follows:

| Theory for the lensing prediction | chi2 variant0 | chi2 optimal |
|---|---|---|
| CMB/BAO mock conventions, mead2020 | 0.018 | 0.030 |
| 2017 `.ini` neutrinos and helium, mead2020 | 0.019 | 0.032 |
| 2017 `.ini` neutrinos and helium, halofit takahashi | 0.005 | 0.010 |

About 70 per cent is the nonlinear prescription (mead2020 against the 2017
halofit default) and the rest the 2017 code and accuracy settings; the
neutrino and helium conventions contribute less than 0.002. In the linear
approximation a data offset with chi2 = X can move any parameter by at most
sqrt(X) of its standard deviation, so the lensing-driven shift in the closure
chains is bounded by 0.14 (variant0) and 0.17 (optimal) sigma in lens-only mode
and by 0.17 and 0.22 sigma in joint mode.

## Reproduce

From the repo root, with the cluster environment loaded
(`module load StdEnv/2023 aocl-lapack/5.1 openblas gsl openmpi fftw cfitsio python`,
`source ~/.bashrc`, cobaya on `PYTHONPATH`):

```bash
python tests/mock/make_mock_lensing.py
COBAYA_NOMPI=1 MPI4PY_RC_INITIALIZE=false python tests/mock/make_mock_cmb_bao.py
COBAYA_NOMPI=1 python -m pytest -p no:flaky tests/test_mock_products.py
```

`make_mock_cmb_bao.py` makes two CAMB calls at the reference accuracy and runs in
about 10 s on a login node. Both generators exit with an error if their checks
fail (reproduction of the packaged lensing fiducial; zero swap-test chi2).
