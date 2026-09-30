# Verification chains on the cosmo2017 mocks

Before running the Legacy extended-model chains, we reproduce a set of reference
chains in which every data vector is a noiseless mock evaluated at one known
cosmology, theta_fid (the `cosmo2017_10K_acc3` fiducial of the DR6 lensing
simulations and of the lensing normalization). A correct likelihood, data package
and theory setup must then recover theta_fid: in our runs every parameter that the
data constrain lies within 0.23 sigma of its input value, and the best-fit sample
reaches chi2 <= 0.08 per likelihood. This page lists what to run, the numbers each
step must reproduce, and the figures.

**Status (2026-09-30).** This README and the figures in `plots/` are on the
branch. The scripts and configs described below (`make_configs.py`,
`tier0_check.py`, `cobaya.sh`, `debug_test.sh`, `recovery.py`,
`plot_triangles.py`, `configs/`) and the `map_frame_cal` option of the lensing
likelihood, which the joint configs use, are not yet committed: the map-frame
calibration term is under review. Our chains and figures were produced from
commit `2f41b62` plus that uncommitted diff, recorded with md5 checksums in
`/scratch/jiaqu/chains/verification/PROVENANCE.txt`.

## What to reproduce

Six reference chains, for the two lensing variants `<v>` = `actbase`, `act_planck`:

| config | likelihoods | theory | run time (one node) |
|---|---|---|---|
| `<v>_lens_classsz` | lensing alone | class_sz | about 4 min |
| `<v>_lensbao_classsz` | lensing + DESI DR2 BAO | class_sz | about 4 min |
| `<v>_lenscmbbao_camb` | lensing + P-ACT lite + DESI DR2 BAO | CAMB | 6.6-7.1 h |

Without primary CMB the class_sz emulator is the default theory; with the CMB,
CAMB. Two further sets exist and need not be reproduced: `<v>_lens_camb` and
`<v>_lensbao_camb` (CAMB cross-check of the class_sz chains, 1.3-2.1 h each), and
`<v>_lenscmbbao_mnu_camb` (the joint chain with Sum m_nu free, which provides the
Sum m_nu posterior; results pending, jobs 2466559 and 2466560).

## Input cosmology and mock data

The mock point, in the parametrization of the chains
(`tests/mock/cosmo2017/theta_fid.yaml`), with the derived quantities the chains
are compared against (CAMB, ACT DR6 reference theory):

| Omega_b h^2 | Omega_c h^2 | 100 theta_MC | ln(10^10 A_s) | n_s | tau | Sum m_nu |
|---|---|---|---|---|---|---|
| 0.02219218 | 0.1203058 | 1.0407345 | 3.068453 | 0.9625356 | 0.06574325 | 0.06 eV, one massive state |

| H_0 [km/s/Mpc] | sigma_8 | Omega_m | S_8^CMBL = sigma_8 (Omega_m/0.3)^0.25 | r_drag [Mpc] |
|---|---|---|---|---|
| 67.024 | 0.82196 | 0.31865 | 0.83444 | 147.22 |

class_sz reports Omega_m without the massive neutrino; in its convention the same
point has Omega_m = 0.31721, sigma_8 = 0.82151 and S_8^CMBL = 0.83304, and
`recovery.py` compares the class_sz chains with these values.

The data are four mocks, each the corresponding likelihood's own prediction at
theta_fid, with the real covariances unchanged (construction and checksums in
`tests/mock/README.md`):

1. DR6+ lensing bandpowers, `tests/mock/cosmo2017/clkk_bandpowers_mock.txt` (the
   binned fiducial C_L^kk, 18 bins; the likelihood keeps 40 < L < 1100). The
   variants `ACTbase` (DR6plus_lensing 9h: ACT-only HILC, 40 patches, plus the two
   daytime legs) and `ACT_Planck` (9a: ACT+Planck HILC plus the daytime legs) are
   blinded, and the likelihood refuses their measured bandpowers unless
   `mock: true`.
2. P-ACT lite high-ell: ACT DR6 CMB-only (TT, TE, EE for ell >= 600) and Planck
   2018 plik-lite (TT to ell = 1000, TE and EE to ell = 600), with A_act = A_planck
   = P_act = 1.
3. DESI DR2 BAO, mock mean with the real DR2 covariance.
4. In place of the real low-ell likelihoods, a Gaussian prior tau ~ N(0.0657,
   0.0065) centred on the input, with the width of the SRoll2 EE likelihood. The
   real SRoll2 EE data give tau = 0.0565 +/- 0.0065 at theta_fid, 1.4 sigma below
   the input, and would pull tau and A_s away from it. (The statement in
   `tests/mock/README.md` that the verification runs use the real low-ell data
   predates this choice.)

## Likelihood settings and priors

Lensing alone and lensing + BAO use the likelihood in lens-only mode with the
variant's analytic CMB-marginalized covariance (`lens_only`, `analytic_marg`,
`cov_file: null`) and the priors of the ACT DR6 lensing runs: ln(10^10 A_s) in
[2, 4], n_s ~ N(0.9625, 0.02), Omega_b h^2 ~ N(0.02219, 0.0005) (both Gaussians
centred on the input), Omega_c h^2 in [0.005, 0.99], H_0 in [40, 100], tau and
Sum m_nu fixed at the input. The class_sz block is the fast-mode (emulator) block
of the ACT DR6 lensing runs (one massive neutrino, N_ur = 2.0328).

The joint configs use the full lensing likelihood, in which the per-variant
normalization and N1 response is applied to the map-frame CMB spectra,
C_th/A_act^2 for TT, /(A_act^2 P_act) for TE and /(A_act^2 P_act^2) for EE
(`map_frame_cal: true`), with the priors of the P-ACT LCDM runs: Omega_b h^2 in
[0.017, 0.027], Omega_c h^2 in [0.09, 0.15], theta_MC in [0.0103, 0.0105],
ln(10^10 A_s) in [2.6, 3.5], n_s in [0.9, 1.1], the tau prior above, A_planck in
[0.5, 1.5] with A_act = A_planck and the calibration prior A_act ~ N(1, 0.003),
and P_act in [0.9, 1.1]. Sum m_nu is fixed at 0.06 eV, or free in [0, 5] eV in
the `_mnu` configs. The CAMB theory is the ACT DR6 reference block the CMB and BAO
mocks were made with (`tests/mock/cosmo2017/inputs/act_dr6_theory_camb.yaml`:
CosmoRec, `lens_potential_accuracy` 8, halofit `mead2020`, one massive neutrino,
N_eff = 3.044); it requires a CAMB build with CosmoRec.

## How to run

From the repository root, with the cluster environment loaded (module block and
`source ~/.bashrc`; on Trillium see `cobaya.sh`):

1. Write the configs for your checkout and output directory (all fourteen are
   written; only the six above need to be run):

       python runs/verification/make_configs.py --out /scratch/$USER/chains/verification

2. Tier 0, one likelihood evaluation of each config at theta_fid (about 1 min on
   one node; run it on a compute or debug allocation):

       COBAYA_NOMPI=1 python runs/verification/tier0_check.py runs/verification/configs/{actbase,act_planck}_{lens_classsz,lensbao_classsz,lenscmbbao_camb}.yaml

   It first checks that every data path of every config lies in `tests/mock/`
   (no real data vector anywhere), then prints the chi2 of each likelihood, the
   derived parameters at theta_fid, and `Tier-0: PASS` or `FAIL`.
3. Optionally, `sbatch runs/verification/debug_test.sh $PWD` repeats Tier 0 for
   all configs and runs a 10-sample chain of one class_sz and one joint config on
   the debug partition.
4. Submit each of the six chains (one node, 8 MPI chains of 24 threads;
   resubmitting the same command resumes it; edit the `#SBATCH` account and log
   path for your allocation):

       for c in {actbase,act_planck}_{lens_classsz,lensbao_classsz,lenscmbbao_camb}; do
           sbatch --job-name=$c runs/verification/cobaya.sh $PWD/runs/verification/configs/$c.yaml
       done

   Each job runs the mock-only check of step 2 before `cobaya-run` and stops if it
   fails.
5. Compare with the input and plot (set the chain directory at the top of each
   script):

       python runs/verification/recovery.py
       python runs/verification/plot_triangles.py

## What to expect

**Tier 0.** At theta_fid the P-ACT lite and BAO mocks give chi2 <= 7e-5 (the real
data give about 219, 176 and 37 at this point, so this also confirms that the
mock files were read). The lensing chi2 is not zero, because the lensing mock is
the 2017 fiducial spectrum while the chains compute C_L^kk with a current theory
code; the values must be reproduced to 1e-3:

| config | ACTbase | ACT_Planck |
|---|---|---|
| lens_classsz, lensbao_classsz | 0.0719 | 0.0981 |
| lenscmbbao_camb (and _mnu) | 0.0349 | 0.0477 |
| lens_camb, lensbao_camb (cross-check) | 0.0224 | 0.0303 |

For the joint configs the check also sets A_planck = 1.003, which must raise the
lensing chi2 by 0.62 (ACTbase) and 0.76 (ACT_Planck); this is the map-frame term.

**Convergence.** All chains reach R-1 = 0.003-0.010. At R-1 = 0.01 the chain means
carry a sampling scatter of about 0.1 sigma, which sets the tolerance for
reproducing the numbers below.

**Posteriors** (getdist, first 30 per cent of each chain removed; shift = (mean
- input)/sigma, with the class_sz chains compared with the class_sz-convention
inputs).

Lensing alone. Only S_8^CMBL is constrained; sigma_8, Omega_m and H_0 follow the
priors along the lensing degeneracy (sigma(H_0) = 17 km/s/Mpc) and are not a test.

| config | S_8^CMBL | input | shift |
|---|---|---|---|
| actbase_lens_classsz | 0.8319 +/- 0.0145 | 0.8330 | -0.08 |
| act_planck_lens_classsz | 0.8321 +/- 0.0131 | 0.8330 | -0.07 |
| actbase_lens_camb (cross-check) | 0.8316 +/- 0.0148 | 0.8344 | -0.19 |
| act_planck_lens_camb (cross-check) | 0.8317 +/- 0.0136 | 0.8344 | -0.20 |

Lensing + BAO:

| config | sigma_8 | Omega_m | H_0 | S_8^CMBL |
|---|---|---|---|---|
| actbase_lensbao_classsz | 0.8232 +/- 0.0104 (+0.16) | 0.3180 +/- 0.0082 (+0.10) | 67.06 +/- 0.55 (+0.06) | 0.8352 +/- 0.0107 (+0.21) |
| act_planck_lensbao_classsz | 0.8232 +/- 0.0096 (+0.18) | 0.3180 +/- 0.0078 (+0.10) | 67.06 +/- 0.55 (+0.07) | 0.8352 +/- 0.0096 (+0.23) |
| actbase_lensbao_camb (cross-check) | 0.8222 +/- 0.0104 (+0.03) | 0.3194 +/- 0.0080 (+0.10) | 67.06 +/- 0.54 (+0.06) | 0.8352 +/- 0.0108 (+0.07) |
| act_planck_lensbao_camb (cross-check) | 0.8225 +/- 0.0098 (+0.06) | 0.3195 +/- 0.0079 (+0.10) | 67.05 +/- 0.54 (+0.05) | 0.8355 +/- 0.0099 (+0.11) |

Lensing + P-ACT lite + BAO (CAMB), Sum m_nu fixed; all sampled parameters lie
within 0.11 sigma of the input:

| parameter | input | actbase_lenscmbbao_camb | act_planck_lenscmbbao_camb |
|---|---|---|---|
| S_8^CMBL | 0.83444 | 0.8350 +/- 0.0056 (+0.11) | 0.8350 +/- 0.0053 (+0.10) |
| sigma_8 | 0.82196 | 0.8224 +/- 0.0047 (+0.10) | 0.8224 +/- 0.0045 (+0.09) |
| Omega_m | 0.31865 | 0.3188 +/- 0.0037 (+0.04) | 0.3188 +/- 0.0036 (+0.05) |
| H_0 | 67.025 | 67.02 +/- 0.26 (-0.04) | 67.01 +/- 0.25 (-0.04) |
| Omega_b h^2 | 0.022192 | 0.022191 +/- 0.000102 (-0.01) | 0.022190 +/- 0.000102 (-0.02) |
| Omega_c h^2 | 0.12031 | 0.12033 +/- 0.00063 (+0.04) | 0.12034 +/- 0.00062 (+0.05) |
| ln(10^10 A_s) | 3.0685 | 3.0694 +/- 0.0114 (+0.08) | 3.0692 +/- 0.0112 (+0.07) |
| n_s | 0.96254 | 0.96258 +/- 0.00308 (+0.01) | 0.96256 +/- 0.00306 (+0.01) |
| tau | 0.0657 | 0.0661 +/- 0.0056 (+0.06) | 0.0660 +/- 0.0056 (+0.04) |
| A_planck | 1 | 1.0001 +/- 0.0029 (+0.02) | 1.0001 +/- 0.0029 (+0.03) |
| P_act | 1 | 1.0000 +/- 0.0024 (+0.01) | 1.0000 +/- 0.0024 (+0.01) |

Three points on interpretation. (i) With omega_nu = 0.06/93.14 added back to the
class_sz Omega_m (as in the figures, so that one set of input lines applies to
every chain), the class_sz and CAMB means agree to 0.07-0.10 sigma in S_8^CMBL and
to 0.01 sigma in Omega_m, and the class_sz S_8^CMBL lies -0.11 sigma (lensing
alone) and +0.16 to +0.18 sigma (with BAO) from the CAMB input 0.83444. (ii) With
Sum m_nu free, the posterior mean of Sum m_nu lies above 0.06 eV because of the
prior boundary at zero, and H_0, Omega_m and sigma_8 shift along their
degeneracies with it; the fixed-Sum m_nu chains are the recovery test. (iii)
Because A_act = 1 in the mocks, these chains test the recovery of the input but
not the calibration frame of the joint correction, which changes the coupling
between A_act and the lensing amplitude and hence the posterior widths, not the
means.

## Figures

Input values (CAMB convention) as dashed lines; filled contours with stroked 68
and 95 per cent edges against unfilled contours (`plot_triangles.py`).

Lensing alone, ACTbase against ACT+Planck (class_sz), `plots/verification_lens.pdf`:

![lensing alone](plots/verification_lens.png)

Lensing + DESI DR2 BAO (class_sz), `plots/verification_lensbao.pdf`:

![lensing + BAO](plots/verification_lensbao.png)

Lensing + P-ACT lite + DESI DR2 BAO (CAMB) with Sum m_nu fixed at 0.06 eV,
`plots/verification_lenscmbbao_mnufixed.pdf`:

![lensing + P-ACT + BAO](plots/verification_lenscmbbao_mnufixed.png)

Cross-check, CAMB against class_sz for ACTbase, lensing alone and with BAO,
`plots/verification_classsz_vs_camb_lens.pdf` and
`plots/verification_classsz_vs_camb_lensbao.pdf`:

![class_sz vs CAMB, lensing](plots/verification_classsz_vs_camb_lens.png)
![class_sz vs CAMB, lensing + BAO](plots/verification_classsz_vs_camb_lensbao.png)

The figure with Sum m_nu free (`plots/verification_lenscmbbao.pdf`) is added once
the `_mnu` chains finish.
