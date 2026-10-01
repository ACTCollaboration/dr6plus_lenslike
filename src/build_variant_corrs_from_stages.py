"""Package the per-variant CMB-correction products of a DR6+ patch-combination
variant for dr6plus_lenslike, from the production mbatch stage outputs
(stage_bpshift, stage_combination) of DR6plus_lensing.

Every leg of the combination is N1-deconvolved and the combination is linear in
the leg bandpowers with fixed per-bin weights, so its binned response to the
lensed CMB spectra is

    M^X_b = w_night,b sum_p w_p,b [M_norm^X + M_n1^X]_{p,b}
          + w_dd,b [M_norm^X + M_n1^X]_{dd,b} + w_dw,b [M_norm^X + M_n1^X]_{dw,b}

with [M_norm^X + M_n1^X] = B A (R_norm^X + R_N1^X) the per-region rows written by
stage_bpshift (profile-hardened normalization derivative, plain-MV N1
derivative of each patch, each region's MCN1-matched deconvolution operator A),
w_p the inverse-variance patch weights (patch_comb.npz weights_invvar) and w
the per-bin GLS-LOO leg weights (patch_comb_nightday_gls-loo.npz weights).

Written to data/v1.0/ (18 bins of binning_matrix_act.txt; the likelihood slices
its band):
  clkk_<v>.txt, covmat_clkk_<v>.txt   uncorrected N1-deconvolved bandpowers
                        (BLINDED, never printed) and their raw sim-bank covariance
  like_corrs_<v>/response_binned.npz  Mb (3, 18, 3000), X = TT, EE, TE, and its
                        norm and N1 parts; acts on C^X_th / cal - C^X_fid in full mode
  like_corrs_<v>/covmat_cmbmarg_analytic.txt  sim-bank covariance + Eq. 34 addition
  like_corrs_<v>/lens_only_delta.txt  Eq. 35 recentering Delta_b = sum_X Mb^X dC^X
                        (= the production combined Delta), 18 bins
  like_corrs_<v>/PROVENANCE.json      sources and md5s, conventions, checks

Hard checks (AssertionError): per region, night leg, day legs and combination,
the rows reproduce the stored Delta; the combined Delta equals the packaged one;
the covariance equals the packaged one.

  python src/build_variant_corrs_from_stages.py --variant actbase
  python src/build_variant_corrs_from_stages.py --variant act_planck
"""
import argparse
import datetime
import hashlib
import json
import os
import sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(REPO, 'src'))
import pact_cmbmarg_actplanck as pm          # noqa: E402
import like_corrs_variant_analytic as lva    # noqa: E402

DDIR = os.path.join(REPO, 'src', 'dr6plus_lenslike', 'data', 'v1.0')
PKG = '/project/rrg-rbond-ac/jiaqu/dr6plus_lensing_variants/bandpowers_bpshift_corrected'
SACC = '/scratch/jiaqu/likelihood_data/data/ACTDR6CMBonly/v1.0/dr6_data_cmbonly.fits'
SPECS = ('TT', 'EE', 'TE')
VARIANTS = {
    'actbase': dict(pkg='dr6plus_variant9h', cmb='act',
                    cmb_data='ACT DR6 CMB-only public v1.0 (dr6_data_cmbonly.fits), TT/TE/EE ell >= 600',
                    note='variant 9h: ACT-only HILC 40 patches (lmin 600) + daydeep + daywide, GLS-LOO'),
    'act_planck': dict(pkg='dr6plus_variant9a', cmb='pact',
                       cmb_data='P-ACT lite: plik_lite v22 PlanckActCut (TT<1000, TE/EE<600; dC below ell 600) '
                                '+ ACT DR6 CMB-only public v1.0 (ell >= 600)',
                       note='variant 9a: ACT+Planck (Planck-inpainted) T+P HILC 40 patches (lmin 100), '
                            'no strip infill, + daydeep + daywide, GLS-LOO'),
}


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def rows(f):
    """(M_norm, M_n1) dicts of (18, 3000) rows, dC dict and delta of a bpshift npz."""
    z = np.load(f, allow_pickle=True)
    Mn = {X: z[f'M_norm_{X}'] for X in SPECS}
    M1 = {X: z[f'M_n1_{X}'] for X in SPECS}
    dC = {X: z[f'dC_{X}'] for X in SPECS}
    return Mn, M1, dC, z['delta']


def close(a, b, tol=1e-12):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--variant', required=True, choices=sorted(VARIANTS))
    ap.add_argument('--ddir', default=DDIR)
    a = ap.parse_args()
    V = VARIANTS[a.variant]
    pkg = os.path.join(PKG, V['pkg'])
    stg = os.path.join(pkg, 'stages')
    checks, inputs = {}, {}

    # ---- combination structure
    fg = os.path.join(stg, 'stage_combination', 'patch_comb_nightday_gls-loo.npz')
    fpc = os.path.join(stg, 'stage_combination', 'patch_comb.npz')
    g, pc = np.load(fg, allow_pickle=True), np.load(fpc, allow_pickle=True)
    W, legs = g['weights'], [str(x) for x in g['names']]
    day_files = [str(f) for f in g['bpshift_files'][1:]]
    assert len(legs) == 3 and W.shape == (3, 18), (legs, W.shape)
    inputs.update({'gls_combination': fg, 'patch_combination': fpc})

    # ---- night leg: inverse-variance sum over the 40 regions
    wp, regions = pc['weights_invvar'], [int(r) for r in pc['regions']]
    Mn_night = {X: np.zeros((18, 3000)) for X in SPECS}
    M1_night = {X: np.zeros((18, 3000)) for X in SPECS}
    d_night, worst_region = np.zeros(18), 0.
    for i, r in enumerate(regions):
        f = os.path.join(stg, 'stage_bpshift', f'bpshift_region{r}.npz')
        Mn, M1, dC, dlt = rows(f)
        worst_region = max(worst_region, close(sum((Mn[X] + M1[X]) @ dC[X] for X in SPECS), dlt))
        for X in SPECS:
            Mn_night[X] += wp[i][:, None] * Mn[X]
            M1_night[X] += wp[i][:, None] * M1[X]
        d_night += wp[i] * dlt
        inputs[f'bpshift_region{r}'] = f
    checks['region_rows_reproduce_delta_max_rel'] = worst_region
    checks['night_delta_vs_patch_comb_max_rel'] = close(d_night, pc['delta_comb_invvar'])
    checks['night_delta_vs_gls_leg_max_rel'] = close(d_night, g['deltas'][0])

    # ---- day legs
    leg_n, leg_1, leg_d = [Mn_night], [M1_night], [d_night]
    for j, f in enumerate(day_files, start=1):
        Mn, M1, dC, dlt = rows(f)
        checks[f'{legs[j]}_rows_reproduce_delta_max_rel'] = close(sum((Mn[X] + M1[X]) @ dC[X] for X in SPECS), dlt)
        checks[f'{legs[j]}_delta_vs_gls_leg_max_rel'] = close(dlt, g['deltas'][j])
        leg_n.append(Mn), leg_1.append(M1), leg_d.append(dlt)
        inputs[f'bpshift_{legs[j]}'] = f

    # ---- GLS-LOO combination
    Mb_norm = {X: sum(W[l][:, None] * leg_n[l][X] for l in range(3)) for X in SPECS}
    Mb_n1 = {X: sum(W[l][:, None] * leg_1[l][X] for l in range(3)) for X in SPECS}
    Mb = {X: Mb_norm[X] + Mb_n1[X] for X in SPECS}
    delta = sum(W[l] * leg_d[l] for l in range(3))
    bp = np.load(os.path.join(pkg, 'bandpowers.npz'), allow_pickle=True)
    checks['combined_delta_vs_gls_delta_comb_max_rel'] = close(delta, g['delta_comb'])
    checks['combined_delta_vs_packaged_delta_max_rel'] = close(delta, bp['delta'])
    cov = np.loadtxt(os.path.join(pkg, 'covmat_clkk.txt'))
    checks['cov_vs_combination_cov_max_rel'] = close(cov, g['cov'])
    nsim = int(bp['nsim'])
    for k, v in checks.items():
        print(f'  {k}: {v:.2e}')
        assert v < 1e-10, f'check failed: {k} = {v}'

    # ---- Eq. 34 addition: ACT sacc (ell >= 600), + plik_lite blocks for P-ACT
    pm.SACC_FILE = SACC
    cov_sacc, sinfo = pm.sacc_blocks()
    add = lva.sacc_addition(Mb, cov_sacc, sinfo)
    inputs['act_sacc'] = SACC
    if V['cmb'] == 'pact':
        cov_plik, keep, lo, hi = pm.plik_bins()
        add = add + lva.plik_addition(Mb, cov_plik, keep, lo, hi)
    sym = float(np.max(np.abs(add - add.T)))
    ev = np.linalg.eigvalsh(add)
    broad = 100 * (np.sqrt(np.diag(cov + add) / np.diag(cov)) - 1)
    print('  sigma broadening % (40<L<1100):', np.round(broad[2:14], 2))

    # ---- write
    out = os.path.join(a.ddir, f'like_corrs_{a.variant}')
    os.makedirs(out, exist_ok=True)
    cents = bp['cents']
    np.savetxt(os.path.join(a.ddir, f'clkk_{a.variant}.txt'), bp['band'])
    np.savetxt(os.path.join(a.ddir, f'covmat_clkk_{a.variant}.txt'), cov)
    np.savez(os.path.join(out, 'response_binned.npz'), Mb=np.stack([Mb[X] for X in SPECS]),
             Mb_norm=np.stack([Mb_norm[X] for X in SPECS]), Mb_n1=np.stack([Mb_n1[X] for X in SPECS]),
             specs=np.array(SPECS), cents=cents, ells=np.arange(3000))
    np.savetxt(os.path.join(out, 'covmat_cmbmarg_analytic.txt'), cov + add,
               header=f'{a.variant}: sim-bank covariance (covmat_clkk_{a.variant}.txt) + Eq. 34 CMB addition; raw')
    np.savetxt(os.path.join(out, 'lens_only_delta.txt'), delta,
               header=f'{a.variant}: Eq. 35 recentering Delta_b (C_L^kk units), 18 bins')
    prov = {
        'variant': a.variant, 'upstream_variant': V['pkg'], 'note': V['note'],
        'written': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
        'script': os.path.relpath(__file__, REPO), 'script_md5': md5(__file__),
        'inputs': {k: {'path': f, 'md5': md5(f)} for k, f in inputs.items()},
        'package_meta': json.load(open(os.path.join(pkg, 'package_meta.json'))),
        'outputs': {os.path.relpath(p, REPO): md5(p) for p in
                    [os.path.join(a.ddir, f'clkk_{a.variant}.txt'), os.path.join(a.ddir, f'covmat_clkk_{a.variant}.txt')]
                    + [os.path.join(out, n) for n in ('response_binned.npz', 'covmat_cmbmarg_analytic.txt',
                                                      'lens_only_delta.txt')]},
        'nsim': nsim, 'legs': legs, 'norm_matrix': 'ph (profile-hardened bhph, stage_normder)',
        'n1_derivative': 'plain MV N1der per patch (stage_n1der), decision 2026-09-28',
        'form': 'N1-deconvolved: response B A R, no dN1/dC_kk term in the theory',
        'cmb_data': V['cmb_data'], 'bpshift_cmb': V['cmb'],
        'checks': {**checks, 'addition_symmetric_max_abs': sym,
                   'addition_min_eigenvalue_over_max': float(ev.min() / ev.max())},
        'peak_sigma_broadening_pct_40_1100': float(np.max(broad[2:14])),
    }
    with open(os.path.join(out, 'PROVENANCE.json'), 'w') as fh:
        json.dump(prov, fh, indent=1)
    print('wrote', out, 'and', f'clkk_{a.variant}.txt / covmat_clkk_{a.variant}.txt')


if __name__ == '__main__':
    main()
