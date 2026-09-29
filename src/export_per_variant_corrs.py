"""Package the per-variant CMB-correction products of a DR6+ likelihood variant
into data/v1.0/like_corrs_<variant>/, from the analytic-forms npz of the
like_corrs builders (built with the canonical profile-hardened, ph, norm
matrices; N1 derivative plain MV).

Written (18 bins of binning_matrix_act.txt; the likelihood slices its band):
  response_binned.npz   Mb (3, 18, 3000): binned response of the N1-deconvolved
                        bandpowers to the lensed CMB spectra, d(B C_kk)_b / dC^X_l,
                        X = TT, EE, TE (BB carries none), kappa units per muK^2,
                        l = 0..2999; acts on C^X_th / cal - C^X_fid in full mode
  covmat_cmbmarg_analytic.txt   base covariance + Eq. 34 addition (raw; Hartlap
                        applied by the likelihood)
  lens_only_delta.txt   Eq. 35 recentering Delta_b = sum_X Mb^X dC^X (data CMB minus
                        fiducial), added to B C_kk in lens_only analytic_marg mode
  PROVENANCE.json       sources, conventions, checks, md5 of inputs

  python src/export_per_variant_corrs.py --variant dr6plus_variant0 \\
      --forms .../variant0_cmbmarg_analytic_forms.npz --add-key add_dec --base-key base_dec
  python src/export_per_variant_corrs.py --variant dr6plus_optimal \\
      --forms .../variant9a_cmbmarg_analytic_forms.npz --add-key add_pact --base-key base
"""
import argparse
import datetime
import hashlib
import json
import os

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DDIR = os.path.join(REPO, 'dr6plus_lenslike', 'data', 'v1.0')
BASE_COV = {'dr6plus_variant0': 'covmat_clkk_dr6plus_variant0.txt',
            'dr6plus_optimal': 'covmat_clkk_hilcTP_nightday_glsloo.txt'}
SPECS = ('TT', 'EE', 'TE')


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--variant', required=True, choices=sorted(BASE_COV))
    ap.add_argument('--forms', required=True)
    ap.add_argument('--add-key', required=True, help='Eq. 34 addition in the forms npz (add_dec, add_pact)')
    ap.add_argument('--base-key', required=True, help='base covariance in the forms npz (base_dec, base)')
    ap.add_argument('--norm', default='ph', help='norm matrix convention the forms were built with (recorded)')
    ap.add_argument('--cmb-data', required=True, help='description of the CMB data behind Eq. 34 and dC (recorded)')
    ap.add_argument('--note', default='')
    ap.add_argument('--ddir', default=DDIR)
    a = ap.parse_args()

    F = np.load(a.forms, allow_pickle=True)
    Mb = np.stack([F[f'Mb_dec_{X}'] for X in SPECS])
    dC = {X: F[f'dC_{X}'] for X in SPECS}
    add, base = F[a.add_key], F[a.base_key]
    delta = sum(Mb[i] @ dC[X] for i, X in enumerate(SPECS))

    bm = np.loadtxt(os.path.join(a.ddir, 'binning_matrix_act.txt'))
    cents = bm @ np.arange(bm.shape[1])
    fbase = os.path.join(a.ddir, BASE_COV[a.variant])
    cov_pkg = np.loadtxt(fbase)
    checks = {
        'bin_centres_max_abs_diff': float(np.max(np.abs(F['cents'] - cents))),
        'base_vs_packaged_cov_max_rel_diff': float(np.max(np.abs(base / cov_pkg - 1))),
        'addition_symmetric_max_abs': float(np.max(np.abs(add - add.T))),
        'addition_min_eigenvalue_over_max': float(np.linalg.eigvalsh(add).min() / np.linalg.eigvalsh(add).max()),
    }
    if 'Delta' in F.files:
        checks['delta_vs_forms_Delta_max_rel_diff'] = float(np.max(np.abs(delta / F['Delta'] - 1)))
    for k, v in checks.items():
        print(f'  {k}: {v:.3e}')
    assert checks['bin_centres_max_abs_diff'] < 1e-8
    assert checks['base_vs_packaged_cov_max_rel_diff'] < 1e-8, 'forms base covariance is not the packaged one'
    assert Mb.shape == (3, 18, 3000)

    out = os.path.join(a.ddir, f'like_corrs_{a.variant}')
    os.makedirs(out, exist_ok=True)
    np.savez(os.path.join(out, 'response_binned.npz'), Mb=Mb, specs=np.array(SPECS), cents=cents,
             ells=np.arange(Mb.shape[2]))
    np.savetxt(os.path.join(out, 'covmat_cmbmarg_analytic.txt'), base + add,
               header=f'{a.variant}: base covariance ({BASE_COV[a.variant]}) + Eq. 34 CMB addition; raw')
    np.savetxt(os.path.join(out, 'lens_only_delta.txt'), delta,
               header=f'{a.variant}: Eq. 35 recentering Delta_b (C_L^kk units), 18 bins')
    prov = {
        'variant': a.variant, 'written': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
        'forms': a.forms, 'forms_md5': md5(a.forms), 'add_key': a.add_key, 'base_key': a.base_key,
        'base_cov_file': BASE_COV[a.variant], 'norm_matrix': a.norm,
        'norm_matrix_files': '/project/rrg-rbond-ac/jiaqu/dr6plus_lensing_variants/like_corrs/<config>/'
                             'norm_correction_matrix_bhph_Lmin0_Lmax4000.npy (see CANONICAL.md there)',
        'n1_derivative': 'plain MV N1der_{TT,EE,BB,TE} (decision 2026-09-25)',
        'form': 'N1-deconvolved: response B A R, no dN1/dC_kk term in the theory',
        'cmb_data': a.cmb_data, 'checks': checks, 'note': a.note,
        'peak_sigma_broadening_pct_40_1100': float(np.max(100 * (np.sqrt(np.diag(base + add) / np.diag(base)) - 1)[2:14])),
    }
    with open(os.path.join(out, 'PROVENANCE.json'), 'w') as fh:
        json.dump(prov, fh, indent=1)
    print('wrote', out)


if __name__ == '__main__':
    main()
