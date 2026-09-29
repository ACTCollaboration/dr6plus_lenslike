"""Plot the C_L^kk bandpowers with error bars and show how adding A_fg*template
shifts the theory curve for A_fg = 1, 2, 10 (extended fiducial variant)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dr6plus_lenslike.dr6plus_lenslike import load_data

# Same settings as runs/fg_marg_test/mcmc_afg_extended*.yaml
d = load_data(
    variant='dr6plus_fiducial_extended',
    lens_only=True,
    trim_lmax=2998,
    apply_hartlap=True,
    like_corrections=False,
    cov_file='/home/jiaqu/dr6plus_lenslike/products/mock_covmat.npy',
    fg_marg=True,
)

L = d['bcents_act']
data = d['data_binned_clkk']
err = np.sqrt(np.diag(d['cov']))[:data.size]
templ = d['fg_template_bandpower']

# For the dr6plus_fiducial variant the data vector is the noiseless binned
# fiducial theory, so the theory bandpowers equal the data bandpowers.
theory = data

s = 1e7  # scale for readability
amps = [(1, '#1f77b4'), (2, '#ff7f0e'), (10, '#2ca02c')]

fig, (ax, axr) = plt.subplots(
    2, 1, figsize=(9, 8), sharex=True,
    gridspec_kw={'height_ratios': [2.2, 1], 'hspace': 0.08})

# Top: bandpowers with error bars + theory and theory+A*template curves.
ax.errorbar(L, s*data, yerr=s*err, fmt='o', ms=5, color='k',
            capsize=2, label='fiducial data', zorder=5)
ax.plot(L, s*theory, '-', color='k', lw=1.5, label=r'theory $C_L^{\kappa\kappa}$')
for A, c in amps:
    ax.plot(L, s*(theory + A*templ), '--', color=c, lw=1.8,
            label=rf'theory $+\,{A}\times$ template')
ax.set_ylabel(r'$10^{7}\,C_L^{\kappa\kappa}$', fontsize=15)
ax.set_title('DR6+ extended lensing bandpowers vs. foreground template', fontsize=13)
ax.legend(fontsize=11, frameon=True)
ax.grid(alpha=0.3)

# Bottom: template shift in units of the per-bin error bar.
axr.axhspan(-1, 1, color='grey', alpha=0.15, label=r'$\pm1\sigma$ band')
axr.axhline(0, color='k', lw=0.8)
for A, c in amps:
    axr.plot(L, A*templ/err, '--o', color=c, ms=4, lw=1.6,
             label=rf'${A}\times$ template')
axr.set_xlabel(r'$L$', fontsize=15)
axr.set_ylabel(r'$A_{\rm fg}\,T(L)\,/\,\sigma_L$', fontsize=14)
axr.legend(fontsize=10, frameon=True, ncol=2, loc='lower left')
axr.grid(alpha=0.3)
fig.tight_layout()
out = '/home/jiaqu/dr6plus_lenslike/notebooks/fg_template_clkk_bandpowers'
fig.savefig(out + '.png', dpi=300)
fig.savefig(out + '.pdf')
print('wrote', out + '.png', 'and', out + '.pdf')

# Diagnostics: how big is the template shift relative to the error bars?
print('\nbin   L      data*1e7   err*1e7   1x/err   2x/err  10x/err')
for i in range(data.size):
    print(f'{i:>3} {L[i]:6.0f} {s*data[i]:9.3f} {s*err[i]:8.3f}'
          f' {abs(templ[i]/err[i]):7.3f} {abs(2*templ[i]/err[i]):7.3f}'
          f' {abs(10*templ[i]/err[i]):7.3f}')
