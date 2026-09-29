"""sigma8 vs A_fg triangle for the lensing + DESI DR2 BAO chains, across widened
A_fg priors. BAO breaks the sigma8-Omega_m degeneracy so sigma8 (not just
S8_lens) is well constrained. Degradation in sigma(sigma8) relative to the
no-fg-marg control is shown in the legend."""
import getdist.plots as gplt
from getdist import loadMCSamples
import matplotlib
matplotlib.use('Agg')

base = '/scratch/jiaqu/chains/fg_marg_test/'
# (root, legend label, contour color). Control has no A_fg (fg_marg=False).
# An A_fg run is listed first so getdist builds the param array from it; the
# control is skipped in the A_fg panels but shown in sigma8.
runs = [
    ('afg_extended_bao',     r'$A_{\rm fg}\sim U[-2,2]$',        '#1f77b4'),
    ('afg_extended_u5_bao',  r'$A_{\rm fg}\sim U[-5,5]$',         '#ff7f0e'),
    ('afg_extended_u10_bao', r'$A_{\rm fg}\sim U[-10,10]$',       '#2ca02c'),
    ('control_extended_bao', r'no fg marg ($A_{\rm fg}\!=\!0$)', '#000000'),
]
samples, errs = [], []
for root, lab, _ in runs:
    s = loadMCSamples(base + root, settings={'ignore_rows': 0.3})
    s.paramNames.parWithName('sigma8').label = r'\sigma_8'
    afg = s.paramNames.parWithName('A_fg')
    if afg is not None:
        afg.label = r'A_{\rm fg}'
    samples.append(s)
    errs.append(s.getMargeStats().parWithName('sigma8').err)

# Degradation in sigma(sigma8) relative to the no-fg-marg control (last entry).
sig_ref = errs[-1]
legend_labels = []
for (root, lab, _), err in zip(runs, errs):
    if root == 'control_extended_bao':
        legend_labels.append(lab + ' [ref]')
    else:
        pct = int(round(100*(err/sig_ref - 1)))  # int(round) collapses -0
        legend_labels.append(lab + f' (${pct:+d}\\%$)')
    print(f'{lab:34s}  sigma(sigma8) = {err:.4f}  {legend_labels[-1]}')

g = gplt.get_subplot_plotter(width_inch=9)
g.settings.legend_fontsize = 12
g.settings.axes_fontsize = 12
g.settings.axes_labelsize = 15
g.settings.alpha_filled_add = 0.55
g.triangle_plot(
    samples, ['sigma8', 'A_fg'],
    filled=True,
    legend_labels=legend_labels,
    legend_loc='upper right',
    contour_colors=[c for _, _, c in runs],
    param_limits={'A_fg': (-10, 10)},
)
out = '/home/jiaqu/dr6plus_lenslike/notebooks/fg_marg_sigma8_Afg_triangle_bao'
g.export(out + '.png', dpi=300)
g.export(out + '.pdf')
print('wrote', out + '.png and .pdf')
