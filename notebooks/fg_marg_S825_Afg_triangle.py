import getdist.plots as gplt
from getdist import loadMCSamples
import matplotlib
matplotlib.use('Agg')

base = '/scratch/jiaqu/chains/fg_marg_test/'
# (root, legend label, contour color). Control has no A_fg (fg_marg=False).
# An A_fg run is listed first so getdist builds the param array from it; the
# control is skipped in the A_fg panels (no such parameter) but shown in S825.
runs = [
    ('afg_extended',     r'$A_{\rm fg}\sim U[-2,2]$',        '#1f77b4'),
    ('afg_extended_u5',  r'$A_{\rm fg}\sim U[-5,5]$',         '#ff7f0e'),
    ('afg_extended_u10', r'$A_{\rm fg}\sim U[-10,10]$',       '#2ca02c'),
    ('control_extended', r'no fg marg ($A_{\rm fg}\!=\!0$)', '#000000'),
]
samples = []
for root, lab, _ in runs:
    s = loadMCSamples(base + root, settings={'ignore_rows': 0.3})
    s.paramNames.parWithName('S825').label = r'S_8^{\rm lens}\equiv\sigma_8(\Omega_m/0.3)^{0.25}'
    afg = s.paramNames.parWithName('A_fg')
    if afg is not None:
        afg.label = r'A_{\rm fg}'
    samples.append(s)
    st = s.getMargeStats().parWithName('S825')
    print(f'{lab:34s}  S8_lens = {st.mean:.4f} +/- {st.err:.4f}')

g = gplt.get_subplot_plotter(width_inch=7)
g.settings.legend_fontsize = 12
g.settings.axes_fontsize = 12
g.settings.axes_labelsize = 15
g.settings.alpha_filled_add = 0.55
g.triangle_plot(
    samples, ['S825', 'A_fg'],
    filled=True,
    legend_labels=[lab for _, lab, _ in runs],
    legend_loc='upper right',
    contour_colors=[c for _, _, c in runs],
    param_limits={'A_fg': (-10, 10), 'S825': (0.74, 0.90)},
)
out = '/home/jiaqu/dr6plus_lenslike/notebooks/fg_marg_S825_Afg_triangle.png'
g.export(out)
print('wrote', out)
