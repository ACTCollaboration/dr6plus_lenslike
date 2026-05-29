import getdist.plots as gplt
from getdist import loadMCSamples
import matplotlib
matplotlib.use('Agg')

base = '/scratch/jiaqu/chains/fg_marg_test/'
runs = [
    ('afg_extended',     r'$A_{\rm fg}\sim U[-2,2]$'),
    ('afg_extended_u5',  r'$A_{\rm fg}\sim U[-5,5]$'),
    ('afg_extended_u10', r'$A_{\rm fg}\sim U[-10,10]$'),
]
samples = []
for root, lab in runs:
    s = loadMCSamples(base + root, settings={'ignore_rows': 0.3})
    s.setParamNames  # noop
    samples.append(s)
    st = s.getMargeStats().parWithName('S825')
    print(f'{lab:32s}  S8_lens = {st.mean:.4f} +/- {st.err:.4f}')

# label S825 nicely
for s in samples:
    s.paramNames.parWithName('S825').label = r'S_8^{\rm lens}\equiv\sigma_8(\Omega_m/0.3)^{0.25}'
    s.paramNames.parWithName('A_fg').label = r'A_{\rm fg}'

g = gplt.get_subplot_plotter(width_inch=7)
g.settings.legend_fontsize = 13
g.settings.axes_fontsize = 12
g.settings.axes_labelsize = 15
g.settings.alpha_filled_add = 0.6
g.triangle_plot(
    samples, ['S825', 'A_fg'],
    filled=True,
    legend_labels=[lab for _, lab in runs],
    legend_loc='upper right',
    contour_colors=['#1f77b4', '#ff7f0e', '#2ca02c'],
    param_limits={'A_fg': (-10, 10), 'S825': (0.74, 0.90)},
)
out = '/home/jiaqu/dr6plus_lenslike/notebooks/fg_marg_S825_Afg_triangle.png'
g.export(out)
print('wrote', out)
