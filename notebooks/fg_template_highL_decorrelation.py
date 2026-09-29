"""Test the hypothesis that high-L bins pin down the foreground template and
decorrelate it from the lensing-amplitude (S8) direction.

For bins up to each L_max we compute, in the data metric C^-1:
  - sigma(A_fg | cosmology fixed) = 1/sqrt(t.Cinv.t)   -> how well A_fg is measured
  - |rho| between the template and the amplitude (~S8) direction (proxy g ~ Clkk)
A falling |rho| with L_max means high-L information rotates the template away
from the S8 direction; a falling sigma(A_fg) means high-L pins down its amplitude.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dr6plus_lenslike.dr6plus_lenslike import load_data

d = load_data(variant='dr6plus_fiducial_extended', lens_only=True, trim_lmax=2998,
              apply_hartlap=True, like_corrections=False,
              cov_file='/home/jiaqu/dr6plus_lenslike/products/mock_covmat.npy',
              fg_marg=True)
L  = d['bcents_act']
cl = d['data_binned_clkk']
t  = d['fg_template_bandpower']
cov = d['cov'][:cl.size, :cl.size]

Lmax, sigA, absrho = [], [], []
for k in range(2, cl.size + 1):
    sl = slice(0, k)
    Ck = np.linalg.inv(cov[sl, sl])
    tt = t[sl] @ Ck @ t[sl]
    ta = t[sl] @ Ck @ cl[sl]
    aa = cl[sl] @ Ck @ cl[sl]
    Lmax.append(L[k-1]); sigA.append(1/np.sqrt(tt))
    absrho.append(abs(ta/np.sqrt(tt*aa)))
Lmax = np.array(Lmax)

fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                             gridspec_kw={'hspace': 0.08})
a1.plot(Lmax, sigA, 'o-', color='#1f77b4', lw=2)
a1.set_ylabel(r'$\sigma(A_{\rm fg}\,|\,{\rm cosmo\ fixed})$', fontsize=13)
a1.set_title('High-$L$ bins pin down the template and rotate it off the $S_8$ direction',
             fontsize=12)
a1.set_yscale('log'); a1.grid(alpha=0.3)
a1.annotate('low-$L$ only:\ntemplate barely measured', (Lmax[1], sigA[1]),
            textcoords='offset points', xytext=(35, 5), fontsize=10,
            arrowprops=dict(arrowstyle='->', color='grey'))

a2.plot(Lmax, absrho, 'o-', color='#d62728', lw=2)
a2.axhline(absrho[-1], ls=':', color='grey',
           label=rf'full survey $|\rho|={absrho[-1]:.2f}$')
a2.set_ylabel(r'$|\rho(\,{\rm template},\ S_8\ {\rm dir}\,)|$', fontsize=13)
a2.set_xlabel(r'$L_{\rm max}$ included', fontsize=14)
a2.set_ylim(0, 1); a2.grid(alpha=0.3); a2.legend(fontsize=11)
a2.annotate('low-$L$ template $\\approx$ pure amplitude shift\n(degenerate with $S_8$)',
            (Lmax[1], absrho[1]), textcoords='offset points', xytext=(30, -35),
            fontsize=10, arrowprops=dict(arrowstyle='->', color='grey'))

fig.tight_layout()
out = '/home/jiaqu/dr6plus_lenslike/notebooks/fg_template_highL_decorrelation'
fig.savefig(out + '.png', dpi=300); fig.savefig(out + '.pdf')
print('wrote', out + '.png and .pdf')
print(f'low-L (L<={Lmax[1]:.0f}):  |rho|={absrho[1]:.2f}, sigma(A_fg|cosmo)={sigA[1]:.1f}')
print(f'full  (L<={Lmax[-1]:.0f}): |rho|={absrho[-1]:.2f}, sigma(A_fg|cosmo)={sigA[-1]:.1f}')
