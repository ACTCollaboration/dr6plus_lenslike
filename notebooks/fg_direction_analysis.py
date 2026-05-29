import numpy as np
from dr6plus_lenslike.dr6plus_lenslike import load_data
from getdist import loadMCSamples

d = load_data(variant='dr6plus_fiducial_extended', lens_only=True, trim_lmax=2998,
              apply_hartlap=True, like_corrections=False,
              cov_file='/home/jiaqu/dr6plus_lenslike/products/mock_covmat.npy',
              fg_marg=True)

L   = d['bcents_act']
cl  = d['data_binned_clkk']          # fiducial bandpowers = amplitude (S8) direction proxy
t   = d['fg_template_bandpower']     # foreground template, bandpower space
cov = d['cov'][:cl.size, :cl.size]
Cinv = np.linalg.inv(cov)

def dot(x, y): return x @ Cinv @ y

# ---- 1. Alignment of template with the S8 / amplitude direction (g ~ Clkk) ----
rho = dot(t, cl) / np.sqrt(dot(t, t) * dot(cl, cl))
print(f"[1] Fisher correlation rho(template, amplitude/S8 dir) = {rho:.3f}")
print(f"    angle between them                                = {np.degrees(np.arccos(rho)):.1f} deg")
print(f"    predicted S8 degradation 1/sqrt(1-rho^2)          = {1/np.sqrt(1-rho**2):.3f}")
print(f"    => marginalizing A_fg should inflate sigma(S8) by ~{100*(1/np.sqrt(1-rho**2)-1):.0f}%\n")

# ---- 2. Posterior correlation measured directly from the chains ----
base = '/scratch/jiaqu/chains/fg_marg_test/'
print("[2] Measured posterior corr rho(S825, A_fg) and sigma(S8) from chains:")
sig_ctrl = loadMCSamples(base+'control_extended', settings={'ignore_rows':0.3}
            ).getMargeStats().parWithName('S825').err
for root in ['afg_extended', 'afg_extended_u5', 'afg_extended_u10']:
    s = loadMCSamples(base+root, settings={'ignore_rows':0.3})
    c = s.cov(['S825', 'A_fg'])
    r = c[0,1]/np.sqrt(c[0,0]*c[1,1])
    sigS8 = np.sqrt(c[0,0])
    print(f"    {root:18s} rho={r:+.3f}  sigma(S8)={sigS8:.4f}  "
          f"infl vs ctrl={sigS8/sig_ctrl:.3f}")
print(f"    (control no-marg sigma(S8) = {sig_ctrl:.4f})\n")

# ---- 3. Collaborator's hypothesis: does high L pin down the template? ----
# Per-bin Fisher information on A_fg (diagonal) and cumulative with full Cinv.
print("[3] Template information vs L (does high L pin down A_fg?):")
# self-SNR of template using full covariance, as a function of Lmax kept:
print("    Lmax_bin   sigma(A_fg|cosmo fixed)   rho(t,amp)   pred.degrad")
for k in range(2, cl.size+1):
    sl = slice(0, k)
    Ck = np.linalg.inv(cov[sl, sl])
    tt = t[sl] @ Ck @ t[sl]
    ta = t[sl] @ Ck @ cl[sl]
    aa = cl[sl] @ Ck @ cl[sl]
    rk = ta/np.sqrt(tt*aa)
    print(f"    L<= {L[k-1]:5.0f}     {1/np.sqrt(tt):8.3f}              "
          f"{rk:+.3f}      {1/np.sqrt(1-rk**2):.3f}")
# total self-SNR
tot = np.sqrt(dot(t, t))
print(f"\n    total template self-SNR sqrt(t.Cinv.t) = {tot:.3f}  "
      f"=> sigma(A_fg|cosmo fixed) = {1/tot:.3f}")
