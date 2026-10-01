#!/usr/bin/env python
"""Eq. 34 covariance addition and Eq. 35 constant of the lens-only CMB
marginalization for one reconstruction configuration, from its like_corrs
derivative products (DR6plus_lensing/like_corrs), in two estimator forms.

Response per unit C_ell (kappa units per muK^2), X = TT, EE, TE:
    R^X[L, ell] = -2 dA_L/dC^X_ell / A_L * C^kk,fid_L + dN1_L/dC^X_ell
Standard estimator (MCN1 subtracted; what dr6plus_lenslike implements today):
    M^X = B R^X                       cov_add = sum_XY M^X N^XY M^Y^T
    const = sum_X R^X dC^X - dN1/dClkk @ clkk_fid        (theory: clkk + dN1/dClkk @ clkk + const)
N1-deconvolved estimator (DR6+ products, requirements note Eq. Mdecon):
    M^X = B A R^X,   A = (I + M_kk)^-1 with the MCN1-matched kernel
    const = A sum_X R^X dC^X                              (theory: plain clkk + const, no dN1/dClkk term)
dC^X = ACT DR6 CMB-only bandpowers minus the window-binned fiducial, ell in
[600, 3000] (src/generate_lens_only_const.py construction); N = the ACT DR6
CMB-only bandpower covariance at ell >= 600 (the covmat_act_cmbmarg_analytic
construction; TT/TE/EE, nine blocks, doubly binned) or, with --pact, the P-ACT
lite covariance (plik_lite PlanckActCut blocks + ACT sacc; src/pact_cmbmarg_actplanck.py).

Writes into the configuration's test data directory (src/like_corrs_test_ddir.py):
    covmat_act_cmbmarg_analytic.txt(+_info.npy)   standard form, base --base-cov
    lens_only_const.npy                           standard form (load_data keys)
    <name>_cmbmarg_analytic_forms.npz             both forms: additions, constants,
                                                  binned responses, bases, sigma tables
Usage (module environment):
  python src/like_corrs_variant_analytic.py --name variant0 \
      --ddir /scratch/jiaqu/like_corrs/lenslike_ddir_variant0 [--sacc FILE] [--pact]
"""
import argparse
import datetime
import json
import os
import sys

import numpy as np

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
import pact_cmbmarg_actplanck as pm                      # noqa: E402
import generate_lens_only_const as glc                   # noqa: E402

D = os.path.join(REPO, "src", "dr6plus_lenslike", "data", "v1.0")
LC_ROOT = "/project/rrg-rbond-ac/jiaqu/dr6plus_lensing_variants/like_corrs"
CONFIGS = "/home/jiaqu/DR6plus_lensing/like_corrs/configs.yaml"
DAL = {"TT": 0, "EE": 1, "BB": 2, "TE": 3}
LMAX = 3000


def load_products(name, lc_root=LC_ROOT, norm="bhph", n1=True):
    """like_corrs products of one configuration. norm: 'plain' (plain-MV
    matrix) or 'bhph' (TT row of the profile-hardened MV of the pipeline,
    norm_derivative_cmb.py --harden); the fiducial A_L is A_MV in both.
    n1=False skips the N1der files (the patch configurations have none)."""
    import yaml
    cfg = yaml.safe_load(open(CONFIGS))["configs"][name]
    lc = os.path.join(lc_root, name)
    meta = json.load(open(os.path.join(lc, "norm_derivative_meta.json")))
    lmin, lmax, kmax = cfg["lmin"], cfg["lmax"], meta["k_ellmax"]
    tag = {"plain": "", "bhph": "bhph_"}[norm]
    dal = np.load(os.path.join(lc, f"norm_correction_matrix_{tag}Lmin0_Lmax{kmax}.npy"))
    fal = np.loadtxt(os.path.join(lc, f"n0mv_fiducial_lmin{lmin}_lmax{lmax}_Lmin0_Lmax{kmax}.txt"))[1]
    dn1 = n1chk = None
    if n1:
        dn1 = {X: np.load(os.path.join(lc, f"N1der_{X}_lmin{lmin}_lmax{lmax}_full.npy")) for X in ("TT", "EE", "BB", "TE")}
        n1chk = np.loadtxt(os.path.join(lc, "n1_derivative_check.txt"))
    return dict(cfg=cfg, dal=dal, fal=fal, dn1=dn1, lmin=lmin, lmax=lmax, n1_check=n1chk, dir=lc, norm=norm)


def response_perCl(dal, fal, dn1, clkk_fid, Lmax=LMAX, lmax_in=LMAX):
    """R^X[L, ell], L, ell < Lmax (norm + N1 terms), X = TT, EE, TE (BB carries none)."""
    R = {}
    good = fal[:Lmax] > 0
    good[:2] = False
    for X in ("TT", "EE", "TE"):
        m = np.zeros((Lmax, lmax_in))
        m[good, :] = (-2 * dal[DAL[X]][:Lmax, :lmax_in][good, :] / fal[:Lmax, None][good]) * clkk_fid[:Lmax, None][good]
        n = min(lmax_in, dn1[X].shape[1])
        m[:min(Lmax, dn1[X].shape[0]), :n] += dn1[X][:min(Lmax, dn1[X].shape[0]), :n]
        R[X] = m
    return R


def sacc_addition(Mb, cov, sinfo):
    """sum_XY Mb^X N^XY Mb^Y^T with the doubly-binned ACT sacc operator."""
    Md = {X: pm.dbin_sacc(Mb[X], sinfo[X]["supports"]) for X in ("TT", "TE", "EE")}
    n = Md["TT"].shape[0]
    add = np.zeros((n, n))
    for X in ("TT", "TE", "EE"):
        for Y in ("TT", "TE", "EE"):
            add += Md[X] @ cov[np.ix_(sinfo[X]["idx"], sinfo[Y]["idx"])] @ Md[Y].T
    return add


def plik_addition(Mb, cov_plik, keep, lo, hi):
    Mp = {X: pm.dbin_plik(Mb[X], lo[X], hi[X]) for X in ("TT", "TE", "EE")}
    n = Mp["TT"].shape[0]
    add = np.zeros((n, n))
    for X in ("TT", "TE", "EE"):
        for Y in ("TT", "TE", "EE"):
            add += Mp[X] @ cov_plik[np.ix_(keep[X], keep[Y])] @ Mp[Y].T
    return add


def sig_pct(base, add):
    return 100 * (np.sqrt(np.diag(base + add) / np.diag(base)) - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--ddir", required=True, help="test data directory (src/like_corrs_test_ddir.py)")
    ap.add_argument("--lc-root", default=LC_ROOT)
    ap.add_argument("--norm", choices=["plain", "bhph"], default="bhph",
                    help="normalization matrix: plain MV or the pipeline's profile-hardened MV")
    ap.add_argument("--operator", default=None,
                    help="npz with A (deconvolution operator) and cov_dec/cov_std "
                         "(default <ddir>/<name>_decon_operator.npz)")
    ap.add_argument("--sacc", default=pm.SACC_FILE, help="ACT CMB-only sacc file for N and dC")
    ap.add_argument("--pact", action="store_true", help="add the plik_lite PlanckActCut blocks (P-ACT lite N)")
    ap.add_argument("--base-cov", default=os.path.join(D, "covmat_act.txt"),
                    help="base lensing covariance written with the standard-form file (stand-in)")
    a = ap.parse_args()

    P = load_products(a.name, a.lc_root, norm=a.norm)
    lp = np.loadtxt(f"{D}/like_corrs/cosmo2017_10K_acc3_lenspotentialCls.dat")
    clkk_fid = np.zeros(int(lp[-1, 0]) + 1)          # indexed by L: the file starts at L = 2
    clkk_fid[lp[:, 0].astype(int)] = lp[:, 5] * 2 * np.pi / 4
    binmat = np.loadtxt(f"{D}/binning_matrix_act.txt")
    Lc = binmat @ np.arange(binmat.shape[1])
    op = np.load(a.operator or os.path.join(a.ddir, f"{a.name}_decon_operator.npz"), allow_pickle=True)
    A = op["A"]
    nL = A.shape[0]

    print(f"[{a.name}] per-ell response (norm + N1)...", flush=True)
    R = response_perCl(P["dal"], P["fal"], P["dn1"], clkk_fid)
    B = binmat[:, :LMAX]
    BA = B[:, :nL] @ A
    Mb_std = {X: B @ R[X] for X in R}                         # (18, 3000)
    Mb_dec = {X: BA @ R[X][:nL, :] for X in R}

    print("CMB covariance blocks...", flush=True)
    pm.SACC_FILE = a.sacc
    cov_sacc, sinfo = pm.sacc_blocks()
    add_std = sacc_addition(Mb_std, cov_sacc, sinfo)
    add_dec = sacc_addition(Mb_dec, cov_sacc, sinfo)
    tag = "ACT DR6 CMB-only sacc, ell >= 600"
    if a.pact:
        cov_plik, keep, lo, hi = pm.plik_bins()
        add_std += plik_addition(Mb_std, cov_plik, keep, lo, hi)
        add_dec += plik_addition(Mb_dec, cov_plik, keep, lo, hi)
        tag = "P-ACT lite (plik_lite PlanckActCut + " + tag + ")"

    # Eq. 35 constant: dC from the sacc data minus window-binned fiducial
    import sacc
    fid = np.loadtxt(f"{D}/like_corrs/cosmo2017_10K_acc3_lensedCls.dat")
    fid_ells = fid[:, 0].astype(int)
    fid_Dl = {"TT": fid[:, 1], "EE": fid[:, 2], "TE": fid[:, 4]}
    s = sacc.Sacc.load_fits(a.sacc)
    ell_range, dCl = glc.build_delta_Cl(s, fid_ells, fid_Dl)
    dC = {X: np.zeros(LMAX) for X in ("TT", "EE", "TE")}
    for X in dC:
        m = ell_range < LMAX
        dC[X][ell_range[m]] = dCl[X][m]
    norm_term = np.zeros(LMAX)
    n1_cmb_term = np.zeros(LMAX)
    good = P["fal"][:LMAX] > 0
    good[:2] = False
    for X in ("TT", "EE", "TE"):
        raw = P["dal"][DAL[X]][:LMAX, :LMAX] @ dC[X]
        norm_term[good] += -2 * raw[good] / P["fal"][:LMAX][good] * clkk_fid[:LMAX][good]
        n1_cmb_term += P["dn1"][X][:LMAX, :LMAX] @ dC[X]
    kern_raw = np.load(str(op["kernel_raw_path"]))
    n1_kk_term = -(kern_raw[:LMAX, :LMAX] @ clkk_fid[:LMAX])
    const_std = norm_term + n1_cmb_term + n1_kk_term
    const_dec = np.zeros(LMAX)
    const_dec[:nL] = A @ (norm_term + n1_cmb_term)[:nL]

    base_std = np.loadtxt(a.base_cov)
    base_dec = op["cov_dec"]
    base_own_std = op["cov_std"]
    t = binmat @ clkk_fid[:binmat.shape[1]]
    np.set_printoptions(precision=3, linewidth=170, suppress=True)
    print(f"\nN = {tag}")
    print("L centres          ", Lc.round(0))
    print("sigma % std form / covmat_act (DR6 night base)   ", sig_pct(base_std, add_std))
    print("sigma % std form / variant standard sim cov      ", sig_pct(base_own_std, add_std))
    print("sigma % deconvolved form / variant deconv sim cov", sig_pct(base_dec, add_dec))
    print("frac % sqrt(add_bb)/t_b: std", (100 * np.sqrt(np.diag(add_std)) / t), "\n                         dec", (100 * np.sqrt(np.diag(add_dec)) / t))
    print("ratio dec/std of the addition diagonal", np.diag(add_dec) / np.diag(add_std))
    bt_std = binmat @ (const_std + kern_raw[:LMAX, :LMAX] @ clkk_fid[:LMAX])   # theory shift at fiducial, std form
    bt_dec = binmat @ const_dec
    print("net recentering at fiducial (% of t_b): std", (100 * bt_std / t), "\n                                       dec", (100 * bt_dec / t))

    info = dict(cov_original=base_std, cov_add=add_std, cov_modified=base_std + add_std,
                diag_increase_pct=sig_pct(base_std, add_std), L_bin_centers=Lc, n_bins=len(Lc),
                date_generated=datetime.datetime.now().isoformat(),
                description=f"Eq. 34 modified covariance, standard-estimator form, like_corrs {a.name}; "
                            f"base {a.base_cov} (stand-in)", source_cmb_data=tag, sacc_file=a.sacc,
                like_corrs_dir=P["dir"], generator="src/like_corrs_variant_analytic.py")
    np.savetxt(os.path.join(a.ddir, "covmat_act_cmbmarg_analytic.txt"), base_std + add_std)
    np.save(os.path.join(a.ddir, "covmat_act_cmbmarg_analytic_info.npy"), info)
    out = dict(L=np.arange(LMAX), eq35_const=const_std, norm_term=norm_term, n1_cmb_term=n1_cmb_term,
               n1_kk_term=n1_kk_term, clkk_fid=clkk_fid[:LMAX],
               source_cmb_data=f"{os.path.basename(a.sacc)} (foreground-marginalized ACT DR6 CMB-only), ell >= 600",
               delta_cl_method="data bandpower minus window-binned fiducial, interpolated to ell in [600, 3000]",
               like_corrs=a.name, like_corrs_dir=P["dir"], kernel_raw_path=str(op["kernel_raw_path"]),
               date_generated=datetime.datetime.now().isoformat(), generator="src/like_corrs_variant_analytic.py")
    np.save(os.path.join(a.ddir, "lens_only_const.npy"), out)
    np.savez(os.path.join(a.ddir, f"{a.name}_cmbmarg_analytic_forms.npz"),
             add_std=add_std, add_dec=add_dec, const_std=const_std, const_dec=const_dec,
             norm_term=norm_term, n1_cmb_term=n1_cmb_term, n1_kk_term=n1_kk_term,
             base_std=base_std, base_own_std=base_own_std, base_dec=base_dec, cents=Lc, t_fid=t,
             **{f"Mb_std_{X}": Mb_std[X] for X in Mb_std}, **{f"Mb_dec_{X}": Mb_dec[X] for X in Mb_dec},
             **{f"dC_{X}": dC[X] for X in dC}, cmb_cov=tag, sacc=a.sacc, pact=a.pact,
             note="standard = B R (MCN1-subtracted estimator, dr6plus_lenslike today); "
                  "dec = B A R, A=(I+M)^-1 MCN1-matched kernel (N1-deconvolved DR6+ product)")
    print(f"\nwrote covmat_act_cmbmarg_analytic.txt, lens_only_const.npy, {a.name}_cmbmarg_analytic_forms.npz in {a.ddir}")


if __name__ == "__main__":
    main()
