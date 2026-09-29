#!/usr/bin/env python
"""CMB-2pt marginalization for JOINT ACT+Planck lensing bandpowers,
using the P-ACT lite CMB covariance (Eq. 34 at the bandpower level).

Both lensing normalizations respond to the same primary CMB, so for the
joint vector [ACT lensing bins ; Planck lensing bins]

    cov_add = sum_D sum_XY  Mtilde_D^X  cov_D^XY  (Mtilde_D^Y)^T ,

with D running over the two independent P-ACT lite CMB blocks
(pact_lite convention, base_config/pact_lite.yaml):
  - plik_lite bins cut per PlanckActCut: TT ell < 1000, TE/EE ell < 600
    (covariance in C_ell^2 units, block order TT/TE/EE);
  - ACT DR6 CMB-only SACC bins at ell >= 600 (covariance in D_ell^2 units,
    all nine TT/TE/EE blocks).
Mtilde stacks the ACT-lensing rows (per-ell M = -2 dAL_dC/fAL * clkk_fid
+ dN1, DR6 night MV matrices, reconstruction 600 < ell < 3000, binned with
binning_matrix_act) on the Planck-lensing rows (P18 per-ell matrices,
reconstruction 100 <= ell <= 2048, binned with binning_matrix_planck).
The ACT x Planck lensing cross block arises because both respond to the
same CMB bins.

Cases:
  A. actplanck: ACT DR6 night MV (18 bins) + Planck (9 bins); compared
     against the upstream chain-based covmat_actplanck_cmbmarg.txt.
  B. hilc + Planck: DR6+ optimal (hilcTP nightday glsloo, kept band
     [2:-4], 12 bins) + Planck (9 bins). CAVEAT: reuses the DR6 night MV
     response matrices for the GLS combination; base covariance assumes no
     ACT x Planck lensing cross term (none is available for this pair).

Outputs: printed sigma-increase profiles per block and added cross-block
correlations; arrays saved to products/pact_cmbmarg_joint.npz.
Run with the cluster module environment (see docs/CODE_README.md).
"""
import os
import numpy as np
from scipy.io import FortranFile

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
D = os.path.join(REPO, "dr6plus_lenslike", "data", "v1.0")
LC = os.path.join(D, "like_corrs")
PL = "/scratch/jiaqu/likelihood_data/data/planck_2018_pliklite_native"
# public ACT DR6 CMB-only release v1.0 (what act_dr6_cmbonly installs; the published chains condition on
# its covariance). Until 2026-09-25 the DR6-ACT-lite repository file act_dr6_cmb_sacc.fits, a different
# data version.
SACC_FILE = "/scratch/jiaqu/likelihood_data/data/ACTDR6CMBonly/v1.0/dr6_data_cmbonly.fits"
SPT12 = "/home/jiaqu/spt_act_likelihood/act_dr6_spt_lenslike/data/v1.2"
OUT = os.path.join(REPO, "products", "pact_cmbmarg_joint.npz")

# PlanckActCut lmax per spectrum. Planck TT and ACT TT overlap over
# 600 < ell < 1000; summing the response against both blocks double-counts it
# there (upper bracket). PLIK_TT_LMAX=600 removes the plik overlap bins
# instead (lower bracket, ACT-only in the overlap).
PACT_CUTS = {"TT": int(os.environ.get("PLIK_TT_LMAX", 1000)), "TE": 600, "EE": 600}
DAL_IDX = {"TT": 0, "EE": 1, "TE": 3}
SACC_TYPE = {"TT": "cl_00", "TE": "cl_0e", "EE": "cl_ee"}


def build_M_perCl(dal_file, fal, dn1_files, clkk_fid, Lmax_out, lmax_in=3000):
    """Per-C_ell response M^X[L, ell] for X in TT, EE, TE."""
    dal = np.load(dal_file)
    M = {}
    for X in ("TT", "EE", "TE"):
        m = np.zeros((Lmax_out, lmax_in))
        good = fal[:Lmax_out] > 0
        good[:2] = False
        m[good, :] = (-2 * dal[DAL_IDX[X]][:Lmax_out, :lmax_in][good, :]
                      / fal[:Lmax_out, None][good]) * clkk_fid[:Lmax_out, None][good]
        print(f"    loading {os.path.basename(dn1_files[X])}")
        dn1 = np.loadtxt(dn1_files[X])
        n = min(lmax_in, dn1.shape[1])
        m[:min(Lmax_out, dn1.shape[0]), :n] += dn1[:min(Lmax_out, dn1.shape[0]), :n]
        M[X] = m
    del dal
    return M


def plik_bins():
    blmin = np.loadtxt(f"{PL}/blmin.dat").astype(int) + 30
    blmax = np.loadtxt(f"{PL}/blmax.dat").astype(int) + 30
    f = FortranFile(f"{PL}/c_matrix_plik_v22.dat", "r")
    cov = f.read_reals(dtype=float).reshape((613, 613))
    cov = np.tril(cov) + np.tril(cov, -1).T
    NTT, NTE = 215, 199
    idx = {"TT": np.arange(0, NTT), "TE": np.arange(NTT, NTT + NTE),
           "EE": np.arange(NTT + NTE, 613)}
    keep = {X: idx[X][blmax[:len(idx[X])] <= PACT_CUTS[X]] for X in idx}
    lo = {X: blmin[keep[X] - idx[X][0]] for X in idx}
    hi = {X: blmax[keep[X] - idx[X][0]] for X in idx}
    print("  plik_lite kept bins:", {X: len(keep[X]) for X in keep})
    return cov, keep, lo, hi


def sacc_blocks():
    import sacc
    s = sacc.Sacc.load_fits(SACC_FILE)
    full = s.covariance.dense
    info = {}
    for X, t in SACC_TYPE.items():
        idx = np.array(s.indices(t))
        ells = np.array([s.data[i].get_tag("ell") for i in idx])
        srt = np.argsort(ells)
        idx, ells = idx[srt], ells[srt]
        m = ells >= 600
        idx = idx[m]
        win = s.get_bandpower_windows(list(idx))
        W = win.weight.T
        wl = win.values.astype(int)
        supports = [wl[W[b] > 1e-10] for b in range(W.shape[0])]
        info[X] = {"idx": idx, "supports": supports}
    print("  ACT sacc kept bins:", {X: len(info[X]["idx"]) for X in info})
    return full, info


def dbin_plik(Mb, lo, hi):
    """Sum per-C_ell response over each kept plik bin."""
    out = np.zeros((Mb.shape[0], len(lo)))
    for B, (a, b) in enumerate(zip(lo, hi)):
        b_eff = min(b, Mb.shape[1] - 1)
        if a <= b_eff:
            out[:, B] = Mb[:, a:b_eff + 1].sum(axis=1)
    return out


def dbin_sacc(Mb, supports):
    """Convert per-C_ell to per-D_ell and sum over each window support."""
    ell = np.arange(Mb.shape[1])
    fac = np.zeros_like(ell, dtype=float)
    fac[1:] = 2 * np.pi / (ell[1:] * (ell[1:] + 1))
    Md = Mb * fac[None, :]
    out = np.zeros((Mb.shape[0], len(supports)))
    for B, sup in enumerate(supports):
        sup = sup[sup < Mb.shape[1]]
        out[:, B] = Md[:, sup].sum(axis=1)
    return out


def main():
    lp = np.loadtxt(f"{LC}/cosmo2017_10K_acc3_lenspotentialCls.dat")
    clkk_fid = np.zeros(int(lp[-1, 0]) + 1)   # indexed by L: the file starts at L = 2
    clkk_fid[lp[:, 0].astype(int)] = lp[:, 5] * 2 * np.pi / 4

    print("Building ACT (DR6 night MV) lensing response...")
    fal_a = np.loadtxt(f"{LC}/n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax4000.txt")[1]
    MA = build_M_perCl(f"{LC}/norm_correction_matrix_Lmin0_Lmax4000.npy", fal_a,
                       {X: f"{LC}/N1der_{X}_lmin600_lmax3000_full.txt"
                        for X in ("TT", "EE", "TE")}, clkk_fid, 3000)
    binmat_a = np.loadtxt(f"{D}/binning_matrix_act.txt")           # (18, 3000)
    MA_b = {X: binmat_a @ MA[X][:binmat_a.shape[1], :] for X in MA}
    del MA

    print("Building Planck lensing response (P18 matrices)...")
    fal_p = np.loadtxt(f"{LC}/PLANCK_n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax3000.txt")[1]
    MP = build_M_perCl(f"{LC}/P18_norm_correction_matrix_Lmin0_Lmax3000.npy", fal_p,
                       {X: f"{LC}/N1_planck_der_{X}_lmin100_lmax2048.txt"
                        for X in ("TT", "EE", "TE")}, clkk_fid, 3000)
    binmat_p = np.loadtxt(f"{D}/binning_matrix_planck.txt")
    nlp = min(binmat_p.shape[1], 3000)
    MP_b = {X: binmat_p[:, :nlp] @ MP[X][:nlp, :] for X in MP}
    del MP

    print("Loading P-ACT lite CMB covariances...")
    cov_plik, keep, lo, hi = plik_bins()
    cov_sacc, sinfo = sacc_blocks()

    def joint_add(rows_act):
        """cov_add for [rows_act of ACT lensing ; all 9 Planck lensing bins]."""
        Mj_plik, Mj_sacc = {}, {}
        for X in ("TT", "TE", "EE"):
            a_pl = dbin_plik(MA_b[X][rows_act, :], lo[X], hi[X])
            p_pl = dbin_plik(MP_b[X], lo[X], hi[X])
            Mj_plik[X] = np.vstack([a_pl, p_pl])
            a_sc = dbin_sacc(MA_b[X][rows_act, :], sinfo[X]["supports"])
            p_sc = dbin_sacc(MP_b[X], sinfo[X]["supports"])
            Mj_sacc[X] = np.vstack([a_sc, p_sc])
        n = Mj_plik["TT"].shape[0]
        add = np.zeros((n, n))
        for X in ("TT", "TE", "EE"):
            for Y in ("TT", "TE", "EE"):
                add += Mj_plik[X] @ cov_plik[np.ix_(keep[X], keep[Y])] @ Mj_plik[Y].T
                add += (Mj_sacc[X] @ cov_sacc[np.ix_(sinfo[X]["idx"], sinfo[Y]["idx"])]
                        @ Mj_sacc[Y].T)
        return add

    def report(tag, base, add, n_act):
        sig = 100 * (np.sqrt(np.diag(base + add) / np.diag(base)) - 1)
        print(f"\n=== {tag} ===")
        print(f"  ACT-side block sigma increase (%):",
              np.array2string(sig[:n_act], precision=2))
        print(f"  Planck block sigma increase (%): ",
              np.array2string(sig[n_act:], precision=2))
        s = np.sqrt(np.diag(base + add))
        xcorr = add[:n_act, n_act:] / np.outer(s[:n_act], s[n_act:])
        print(f"  added ACTxPlanck cross correlation: max {np.max(np.abs(xcorr)):.3f}"
              f" (mean |.| {np.mean(np.abs(xcorr)):.3f})")
        return sig, xcorr

    # --- Case A: actplanck (18 + 9), vs upstream chain-based marged ---
    base_ap = np.loadtxt(f"{D}/covmat_actplanck.txt")
    add_ap = joint_add(np.arange(18))
    sigA, xA = report("A: ACT DR6 night MV + Planck (P-ACT lite analytic)",
                      base_ap, add_ap, 18)
    apm = np.loadtxt(f"{SPT12}/covmat_actplanck_cmbmarg.txt")
    sig_up = 100 * (np.sqrt(np.diag(apm) / np.diag(base_ap)) - 1)
    print("  upstream chain-based, ACT block:   ",
          np.array2string(sig_up[:18], precision=2))
    print("  upstream chain-based, Planck block:",
          np.array2string(sig_up[18:], precision=2))
    s0 = np.sqrt(np.diag(base_ap))
    x_up = (apm - base_ap)[:18, 18:] / np.outer(np.sqrt(np.diag(apm))[:18],
                                                np.sqrt(np.diag(apm))[18:])
    print(f"  upstream added cross correlation: max {np.max(np.abs(x_up)):.3f}"
          f" (mean |.| {np.mean(np.abs(x_up)):.3f})")

    # --- Case B: DR6+ hilc glsloo kept band (12) + Planck (9) ---
    g = np.loadtxt(f"{D}/covmat_clkk_hilcTP_nightday_glsloo.txt")[2:-4, 2:-4]
    base_hp = np.zeros((21, 21))
    base_hp[:12, :12] = g
    base_hp[12:, 12:] = base_ap[18:, 18:]
    add_hp = joint_add(np.arange(2, 14))
    sigB, xB = report("B: DR6+ hilc (glsloo, 12 bins 40<L<1100) + Planck "
                      "(DR6-night M approximation)", base_hp, add_hp, 12)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, add_actplanck=add_ap, base_actplanck=base_ap,
             add_hilcplanck=add_hp, base_hilcplanck=base_hp,
             note="P-ACT lite analytic Eq.34 additions; see src/pact_cmbmarg_actplanck.py")
    print(f"\nsaved {OUT}")


if __name__ == "__main__":
    main()
