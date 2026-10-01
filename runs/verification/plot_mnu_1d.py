"""1D posterior of Sum m_nu from the free-neutrino-mass verification chains.

Lensing + P-ACT lite + DESI DR2 BAO on the cosmo2017 mocks (input 0.06 eV),
ACTbase (crimson, filled) and ACT+Planck (navy line). Style and chain loading
follow plot_triangles.py. Prints mean, standard deviation and the 95 per cent
upper limit of each chain.

Output: runs/verification/plots/verification_mnu_1d.{pdf,png}

Usage: python runs/verification/plot_mnu_1d.py   (cluster env loaded)
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from plot_triangles import CRIMSON, CRIMSON_EDGE, NAVY, LW, MARK, INPUT, OUT, load

ROOTS = {"ACTbase": "actbase_lenscmbbao_mnu_camb", "ACT+Planck": "act_planck_lenscmbbao_mnu_camb"}


def main():
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    handles, labels = [], []
    for (name, root), color in zip(ROOTS.items(), (CRIMSON, NAVY)):
        s = load(root)
        d = s.get1DDensity("mnu")
        x, p = d.x, d.P / d.P.max()
        m = s.mean("mnu")
        sd = s.std("mnu")
        ul = s.confidence("mnu", 0.05, upper=True)
        if name == "ACTbase":
            ax.fill_between(x, 0, p, color=CRIMSON, alpha=0.85, lw=0, zorder=2)
            ax.plot(x, p, color=CRIMSON_EDGE, lw=LW, zorder=3)
            handles.append(Patch(facecolor=CRIMSON, edgecolor=CRIMSON_EDGE, lw=LW))
        else:
            ax.plot(x, p, color=NAVY, lw=LW, zorder=4)
            handles.append(Line2D([], [], color=NAVY, lw=LW))
        labels.append(rf"{name}: $\Sigma m_\nu < {ul:.3f}$ eV (95%)")
        print(f"{name:11s} {root}: mean {m:.4f}  std {sd:.4f}  95% UL {ul:.4f} eV  (input {INPUT['mnu']})")
    ax.axvline(INPUT["mnu"], zorder=5, **MARK)
    handles.append(Line2D([], [], **MARK))
    labels.append("input 0.06 eV")
    ax.set_xlabel(r"$\Sigma m_\nu\,[{\rm eV}]$", fontsize=15)
    ax.set_ylabel(r"$P/P_{\max}$", fontsize=15)
    ax.set_xlim(0, 0.3)
    ax.set_ylim(0, 1.08)
    ax.tick_params(labelsize=11)
    ax.set_title("Lensing + P-ACT + DESI DR2 BAO (mock)", fontsize=11)
    ax.legend(handles, labels, frameon=False, fontsize=10, loc="upper right",
              handlelength=2.6, alignment="left")
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"verification_mnu_1d.{ext}"), dpi=200 if ext == "png" else None)
    print(os.path.join(OUT, "verification_mnu_1d"))


if __name__ == "__main__":
    main()
