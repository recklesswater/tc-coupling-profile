"""Figure for the B-factor check: figures/fig3_bfactor_validation.png.

Left panel: per-protein correlation of BD and of the GNM mean-square fluctuation
against experimental B-factors. Points above the diagonal are proteins where the
mean-field diagonal wins.

Right panel: the paired difference, so the direction of the result is readable
without eyeballing the scatter.
"""

from __future__ import annotations

import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 140
plt.rcParams["savefig.bbox"] = "tight"

BD_COLOR = "#2b6cb0"
MSF_COLOR = "#2f855a"


def load():
    path = os.path.join(ROOT, "results", "bfactor_validation.csv")
    ids, bd, msf = [], [], []
    with open(path, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            ids.append(row["pdb"])
            bd.append(float(row["r_bd_bfactor"]))
            msf.append(float(row["r_msf_bfactor"]))
    return ids, np.array(bd), np.array(msf)


def main() -> None:
    ids, bd, msf = load()
    diff = msf - bd

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))

    ax = axes[0]
    lim = (-0.35, 1.0)
    ax.plot(lim, lim, color="#a0aec0", lw=1.2, ls="--", zorder=1)
    ax.scatter(bd, msf, s=46, color=MSF_COLOR, alpha=0.85,
               edgecolor="white", linewidth=0.8, zorder=2)
    ax.axhline(0, color="k", lw=0.8, alpha=0.5)
    ax.axvline(0, color="k", lw=0.8, alpha=0.5)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("r(BD, B-factor), per protein")
    ax.set_ylabel("r(log MSF, B-factor), per protein")
    ax.set_title("above the diagonal: the mean-field diagonal wins", fontsize=11)
    ax.text(0.04, 0.94, "n = %d proteins" % len(ids), transform=ax.transAxes,
            va="top", fontsize=9, color="#4a5568")
    ax.grid(alpha=0.25)

    ax = axes[1]
    order = np.argsort(diff)
    ypos = np.arange(len(diff))
    colors = [MSF_COLOR if d > 0 else BD_COLOR for d in diff[order]]
    ax.barh(ypos, diff[order], color=colors, alpha=0.9, height=0.72)
    ax.axvline(0, color="k", lw=1)
    ax.axvline(np.median(diff), color="#e53e3e", lw=1.5, ls="--",
               label="median %+.3f" % np.median(diff))
    ax.set_yticks(ypos)
    ax.set_yticklabels([ids[i] for i in order], fontsize=6.5)
    ax.set_xlabel("r(log MSF, B) - r(BD, B)")
    ax.set_title("MSF ahead in %d / %d proteins" % (int((diff > 0).sum()), len(diff)),
                 fontsize=11)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.25, axis="x")

    fig.suptitle("BD against an independent axis: experimental B-factors", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    outdir = os.path.join(ROOT, "figures")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "fig3_bfactor_validation.png")
    fig.savefig(path)
    plt.close(fig)
    print("wrote %s" % path)


if __name__ == "__main__":
    main()
