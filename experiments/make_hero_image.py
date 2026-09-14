"""Build the repository hero image: 1UBQ coloured by TC, next to the pooled TC-flexibility scatter.

Writes figures/fig0_hero.png
"""

from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401,E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.pdb_io import ca_trace  # noqa: E402
from src.tc_profile import contact_map, fetch_pdb, gnm_correlation, tc_profile  # noqa: E402

PDB_LIST = [
    "1CRN", "1UBQ", "1L2Y", "2CI2", "1PGB", "1BDD", "1SHG", "1TIT",
    "1ENH", "1RIS", "1APS", "1STN", "3CHY", "1TEN", "1AJ3", "1PHT",
    "1IGD", "1PIN", "2PTL", "1CTF", "1FKB", "1CSP", "1SHF", "1PBA",
    "1HZ6", "1UBI", "2GB1", "1DIV", "1BEB", "1NLS",
]
MIN_LEN, MAX_LEN, WINDOW = 40, 160, 7

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.linewidth": 0.9,
    "figure.dpi": 300,
})


def tc_and_flexibility(pdb_id):
    ca, names, _ = ca_trace(fetch_pdb(pdb_id))
    contacts = contact_map(ca)
    gamma = np.diag(contacts.sum(axis=1)) - contacts
    flexibility = np.diag(np.linalg.pinv(gamma))
    prof = tc_profile(gnm_correlation(contacts), WINDOW)
    return ca, prof, flexibility


def main():
    # ---- left panel: ubiquitin coloured by TC ----
    ca, prof, _ = tc_and_flexibility("1UBQ")
    tail = slice(69, 76)  # residues 70-76, the flexible C-terminal tail

    fig = plt.figure(figsize=(11.5, 5.2))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.plot(ca[:, 0], ca[:, 1], ca[:, 2], color="#B8BFC6", linewidth=1.6, zorder=1)
    sc = ax.scatter(ca[:, 0], ca[:, 1], ca[:, 2], c=prof, cmap="viridis",
                    s=42, edgecolor="white", linewidth=0.4, depthshade=False, zorder=2)
    ax.plot(ca[tail, 0], ca[tail, 1], ca[tail, 2], color="#D62728", linewidth=3.4,
            zorder=3, solid_capstyle="round")
    ax.scatter(ca[tail, 0], ca[tail, 1], ca[tail, 2], color="#D62728", s=58,
               edgecolor="white", linewidth=0.5, depthshade=False, zorder=4)
    ax.text(ca[73, 0], ca[73, 1], ca[73, 2] + 6, "C-terminal tail\n(residues 70-76)",
            color="#D62728", fontsize=10, ha="center", va="bottom")
    ax.set_axis_off()
    ax.view_init(elev=20, azim=-58)
    ax.set_box_aspect((1, 1, 0.85))
    margin = 4.0
    for setter, values in (
        (ax.set_xlim, ca[:, 0]),
        (ax.set_ylim, ca[:, 1]),
        (ax.set_zlim, ca[:, 2]),
    ):
        setter(values.min() - margin, values.max() + margin)
    cbar = fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.04, pad=0.0, shrink=0.7)
    cbar.set_label("TC (nats)", fontsize=10)

    # ---- right panel: pooled TC vs flexibility ----
    xs, ys = [], []
    for pdb_id in PDB_LIST:
        try:
            _, prof_i, flex_i = tc_and_flexibility(pdb_id)
        except Exception:
            continue
        if not (MIN_LEN <= len(prof_i) <= MAX_LEN):
            continue
        mask = np.isfinite(prof_i) & np.isfinite(flex_i)
        if mask.sum() < 30:
            continue
        tc_z = (prof_i[mask] - prof_i[mask].mean()) / prof_i[mask].std()
        fx_z = (flex_i[mask] - flex_i[mask].mean()) / flex_i[mask].std()
        xs.append(tc_z)
        ys.append(fx_z)
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    r = np.corrcoef(x, y)[0, 1]

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.scatter(x, y, s=6, alpha=0.28, color="#2F6DB5", edgecolor="none")
    slope, intercept = np.polyfit(x, y, 1)
    grid = np.linspace(x.min(), x.max(), 50)
    ax2.plot(grid, slope * grid + intercept, color="#D62728", linewidth=1.8)
    ax2.set_xlabel("TC, z-scored within each protein", fontsize=11)
    ax2.set_ylabel("GNM flexibility, z-scored", fontsize=11)
    ax2.text(0.03, 0.95, f"n = {len(x):,} residues\n{len(xs)} proteins\nr = {r:+.3f}",
             transform=ax2.transAxes, va="top", fontsize=11)
    ax2.grid(False)
    ax2.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "figures", "fig0_hero.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print("wrote %s  (pooled r = %+.3f, n = %d residues from %d proteins)"
          % (out, r, len(x), len(xs)))


if __name__ == "__main__":
    main()
