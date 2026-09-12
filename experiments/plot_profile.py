"""Plot the TC profile for one structure, next to its contact map.

Writes ``figures/fig2_profile_ubiquitin.png``.
"""

from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.pdb_io import burial_proxy, ca_trace  # noqa: E402
from src.tc_profile import (  # noqa: E402
    contact_map,
    fetch_pdb,
    gnm_correlation,
    tc_profile,
)

plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 140
plt.rcParams["savefig.bbox"] = "tight"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDB_ID = "1UBQ"
WINDOW = 7


def main():
    path = fetch_pdb(PDB_ID)
    ca, names, _ = ca_trace(path)
    contacts = contact_map(ca)
    degree = contacts.sum(axis=1)
    corr = gnm_correlation(contacts)
    prof = tc_profile(corr, WINDOW)

    n = len(ca)
    resids = np.arange(1, n + 1)
    burial = burial_proxy(ca, cutoff=10.0)
    burial_n = (burial - burial.min()) / max(1e-9, burial.max() - burial.min())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))

    ax = axes[0]
    im = ax.imshow(contacts, cmap="Greys", interpolation="nearest")
    ax.set_title("(a) contact map, C-alpha < 8 A\nmean degree %.1f" % degree.mean())
    ax.set_xlabel("residue")
    ax.set_ylabel("residue")
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[1]
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
    ax.set_title("(b) GNM correlation matrix")
    ax.set_xlabel("residue")
    ax.set_ylabel("residue")
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[2]
    ax.plot(resids, prof, "-o", ms=3, color="#2f855a", label="TC profile")
    ax.plot(resids, burial_n * np.nanmax(prof), "-", lw=1.2, color="#2b6cb0",
            alpha=0.7, label="burial proxy, scaled")
    peak = int(np.nanargmax(prof))
    ax.axvspan(max(1, peak - 3), min(n, peak + 4), color="#e53e3e", alpha=0.13)
    ax.annotate("peak TC at residue %d" % resids[peak],
                xy=(resids[peak], prof[peak]), xytext=(resids[peak] - 30, prof[peak]),
                fontsize=9, arrowprops=dict(arrowstyle="->", lw=1))
    ax.set_xlabel("residue number")
    ax.set_ylabel("TC (nats)")
    ax.set_title("(c) TC profile, window = %d" % WINDOW)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.25)

    fig.suptitle("%s: coupling indicator computed from the structure alone"
                 % PDB_ID, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(ROOT, "figures", "fig2_profile_ubiquitin.png")
    fig.savefig(out)
    plt.close(fig)

    print("%s: %d residues" % (PDB_ID, n))
    print("TC range %.3f - %.3f, mean %.3f"
          % (np.nanmin(prof), np.nanmax(prof), np.nanmean(prof)))
    print("highest TC at residue %d (%s)" % (resids[peak], names[peak]))
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
