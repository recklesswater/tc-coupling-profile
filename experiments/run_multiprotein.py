"""Validate the BD profile across a set of single-chain proteins.

For each structure we compute the BD profile and correlate it against four
per-residue properties: GNM flexibility (mean-square fluctuation), contact
degree, burial (1 - normalised SASA) and hydrophobicity (Kyte-Doolittle).

Writes ``results/multiprotein_correlations.csv`` and
``figures/fig1_multiprotein.png``.
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
from src.bd_profile import (  # noqa: E402
    contact_map,
    fetch_pdb,
    gnm_correlation,
    bd_profile,
)

plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 140
plt.rcParams["savefig.bbox"] = "tight"

PDB_LIST = [
    "1CRN", "1UBQ", "1L2Y", "2CI2", "1PGB", "1BDD", "1SHG", "1TIT",
    "1ENH", "1RIS", "1APS", "1STN", "3CHY", "1TEN", "1AJ3", "1PHT",
    "1IGD", "1PIN", "2PTL", "1CTF", "1FKB", "1CSP", "1SHF", "1PBA",
    "1HZ6", "1UBI", "2GB1", "1DIV", "1BEB", "1NLS",
]
MIN_LEN, MAX_LEN = 40, 160
WINDOW = 7

KD = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    rows = []
    for pdb_id in PDB_LIST:
        try:
            path = fetch_pdb(pdb_id)
            ca, names, _ = ca_trace(path)
        except Exception as exc:  # noqa: BLE001
            print("  %-6s skipped (%s)" % (pdb_id, repr(exc)[:40]))
            continue
        n = len(ca)
        if not (MIN_LEN <= n <= MAX_LEN):
            print("  %-6s skipped (length %d)" % (pdb_id, n))
            continue

        contacts = contact_map(ca)
        degree = contacts.sum(axis=1)
        gamma = np.diag(degree) - contacts
        flexibility = np.diag(np.linalg.pinv(gamma))
        corr = gnm_correlation(contacts)
        prof = bd_profile(corr, WINDOW)

        hydro = np.array([KD.get(nm, np.nan) for nm in names])
        # Burial proxy: number of C-alpha neighbours within 10 A. High = buried.
        # (A real SASA would need a Biopython compiled extension; see
        #  src/pdb_io.py for why that dependency was removed.)
        burial = burial_proxy(ca, cutoff=10.0)

        m = np.isfinite(prof) & np.isfinite(hydro) & np.isfinite(burial)
        if m.sum() < 30:
            continue
        rows.append({
            "pdb": pdb_id,
            "n_res": n,
            "mean_contact_degree": float(degree.mean()),
            "r_flexibility": float(np.corrcoef(prof[m], flexibility[m])[0, 1]),
            "r_contact_degree": float(np.corrcoef(prof[m], degree[m])[0, 1]),
            "r_burial": float(np.corrcoef(prof[m], burial[m])[0, 1]),
            "r_hydrophobicity": float(np.corrcoef(prof[m], hydro[m])[0, 1]),
        })
        r = rows[-1]
        print("  %-6s n=%3d  flex %+.3f  degree %+.3f  burial %+.3f  hydro %+.3f"
              % (pdb_id, n, r["r_flexibility"], r["r_contact_degree"],
                 r["r_burial"], r["r_hydrophobicity"]))

    if not rows:
        print("no structures parsed")
        return

    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    keys = list(rows[0].keys())
    out = os.path.join(ROOT, "results", "multiprotein_correlations.csv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(
                r[k] if isinstance(r[k], str) else "%.6f" % r[k] for k in keys
            ) + "\n")

    print()
    print("  %d proteins" % len(rows))
    summary = {}
    for key, label in (("r_flexibility", "TC vs flexibility"),
                       ("r_contact_degree", "TC vs contact degree"),
                       ("r_burial", "TC vs burial"),
                       ("r_hydrophobicity", "TC vs hydrophobicity")):
        v = np.array([r[key] for r in rows])
        summary[key] = v
        print("    %-24s mean %+.3f  median %+.3f  negative %2d/%d"
              % (label, v.mean(), np.median(v), int((v < 0).sum()), len(v)))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    colors = {"r_flexibility": "#2f855a", "r_contact_degree": "#2b6cb0",
              "r_burial": "#805ad5", "r_hydrophobicity": "#dd6b20"}
    panels = [("r_flexibility", "correlation with flexibility"),
              ("r_hydrophobicity", "correlation with hydrophobicity")]
    for ax, (key, label) in zip(axes, panels):
        v = summary[key]
        ax.hist(v, bins=14, color=colors[key], alpha=0.85)
        ax.axvline(0, color="k", lw=1)
        ax.axvline(v.mean(), color="#e53e3e", lw=1.6, ls="--",
                   label="mean %+.3f" % v.mean())
        ax.set_xlabel("Pearson r (per protein)")
        ax.set_ylabel("number of proteins")
        ax.set_title(label)
        ax.legend()
        ax.grid(alpha=0.25)
    fig.suptitle("BD profile across %d single-chain proteins" % len(rows), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    figpath = os.path.join(ROOT, "figures", "fig1_multiprotein.png")
    fig.savefig(figpath)
    plt.close(fig)
    print()
    print("  wrote %s" % out)
    print("  wrote %s" % figpath)


if __name__ == "__main__":
    main()
