"""Does TC concentrate on ligand-binding sites?

The claim worth testing before making it in public: computing the BD profile
from a structure can point at functionally important regions without wet-lab
work.

Test: take structures that have a bound ligand (HETATM, excluding water and
common ions), find the residues within 5 A of any ligand heavy atom, and
compare the mean TC there against the mean over the whole chain. Each protein
is its own control, and we report the difference in units of the profile's
own standard deviation.

This is a weak proxy -- a bound ligand marks a binding site, not necessarily a
catalytic one -- but it is fully automated and does not depend on curated
annotation being present in the file.
"""

from __future__ import annotations

import os
import statistics
import sys
import urllib.request

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.pdb_io import ca_trace, ligand_atoms  # noqa: E402
from src.bd_profile import (  # noqa: E402
    contact_map,
    fetch_pdb,
    gnm_correlation,
    bd_profile,
)

# Classic holo structures across several fold / function classes.
HOLO_LIST = [
    "3PTB",   # trypsin + benzamidine
    "1HSG",   # HIV-1 protease + inhibitor
    "4DFR",   # DHFR + methotrexate
    "1STP",   # streptavidin + biotin
    "3ERT",   # oestrogen receptor + tamoxifen
    "1M17",   # EGFR kinase + erlotinib
    "1LYZ",   # lysozyme
    "2CBA",   # carbonic anhydrase
]

WINDOW = 7
CUTOFF = 5.0
SKIP_HET = {"HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "ACT",
            "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "CU", "NI", "CD",
            "IOD", "BR", "FMT", "TRS", "MES", "DMS", "NH4", "NO3", "CIT"}


def analyse(pdb_id: str):
    path = fetch_pdb(pdb_id)
    ca, names, resseqs = ca_trace(path)
    if len(ca) < 40:
        return None
    lig = ligand_atoms(path, skip_names=SKIP_HET)
    if lig.size == 0:
        return None

    # residues whose CA is near any ligand atom
    d = np.linalg.norm(ca[:, None, :] - lig[None, :, :], axis=-1).min(axis=1)
    site_idx = np.where(d < CUTOFF)[0]
    if len(site_idx) < 3:
        return None

    contacts = contact_map(ca)
    prof = bd_profile(gnm_correlation(contacts), WINDOW)
    site_tc = float(np.nanmean(prof[site_idx]))
    bg_tc = float(np.nanmean(prof))
    spread = float(np.nanstd(prof))
    z = (site_tc - bg_tc) / spread if spread > 0 else float("nan")
    return dict(pdb=pdb_id, n_res=len(ca), n_site=len(site_idx),
                site_tc=site_tc, bg_tc=bg_tc, z=z)


def main():
    rows = []
    for pdb_id in HOLO_LIST:
        try:
            r = analyse(pdb_id)
        except Exception as exc:  # noqa: BLE001
            print("  %-6s skipped (%s)" % (pdb_id, repr(exc)[:45]))
            continue
        if r is None:
            print("  %-6s skipped (no usable ligand / too small)" % pdb_id)
            continue
        rows.append(r)
        print("  %-6s n=%3d  site residues=%2d  site TC %.3f  background %.3f  z=%+.2f"
              % (r["pdb"], r["n_res"], r["n_site"], r["site_tc"], r["bg_tc"], r["z"]))

    print()
    if not rows:
        print("  nothing usable")
        return

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(os.path.join(root, "results"), exist_ok=True)
    out = os.path.join(root, "results", "ligand_site_check.csv")
    keys = list(rows[0].keys())
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(
                r[k] if isinstance(r[k], str) else "%.6f" % r[k] for k in keys
            ) + "\n")
    print("  wrote %s" % out)
    print()

    zs = np.array([r["z"] for r in rows])
    print("  %d structures" % len(rows))
    print("  mean z  = %+.3f" % zs.mean())
    print("  median z = %+.3f" % statistics.median(zs))
    print("  z > 0: %d / %d" % (int((zs > 0).sum()), len(zs)))
    print()
    if zs.mean() > 0.2 and (zs > 0).mean() > 0.6:
        print("  => weak positive: TC tends to be higher at ligand-binding residues")
    else:
        print("  => no clear enrichment; do NOT claim that TC identifies binding sites")


if __name__ == "__main__":
    main()
