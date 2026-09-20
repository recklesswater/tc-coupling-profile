"""Does BD concentrate on annotated functional sites?

The claim worth testing before making it in public: computing BD from a
structure can point at functionally important regions without wet-lab work.

We use PDB SITE records as a proxy for annotated functional sites (catalytic
residues, ligand-binding residues, and similar). For every structure with at
least one SITE record we compare the mean BD at those residues against the
mean over all residues of the same protein, so each protein acts as its own
control.

This is a deliberately weak test: SITE records are sparse and inconsistently
curated. A negative result here does not prove there is no effect; a positive
result would be encouraging but would still need a better annotation source.
"""

from __future__ import annotations

import os
import statistics
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.bd_profile import (  # noqa: E402
    contact_map,
    fetch_pdb,
    gnm_correlation,
    parse_ca,
    bd_profile,
)

PDB_LIST = [
    "1CRN", "1UBQ", "1L2Y", "2CI2", "1PGB", "1BDD", "1SHG", "1TIT",
    "1ENH", "1RIS", "1APS", "1STN", "3CHY", "1TEN", "1AJ3", "1PHT",
    "1IGD", "1PIN", "2PTL", "1CTF", "1FKB", "1CSP", "1SHF", "1PBA",
    "1HZ6", "1UBI", "2GB1", "1DIV", "1BEB", "1NLS",
]
MIN_LEN, MAX_LEN = 40, 160
WINDOW = 7


def parse_site_residues(path):
    """Residue sequence numbers listed in SITE records (first chain only)."""
    sites = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("SITE"):
                continue
            n = line[15:17].strip()
            try:
                n = int(n)
            except ValueError:
                continue
            for k in range(n):
                seq = line[22 + 11 * k:26 + 11 * k].strip()
                if seq.lstrip("-").isdigit():
                    sites.append(int(seq))
    return sorted(set(sites))


def main():
    hits = []
    for pdb_id in PDB_LIST:
        try:
            path = fetch_pdb(pdb_id)
            ca, names, _ = parse_ca(path, pdb_id)
        except Exception:  # noqa: BLE001
            continue
        n = len(ca)
        if not (MIN_LEN <= n <= MAX_LEN):
            continue

        sites = [s for s in parse_site_residues(path) if 1 <= s <= n]
        if not sites:
            continue

        corr = gnm_correlation(contact_map(ca))
        prof = bd_profile(corr, WINDOW)
        idx = [s - 1 for s in sites]
        site_tc = np.nanmean(prof[idx])
        bg_tc = np.nanmean(prof)
        # normalise by the spread of the profile so proteins are comparable
        spread = np.nanstd(prof)
        z = (site_tc - bg_tc) / spread if spread > 0 else float("nan")
        hits.append((pdb_id, len(sites), site_tc, bg_tc, z))
        print("  %-6s n_sites=%2d  site BD %.3f  background %.3f  z = %+.2f"
              % (pdb_id, len(sites), site_tc, bg_tc, z))

    print()
    if not hits:
        print("  no usable SITE annotations in this set")
        return

    zs = np.array([h[4] for h in hits], dtype=float)
    zs = zs[np.isfinite(zs)]
    print("  %d proteins with SITE records" % len(hits))
    print("  mean z (site BD vs protein background) = %+.3f" % zs.mean())
    print("  median z = %+.3f" % statistics.median(zs))
    print("  proteins with z > 0 = %d / %d" % (int((zs > 0).sum()), len(zs)))
    print()
    if zs.mean() > 0.2 and (zs > 0).mean() > 0.6:
        print("  => weak evidence that BD is enriched at annotated sites")
    else:
        print("  => no clear enrichment; do NOT claim site identification")


if __name__ == "__main__":
    main()
