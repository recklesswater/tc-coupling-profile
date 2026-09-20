"""Independent check: does the BD profile track experimental B-factors?

Everything else in this repository validates BD against quantities that come
out of the same elastic network (mean-square fluctuation, contact degree,
burial computed from the same C-alpha trace). Those tests cannot establish
that BD tracks real motion, because they are correlations between functionals
of one matrix.

B-factors are different: they ship with the PDB file and owe nothing to the
GNM. They are a weak proxy -- crystal packing, resolution and refinement all
leak into them -- but they are *independent*, which is the property that was
missing.

Per protein we correlate BD, the GNM mean-square fluctuation (MSF) and contact
degree against the C-alpha B-factor of the first chain, and we compute the
partial correlation of BD with B-factor controlling for contact degree.

Two pooling conventions are reported, because they disagree and the
disagreement is the point:

* ``pooled_raw``   -- concatenate residues from all proteins as they are.
  Between-protein variation (some structures are simply more mobile than
  others) enters the correlation.
* ``pooled_z``     -- z-score each quantity within its own protein first, then
  pool. This removes between-protein variation and is the honest analogue of
  the per-protein correlations.

Averaging Pearson r over proteins is also skewed, so both the arithmetic mean
and the Fisher-z mean (mean of artanh(r), back-transformed) are reported.

Writes ``results/bfactor_validation.csv`` (per protein) and
``results/bfactor_validation_summary.csv`` (pooled and averaged values).
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.pdb_io import burial_proxy, ca_trace  # noqa: E402
from src.bd_profile import (  # noqa: E402
    bd_profile,
    contact_map,
    fetch_pdb,
    gnm_correlation,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDB_DIR = os.path.join(ROOT, "data", "pdb")
RESULTS = os.path.join(ROOT, "results")

CUTOFF = 8.0
WINDOW = 7
MIN_LEN = 20

# The 30 structures used in the multi-protein validation; downloaded on first run.
PDB_LIST = [
    "1CRN", "1UBQ", "1L2Y", "2CI2", "1PGB", "1BDD", "1SHG", "1TIT",
    "1ENH", "1RIS", "1APS", "1STN", "3CHY", "1TEN", "1AJ3", "1PHT",
    "1IGD", "1PIN", "2PTL", "1CTF", "1FKB", "1CSP", "1SHF", "1PBA",
    "1HZ6", "1UBI", "2GB1", "1DIV", "1BEB", "1NLS",
]


def ca_bfactors(path: str) -> np.ndarray:
    """C-alpha B-factors of the first chain, in file order.

    ``ca_trace`` reads the first chain of the first model, so mirror that here
    to keep the two arrays aligned.
    """
    out = []
    chain_id = None
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            record = line[0:6].strip()
            if record == "ENDMDL":
                break
            if record != "ATOM":
                continue
            if chain_id is None:
                chain_id = line[21:22].strip()
            if line[21:22].strip() != chain_id:
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                out.append(float(line[60:66]))
            except ValueError:
                out.append(np.nan)
    return np.asarray(out, dtype=float)


def partial_r(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """r(x, y | z), by linear residuals."""
    x, y, z = (np.asarray(v, dtype=float) for v in (x, y, z))
    if z.std() == 0:
        return float("nan")
    rx = x - np.polyval(np.polyfit(z, x, 1), z)
    ry = y - np.polyval(np.polyfit(z, y, 1), z)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def fisher_z(values: np.ndarray) -> float:
    """Mean of artanh(r), back-transformed."""
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    v = np.clip(v, -0.999999, 0.999999)
    return float(np.tanh(np.arctanh(v).mean()))


def zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    s = v.std()
    return (v - v.mean()) / (s if s > 0 else 1.0)


def main() -> None:
    os.makedirs(PDB_DIR, exist_ok=True)
    entries = [(pid, fetch_pdb(pid, PDB_DIR)) for pid in PDB_LIST]
    rows, pooled = [], []
    dropped = []

    for pdb_id, path in entries:
        ca, _names, _resseq = ca_trace(path)
        if len(ca) < MIN_LEN:
            continue
        bfac = ca_bfactors(path)
        if len(bfac) != len(ca):
            bfac = np.full(len(ca), np.nan)
        if not np.isfinite(bfac).all() or np.nanstd(bfac) < 1e-6:
            # Uniform B-factors mean a predicted model or fixed-B refinement;
            # there is nothing to correlate against.
            dropped.append(pdb_id)
            continue

        contacts = contact_map(ca, CUTOFF)
        degree = contacts.sum(axis=1)
        gamma = np.diag(degree) - contacts
        cov = np.linalg.pinv(gamma)
        msf = np.clip(np.diag(cov), 1e-12, None)
        corr = gnm_correlation(contacts)
        prof = bd_profile(corr, WINDOW)

        rows.append({
            "pdb": pdb_id,
            "n_res": len(ca),
            "r_bd_bfactor": float(np.corrcoef(prof, bfac)[0, 1]),
            "r_msf_bfactor": float(np.corrcoef(np.log(msf), bfac)[0, 1]),
            "r_degree_bfactor": float(np.corrcoef(degree, bfac)[0, 1]),
            "partial_r_bd_bfactor_given_degree":
                partial_r(prof, bfac, degree),
        })
        pooled.append({
            "pdb": pdb_id,
            "BD": prof,
            "logMSF": np.log(msf),
            "degree": degree.astype(float),
            "burial": burial_proxy(ca, 10.0),
            "bfactor": bfac,
        })
        r = rows[-1]
        print("  %-6s n=%3d  BD %+.3f  MSF %+.3f  degree %+.3f  partial(BD|deg) %+.3f"
              % (pdb_id, r["n_res"], r["r_bd_bfactor"], r["r_msf_bfactor"],
                 r["r_degree_bfactor"],
                 r["partial_r_bd_bfactor_given_degree"]))

    if not rows:
        print("no structures with usable B-factors")
        return

    os.makedirs(RESULTS, exist_ok=True)
    keys = list(rows[0].keys())
    per_path = os.path.join(RESULTS, "bfactor_validation.csv")
    with open(per_path, "w", encoding="utf-8") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(
                r[k] if isinstance(r[k], str) else "%.6f" % r[k] for k in keys
            ) + "\n")

    raw = {k: np.concatenate([p[k] for p in pooled]) for k in
           ("BD", "logMSF", "degree", "bfactor")}
    zsc = {k: np.concatenate([zscore(p[k]) for p in pooled]) for k in
           ("BD", "logMSF", "degree", "bfactor")}

    def pooled_stats(d):
        return {
            "r_bd_bfactor": float(np.corrcoef(d["BD"], d["bfactor"])[0, 1]),
            "r_msf_bfactor": float(np.corrcoef(d["logMSF"], d["bfactor"])[0, 1]),
            "r_degree_bfactor": float(np.corrcoef(d["degree"], d["bfactor"])[0, 1]),
            "partial_r_bd_bfactor_given_degree":
                partial_r(d["BD"], d["bfactor"], d["degree"]),
            "r_bd_logmsf": float(np.corrcoef(d["BD"], d["logMSF"])[0, 1]),
        }

    prm, praw = pooled_stats(zsc), pooled_stats(raw)
    arr = {k: np.array([r[k] for r in rows]) for k in keys if k != "pdb"}

    print()
    print("  structures with usable B-factors: %d (%d residues)"
          % (len(rows), sum(r["n_res"] for r in rows)))
    if dropped:
        print("  excluded, uniform or missing B-factors: %s" % ", ".join(dropped))
    print()
    print("  per-protein correlations")
    for key, label in (("r_bd_bfactor", "BD vs B-factor"),
                       ("r_msf_bfactor", "log MSF vs B-factor"),
                       ("r_degree_bfactor", "contact degree vs B-factor"),
                       ("partial_r_bd_bfactor_given_degree",
                        "partial: BD vs B-factor | degree")):
        v = arr[key]
        print("    %-34s median %+.3f  fisher-z %+.3f  positive %2d/%d  range [%+.2f, %+.2f]"
              % (label, np.median(v), fisher_z(v), int((v > 0).sum()), len(v),
                 v.min(), v.max()))
    print()
    print("  pooled, residues concatenated as they are (between-protein variation kept)")
    for key in praw:
        print("    %-34s %+.3f" % (key, praw[key]))
    print("  pooled, z-scored within each protein first (between-protein variation removed)")
    for key in prm:
        print("    %-34s %+.3f" % (key, prm[key]))

    summary = []
    for key in praw:
        summary.append({
            "quantity": key,
            "per_protein_median": float(np.median(arr[key])) if key in arr else float("nan"),
            "per_protein_fisher_z": fisher_z(arr[key]) if key in arr else float("nan"),
            "n_positive": int((arr[key] > 0).sum()) if key in arr else 0,
            "n_proteins": len(rows) if key in arr else 0,
            "pooled_raw": praw[key],
            "pooled_z": prm[key],
        })
    summary_path = os.path.join(RESULTS, "bfactor_validation_summary.csv")
    skeys = list(summary[0].keys())
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(",".join(skeys) + "\n")
        for s in summary:
            fh.write(",".join(
                s[k] if isinstance(s[k], str) else "%.6f" % s[k] for k in skeys
            ) + "\n")

    print()
    print("  wrote %s" % per_path)
    print("  wrote %s" % summary_path)


if __name__ == "__main__":
    main()
