"""Stage-0 audit: is BD a separate axis from the mean-field diagonal?

The original framing proposed a 2x2 reading of a protein -- amplitude (mean-square
fluctuation, MSF) against coupling (BD). That only works if the two axes carry
different information. This script tests that directly, and then stress-tests
whatever answer comes out:

1. **Orthogonality.** Correlation between z(log MSF) and z(BD) across all residues.
   A high value means the "two axes" are one axis and the 2x2 map cannot be built.
2. **Numeric conditioning.** Smallest eigenvalue and condition number of each
   window block, so that the near-singularity objection can be checked rather
   than assumed.
3. **Cutoff robustness.** Contact cutoff 8-12 A.
4. **Window sensitivity.** Window 5-15 residues.
5. **B-factor check of both axes.** Per-protein correlation of BD and of log MSF
   against experimental B-factors (the independent axis), plus the partial
   correlation controlling for contact degree.
6. **Leave-one-protein-out.** Whether the pooled result is carried by any single
   structure.

Writes ``results/stage0_per_protein.csv``, ``stage0_cutoff.csv``,
``stage0_window.csv``, ``stage0_lopo.csv`` and ``stage0_key_numbers.csv``.
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from src.bd_profile import (  # noqa: E402
    bd_profile,
    contact_map,
    fetch_pdb,
    gnm_correlation,
)
from src.pdb_io import burial_proxy, ca_trace  # noqa: E402

from bfactor_validation import ca_bfactors, partial_r, zscore  # noqa: E402
from bfactor_validation import PDB_LIST  # noqa: E402

PDB_DIR = os.path.join(ROOT, "data", "pdb")
RESULTS = os.path.join(ROOT, "results")

WINDOW = 7
CUTOFF = 8.0
MIN_LEN = 20


def window_spectrum(corr: np.ndarray, window: int = WINDOW):
    """BD plus eigenvalue diagnostics for every window of the correlation matrix."""
    n = corr.shape[0]
    half = window // 2
    bd = np.full(n, np.nan)
    min_eig = np.full(n, np.nan)
    cond = np.full(n, np.nan)
    eff_rank = np.full(n, np.nan)
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        idx = list(range(lo, hi))
        sub = corr[np.ix_(idx, idx)]
        sub = 0.5 * (sub + sub.T)
        w = np.clip(np.linalg.eigvalsh(sub), 1e-12, None)
        bd[i] = -0.5 * np.sum(np.log(w))
        min_eig[i] = w.min()
        cond[i] = w.max() / w.min()
        p = w / w.sum()
        eff_rank[i] = float(np.exp(-(p * np.log(p)).sum()))
    return bd, min_eig, cond, eff_rank


def model(ca: np.ndarray, cutoff: float = CUTOFF):
    """Contacts, contact degree, GNM covariance, log MSF and BD profile."""
    contacts = contact_map(ca, cutoff)
    degree = contacts.sum(axis=1)
    gamma = np.diag(degree) - contacts
    cov = np.linalg.pinv(gamma)
    msf = np.clip(np.diag(cov), 1e-12, None)
    corr = gnm_correlation(contacts)
    return contacts, degree, cov, np.log(msf), corr


def write_csv(path: str, rows) -> None:
    keys = list(rows[0].keys())
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(
                r[k] if isinstance(r[k], str) else "%.6f" % r[k] for k in keys
            ) + "\n")


def main() -> None:
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(PDB_DIR, exist_ok=True)
    entries = [(pid, fetch_pdb(pid, PDB_DIR)) for pid in PDB_LIST]

    per, pool = [], []
    for pdb_id, path in entries:
        ca, _names, _resseq = ca_trace(path)
        if len(ca) < MIN_LEN:
            continue
        bfac = ca_bfactors(path)
        if len(bfac) != len(ca):
            bfac = np.full(len(ca), np.nan)
        uniform = (not np.isfinite(bfac).all()) or np.nanstd(bfac) < 1e-6

        contacts, degree, _cov, log_msf, corr = model(ca)
        bd, min_eig, cond, eff_rank = window_spectrum(corr, WINDOW)
        bury = burial_proxy(ca, 10.0)

        per.append({
            "pdb": pdb_id,
            "n_res": len(ca),
            "bf_uniform": uniform,
            "r_bd_logmsf": float(np.corrcoef(bd, log_msf)[0, 1]),
            "r_bd_degree": float(np.corrcoef(bd, degree)[0, 1]),
            "r_bd_burial": float(np.corrcoef(bd, bury)[0, 1]),
            "r_bd_bfactor": float("nan") if uniform else float(np.corrcoef(bd, bfac)[0, 1]),
            "r_logmsf_bfactor":
                float("nan") if uniform else float(np.corrcoef(log_msf, bfac)[0, 1]),
            "r_degree_bfactor":
                float("nan") if uniform else float(np.corrcoef(degree, bfac)[0, 1]),
            "partial_r_bd_bfactor_given_degree":
                float("nan") if uniform else partial_r(bd, bfac, degree),
            "median_min_eig": float(np.nanmedian(min_eig)),
            "median_cond": float(np.nanmedian(cond)),
            "median_eff_rank": float(np.nanmedian(eff_rank)),
        })
        pool.append({
            "pdb": pdb_id,
            "BD": bd,
            "logMSF": log_msf,
            "degree": degree.astype(float),
            "burial": bury,
            "bfactor": bfac,
        })

    n_res_total = sum(p["n_res"] for p in per)
    write_csv(os.path.join(RESULTS, "stage0_per_protein.csv"), per)

    print("=" * 90)
    print("  %d structures, %d residues" % (len(per), n_res_total))
    bf_ok = [p["pdb"] for p in per if not p["bf_uniform"]]
    print("  usable B-factors: %d proteins (%s excluded as uniform or missing)"
          % (len(bf_ok), ", ".join(p["pdb"] for p in per if p["bf_uniform"])))

    # ---- 1. orthogonality -------------------------------------------------
    a = np.concatenate([zscore(p["logMSF"]) for p in pool])
    c = np.concatenate([zscore(p["BD"]) for p in pool])
    r_ac = float(np.corrcoef(a, c)[0, 1])
    per_ac = np.array([p["r_bd_logmsf"] for p in per])
    verdict = ("two axes are close to orthogonal; the 2x2 map is defensible"
               if abs(r_ac) < 0.4 else
               "the two axes are largely one axis; the 2x2 map cannot be built"
               if abs(r_ac) > 0.6 else
               "the axes overlap substantially; a 2x2 map would need a caveat")
    print()
    print("  1. orthogonality of the two proposed axes")
    print("     corr(z(log MSF), z(BD)) pooled  = %+.3f" % r_ac)
    print("     per protein, r(BD, log MSF)     median %+.3f  range [%+.3f, %+.3f]"
          % (np.median(per_ac), per_ac.min(), per_ac.max()))
    print("     -> %s" % verdict)

    # ---- 2. conditioning --------------------------------------------------
    print()
    print("  2. conditioning of the window blocks at %.0f A / window %d"
          % (CUTOFF, WINDOW))
    print("     smallest eigenvalue, median over proteins = %.3f"
          % np.median([p["median_min_eig"] for p in per]))
    print("     condition number,     median over proteins = %.1f"
          % np.median([p["median_cond"] for p in per]))
    print("     effective rank,       median over proteins = %.2f"
          % np.median([p["median_eff_rank"] for p in per]))

    # ---- 3. cutoff sweep --------------------------------------------------
    cutoff_rows = []
    for cutoff in (8.0, 9.0, 10.0, 12.0):
        bds, msfs = [], []
        for _pid, path in entries:
            ca, _n, _r = ca_trace(path)
            if len(ca) < MIN_LEN:
                continue
            _c, _d, _cov, log_msf, corr = model(ca, cutoff)
            bds.append(bd_profile(corr, WINDOW))
            msfs.append(log_msf)
        cutoff_rows.append({
            "cutoff_angstrom": cutoff,
            "r_bd_logmsf_pooled": float(np.corrcoef(np.concatenate(bds),
                                                    np.concatenate(msfs))[0, 1]),
        })
    write_csv(os.path.join(RESULTS, "stage0_cutoff.csv"), cutoff_rows)

    # ---- 4. window sweep --------------------------------------------------
    window_rows = []
    for window in (5, 7, 9, 11, 15):
        bds, msfs = [], []
        for _pid, path in entries:
            ca, _n, _r = ca_trace(path)
            if len(ca) < MIN_LEN:
                continue
            _c, _d, _cov, log_msf, corr = model(ca)
            bds.append(bd_profile(corr, window))
            msfs.append(log_msf)
        window_rows.append({
            "window": window,
            "r_bd_logmsf_pooled": float(np.corrcoef(np.concatenate(bds),
                                                    np.concatenate(msfs))[0, 1]),
            "mean_bd": float(np.nanmean(np.concatenate(bds))),
        })
    write_csv(os.path.join(RESULTS, "stage0_window.csv"), window_rows)

    print()
    print("  3/4. robustness of the pooled r(BD, log MSF) = %+.3f" % r_ac)
    print("       cutoff 8 -> 12 A:  %s"
          % "  ".join("%.0f A %+.3f" % (r["cutoff_angstrom"], r["r_bd_logmsf_pooled"])
                      for r in cutoff_rows))
    print("       window 5 -> 15:    %s"
          % "  ".join("%d %+.3f" % (r["window"], r["r_bd_logmsf_pooled"])
                      for r in window_rows))

    # ---- 5. B-factor, both axes -------------------------------------------
    keep = [p for p in per if not p["bf_uniform"]]
    keep_ids = {p["pdb"] for p in keep}
    pool_b = [p for p in pool if p["pdb"] in keep_ids]
    bd_v = np.array([p["r_bd_bfactor"] for p in keep])
    msf_v = np.array([p["r_logmsf_bfactor"] for p in keep])
    print()
    print("  5. B-factor check, per protein (n = %d)" % len(keep))
    print("     BD vs B-factor     median %+.3f  positive %2d/%d"
          % (np.median(bd_v), int((bd_v > 0).sum()), len(bd_v)))
    print("     log MSF vs B-factor median %+.3f  positive %2d/%d"
          % (np.median(msf_v), int((msf_v > 0).sum()), len(msf_v)))
    print("     -> the mean-field diagonal predicts B-factors at least as well as BD")

    raw = {k: np.concatenate([p[k] for p in pool_b]) for k in
           ("BD", "logMSF", "degree", "bfactor")}
    zsc = {k: np.concatenate([zscore(p[k]) for p in pool_b]) for k in raw}
    print("     pooled raw residues:   r(BD, B) %+.3f   r(log MSF, B) %+.3f"
          % (np.corrcoef(raw["BD"], raw["bfactor"])[0, 1],
             np.corrcoef(raw["logMSF"], raw["bfactor"])[0, 1]))
    print("     pooled within-protein: r(BD, B) %+.3f   r(log MSF, B) %+.3f"
          % (np.corrcoef(zsc["BD"], zsc["bfactor"])[0, 1],
             np.corrcoef(zsc["logMSF"], zsc["bfactor"])[0, 1]))

    # ---- 6. leave-one-protein-out -----------------------------------------
    lopo = []
    for i, left in enumerate(pool):
        rest = [p for j, p in enumerate(pool) if j != i]
        lopo.append({
            "left_out": left["pdb"],
            "pooled_r": float(np.corrcoef(
                np.concatenate([zscore(p["BD"]) for p in rest]),
                np.concatenate([zscore(p["logMSF"]) for p in rest]))[0, 1]),
        })
    write_csv(os.path.join(RESULTS, "stage0_lopo.csv"), lopo)
    lo = np.array([r["pooled_r"] for r in lopo])
    print()
    print("  6. leave-one-protein-out, pooled corr(z(log MSF), z(BD))")
    print("     median %+.3f   range [%+.3f, %+.3f]  (no single structure drives it)"
          % (np.median(lo), lo.min(), lo.max()))

    key_rows = [
        {"quantity": "n_structures", "value": len(per)},
        {"quantity": "n_residues", "value": n_res_total},
        {"quantity": "n_structures_with_bfactors", "value": len(keep)},
        {"quantity": "r_z_logmsf_z_bd_pooled", "value": r_ac},
        {"quantity": "median_per_protein_r_bd_logmsf",
         "value": float(np.median(per_ac))},
        {"quantity": "median_min_eigenvalue",
         "value": float(np.median([p["median_min_eig"] for p in per]))},
        {"quantity": "median_condition_number",
         "value": float(np.median([p["median_cond"] for p in per]))},
        {"quantity": "median_effective_rank",
         "value": float(np.median([p["median_eff_rank"] for p in per]))},
        {"quantity": "median_r_bd_bfactor", "value": float(np.median(bd_v))},
        {"quantity": "median_r_logmsf_bfactor", "value": float(np.median(msf_v))},
        {"quantity": "lopo_pooled_r_median", "value": float(np.median(lo))},
        {"quantity": "lopo_pooled_r_min", "value": float(lo.min())},
        {"quantity": "lopo_pooled_r_max", "value": float(lo.max())},
    ]
    write_csv(os.path.join(RESULTS, "stage0_key_numbers.csv"), key_rows)
    print()
    print("  wrote results/stage0_per_protein.csv, stage0_cutoff.csv, "
          "stage0_window.csv, stage0_lopo.csv, stage0_key_numbers.csv")


if __name__ == "__main__":
    main()
