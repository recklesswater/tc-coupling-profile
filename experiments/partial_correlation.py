"""Is TC more than a repackaging of contact density?

This is the sharpest objection the work faces, and it is a fair one. TC and
contact degree are both functions of the same contact graph, so a correlation
between them is guaranteed and tells us nothing. What matters is whether TC
carries information that contact density does not.

Two things are done here.

1. **Fisher z averaging.** Pearson r is skewed, so averaging r values across
   proteins is not the right thing to do; the standard procedure is to average
   ``z = artanh(r)`` and transform back. Both are reported so the size of the
   correction is visible.

2. **Partial correlation.** Within each protein we z-score TC, contact degree
   and flexibility, pool the residues, and ask whether TC still predicts
   flexibility once contact degree is controlled for:

       r_partial(TC, flex | degree)
         = (r_TC,flex - r_TC,deg * r_flex,deg)
           / sqrt((1 - r_TC,deg^2) * (1 - r_flex,deg^2))

   A partial correlation near zero means TC adds nothing beyond contact
   density. A clearly non-zero value means it does.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.pdb_io import burial_proxy, ca_trace  # noqa: E402
from src.tc_profile import (  # noqa: E402
    contact_map,
    fetch_pdb,
    gnm_correlation,
    tc_profile,
)

PDB_LIST = [
    "1CRN", "1UBQ", "1L2Y", "2CI2", "1PGB", "1BDD", "1SHG", "1TIT",
    "1ENH", "1RIS", "1APS", "1STN", "3CHY", "1TEN", "1AJ3", "1PHT",
    "1IGD", "1PIN", "2PTL", "1CTF", "1FKB", "1CSP", "1SHF", "1PBA",
    "1HZ6", "1UBI", "2GB1", "1DIV", "1BEB", "1NLS",
]
MIN_LEN, MAX_LEN = 40, 160
WINDOW = 7


def fisher_mean(rs):
    """Average correlations the correct way: mean of artanh(r), back-transformed."""
    r = np.clip(np.asarray(rs, float), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(r))))


def zscore(x):
    s = np.std(x)
    return (x - np.mean(x)) / s if s > 1e-12 else np.zeros_like(x)


def main():
    rows = []
    pool_tc, pool_deg, pool_flex, pool_bur = [], [], [], []

    for pdb_id in PDB_LIST:
        try:
            path = fetch_pdb(pdb_id)
            ca, names, _ = ca_trace(path)
        except Exception:  # noqa: BLE001
            continue
        n = len(ca)
        if not (MIN_LEN <= n <= MAX_LEN):
            continue

        contacts = contact_map(ca)
        degree = contacts.sum(axis=1).astype(float)
        gamma = np.diag(degree) - contacts
        flexibility = np.diag(np.linalg.pinv(gamma))
        corr = gnm_correlation(contacts)
        prof = tc_profile(corr, WINDOW)
        burial = burial_proxy(ca, 10.0)

        m = np.isfinite(prof) & np.isfinite(flexibility)
        if m.sum() < 30:
            continue
        p, d, f, b = prof[m], degree[m], flexibility[m], burial[m]
        rows.append({
            "pdb": pdb_id,
            "r_flex": np.corrcoef(p, f)[0, 1],
            "r_deg": np.corrcoef(p, d)[0, 1],
            "r_bur": np.corrcoef(p, b)[0, 1],
            # within-protein partial correlation of TC and flexibility,
            # controlling contact degree
            "r_partial": partial(p, f, d),
        })
        pool_tc.append(zscore(p))
        pool_deg.append(zscore(d))
        pool_flex.append(zscore(f))
        pool_bur.append(zscore(b))

    print("=" * 78)
    print("Part 1. Fisher z averaging vs naive arithmetic mean  (n = %d proteins)"
          % len(rows))
    print("=" * 78)
    for key, label in (("r_flex", "TC vs flexibility"),
                       ("r_deg", "TC vs contact degree"),
                       ("r_bur", "TC vs burial")):
        vals = np.array([r[key] for r in rows])
        print("  %-24s  mean %.4f   Fisher-z mean %.4f   (difference %+.4f)"
              % (label, vals.mean(), fisher_mean(vals),
                 fisher_mean(vals) - vals.mean()))

    tc = np.concatenate(pool_tc)
    deg = np.concatenate(pool_deg)
    flex = np.concatenate(pool_flex)
    bur = np.concatenate(pool_bur)
    print()
    print("  pooled residues: %d" % len(tc))

    print()
    print("=" * 78)
    print("Part 2. Does TC add anything beyond contact density?")
    print("=" * 78)
    r_tf = np.corrcoef(tc, flex)[0, 1]
    r_td = np.corrcoef(tc, deg)[0, 1]
    r_fd = np.corrcoef(flex, deg)[0, 1]
    print("  r(TC, flexibility)        = %+.3f" % r_tf)
    print("  r(TC, contact degree)     = %+.3f" % r_td)
    print("  r(flexibility, degree)    = %+.3f" % r_fd)
    print()
    pc = partial(tc, flex, deg)
    print("  PARTIAL r(TC, flexibility | contact degree) = %+.3f" % pc)
    print()
    n = len(tc)
    # Fisher z standard error for a partial correlation with one controlled var
    se = 1.0 / np.sqrt(max(n - 3, 1))
    z = np.arctanh(np.clip(pc, -0.999999, 0.999999))
    print("  Fisher z = %+.2f, SE = %.4f, z/SE = %+.1f" % (z, se, z / se))
    if abs(z / se) > 3:
        print("  => TC retains a substantial, statistically clear association")
        print("     with flexibility after controlling for contact density.")
    else:
        print("  => TC does NOT add much beyond contact density on this measure.")

    print()
    print("-" * 78)
    print("Per-protein partial correlation r(TC, flexibility | contact degree)")
    vals = np.array([r["r_partial"] for r in rows])
    print("  mean %+.3f   median %+.3f   positive %d/%d"
          % (vals.mean(), np.median(vals), int((vals > 0).sum()), len(vals)))
    print()
    worst = min(rows, key=lambda r: r["r_partial"])
    best = max(rows, key=lambda r: r["r_partial"])
    print("  weakest : %s  %+.3f" % (worst["pdb"], worst["r_partial"]))
    print("  strongest: %s  %+.3f" % (best["pdb"], best["r_partial"]))

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "results", "partial_correlation.csv")
    keys = list(rows[0].keys())
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(r[k] if isinstance(r[k], str) else "%.6f" % r[k]
                              for k in keys) + "\n")
    print()
    print("  wrote %s" % out)


def partial(x, y, z):
    rxy = np.corrcoef(x, y)[0, 1]
    rxz = np.corrcoef(x, z)[0, 1]
    ryz = np.corrcoef(y, z)[0, 1]
    denom = np.sqrt(max(1e-12, (1 - rxz ** 2) * (1 - ryz ** 2)))
    return float((rxy - rxz * ryz) / denom)


if __name__ == "__main__":
    main()
