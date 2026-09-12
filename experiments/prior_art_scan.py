"""Prior-art check for the specific claim: block total correlation of the GNM
correlation matrix, used as a per-residue structural indicator.

Phrase-level search over OpenAlex titles + abstracts. This is not a systematic
review -- it is a first pass intended to catch obvious collisions.
"""

from __future__ import annotations

import json
import ssl
import sys
import time
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
PROXY = "http://127.0.0.1:7890"
OPENER = urllib.request.build_opener(
    urllib.request.ProxyHandler({"http": PROXY, "https": PROXY}),
    urllib.request.HTTPSHandler(context=CTX),
)

QUERIES = [
    ("GNM + mutual information",
     '"Gaussian network model" AND "mutual information"'),
    ("GNM + total correlation",
     '"Gaussian network model" AND "total correlation"'),
    ("elastic network + mutual information",
     '"elastic network model" AND "mutual information"'),
    ("elastic network + total correlation",
     '"elastic network model" AND "total correlation"'),
    ("protein + total correlation + residues",
     'protein AND "total correlation" AND residue'),
    ("protein + multi-information",
     'protein AND "multi-information"'),
    ("residue coupling + covariance + structure",
     '"residue coupling" AND covariance AND structure'),
    ("dynamical cross-correlation + protein",
     '"dynamical cross-correlation" AND protein'),
    ("local coupling + protein + elastic network",
     '"local coupling" AND protein AND "elastic network"'),
    ("protein + independence assumption + fluctuations",
     'protein AND "independence assumption" AND fluctuation'),
    ("block covariance + protein structure",
     '"block covariance" AND "protein structure"'),
    ("entropy + residue + normal mode",
     'entropy AND residue AND "normal mode" AND protein'),
    ("protein + information theory + flexibility",
     'protein AND "information theory" AND flexibility'),
    ("allostery + mutual information + structure",
     'allostery AND "mutual information" AND structure'),
    ("protein + correlation matrix + sparsity",
     'protein AND "correlation matrix" AND sparse'),
]


def count(expr):
    url = ("https://api.openalex.org/works?filter=title_and_abstract.search:%s"
           "&per_page=1&mailto=research@example.org" % urllib.parse.quote(expr))
    req = urllib.request.Request(url, headers={"User-Agent": "prior-art"})
    with OPENER.open(req, timeout=45) as r:
        return json.load(r).get("meta", {}).get("count", -1)


def top_titles(expr, n=3):
    url = ("https://api.openalex.org/works?filter=title_and_abstract.search:%s"
           "&per_page=%d&sort=cited_by_count:desc&mailto=research@example.org"
           % (urllib.parse.quote(expr), n))
    req = urllib.request.Request(url, headers={"User-Agent": "prior-art"})
    with OPENER.open(req, timeout=45) as r:
        data = json.load(r)
    out = []
    for w in data.get("results", []):
        year = w.get("publication_year", "?")
        out.append("      %s (%s)" % (w.get("title", "?"), year))
    return out


def main():
    print("Prior-art scan: OpenAlex, title + abstract field")
    print("=" * 78)
    hits = []
    for label, expr in QUERIES:
        try:
            n = count(expr)
        except Exception as exc:  # noqa: BLE001
            n = "error:%s" % repr(exc)[:30]
        hits.append((label, expr, n))
        print("  %-42s %s" % (label, n))
        time.sleep(1.2)

    print()
    print("=" * 78)
    print("The closest queries, and their most-cited hits:")
    for label, expr, n in hits:
        if isinstance(n, int) and 0 < n <= 40:
            print()
            print("  [%d hits] %s" % (n, label))
            try:
                for t in top_titles(expr):
                    print(t)
            except Exception:  # noqa: BLE001
                pass
            time.sleep(1.2)


if __name__ == "__main__":
    main()
