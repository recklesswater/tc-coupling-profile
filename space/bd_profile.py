"""BD coupling profile: a structure-derived indicator of local coupling in proteins.

The indicator is the **block total correlation** of a window of residues,

    BD(B) = KL( N(0, R_BB) || N(0, diag(R_BB)) ) = -1/2 * ln det R_BB

i.e. the amount of information, in nats, that is lost if the motions of the
residues inside that window are assumed to be independent. Large BD means the
residues in that window move in a correlated way and are poorly described by
independent (diagonal) fluctuations; small BD means they are nearly
independent.

`R` is the correlation matrix of residue fluctuations under the **Gaussian
Network Model** (GNM), a coarse elastic-network model derived from the contact
map. Nothing here needs a simulation or experimental data -- only a structure.

What the indicator is good for
------------------------------
It gives a per-residue scalar that can be computed for any structure and used
to rank regions by how correlated their motion is. On 28 single-chain proteins
it correlates with the GNM flexibility profile (+0.56) and negatively with
contact density (-0.41); see ``docs/NOTES.md`` for the full validation and for
what it does *not* mean.

Reference implementation only -- see ``docs/NOTES.md`` for caveats.
"""

from __future__ import annotations

import os
import urllib.request

import numpy as np

from pdb_io import burial_proxy, ca_trace, ligand_atoms

__all__ = [
    "fetch_pdb",
    "ca_trace",
    "burial_proxy",
    "ligand_atoms",
    "contact_map",
    "gnm_correlation",
    "block_dependence",
    "bd_profile",
]

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


def fetch_pdb(pdb_id: str, cache_dir: str = "data/pdb") -> str:
    """Download a PDB entry into a local cache and return the file path."""
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{pdb_id}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    req = urllib.request.Request(
        RCSB_URL.format(pdb_id=pdb_id), headers={"User-Agent": "tc-profile"}
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = resp.read()
    with open(path, "wb") as fh:
        fh.write(data)
    return path


def contact_map(ca: np.ndarray, cutoff: float = 8.0) -> np.ndarray:
    """GNM contact map: C-alpha pairs closer than ``cutoff`` angstrom."""
    d = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=-1)
    out = (d < cutoff).astype(float)
    np.fill_diagonal(out, 0.0)
    return out


def gnm_correlation(contacts: np.ndarray) -> np.ndarray:
    """Correlation matrix of residue fluctuations under the GNM.

    ``Gamma`` is the Kirchhoff matrix; it has one zero eigenvalue
    (global translation), hence the pseudo-inverse.
    """
    degree = contacts.sum(axis=1)
    gamma = np.diag(degree) - contacts
    cov = np.linalg.pinv(gamma)
    sd = np.sqrt(np.clip(np.diag(cov), 1e-12, None))
    corr = cov / np.outer(sd, sd)
    np.fill_diagonal(corr, 1.0)
    return corr


def block_dependence(corr: np.ndarray, idx) -> float:
    """BD(B) = -1/2 ln det R_BB, in nats.

    Equals the KL divergence between the true zero-mean Gaussian on that block
    and its diagonal approximation, so it is the exact information lost by
    treating the block as independent.
    """
    sub = corr[np.ix_(list(idx), list(idx))]
    sub = 0.5 * (sub + sub.T)
    sub = sub + np.eye(sub.shape[0]) * 1e-8  # ridge for near-singular blocks
    sign, logdet = np.linalg.slogdet(sub)
    if sign <= 0:
        return float("nan")
    return float(-0.5 * logdet)


def bd_profile(corr: np.ndarray, window: int = 7) -> np.ndarray:
    """Per-residue BD, using a centred window of ``window`` residues."""
    n = corr.shape[0]
    half = window // 2
    prof = np.full(n, np.nan)
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        prof[i] = block_dependence(corr, range(lo, hi))
    return prof


def profile_from_pdb(pdb_id: str, window: int = 7, cutoff: float = 8.0,
                     cache_dir: str = "data/pdb"):
    """Convenience wrapper: PDB id -> (BD profile, correlation matrix, names)."""
    path = fetch_pdb(pdb_id, cache_dir)
    ca, names, _ = ca_trace(path)
    corr = gnm_correlation(contact_map(ca, cutoff))
    return bd_profile(corr, window), corr, names
