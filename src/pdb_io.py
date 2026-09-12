"""Minimal PDB reader -- only what this project needs.

Deliberately dependency-free. The published version of this code used
Biopython, but Biopython pulls in a compiled extension that some Windows
application-control policies block outright, which makes the tool fail for
reasons that have nothing to do with the science. Everything required here is
a handful of fixed-column slices.

Only the fields we use are parsed: record type, residue name, chain id,
residue number, atom name, element and coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

WATER = {"HOH", "DOD", "WAT"}


@dataclass
class Atom:
    record: str
    resname: str
    chain: str
    resseq: int
    icode: str
    name: str
    element: str
    coord: np.ndarray


def read_atoms(path: str):
    """Parse ATOM/HETATM records from a PDB file.

    Fixed-column format, per the PDB specification. Lines that are too short
    or malformed are skipped rather than raising, because real PDB files
    routinely contain both.
    """
    atoms = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            record = line[0:6].strip()
            # Multi-model files (NMR ensembles) repeat the whole coordinate
            # block; only the first model is wanted, otherwise residue counts
            # are multiplied by the number of models.
            if record == "ENDMDL":
                break
            if record not in ("ATOM", "HETATM"):
                continue
            try:
                coord = np.array([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ])
            except ValueError:
                continue
            element = line[76:78].strip().upper()
            name = line[12:16].strip()
            if not element:  # older files omit the element column
                element = name[0].upper() if name else ""
            atoms.append(Atom(
                record=record,
                resname=line[17:20].strip().upper(),
                chain=line[21:22].strip(),
                resseq=int(line[22:26]) if line[22:26].strip().lstrip("-").isdigit() else 0,
                icode=line[26:27].strip(),
                name=name,
                element=element,
                coord=coord,
            ))
    return atoms


def first_chain(atoms):
    """Keep the first chain that contains at least one ATOM record."""
    for atom in atoms:
        if atom.record == "ATOM" and atom.chain:
            chain_id = atom.chain
            break
    else:
        return []
    return [a for a in atoms if a.chain == chain_id]


def ca_trace(path: str):
    """Return ``(coords, resnames, resseqs)`` for the first chain's C-alphas."""
    chain = [
        a for a in first_chain(read_atoms(path))
        if a.name == "CA" and a.resname not in WATER
    ]
    coords = np.asarray([a.coord for a in chain], dtype=float)
    names = [a.resname for a in chain]
    resseqs = [a.resseq for a in chain]
    return coords, names, resseqs


def ligand_atoms(path: str, skip_names=None):
    """Heavy atoms of every non-water HETATM group in the first chain."""
    skip = {"HOH", "DOD"} | set(skip_names or ())
    out = []
    for a in first_chain(read_atoms(path)):
        if a.record != "HETATM" or a.resname in skip:
            continue
        if a.element == "H":
            continue
        out.append(a.coord)
    return np.asarray(out, dtype=float).reshape(-1, 3)


def burial_proxy(coords: np.ndarray, cutoff: float = 10.0) -> np.ndarray:
    """Per-residue burial proxy: number of C-alphas within ``cutoff``.

    A solvent-exposed residue has few C-alpha neighbours; a buried one has
    many. This is cruder than a real solvent-accessible surface area, but it
    needs nothing beyond the C-alpha trace and it is monotone in burial, which
    is all the correlation analysis requires.
    """
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    counts = (d < cutoff).sum(axis=1) - 1  # exclude self
    return counts.astype(float)
