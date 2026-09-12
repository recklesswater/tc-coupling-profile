# tc-coupling-profile

**Mean field is almost always good enough. This tells you where "almost" fails.**

Independence is the default assumption when modelling the flexibility of a
protein: each residue is treated as moving on its own. That assumption is
cheap, and almost everywhere it is fine. But "almost everywhere" is not
"everywhere", and the places where it breaks down are not random.

This repository computes a single number per residue -- **TC**, the block total
correlation -- measuring how much is lost by assuming that the residues in a
local window move independently:

```
TC(B) = KL( N(0, R_BB) || N(0, diag(R_BB)) ) = -1/2 * ln det R_BB
```

Large TC means the residues in that window are strongly coupled and an
independent description of them is a poor one. Small TC means independence is
a good description.

The input is a PDB file. There is no simulation, no training, no experimental
data and no fitting.

## What it is for

It is a **prior on where to look**, not an answer. If effort is going to be
spent modelling local flexibility -- a normal-mode calculation, a covariance
estimate, a simulation at coarse-grained resolution -- TC says which windows
are the ones where the independent-motion assumption costs the most, and gives
that ranking in under a second from the structure alone.

## Use it

```bash
pip install -r requirements.txt

python experiments/plot_profile.py          # profile for one structure
python experiments/run_multiprotein.py      # validation over 28 proteins
python experiments/check_ligand_sites.py    # does it point at binding sites?
```

```python
from src.tc_profile import profile_from_pdb

profile, corr, names = profile_from_pdb("1UBQ", window=7)
```

`profile[i]` is the TC of a 7-residue window centred on residue `i+1`.

## What it looks like

![profile](figures/fig2_profile_ubiquitin.png)

For ubiquitin the largest value sits at the C-terminal tail (residues 70-76),
the segment that is known to be flexible. Below it, the burial proxy for the
same residues on a comparable scale.

## What it correlates with

28 single-chain proteins, 40-160 residues, diverse folds, window 7. Each
protein is correlated independently; the table reports the mean over proteins.

| property | mean r | median r | same sign |
| --- | --- | --- | --- |
| flexibility (GNM mean-square fluctuation) | **+0.562** | +0.569 | 28 / 28 positive |
| contact degree | -0.405 | -0.434 | 27 / 28 negative |
| burial (C-alpha coordination number) | -0.469 | -0.529 | 27 / 28 negative |
| hydrophobicity (Kyte-Doolittle) | -0.135 | -0.115 | 25 / 28 negative |

![multiprotein](figures/fig1_multiprotein.png)

The sign pattern is the interesting part, and it is mildly counter-intuitive:

> **TC is highest where the protein is soft and exposed, not where it is
> densely packed.**

A flexible segment swings as a unit, so the motions of its residues are almost
perfectly correlated with each other. A residue in a tightly packed core is
held by many neighbours at once, barely moves, and correlates less with its
immediate neighbours. Contact density increases the *number* of constraints
without increasing the *correlation* between neighbouring displacements, which
is why the contact-degree row is negative.

The hydrophobicity row has the sign one might expect from the intuition that
polar surfaces couple more strongly, but the effect is weak and largely
explained by burial; do not lean on it.

## What it is **not**

**It does not identify binding sites or active pockets.** We tested this
directly, because it is the natural thing to hope for: take structures with a
bound ligand, find the residues within 5 A of the ligand, and compare their TC
against the chain average. The result is the opposite of the hope.

| structure | ligand-site residues | site TC | chain mean TC | z |
| --- | --- | --- | --- | --- |
| 3PTB (trypsin) | 8 | 0.301 | 0.527 | -0.73 |
| 4DFR (DHFR) | 8 | 0.473 | 0.640 | -0.68 |
| 1STP (streptavidin) | 5 | 0.546 | 0.578 | -0.08 |
| 3ERT (oestrogen receptor) | 8 | 0.883 | 1.001 | -0.28 |
| 1M17 (EGFR kinase) | 14 | 0.856 | 1.248 | -0.34 |

Mean z = -0.42, and 0 of 5 structures show enrichment. Ligand-binding sites sit
at *lower* TC than the chain average, consistent with the main result: binding
sites tend to be the rigid part of a protein, while TC marks the soft part. So
TC is, if anything, mildly **anti**-correlated with binding sites. The honest
statement of the use case is "where independence costs the most", not "where
the active site is".

This test uses a small set and only five structures produced usable ligand
geometry, so treat it as a caution rather than a result. It is included because
a negative check is more useful than an unexamined claim.

Other caveats:

* the underlying model is the GNM, a single-parameter elastic network. It
  captures contact topology and nothing else -- no side chains, no chemistry,
  no solvent, no sequence conservation;
* the validation set is 28 small single-chain proteins; nothing here says the
  numbers transfer to large multi-domain proteins or complexes;
* window length is a free parameter. 7 is a reasonable default, but TC grows
  with window size, so profiles are comparable in *shape* across window
  lengths, not in magnitude;
* zero hits in the literature scan in `docs/NOTES.md` are not evidence of
  novelty, only that a phrase search did not find them.

## Files

```
src/tc_profile.py            the indicator: GNM covariance and block TC
src/pdb_io.py                minimal dependency-free PDB reader
experiments/plot_profile.py  contact map, GNM correlation, TC profile
experiments/run_multiprotein.py     the 28-protein validation
experiments/check_ligand_sites.py   the binding-site check above
results/                     CSV output
figures/                     figures used above
docs/NOTES.md                derivation, worked example, literature scan
```

## Note

A separate piece of work uses this indicator to decide where a non-diagonal
covariance block should be placed in a variational model of conformational
ensembles. That is not part of this repository.

## Attribution and AI assistance

Every mathematical ingredient used here is published work and is cited in
[`docs/NOTES.md`](docs/NOTES.md): the Gaussian Network Model, total
correlation, the modified Cholesky decomposition and the Kyte-Doolittle scale.
What this repository claims is only the specific formulation and its
validation, not the underlying methods -- and a phrase search is not a
systematic review, so please read the prior-art note before treating anything
here as new.

The code was written with AI coding assistance under the author's direction.
The research question, the design decisions and the interpretation of results
are the author's, and the author is responsible for their correctness.
