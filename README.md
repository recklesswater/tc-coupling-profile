# tc-coupling-profile

A small indicator, computed from a protein structure alone, that says how
strongly the motions of a local group of residues are coupled.

For a window of residues `B`:

```
TC(B) = KL( N(0, R_BB) || N(0, diag(R_BB)) ) = -1/2 * ln det R_BB
```

This is the **total correlation** of that block: how much information (in
nats) is lost if you assume the residues inside it fluctuate independently.
Large TC means the residues move together and a diagonal (independent)
description of their motion is a poor one; small TC means they are nearly
independent.

`R` is the correlation matrix of residue fluctuations from the **Gaussian
Network Model**, which needs nothing but a contact map. So the whole thing
costs a fraction of a second and needs no simulation, no experimental data and
no training.

## Install and run

```bash
pip install -r requirements.txt

# per-residue TC profile for a structure
python experiments/plot_profile.py

# correlation against four per-residue properties over 28 proteins
python experiments/run_multiprotein.py
```

Two lines for your own structure:

```python
from src.tc_profile import profile_from_pdb

profile, corr, names = profile_from_pdb("1UBQ", window=7)
```

`profile[i]` is the TC of a 7-residue window centred on residue `i+1`.

## What it looks like

![profile](figures/fig2_profile_ubiquitin.png)

For ubiquitin the highest value sits at the C-terminal tail (residues 70-76),
which is the flexible segment. Below the profile is the solvent accessibility
of the same residues, plotted on a comparable scale, so the two can be
compared directly.

## What it correlates with

Across 28 single-chain proteins (40-160 residues, diverse folds) we correlated
the profile with four per-residue properties:

| property | mean r | same sign across proteins |
| --- | --- | --- |
| flexibility (GNM mean-square fluctuation) | **+0.563** | 28 / 28 positive |
| contact degree | -0.406 | 27 / 28 negative |
| burial (1 - normalised SASA) | -0.254 | 26 / 28 negative |
| hydrophobicity (Kyte-Doolittle) | -0.135 | 25 / 28 negative |

![multiprotein](figures/fig1_multiprotein.png)

The first two rows are the interesting part and are mildly counter-intuitive:
**TC is highest where the protein is soft, not where it is densely packed.**
A flexible segment swings as a unit, so its internal correlations are close to
one, while a residue in a tightly packed core is held by many neighbours at
once, fluctuates little, and correlates less with its immediate neighbours.

The negative correlation with hydrophobicity is consistent with the intuition
that polar, solvated surfaces couple more strongly, but it is weak and largely
explained by burial: hydrophobicity and burial themselves correlate at
`+0.524`, and the partial correlation between TC and hydrophobicity after
controlling for burial is `+0.116`.

## What this is not

It is not a predictor of function, binding sites, or allostery. It is a
restatement of a correlation structure that comes out of a coarse elastic
network model, expressed as a single number per residue. Whether that number
is *useful* for anything is a separate question, and the honest answer on the
evidence collected so far is: it is a reasonable way to rank regions by how
correlated their motion is, and nothing more has been established.

Concretely, the caveats:

* the underlying model is the GNM, a single-parameter elastic network. It
  captures contact topology and nothing else -- no side chains, no chemistry,
  no solvent, no sequence conservation;
* the validation set is 28 small single-chain proteins. Nothing here says the
  numbers transfer to large multi-domain proteins or to complexes;
* zero hits in the literature scan in `docs/NOTES.md` are not evidence of
  novelty, only that a phrase search did not find them;
* the window length is a free parameter. 7 is a reasonable default but the
  profile changes with it, so results should be reported for more than one
  value.

## Files

```
src/tc_profile.py            the indicator: GNM covariance and block TC
experiments/plot_profile.py  contact map, GNM correlation, TC profile
experiments/run_multiprotein.py   the 28-protein validation
results/                     CSV output
figures/                     figures used above
docs/NOTES.md                definition, a worked derivation, literature scan
```

## Note

A separate piece of work uses this indicator to decide where a non-diagonal
covariance block should be placed in a variational model of conformational
ensembles. That is not part of this repository.
