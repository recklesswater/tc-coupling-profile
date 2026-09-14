---
title: TC Coupling Profile
emoji: 🧬
colorFrom: green
colorTo: blue
sdk: gradio
sdk_version: 5.50.0
app_file: app.py
pinned: false
license: mit
short_description: Where does assuming independent residue motion cost the most?
---

Type a PDB ID and see, in about a second, which parts of a protein chain are
the ones where treating residues as independently moving is a poor
approximation.

The indicator is the **block total correlation** of the Gaussian Network Model
correlation matrix:

```
TC(B) = KL( N(0, R_BB) || N(0, diag(R_BB)) ) = -1/2 * ln det R_BB
```

computed over a sliding window of residues. Large TC means the residues in that
window move together and an independent (diagonal) description of their motion
costs the most; small TC means independence is a good approximation.

Everything is derived from the C-alpha trace. No simulation, no training, no
experimental data, no fitting.

**What it is for:** it is a prior on where to look, not an answer. It ranks
local windows by how much the independent-motion assumption costs there, which
is useful before spending effort on a normal-mode calculation, a covariance
estimate, or a coarse-grained simulation.

**What it is not:** it does not identify binding sites or active pockets. On
five holo structures, ligand-binding residues had *lower* TC than the chain
average, not higher -- binding sites tend to be the rigid part of a protein,
while TC marks the soft part.

Source code, the 28-protein validation and the caveats are on GitHub:
https://github.com/recklesswater/tc-coupling-profile
