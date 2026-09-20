# Validation plan: a BD-based druggability filter

Status: proposed, not yet run. This document specifies the benchmark that would turn the
filter described in the README from a hypothesis derived from the correlation statistics into
a measured enrichment claim. Until it is run, the README states the criterion as proposed.

## 1. What is being tested

**Hypothesis.** Among surface-exposed pockets of comparable geometry, those with **low BD**
(rigid, well-coupled neighbourhood) are more likely to coincide with a true ligand-binding
site than those with **high BD** (soft, flexible neighbourhood).

Two specific claims have to be separated, because they can fail independently:

1. **BD adds signal.** Within a stratum of matched SASA and volume, pockets containing the
   true ligand have lower BD than pockets that do not.
2. **BD adds signal beyond SASA.** The BD term improves ranking over a geometry-only score
   (SASA, volume, depth) on the same pocket set.

Claim 1 is the mechanism; claim 2 is the product value. A result of "claim 1 true, claim 2
false" means BD explains part of pocket quality but is redundant with geometry, and the README
should say so.

## 2. Data

| source | content | usable structures | notes |
| --- | --- | --- | --- |
| COACH420 / HOLO4K | holo structures with a curated primary ligand | ~400 / ~4000 | single-chain, moderate resolution; the standard pocket-detection benchmark |
| sc-PDB | drug-like ligands in PDB with UniProt mapping | ~3000 | filter to one ligand per chain, no covalent binders |
| PDBbind refined set | protein-ligand complexes with affinity | ~5000 | use if a link to measured affinity is wanted later |

Start with **COACH420** (about 400 structures) as the primary set, because it is small enough
to iterate on and large enough for a first enrichment estimate. Exclusions: multi-chain
interfaces, structures with ligands shorter than 8 heavy atoms, metal-only or buffer
components, and any chain below 40 residues (the BD profile needs a window).

Sample size: with 400 structures and roughly 3-8 pockets per structure, the pocket-level
analysis has on the order of 2000 observations, which is adequate to detect an AUC difference
of about 0.05 with acceptable power.

## 3. Pocket definition

Geometry comes from an existing detector, not from BD. Two options, in order of preference:

1. **fpocket** (open source, scriptable): `fpocket -f structure.pdb`, keep pockets with a
   score above the default cut-off. Gives pocket residues, volume and an alpha-sphere
   description.
2. **DoGSiteScorer / SiteMap** if a commercial licence is available: gives surface area and
   volume per pocket directly, which is closer to the "high SASA" half of the criterion.

Each pocket is reduced to a residue list (pocket-lining residues within 5 Å of any pocket
alpha sphere or of the ligand, depending on which tool is used).

The **reference site** is the set of residues within 5 Å of the primary ligand in the same
structure. A pocket counts as a hit if its residue set overlaps the reference site by at
least 50% of the reference residues (this threshold should be reported; a stricter 80% and a
looser 30% version are worth adding as sensitivity analyses).

## 4. Features per pocket

| feature | source | role |
| --- | --- | --- |
| mean BD over pocket residues | BD profile, window 7 | the new term |
| max BD over pocket residues | BD profile | captures a locally soft patch inside an otherwise rigid pocket |
| pocket SASA / volume / depth | fpocket or DoGSiteScorer | the geometry baseline |
| chain-mean BD | BD profile | used to normalise BD within each structure |

Because BD magnitudes are not comparable across proteins (it grows with window length and
with chain length), every BD feature must be **z-scored within its own protein** before being
pooled, exactly as in the contact-density analysis in `docs/NOTES.md`.

## 5. Analysis

**Step 1 — stratified comparison (claim 1).** Bin pockets by SASA into terciles. Within each
tercile, compare the z-scored mean BD of hit versus non-hit pockets (Mann-Whitney, effect
size). The hypothesis predicts hits have systematically lower BD.

**Step 2 — incremental value (claim 2).** Fit a logistic model per pocket:

```
hit ~ SASA + volume                     (geometry only)
hit ~ SASA + volume + bd_z              (geometry + BD)
```

Report the AUC of each model with a grouped cross-validation split **by protein** (never by
pocket, otherwise the same structure leaks across folds) and compare with a paired test on
the out-of-fold predictions. The deliverable is one number: the AUC gained by adding BD.

**Step 3 — triage curve.** Rank pockets by the geometry score, then by the combined score,
and plot precision (fraction of hits) at the top 1, 2, 5 and 10% of the ranking. This is the
form the result takes in a screening context, and it is the figure to put in the README.

**Step 4 — negative control.** Repeat step 2 with a shuffled BD profile (permute BD values
within each protein) to confirm that any gain is not an artefact of the pipeline.

## 6. Pre-registered outcome interpretation

| outcome | reading | action |
| --- | --- | --- |
| hits show lower BD within SASA terciles, and adding BD raises grouped AUC by ≥ 0.03 with a non-significant shuffled control | filter validated | move the module in the README from "proposed" to "validated", add the triage figure |
| BD difference present but AUC gain < 0.03 | mechanism real, product value not demonstrated | keep the filter as a secondary annotation, state the redundancy with geometry |
| no BD difference within SASA terciles | hypothesis fails | remove the druggability-filter section, keep the "softness indicator" framing of the original README |

The third row is a realistic outcome and should be written down before the run, not after.

## 7. Cost

Structure acquisition and pocket detection dominate: fpocket runs in a few seconds per
structure, BD in under a second. The whole COACH420 benchmark is a few hours of compute and no
GPU. The analysis scripts can follow the layout of `experiments/check_ligand_sites.py`, which
already implements the ligand-site comparison on five structures.

## 8. Known confounds to check before publishing

* **Ligand size.** Larger ligands contact more residues and mechanically raise the number of
  rigid contacts. Check that the BD effect survives a ligand-heavy-atom-count covariate.
* **Flexibility of the holo structure itself.** Some holo structures are of the same protein
  in a different conformational state than the apo form. Report which PDB entries differ.
* **Resolution.** Low-resolution structures have noisier C-alpha positions; check that BD
  differences are not concentrated in the low-resolution tail of the set.
* **Pocket count per structure.** Structures with many detected pockets contribute more
  observations; verify that grouped cross-validation by protein keeps the weighting honest.
