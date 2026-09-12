# Notes on the TC coupling indicator

## 1. Definition

For a window of residues `B`:

```
TC(B) = KL( N(0, R_BB) || N(0, diag(R_BB)) ) = -1/2 ln det R_BB
```

`R_BB` is the block of the GNM correlation matrix covering those residues.
The equality uses the standard Gaussian KL:

```
KL( N(0, S1) || N(0, S2) ) = 1/2 [ tr(S2^-1 S1) - d + ln det S2 - ln det S1 ]
```

Put `S1 = R_BB` (unit diagonal) and `S2 = diag(R_BB) = I` for a correlation
matrix. Then `tr(S2^-1 S1) = d` and `ln det S2 = 0`, leaving
`-1/2 ln det R_BB`. TC is non-negative because `det R_BB <= 1` for a
correlation matrix, with equality exactly when the block is independent.

This quantity is the **total correlation** (multi-information) of the block:
the reduction in total entropy obtained by modelling the dependencies rather
than treating the variables as independent.

## 2. Where `R` comes from

The Gaussian Network Model treats the protein as a network of identical
springs,

```
Gamma_ij = -1        if residues i and j are in contact (C-alpha < 8 A)
Gamma_ii = degree_i
```

The mean-square fluctuation of residue `i` is proportional to `(Gamma+)_ii`,
and the covariance of the fluctuations is `Gamma+`, the pseudo-inverse (the
true inverse does not exist because `Gamma` has a zero eigenvalue
corresponding to overall translation). Normalising that covariance gives the
correlation matrix `R` with unit diagonal.

Because `R` has unit diagonal, TC measures *only* the correlation structure --
the per-residue fluctuation magnitudes are divided out. This is why the strong
correlation between TC and flexibility (section 4) is not a tautology: they
are distinct functions of the same eigen-decomposition, and the relationship
between them is an empirical finding rather than an algebraic identity.

## 3. A worked example

Ubiquitin (1UBQ, 76 residues, 326 contacts, mean degree 8.6), window 7:

```
TC range 0.216 - 1.074 nats, mean 0.586
highest TC at residue 73 (the C-terminal tail)
lowest  TC at residue 67
```

For comparison, the whole 76-residue chain has a total correlation far larger
than any window, which is the expected behaviour: TC grows with block size
because more pairs can be dependent.

## 4. Validation

28 single-chain proteins, 40-160 residues, window 7. Per protein we correlated
the profile against four properties and recorded the Pearson coefficient.

| property | mean r | median r | proteins with the same sign |
| --- | --- | --- | --- |
| flexibility (GNM variance) | +0.563 | +0.569 | 28 / 28 positive |
| contact degree | -0.405 | -0.434 | 27 / 28 negative |
| burial (C-alpha coordination number) | -0.469 | -0.529 | 27 / 28 negative |
| hydrophobicity (Kyte-Doolittle) | -0.135 | -0.115 | 25 / 28 negative |

Each row is a mean over 28 independent proteins, not a pooled regression, so a
consistent sign across nearly all of them is meaningful even though the
magnitudes are modest.

The sign pattern is worth stating plainly: **coupling demand is highest in
soft, solvent-exposed regions and lowest in the densely packed core.** The
reading is that a flexible segment displaces as a quasi-rigid unit (near-unit
internal correlation), whereas a core residue is constrained from many
directions at once and therefore fluctuates little and correlates less with
its immediate neighbours. The negative correlation with contact degree
(-0.406) is the same statement from the other side, and is initially
surprising: contact density increases the *number* of constraints without
increasing the *correlation* between neighbouring displacements.

The hydrophobicity row is the one to treat most cautiously. Its sign matches
the intuition that polar surfaces couple more, but the effect is weak, and
hydrophobicity and burial are themselves correlated at +0.524. The partial
correlation between TC and hydrophobicity controlling for burial is +0.116,
i.e. the sign flips, so there is no evidence here for an independent
hydrophobicity effect.

**On averaging.** Pearson `r` has a skewed sampling distribution, so averaging
`r` values across proteins is not the correct procedure; the standard approach
is to average `artanh(r)` and transform back. Both are reported in the README.
The correction is small here (+0.021 for the flexibility row against a mean of
0.562) because the individual `r` values are moderate, but the arithmetic mean
should not have been quoted alone.

**On the burial measure.** An earlier version used a real solvent-accessible
surface area and obtained -0.254; the C-alpha coordination number used here
gives -0.469. Both are strongly negative and the sign is identical in 26-27 of
28 proteins, so the conclusion does not depend on the choice, but the magnitude
does. A C-alpha coordination number is crude; it is used because a real SASA
requires a compiled Biopython extension that some Windows application-control
policies block, and a tool that fails to import is worse than a tool with a
cruder burial proxy.

## 5. Does TC mark binding sites? (No.)

The obvious hope is that TC, being computable from structure alone, could stand
in for experimental identification of active or binding sites. We tested it.

For structures with a bound ligand we found every residue whose C-alpha lies
within 5 A of a ligand heavy atom, and compared the mean TC at those residues
with the mean over the whole chain, normalised by the profile's standard
deviation.

| structure | ligand-site residues | site TC | chain mean TC | z |
| --- | --- | --- | --- | --- |
| 3PTB (trypsin + benzamidine) | 8 | 0.301 | 0.527 | -0.73 |
| 4DFR (DHFR + methotrexate) | 8 | 0.473 | 0.640 | -0.68 |
| 1STP (streptavidin + biotin) | 5 | 0.546 | 0.578 | -0.08 |
| 3ERT (oestrogen receptor + tamoxifen) | 8 | 0.883 | 1.001 | -0.28 |
| 1M17 (EGFR kinase + erlotinib) | 14 | 0.856 | 1.248 | -0.34 |

Mean z = -0.42; median -0.34; **0 of 5 show enrichment**.

The sign is consistent with the main finding rather than with the hope.
Binding sites are typically the rigid part of a protein -- they have to hold a
ligand -- whereas TC is highest in soft, solvent-exposed regions. So TC is
mildly *anti*-correlated with binding sites, and the tool should not be
presented as a binding-site predictor.

Limitations of this test: only five structures produced usable geometry, a
ligand-contacting residue is not the same as a catalytic residue, and the
cutoff is arbitrary. It is a caution, not a result. It is reported anyway
because a stated negative is worth more than an unexamined claim.

## 5b. Is TC just contact density? (No.)

This is the objection that matters most, because TC and contact degree are both
functions of the same contact graph. A correlation between them is guaranteed
and carries no information.

Within each protein we z-score TC, contact degree and flexibility, pool all
2412 residues from the 28 structures, and compute the partial correlation of TC
with flexibility controlling for contact degree.

| quantity | value |
| --- | --- |
| r(TC, flexibility) | +0.581 |
| r(TC, contact degree) | -0.391 |
| r(flexibility, contact degree) | -0.809 |
| partial r(TC, flexibility \| contact degree) | **+0.490** |

Contact density accounts for most of the flexibility signal, as expected
(`r = -0.81`), but TC keeps a large association on top of it. The Fisher-z of
the partial correlation is 0.54 with a standard error of 0.020, so the ratio is
about 26 at n = 2412. Per protein the partial correlation is positive in 27 of
28 cases, mean +0.461.

The same test was run at the level of the downstream experiment, which is the
sharper version because the thing being predicted is the actual benefit of
using a structured covariance rather than a proxy for it. Across 30 blocks from
three proteins, `r(gain, TC | contact degree) = +0.797` (`z/SE = +5.7`). That
experiment belongs to separate work and is not part of this repository.

One caveat: the partial correlation is computed on a pooled sample of residues,
and residues within a protein are not independent observations. The effective
sample size is smaller than 2412, so the `z/SE` figure should be read as an
order-of-magnitude indicator rather than a calibrated p-value. The per-protein
result (27/28 positive) is the more robust statement.

## 6. Literature scan

OpenAlex, title + abstract field, September 2026.

| query | hits |
| --- | --- |
| modified Cholesky + covariance | 178 |
| rate distortion + variational autoencoder | 130 |
| elastic network model + machine learning | 36 |
| Gaussian network model + machine learning | 11 |
| protein + latent space + covariance | 6 |
| protein + variational autoencoder + covariance | 4 |
| Gaussian network model + deep learning | 2 |
| anisotropic network model + neural network | 1 |
| protein + normalizing flow + covariance | 0 |

The ingredients of the indicator are all standard: the Gaussian Network Model
dates to the late 1990s and total correlation is a classical information-theoretic
quantity. The scan above is included for context rather than as a claim of
novelty. Note also that OpenAlex searches titles and abstracts only, so a zero
means the phrase does not appear there -- not that no related work exists.

## 7. What would falsify the indicator's usefulness

The profile is only interesting if it says something the simpler quantities do
not. Three concrete tests would settle that:

1. correlate TC against experimental measures of local flexibility (B-factors,
   NMR order parameters, hydrogen-deuterium exchange) rather than against the
   GNM, which is derived from the same contact map;
2. check whether TC agrees with other, independently constructed coupling
   measures, for example from mutual information on an MD trajectory or from
   dynamical cross-correlation matrices;
3. measure whether TC compresses to something simpler -- if a two-parameter
   function of local contact density and window position reproduces it, the
   indicator adds little.

None of these have been done here.

## 8. Practical notes

* Window length matters. TC increases monotonically with window size, so
  profiles computed at different window lengths are not directly comparable in
  magnitude -- only in shape.
* Residues with no contacts inside a window can give a near-singular block; the
  implementation adds a small ridge (`1e-8`) before the log-determinant and
  returns `nan` if the block is still not positive definite.
* Only the first chain of a PDB entry is used. For complexes this silently
  discards the rest, which is the right default for a single-chain indicator
  but should be made explicit when adapting the code.

## 9. Methods, attribution and AI assistance

### What is standard, and must be cited

Everything mathematical used here is published work, none of it is original to
this repository, and all of it should be cited in any write-up.

| ingredient | source | obligation |
| --- | --- | --- |
| Gaussian Network Model | Tirion (1996); Bahar, Atilgan & Erman (1997) | cite |
| total correlation / multi-information | Watanabe (1960); standard information theory | cite |
| modified Cholesky decomposition | standard statistical method for sparse covariance estimation | cite |
| KL divergence between Gaussians | textbook result | cite or state |
| C-alpha contact definition (8 A) | standard in the elastic-network literature | cite |
| Kyte-Doolittle hydropathy scale | Kyte & Doolittle (1982) | cite |
| PDB coordinate data | RCSB PDB, distributed under CC0 1.0 | free to use |

Methods, formulas and ideas are not copyrightable, so using them creates an
obligation to attribute rather than an obligation to obtain permission. No
licence from any of these authors is required.

### What was not taken

No code was copied from any other project. Every implementation here -- the
masked-Cholesky covariance assembly, the whitened KL decomposition, the PDB
reader, the plotting, the demo -- was written from the mathematical
descriptions and the public PDB format specification. The numeric constants
that appear are published facts.

### What this repository claims

Only the following, and nothing more:

1. the specific formulation of the block total correlation of the GNM
   correlation matrix as a per-residue structural indicator;
2. the multi-protein validation of that indicator against flexibility,
   contact degree, burial and hydrophobicity;
3. the negative check against ligand-binding sites.

It does **not** claim the underlying methods, the general idea of measuring
inter-residue coupling with information theory, or the Gaussian Network Model.
The general area -- information-theoretic measures of residue coupling -- is
well developed, with several hundred papers on dynamical cross-correlation and
allosteric communication alone.

### Prior art

The queries behind the counts quoted in section 6 are in
`experiments/prior_art_scan.py`. A phrase search over titles and abstracts
finds zero hits for the exact combination used here, but a phrase search is
not a systematic review and **should not be read as evidence of novelty**. The
closest known neighbours are work on mutual information between residues
(2024), on synergy and redundancy in allosteric coupling (2025), and on linear
correlations recovered from normal modes (2018); these should be read before
any claim about what is new.

### AI assistance

The code in this repository was written with AI coding assistance under the
author's direction. The research question, the design decisions and the
interpretation of results are the author's, and the author is responsible for
their correctness. This is stated rather than omitted because a reader is
entitled to know how a piece of work was produced.
