# Per-position contact profiles — batch 1 (2026-11-17)

Derived analysis of the batch-1 BoltzGen designs: which target residues the
designed binders actually touch, and how often.

## Input

`refold_cif/<target>_<id>.cif` from each batch-1 run (Boltz-2 refolded complex;
chain A = designed binder, chain B = target peptide). **15 complexes per target**
out of the 60 available, single seed. Raw structures are not in this repository.

## Definition

A target residue counts as *in contact* in a given design when any of its heavy
atoms lies within the cutoff of any binder heavy atom. Contact frequency of a
position = fraction of designs in which it is in contact. Both 4.5 Å (close
contact) and 5.0 Å (the usual interface definition) are reported; the two agree.
Heavy atoms rather than CB are used because long side chains (Arg, Lys) make up a
large part of these peptides.

Script: `analysis/contact_profile.py` (`python contact_profile.py <cif_dir> <out_dir>`,
files grouped by the filename prefix before the first `_`).

## Result

| | bnp | anp | scr (shuffled control) |
|---|---|---|---|
| designs | 15 | 15 | 15 |
| target length | 32 | 28 | 32 |
| mean contact freq (5.0 Å) | 0.69 | 0.63 | 0.65 |
| SD across positions | 0.295 | 0.271 | 0.193 |
| lowest position | 0.07 (Ser8) | 0.13 (Ser1 / Arg3) | 0.33 |
| residues contacted per design | 22.0 / 32 | 17.6 / 28 | 20.9 / 32 |

The profiles are **structured, not flat, and not wrap-around**: on average a
design touches about two thirds of the target, and some positions are contacted
in 1–2 designs out of 15. Both real targets concentrate contacts on the
disulfide ring and the C-terminal tail and avoid the segment preceding the first
Cys (BNP 5–9, ANP 1–6); the shuffled control has no such avoided region and is
the flattest of the three.

Aligning the two rings Cys-to-Cys (BNP 10–26 vs ANP 7–23, 17 positions,
11/17 identical residues), the two contact profiles correlate at **r = 0.56**;
bnp vs scr at identical numbering gives r = 0.07. Suggestive, not conclusive
(17 autocorrelated points, one seed): the pipeline appears to pick the *same*,
conserved epitope on both homologues, which would explain why its confidence
scores do not separate them.

## Caveats

- n = 15 per target, single seed; per-position standard error ≤ 0.13, so only
  the overall shape, not individual positions, should be read.
- Contacts are computed on Boltz-2 refolded model structures, not experimental
  ones. The "epitope" is the model's, with no wet-lab support.

## Files

- `contact_profile.png` — three-panel figure, disulfide ring shaded
- `{bnp,anp,scr}_perposition.csv` — position, residue, contact frequency at both
  cutoffs with standard errors, mean and SD of the minimum heavy-atom distance
- `contact_profile_summary.csv` — per-target summary
