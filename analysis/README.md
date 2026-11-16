# analysis/

## first_batch_stats.py

Descriptive comparison of the three first-batch targets (`bnp`, `anp`, and the
shuffled-sequence negative control `scr`), 60 designs each, BoltzGen 0.3.2 on a
Colab L4, single seed.

```
python first_batch_stats.py results analysis
```

It expects `<results_dir>/{bnp,anp,scr}_all_designs_metrics.csv` (the raw
BoltzGen metric tables) and writes `first_batch_descriptives.csv` and
`first_batch_distributions.png` into `<out_dir>`.

### Rules this script encodes (do not change without a reason)

- Ranking score is `design_to_target_iptm`, with `min_design_to_target_pae`
  (lower is better) as the secondary. `design_to_target_ipsae` is reported to
  replicate Junker & Schoeder.
- `quality_score` and every `rank_*` column are **within-batch rank-normalised**
  and must never be compared across targets.
- `design_ptm` is not used: its batch-to-batch drift (0.027) exceeds its own
  per-design SD (0.024).
- All of these are model **confidence, not affinity**. Boltz-2's affinity module
  is small-molecule only.
- Noise floor, measured on the same target with a new seed (n=20 per batch):
  batch-mean iptm difference 0.022, per-design SD 0.107, SE of the mean ~0.024.
  Means therefore reproduce only to about **+/-0.05**, and any between-target
  difference has to clearly exceed that to mean anything.

Findings are written up in `docs/first_batch_analysis_2026-11-16.md`.
