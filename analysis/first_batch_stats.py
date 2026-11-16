"""First-batch descriptive comparison of BoltzGen confidence metrics across
three targets (bnp, anp, scr=shuffled negative control). 60 designs each.

Metrics used (fixed by the project plan):
  design_to_target_iptm      primary ranking score
  min_design_to_target_pae   secondary (lower is better)
  design_to_target_ipsae     ipSAE replication (Junker & Schoeder)
NEVER use quality_score across targets (within-batch rank-normalised),
and do not use design_ptm (batch-to-batch drift exceeds its own SD).
All of these are model CONFIDENCE, not affinity.

Usage:  python first_batch_stats.py <results_dir> <out_dir>
"""
import sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = sys.argv[1] if len(sys.argv) > 1 else "results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "."
TARGETS = ["bnp", "anp", "scr"]
METRICS = ["design_to_target_iptm", "min_design_to_target_pae", "design_to_target_ipsae"]
NOISE_FLOOR_IPTM = 0.05  # mean reproducibility, same target / new seed, n=20

d = {t: pd.read_csv(f"{RESULTS}/{t}_all_designs_metrics.csv") for t in TARGETS}

print("== data health ==")
for t, df in d.items():
    print(f"{t}: rows={len(df)} num_design={sorted(df.num_design.unique())} "
          f"dup_seq={df.designed_sequence.duplicated().sum()} "
          f"missing={ {m: int(df[m].isna().sum()) for m in METRICS} }")

print("\n== descriptives ==")
rows = []
for m in METRICS:
    for t in TARGETS:
        x = d[t][m].dropna()
        q = np.percentile(x, [25, 50, 75])
        rows.append(dict(metric=m, target=t, n=len(x), mean=x.mean(), sd=x.std(ddof=1),
                         min=x.min(), q25=q[0], median=q[1], q75=q[2], max=x.max()))
tab = pd.DataFrame(rows)
print(tab.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
tab.to_csv(f"{OUT}/first_batch_descriptives.csv", index=False)

def cohen_d(a, b):
    na, nb = len(a), len(b)
    s = np.sqrt(((na-1)*a.var(ddof=1) + (nb-1)*b.var(ddof=1)) / (na+nb-2))
    return (a.mean() - b.mean()) / s

print("\n== pairwise (Mann-Whitney U, two-sided; Cohen's d) ==")
for m in METRICS:
    for a, b in [("bnp", "anp"), ("bnp", "scr"), ("anp", "scr")]:
        x, y = d[a][m].dropna(), d[b][m].dropna()
        u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        diff = x.mean() - y.mean()
        note = ""
        if m == "design_to_target_iptm":
            note = ("  [|diff| <= noise floor 0.05 -> NOT interpretable]"
                    if abs(diff) <= NOISE_FLOOR_IPTM else
                    "  [|diff| > noise floor 0.05]")
        print(f"{m:26s} {a} vs {b}: diff={diff:+.4f} d={cohen_d(x,y):+.3f} "
              f"U={u:.0f} p={p:.3g}{note}")

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
colors = {"bnp": "#1f77b4", "anp": "#ff7f0e", "scr": "#7f7f7f"}
for ax, m in zip(axes, METRICS):
    lo = min(d[t][m].min() for t in TARGETS)
    hi = max(d[t][m].max() for t in TARGETS)
    bins = np.linspace(lo, hi, 18)
    for t in TARGETS:
        x = d[t][m].dropna()
        ax.hist(x, bins=bins, alpha=0.45, color=colors[t], label=f"{t} (n={len(x)})")
        ax.axvline(x.mean(), color=colors[t], ls="--", lw=1.5)
    ax.set_xlabel(m); ax.set_ylabel("designs")
    ax.set_title(m, fontsize=10)
    ax.legend(fontsize=8, frameon=False)
fig.suptitle("BoltzGen confidence metrics, 60 designs per target "
             "(dashed = mean; scr = shuffled-sequence negative control)", fontsize=11)
fig.tight_layout()
fig.savefig(f"{OUT}/first_batch_distributions.png", dpi=160)
print(f"\nwrote {OUT}/first_batch_distributions.png and {OUT}/first_batch_descriptives.csv")
