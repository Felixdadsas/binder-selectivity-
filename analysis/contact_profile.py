#!/usr/bin/env python3
"""Per-position target contact frequency for BoltzGen binder designs.

For each refolded complex (mmCIF with two polymer chains: A = designed binder,
B = target peptide) we compute, for every target residue, the minimum
heavy-atom distance to any binder heavy atom. A residue counts as "in contact"
when that distance is <= the cutoff. The contact frequency of a target position
is the fraction of designs in which it is in contact.

Heavy atoms (all non-hydrogen atoms, side chains included) rather than CB are
used because the question is which parts of the target surface the binder
actually touches; CB-only would blur long side chains (Arg, Lys) that dominate
these peptides. Two cutoffs are reported: 4.5 A (close, roughly a vdW/H-bond
contact) and 5.0 A (loose, the usual "interface residue" definition).

Usage: contact_profile.py <cif_dir> <out_dir>
Files are grouped by the prefix of their filename before the first "_".
"""
import sys, os, glob, csv
import numpy as np

AA3TO1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
}
CUTOFFS = (4.5, 5.0)


def parse_cif_atoms(path):
    """Return dict of columns -> list for the _atom_site loop (ATOM records)."""
    cols, rows, in_loop, header_done = [], [], False, False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('loop_'):
                if cols and rows:
                    break
                cols, in_loop, header_done = [], False, False
                continue
            if s.startswith('_atom_site.'):
                cols.append(s.split('.', 1)[1].split()[0])
                in_loop = True
                continue
            if in_loop and cols:
                if not s or s.startswith('#'):
                    if header_done:
                        break
                    continue
                if s.startswith('_'):
                    break
                parts = s.split()
                if len(parts) < len(cols):
                    continue
                header_done = True
                rows.append(parts[:len(cols)])
    if not cols or not rows:
        raise ValueError(f'no atom_site loop in {path}')
    idx = {c: i for i, c in enumerate(cols)}
    return idx, rows


def load_complex(path):
    """-> (binder_xyz, target_xyz_by_residue, target_residues)"""
    idx, rows = parse_cif_atoms(path)
    ci = lambda name: idx[name]
    binder, target = [], {}
    for r in rows:
        if r[0] not in ('ATOM', 'HETATM'):
            continue
        if r[ci('type_symbol')] == 'H':
            continue
        xyz = (float(r[ci('Cartn_x')]), float(r[ci('Cartn_y')]), float(r[ci('Cartn_z')]))
        chain = r[ci('label_asym_id')]
        if chain == 'A':
            binder.append(xyz)
        elif chain == 'B':
            pos = int(r[ci('label_seq_id')])
            comp = r[ci('label_comp_id')]
            target.setdefault((pos, comp), []).append(xyz)
    if not binder or not target:
        raise ValueError(f'missing chain A or B in {path}')
    keys = sorted(target)
    return np.asarray(binder), [np.asarray(target[k]) for k in keys], keys


def min_dists(path):
    binder, tgt_atoms, keys = load_complex(path)
    out = []
    for atoms in tgt_atoms:
        d = np.sqrt(((atoms[:, None, :] - binder[None, :, :]) ** 2).sum(-1))
        out.append(d.min())
    return np.asarray(out), keys


def main():
    cif_dir, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    groups = {}
    for p in sorted(glob.glob(os.path.join(cif_dir, '*.cif'))):
        groups.setdefault(os.path.basename(p).split('_')[0], []).append(p)

    summary, profiles = [], {}
    for target, paths in sorted(groups.items()):
        mats, keys_ref = [], None
        for p in paths:
            d, keys = min_dists(p)
            if keys_ref is None:
                keys_ref = keys
            elif keys != keys_ref:
                raise ValueError(f'target residue numbering differs in {p}')
            mats.append(d)
        M = np.vstack(mats)                       # designs x positions
        n = M.shape[0]
        seq = ''.join(AA3TO1.get(c, 'X') for _, c in keys_ref)
        rows = []
        for j, (pos, comp) in enumerate(keys_ref):
            row = {'position': pos, 'resname': comp, 'aa': AA3TO1.get(comp, 'X'),
                   'n_designs': n,
                   'mean_min_dist_A': round(float(M[:, j].mean()), 3),
                   'sd_min_dist_A': round(float(M[:, j].std(ddof=1)), 3)}
            for c in CUTOFFS:
                f = float((M[:, j] <= c).mean())
                row[f'contact_freq_{c}A'] = round(f, 4)
                row[f'se_{c}A'] = round((f * (1 - f) / n) ** 0.5, 4)
            rows.append(row)
        csv_path = os.path.join(out_dir, f'{target}_perposition.csv')
        with open(csv_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        f5 = np.array([r['contact_freq_5.0A'] for r in rows])
        cys = [pos for (pos, comp) in keys_ref if comp == 'CYS']
        summary.append({
            'target': target, 'n_designs': n, 'n_target_residues': len(rows),
            'sequence': seq, 'cys_positions': '|'.join(map(str, cys)),
            'mean_contact_freq_5.0A': round(float(f5.mean()), 4),
            'min_contact_freq_5.0A': round(float(f5.min()), 4),
            'max_contact_freq_5.0A': round(float(f5.max()), 4),
            'sd_across_positions': round(float(f5.std(ddof=1)), 4),
            'n_positions_freq_ge_0.9': int((f5 >= 0.9).sum()),
            'n_positions_freq_le_0.5': int((f5 <= 0.5).sum()),
            'mean_residues_in_contact_per_design':
                round(float((np.vstack(mats) <= 5.0).sum(1).mean()), 2),
        })
        profiles[target] = (rows, f5, cys, seq)

    with open(os.path.join(out_dir, 'contact_profile_summary.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)

    # ---- figure ----
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    order = [t for t in ('bnp', 'anp', 'scr') if t in profiles] + \
            [t for t in profiles if t not in ('bnp', 'anp', 'scr')]
    fig, axes = plt.subplots(len(order), 1, figsize=(11, 2.6 * len(order)), sharey=True)
    if len(order) == 1:
        axes = [axes]
    for ax, t in zip(axes, order):
        rows, f5, cys, seq = profiles[t]
        x = [r['position'] for r in rows]
        f45 = [r['contact_freq_4.5A'] for r in rows]
        ax.bar(x, f5, color='#4C78A8', label='contact freq (5.0 Å)')
        ax.plot(x, f45, color='#E45756', marker='o', ms=3, lw=1,
                label='contact freq (4.5 Å)')
        if len(cys) >= 2:
            ax.axvspan(min(cys), max(cys), color='#FFD166', alpha=0.25, zorder=0,
                       label=f'disulfide ring Cys{min(cys)}–Cys{max(cys)}')
        ax.set_xticks(x)
        ax.set_xticklabels([f'{a}\n{p}' for (a, p) in zip(seq, x)], fontsize=6)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel('frequency')
        ax.set_title(f'{t}  (n = {rows[0]["n_designs"]} designs)', fontsize=10, loc='left')
        if ax is axes[0]:
            ax.legend(fontsize=8, ncol=3, loc='lower center',
                      bbox_to_anchor=(0.5, 1.18), frameon=False)
    axes[-1].set_xlabel('target residue (one-letter code / position)')
    fig.suptitle('Per-position target contact frequency, BoltzGen batch 1 (single seed)',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(out_dir, 'contact_profile.png'), dpi=180)

    for s in summary:
        print(s)


if __name__ == '__main__':
    main()
