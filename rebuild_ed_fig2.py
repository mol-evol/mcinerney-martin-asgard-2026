#!/usr/bin/env python3
"""Rebuild Extended Data Fig. 2.

Revision 2: panel a now shows Asgard:bacterial ratios of raw stem, normaliser and
normalised stem separately for broad (LECA-like) and narrow eukaryotic clades; the
normaliser explanation holds only in narrow clades. Panel b no longer invokes
"adaptive distance". Needs EPOC_data.tsv in DATA_DIR for clade breadth.

Earlier (v15) framing:

Panel content and statistics are unchanged from v13. The change is the target of the
argument: Tobiasson et al. explicitly reject the reading of stem length as acquisition
time and propose post-acquisition adaptation instead, so panel b's title no longer
attacks a timing interpretation they do not hold.

Reproduces scripts/asgard_reanalysis.py part9_stemlength / _stem_figure exactly, but
reads the reduced inputs in this directory (EPOC_slim_columns.tsv, KEGG_annotation_3col.tsv.gz)
so the 1.2 GB tree archive is not required.

Usage: python rebuild_ed_fig2.py [DATA_DIR] [OUT_DIR]   (DATA_DIR supplies KEGG_category_mapping.tsv)
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))          # this scripts/ directory
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(HERE)   # Zenodo deposit dir
OUT  = sys.argv[2] if len(sys.argv) > 2 else HERE

CORE_THRESHOLD = 0.40
ARCH = {"Euryarchaeota", "TACK group", "DPANN group", "Thermoplasmata"}
CS, CC, CG = '#1f4e79', '#c25e00', '#999999'

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'], 'font.size': 8,
                     'axes.labelsize': 8.5, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
                     'legend.fontsize': 7, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'svg.fonttype': 'path'})

# ---------------------------------------------------------------- data
x = pd.read_csv(os.path.join(HERE, 'EPOC_slim_columns.tsv'), sep='\t', low_memory=False)
x['celw'] = pd.to_numeric(x['c-ELW'], errors='coerce').fillna(0.0)
for c in ['stem_length', 'raw_stem_length', 'median_euk_leaf_dist']:
    x[c] = pd.to_numeric(x[c], errors='coerce')
x['key'] = list(zip(x.tree_name, x.euk_clade_rep))
mx = x.groupby('key')['celw'].max()
core = set(mx[mx > CORE_THRESHOLD].index)
d = x[x.key.isin(core)]
win = d.loc[d.groupby('key')['celw'].idxmax()].copy()
dom = lambda t: "Asgard" if t == "Asgard" else ("Archaea" if t in ARCH else "Bacteria")
win['domn'] = win['prok_taxa'].map(dom)
win = win[np.isfinite(win['stem_length'])]
# revision 2: eukaryotic clade breadth, from the full deposited file
br = pd.read_csv(os.path.join(DATA, 'EPOC_data.tsv'), sep='\t', usecols=['tree_name', 'euk_clade_rep', 'euk_LCA', 'euk_scope_len'])
win = win.merge(br.drop_duplicates(['tree_name', 'euk_clade_rep']), on=['tree_name', 'euk_clade_rep'], how='left')
win['breadth'] = np.where((win.euk_LCA == 'Eukaryota') & (win.euk_scope_len > 5), 'broad',
                          np.where(win.euk_scope_len <= 5, 'narrow', 'intermediate'))

# reproduce part9's mapping exactly: pandas sort_values('Prob') + drop_duplicates(keep='last')
a = pd.read_csv(os.path.join(HERE, 'KEGG_annotation_3col.tsv.gz'), sep='\t',
                usecols=['Query', 'Target', 'Prob'], dtype=str, low_memory=False)
a['Prob'] = pd.to_numeric(a['Prob'], errors='coerce')
# NOTE: the original uses sort_values('Prob') with the default (unstable) quicksort.
# 24.5% of queries have their maximum Prob tied across >1 distinct target, so that
# assignment is not reproducible. kind='mergesort' makes the tie-break deterministic
# (last occurrence in file order among the maximum).
b = a.sort_values('Prob', kind='mergesort').drop_duplicates('Query', keep='last')
t2k = dict(zip(b['Query'], b['Target']))
cm = pd.read_csv(os.path.join(DATA, 'KEGG_category_mapping.tsv'), sep='\t', dtype=str)
k2c = cm.groupby('kogid')['category_id'].apply(set).to_dict()
win['cats'] = win['tree_name'].map(t2k).map(lambda k: k2c.get(k, set()))

donors = ["Asgard", "Cyanobacteriota", "Alphaproteobacteria", "Gammaproteobacteria", "Actinomycetota"]
perdon = {t: dict(raw=round(win[win.prok_taxa == t].raw_stem_length.median(), 4),
                  norm=round(win[win.prok_taxa == t].stem_length.median(), 4),
                  denom=round(win[win.prok_taxa == t].median_euk_leaf_dist.median(), 3)) for t in donors}

A_ = win[win.domn == 'Asgard']; B_ = win[win.domn == 'Bacteria']
print('raw   Asgard %.4f | Bacteria %.4f' % (A_.raw_stem_length.median(), B_.raw_stem_length.median()))
print('norm  Asgard %.4f | Bacteria %.4f' % (A_.stem_length.median(), B_.stem_length.median()))
print('denom Asgard %.3f | Bacteria %.3f  (%.0f%% larger)' % (
    A_.median_euk_leaf_dist.median(), B_.median_euk_leaf_dist.median(),
    100 * (A_.median_euk_leaf_dist.median() / B_.median_euk_leaf_dist.median() - 1)))
for t in donors:
    print(f'  {t:22s} raw {perdon[t]["raw"]:.4f}  norm {perdon[t]["norm"]:.4f}  denom {perdon[t]["denom"]:.3f}')

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(8.0, 3.5))
gs = GridSpec(1, 2, figure=fig, wspace=0.95, left=0.075, right=0.95, top=0.84, bottom=0.17, width_ratios=[1.15, 1])
A = fig.add_subplot(gs[0, 0]); B = fig.add_subplot(gs[0, 1])
nospine = lambda ax: [ax.spines[s].set_visible(False) for s in ('top', 'right')]
letter = lambda ax, l: ax.text(-0.17, 1.10, l, transform=ax.transAxes, fontsize=12, fontweight='bold', va='top')

METRICS = [('raw_stem_length', 'Raw stem'), ('median_euk_leaf_dist', 'Normaliser'), ('stem_length', 'Normalised stem')]
GROUPS = [('broad', 'Broad clades (LECA-like)', CS), ('narrow', 'Narrow clades (≤5 taxa)', CC)]
xx = np.arange(len(METRICS)); w = 0.38
for j, (g, glab, col) in enumerate(GROUPS):
    sub = win[win.breadth == g]
    ratios = [sub[sub.domn == 'Asgard'][m].median() / sub[sub.domn == 'Bacteria'][m].median() for m, _ in METRICS]
    na, nb = (sub.domn == 'Asgard').sum(), (sub.domn == 'Bacteria').sum()
    A.bar(xx + (j - 0.5) * w, ratios, w, color=col, edgecolor='white', linewidth=0.5, label=f'{glab}; n = {na:,}/{nb:,}')
    for i, r in enumerate(ratios):
        A.text(i + (j - 0.5) * w, r + 0.015, f'{r:.2f}', ha='center', fontsize=6.5, color=col)
    print(f'  {g:7s} Asgard/Bacteria medians: ' + ', '.join(f'{l} {r:.3f}' for (_, l), r in zip(METRICS, ratios)))
from scipy.stats import mannwhitneyu
for g in ['all', 'broad', 'narrow']:
    sub = win if g == 'all' else win[win.breadth == g]
    x, y = sub[sub.domn == 'Asgard'].stem_length.dropna(), sub[sub.domn == 'Bacteria'].stem_length.dropna()
    U, p = mannwhitneyu(x, y)
    print(f'  {g:7s} normalised stem, Asgard vs Bacteria: rank-biserial {2 * U / (len(x) * len(y)) - 1:+.3f} (p = {p:.1e})')
A.axhline(1, color='black', lw=0.5, alpha=0.6)
A.set_xticks(xx); A.set_xticklabels([l for _, l in METRICS], fontsize=7.2)
A.set_ylabel('Asgard : bacterial median'); A.set_ylim(0.6, 1.4)
A.legend(loc='upper left', frameon=False, fontsize=6.3); nospine(A)
A.set_title("Source of the Asgard–bacterial difference\ndiffers with eukaryotic clade breadth", fontsize=8, loc='left', pad=3)
letter(A, 'a')

sets = [("Ribosome\n(Asgard, 'vertical')", win.cats.map(lambda s: 'map03010' in s) & (win.domn == 'Asgard'), CS),
        ("Oxidative phosph.\n(Alpha, one symbiosis)", win.cats.map(lambda s: 'map00190' in s) & (win.prok_taxa == 'Alphaproteobacteria'), CC)]
for i, (lab, m, c) in enumerate(sets):
    v = win.loc[m, "stem_length"]; v = v[v > 1e-6].values
    y = np.full(len(v), i) + (np.random.RandomState(0).rand(len(v)) - 0.5) * 0.28
    B.scatter(v, y, s=9, color=c, alpha=0.55, edgecolor='none', zorder=2)
    lo, hi = np.quantile(v, .025), np.quantile(v, .975)
    B.plot([lo, hi], [i, i], color=c, lw=2.4, alpha=0.9, zorder=3, solid_capstyle='round')
    B.text(hi * 1.15, i, f"{hi/lo:.0f}x spread", va='center', fontsize=7, color=c, fontweight='bold')
    print(f'  set {lab.splitlines()[0]:20s} n={len(v):4d}  {lo:.3f}-{hi:.3f}  {hi/lo:.0f}x')
B.set_xscale('log'); B.set_yticks([0, 1]); B.set_yticklabels([s[0] for s in sets], fontsize=7.2)
B.set_xlim(0.01, 6); B.set_ylim(-0.6, 1.6)
B.set_xlabel('Normalised stem length (log scale)')
B.set_title("One donor, one event:\nnormalised stems span 50–70-fold", fontsize=8, loc='left', pad=3)
nospine(B); letter(B, 'b')

for ext in ('svg', 'pdf', 'png'):
    fig.savefig(os.path.join(OUT, f"ExtendedData_Fig2_v15.{ext}"), dpi=300 if ext == 'png' else None, bbox_inches='tight')
plt.close(fig)
print('written ExtendedData_Fig2_v15.*')
