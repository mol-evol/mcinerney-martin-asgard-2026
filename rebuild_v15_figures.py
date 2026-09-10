#!/usr/bin/env python3
"""Rebuild main Fig. 1 and Extended Data Fig. 1 in the v15 layout.

Panel content is unchanged from v13; the changes are the ones required at revision:
  main Fig. 1  -- panels swapped (filter effect = a, occupancy diagnostic = b),
                  legend label "Head-to-head" -> "Conditional"
  ED Fig. 1    -- 8 panels -> 5: the Asgard:Alpha stratification and the filter sweep
                  are deleted (they are now main Fig. 1a,b) along with the Asgard
                  monophyly panel; d->c, f->d, g->e; full clade names in panel a.

Styling follows scripts/asgard_reanalysis.py and scripts/make_main_figure.py exactly.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))          # this scripts/ directory
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(HERE)   # Zenodo deposit dir
OUT  = sys.argv[2] if len(sys.argv) > 2 else HERE

CORE = 0.40
CS, CC, CI = '#1f4e79', '#c25e00', '#666666'

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'], 'font.size': 8,
                     'axes.labelsize': 8.5, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
                     'legend.fontsize': 7, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'svg.fonttype': 'path'})

nospine = lambda ax: [ax.spines[s].set_visible(False) for s in ('top', 'right')]

# ---------------------------------------------------------------- data
d = pd.read_csv(os.path.join(HERE, 'EPOC_slim_columns.tsv'), sep='\t', low_memory=False)
d['celw'] = pd.to_numeric(d['c-ELW'], errors='coerce').fillna(0.0)
d['key'] = list(zip(d.tree_name, d.euk_clade_rep))
mx = d.groupby('key')['celw'].max()
core = set(mx[mx > CORE].index)
D = d[d.key.isin(core)].copy()
NC = len(core)

TW = D.groupby(['key', 'prok_taxa'])['celw'].sum().groupby('prok_taxa').sum() / NC
OCC = D.groupby('prok_taxa')['key'].nunique() / NC
PIV = D.pivot_table(index='key', columns='prok_taxa', values='celw', aggfunc='sum')

scatter = pd.DataFrame({'taxon': TW.index, 'aelw': TW.values,
                        'in_set_pct': 100 * OCC.reindex(TW.index).values,
                        'aelw_cond': (TW / OCC.reindex(TW.index)).values})

LAB6 = ['Alphaproteobacteria', 'Cyanobacteriota', 'Actinomycetota',
        'Betaproteobacteria', 'Myxococcota', 'Gammaproteobacteria']
rows = []
for t in LAB6:
    sub = PIV[['Asgard', t]].dropna()
    rows.append(dict(label=t.replace('proteobacteria', '').replace('bacteriota', '').replace('mycetota', '')
                     .replace('coccota', ''),
                     standard_ratio=TW['Asgard'] / TW[t],
                     h2h_ratio=sub['Asgard'].mean() / sub[t].mean(), n=len(sub)))
h2h6 = pd.DataFrame(rows)
h2h6['label'] = ['Alpha', 'Cyano', 'Actino', 'Beta', 'Myxo', 'Gamma']

CCOUNT = D.groupby('key').size()
def strat(target):
    out = []
    for m in range(1, 13):
        ks = set(CCOUNT[CCOUNT >= m].index); dd = D[D.key.isin(ks)]; n = len(ks)
        tw = dd.groupby(['key', 'prok_taxa'])['celw'].sum().groupby('prok_taxa').sum() / n
        oc = dd.groupby('prok_taxa')['key'].nunique() / n
        t = dd.pivot_table(index='key', columns='prok_taxa', values='celw', aggfunc='sum')
        sub = t[['Asgard', target]].dropna()
        out.append(dict(min_c=m, std_ratio=tw['Asgard'] / tw[target],
                        h2h_ratio=sub['Asgard'].mean() / sub[target].mean(),
                        inc_ratio=oc['Asgard'] / oc[target]))
    return pd.DataFrame(out)

SA = strat('Alphaproteobacteria')
SC = strat('Cyanobacteriota')
sweep = pd.read_csv(os.path.join(HERE, 'filter_sweep.csv')).set_index('cutoff').reindex(['s0', 's10', 's25', 's67'])
flip = pd.read_csv(os.path.join(HERE, 'winner_flip_s10_vs_s0.csv')).set_index('s10_call')
kk = pd.read_csv(os.path.join(HERE, 'kegg_stratification.csv')).set_index('kegg_class')
trees = pd.read_csv(os.path.join(HERE, 'asgard_winner_tree_analysis.csv'))
NOBACT = 100 * (~trees['has_bact'].astype(bool)).mean()
cats = [c for c in ['informational', 'metabolic', 'other'] if c in kk.index]

print(f'core {NC} | aELW ratio {TW["Asgard"]/TW["Alphaproteobacteria"]:.2f} | '
      f'h2h {SA.h2h_ratio.iloc[0]:.2f} | strat {SA.std_ratio.iloc[0]:.2f}->{SA.std_ratio.iloc[-1]:.2f} | '
      f'no-bact {NOBACT:.1f}% | trees {len(trees)}')

def stratplot(ax, data, tgt, ylim, letter_):
    ax.plot(data.min_c, data.std_ratio, 'o-', color=CS, lw=1.3, ms=4.5, mec='white', mew=0.6,
            label='Standard aELW ratio', zorder=3)
    ax.plot(data.min_c, data.h2h_ratio, 's-', color=CC, lw=1.3, ms=4, mec='white', mew=0.6,
            label='Conditional c-ELW ratio', zorder=3)
    ax.plot(data.min_c, data.inc_ratio, '^--', color=CI, lw=1.1, ms=4, mec='white', mew=0.6,
            label='Inclusion ratio', alpha=0.85, zorder=2)
    ax.axhline(1, color='black', lw=0.5, alpha=0.5)
    ax.text(12.4, 1 + (ylim[1] - ylim[0]) * 0.012, 'parity', fontsize=6.5, ha='right', alpha=0.7, style='italic')
    ax.annotate(f"{data.std_ratio.iloc[0]:.2f}", xy=(1, data.std_ratio.iloc[0]), xytext=(8, 0),
                textcoords='offset points', fontsize=7, color=CS, fontweight='bold', va='center')
    ax.annotate(f"{data.std_ratio.iloc[-1]:.2f}", xy=(12, data.std_ratio.iloc[-1]), xytext=(-8, 6),
                textcoords='offset points', fontsize=7, color=CS, ha='right')
    ax.set_xlabel('Minimum candidates per EPOC tuple'); ax.set_ylabel(f'Asgard : {tgt} ratio')
    ax.set_xticks(range(1, 13)); ax.set_xlim(0.5, 12.5); ax.set_ylim(*ylim)
    nospine(ax); ax.legend(loc='upper right', frameon=False, fontsize=7)
    ax.text(-0.16, 1.08, letter_, transform=ax.transAxes, fontsize=12, fontweight='bold', va='top')

def sweepplot(ax, letter_):
    ax.bar(range(4), sweep['ratio'], color=['#c0392b', CS, CS, CS])
    ax.axhline(1, ls='--', c='grey', lw=0.8); ax.set_yscale('log'); ax.set_ylim(0.3, 25)
    for i, v in enumerate(sweep['ratio']):
        ax.text(i, v * 1.07, f"{v:.2f}", ha='center', va='bottom', fontsize=7.5)
    ax.set_xticks(range(4)); ax.set_xticklabels(["s0\n(no filter)", "s10\n(headline)", "s25", "s67"])
    ax.set_ylabel('Asgard : Alphaproteobacteria aELW ratio'); nospine(ax)
    ax.set_title('Remove the soft-core filter (s0)', fontsize=8.5, loc='left', pad=2)
    ax.text(-0.16, 1.08, letter_, transform=ax.transAxes, fontsize=12, fontweight='bold', va='top')

# ---------------------------------------------------------------- main Fig. 1
fig = plt.figure(figsize=(7.2, 3.2))
gs = GridSpec(1, 2, figure=fig, wspace=0.34, left=0.09, right=0.97, top=0.9, bottom=0.16)
sweepplot(fig.add_subplot(gs[0, 0]), 'a')
stratplot(fig.add_subplot(gs[0, 1]), SA, 'Alphaproteobacteria', (0.5, 9), 'b')
for ext in ('svg', 'pdf', 'png'):
    fig.savefig(os.path.join(OUT, f"Figure1_v15_regen.{ext}"), dpi=300 if ext == 'png' else None, bbox_inches='tight')
plt.close(fig)

# ---------------------------------------------------------------- ED Fig. 1
fig = plt.figure(figsize=(7.5, 8.9))
gs = GridSpec(3, 4, figure=fig, hspace=0.46, wspace=0.95, left=0.085, right=0.975, top=0.97, bottom=0.06)
A = fig.add_subplot(gs[0, 0:2]); B = fig.add_subplot(gs[0, 2:4]); C = fig.add_subplot(gs[1, 0:2])
Dx = fig.add_subplot(gs[1, 2:4]); E = fig.add_subplot(gs[2, 0:2]); F = fig.add_subplot(gs[2, 2:4])
letter = lambda ax, l: ax.text(-0.16, 1.06, l, transform=ax.transAxes, fontsize=12, fontweight='bold', va='top')

# a
A.scatter(scatter.in_set_pct, scatter.aelw, s=30, c=CS, edgecolor='white', linewidth=0.5, label='Standard aELW', zorder=3)
A.scatter(scatter.in_set_pct, scatter.aelw_cond, s=30, c=CC, edgecolor='white', linewidth=0.5, label='Conditional aELW', zorder=3)
xr = np.array([0, 70])
def fit(x, y):
    b1, b0 = np.polyfit(x, y, 1); return b0, b1, float(np.corrcoef(x, y)[0, 1] ** 2)
i_s, ss, rs = fit(scatter.in_set_pct, scatter.aelw); i_c, sc, rc = fit(scatter.in_set_pct, scatter.aelw_cond)
A.plot(xr, i_s + ss * xr, '--', c=CS, lw=1.1, alpha=0.7); A.plot(xr, i_c + sc * xr, '--', c=CC, lw=1.1, alpha=0.7)
# labels placed in the empty right-hand region with leader lines (full names, referee request)
STACK = {'Cyanobacteriota': 0.660, 'Actinomycetota': 0.560, 'Alphaproteobacteria': 0.460,
         'Betaproteobacteria': 0.360, 'Myxococcota': 0.260, 'Gammaproteobacteria': 0.160}
FULL = {'Asgard': 'Asgardarchaeota'}
sc = scatter.set_index('taxon')
for name, ly in STACK.items():
    px, py = sc.loc[name, 'in_set_pct'], sc.loc[name, 'aelw_cond']
    A.annotate(name, xy=(px, py), xytext=(38.0, ly), textcoords='data', fontsize=6.6,
               ha='left', va='center',
               arrowprops=dict(arrowstyle='-', lw=0.6, color='#555555', shrinkA=1, shrinkB=2,
                               connectionstyle='arc3,rad=0.0'))
A.annotate(FULL['Asgard'], (sc.loc['Asgard', 'in_set_pct'], sc.loc['Asgard', 'aelw_cond']),
           xytext=(-8, 9), textcoords='offset points', fontsize=6.8, ha='right', fontweight='bold')
A.text(0.04, 0.96, f'Standard: $R^2$ = {rs:.2f}', transform=A.transAxes, fontsize=7, color=CS, va='top')
A.text(0.04, 0.88, f'Conditional: $R^2$ = {rc:.2f}', transform=A.transAxes, fontsize=7, color=CC, va='top')
A.set_xlabel('Candidate-set inclusion (% of core EPOC tuples)'); A.set_ylabel('aELW')
A.set_xlim(0, 70); A.set_ylim(-0.02, 0.78)
A.legend(loc='lower right', frameon=False, fontsize=6.8); nospine(A); letter(A, 'a')

# b
xp = np.arange(len(h2h6)); w = 0.38
B.bar(xp - w / 2, h2h6.standard_ratio, w, color=CS, edgecolor='white', linewidth=0.5, label='Standard aELW ratio')
B.bar(xp + w / 2, h2h6.h2h_ratio, w, color=CC, edgecolor='white', linewidth=0.5, label='Conditional c-ELW ratio')
B.axhline(1, color='black', lw=0.5, alpha=0.6)
B.text(len(h2h6) - 0.5, 1.2, 'parity', fontsize=6.5, ha='right', alpha=0.7, style='italic')
for i, (sv, hv) in enumerate(zip(h2h6.standard_ratio, h2h6.h2h_ratio)):
    B.text(i - w / 2, sv + 0.25, f'{sv:.1f}', ha='center', fontsize=6.3, color=CS)
    B.text(i + w / 2, hv + 0.25, f'{hv:.2f}', ha='center', fontsize=6.3, color=CC)
B.set_xticks(xp); B.set_xticklabels(h2h6.label, fontsize=7.5)
B.set_xlabel('Bacterial candidate (vs Asgard)'); B.set_ylabel('Asgard : taxon ratio'); B.set_ylim(0, 17)
B.legend(loc='upper left', frameon=False, fontsize=7); nospine(B); letter(B, 'b')

# c
stratplot(C, SC, 'Cyano', (0.5, 5), 'c')

# d -- within-gene flip
gB = [g for g in ['Asgard', 'Bacteria'] if g in flip.index]; bott = np.zeros(len(gB))
for col, c, lab in [('to_Bacteria', CC, 'to Bacteria'), ('to_Asgard', CS, 'to Asgard'), ('to_Archaea', '#7f8c8d', 'to Archaea')]:
    vals = np.array([flip.loc[g, col] for g in gB])
    Dx.bar(range(len(gB)), vals, bottom=bott, color=c, label=lab)
    for i, (v, b) in enumerate(zip(vals, bott)):
        if v >= 8: Dx.text(i, b + v / 2, f"{v:.0f}%", ha='center', va='center', fontsize=7.5, color='white', fontweight='bold')
    bott = bott + vals
Dx.set_xticks(range(len(gB))); Dx.set_xticklabels([g + " call\n(s10)" for g in gB])
Dx.set_ylabel('re-assignment, filter OFF (%)'); Dx.set_ylim(0, 100)
Dx.legend(frameon=False, fontsize=6.5, loc='lower center'); nospine(Dx); letter(Dx, 'd')
Dx.set_title('Same genes, filter off', fontsize=8.5, loc='left', pad=2)

# e -- tree bacterial depletion
vals = [kk.loc[c, 'pct_no_bacteria'] for c in cats]
E.bar(range(len(cats)), vals, color='#16a085'); E.axhline(NOBACT, ls='--', c='grey', lw=0.8)
E.text(len(cats) - 0.5, NOBACT + 1.5, f"all {NOBACT:.0f}%", ha='right', fontsize=6.5, color='grey')
for i, v in enumerate(vals): E.text(i, v + 1, f"{v:.0f}%", ha='center', fontsize=7.5)
E.set_xticks(range(len(cats))); E.set_xticklabels(cats, rotation=12, fontsize=7); E.set_ylim(0, 100)
E.set_ylabel('% Asgard-winner trees\nwith NO bacterium'); nospine(E); letter(E, 'e')
E.set_title('Filter empties trees of bacteria', fontsize=8.5, loc='left', pad=2)

# f -- same gene families with and without the filter, overall and by function (revision 2)
mf = pd.read_csv(os.path.join(HERE, 'matched_family_asgard_share_by_class.csv')).set_index('kegg_class')
rows = ['all', 'informational', 'metabolic', 'other']
xf = np.arange(len(rows)); wf = 0.38
F.bar(xf - wf / 2, mf.loc[rows, 'Asgard%_shared_s0'], wf, color=CI, edgecolor='white', linewidth=0.5, label='No filter (s0)')
F.bar(xf + wf / 2, mf.loc[rows, 'Asgard%_shared_s10'], wf, color=CS, edgecolor='white', linewidth=0.5, label='Filter (s10)')
for i, r in enumerate(rows):
    a0, a1 = mf.loc[r, 'Asgard%_shared_s0'], mf.loc[r, 'Asgard%_shared_s10']
    F.text(i - wf / 2, a0 + 1, f"{a0:.0f}", ha='center', fontsize=6.5, color=CI)
    F.text(i + wf / 2, a1 + 1, f"{a1:.0f}", ha='center', fontsize=6.5, color=CS)
F.set_xticks(xf); F.set_xticklabels([f"{r}\n(n = {int(mf.loc[r, 'fam_shared']):,})" for r in rows], fontsize=6.5)
F.set_ylabel('Asgard share of summed\ncore c-ELW (%)'); F.set_ylim(0, 45)
F.legend(frameon=False, fontsize=6.5, loc='upper left'); nospine(F); letter(F, 'f')
F.set_title('Same gene families, filter off vs on', fontsize=8.5, loc='left', pad=2)

for ext in ('svg', 'pdf', 'png'):
    fig.savefig(os.path.join(OUT, f"ExtendedData_Fig1_v15.{ext}"), dpi=300 if ext == 'png' else None, bbox_inches='tight')
plt.close(fig)
print('written: Figure1_v15_regen.* and ExtendedData_Fig1_v15.*')
