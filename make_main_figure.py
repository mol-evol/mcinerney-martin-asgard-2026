#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main-text Figure 1 (two panels) for the Matters Arising.
  a: Asgard:Alphaproteobacteria aELW ratio across the four deposited soft-core
     cut-offs (s0=no filter, s10, s25, s67) -- the filter effect; the ratio inverts
     to 0.44 when the filter is removed.
  b: Asgard:Alphaproteobacteria ratio vs minimum candidate count per tuple -- the
     standard (occupancy-weighted) ratio collapses onto the stable conditional value
     (the occupancy diagnostic).
Panel order follows the v15 revision: the filter effect leads, the occupancy
decomposition is the diagnostic. Neither panel is repeated in Extended Data Fig. 1.
Drawn from the deposited data (Zenodo 15048010).
Usage: python make_main_figure.py [DATA_DIR] [OUT_DIR]
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(HERE)
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(HERE), "figures")
os.makedirs(OUT, exist_ok=True)
CORE = 0.40
SF = {"s0": "EPOC_data.pangenome_s0.tsv", "s10": "EPOC_data.tsv",
      "s25": "EPOC_data.pangenome_s25.tsv", "s67": "EPOC_data.pangenome_s67.tsv"}

def load(fn):
    df = pd.read_csv(fn, sep="\t", usecols=["tree_name", "euk_clade_rep", "prok_taxa", "c-ELW"], low_memory=False)
    df["celw"] = pd.to_numeric(df["c-ELW"], errors="coerce").fillna(0.0)
    df["key"] = list(zip(df["tree_name"], df["euk_clade_rep"]))
    return df

def core_rows(df):
    mx = df.groupby("key")["celw"].max(); core = set(mx[mx > CORE].index)
    return df[df["key"].isin(core)], len(core)

# panel a: stratification on s10
d, ncore = core_rows(load(os.path.join(DATA, SF["s10"])))
cc = d.groupby("key").size()
mc, std, h2h, inc = [], [], [], []
for m in range(1, 13):
    ks = set(cc[cc >= m].index); dd = d[d["key"].isin(ks)]; nc = len(ks)
    tw = dd.groupby(["key", "prok_taxa"])["celw"].sum().groupby("prok_taxa").sum() / nc
    occ = dd.groupby("prok_taxa")["key"].nunique() / nc
    t = dd.pivot_table(index="key", columns="prok_taxa", values="celw", aggfunc="sum")
    sub = t[["Asgard", "Alphaproteobacteria"]].dropna()
    mc.append(m); std.append(tw["Asgard"] / tw["Alphaproteobacteria"])
    h2h.append(sub["Asgard"].mean() / sub["Alphaproteobacteria"].mean())
    inc.append(occ["Asgard"] / occ["Alphaproteobacteria"])

# panel b: filter sweep
sweep = {}
for tag, fn in SF.items():
    dd, nc = core_rows(load(os.path.join(DATA, fn)))
    tw = dd.groupby(["key", "prok_taxa"])["celw"].sum().groupby("prok_taxa").sum() / nc
    sweep[tag] = tw["Asgard"] / tw["Alphaproteobacteria"]
order = ["s0", "s10", "s25", "s67"]; ratios = [sweep[t] for t in order]

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                     'font.size': 8, 'axes.labelsize': 8.5, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
                     'legend.fontsize': 7, 'axes.linewidth': 0.7})
CS, CC, CI = '#1f4e79', '#c25e00', '#666666'
fig = plt.figure(figsize=(7.2, 3.2)); gs = GridSpec(1, 2, figure=fig, wspace=0.34, left=0.09, right=0.97, top=0.9, bottom=0.16)
# v15: filter sweep (B) is panel a on the left; stratification (A) is panel b on the right
B = fig.add_subplot(gs[0, 0]); A = fig.add_subplot(gs[0, 1])
def nospine(ax): [ax.spines[s].set_visible(False) for s in ('top', 'right')]
def letter(ax, l): ax.text(-0.16, 1.08, l, transform=ax.transAxes, fontsize=12, fontweight='bold', va='top')

A.plot(mc, std, 'o-', color=CS, lw=1.3, ms=4.5, mec='white', mew=0.6, label='Standard aELW ratio', zorder=3)
A.plot(mc, h2h, 's-', color=CC, lw=1.3, ms=4, mec='white', mew=0.6, label='Conditional c-ELW ratio', zorder=3)
A.plot(mc, inc, '^--', color=CI, lw=1.1, ms=4, mec='white', mew=0.6, label='Inclusion ratio', alpha=0.85, zorder=2)
A.axhline(1, color='black', lw=0.5, alpha=0.5); A.text(12.3, 1.15, 'parity', fontsize=6.5, ha='right', alpha=0.7, style='italic')
A.annotate(f"{std[0]:.2f}", xy=(1, std[0]), xytext=(8, 0), textcoords='offset points', fontsize=7.5, color=CS, fontweight='bold', va='center')
A.annotate(f"{std[-1]:.2f}", xy=(12, std[-1]), xytext=(-6, 7), textcoords='offset points', fontsize=7.5, color=CS, ha='right')
A.set_xlabel('Minimum candidates per EPOC tuple'); A.set_ylabel('Asgard : Alphaproteobacteria ratio')
A.set_xticks(range(1, 13)); A.set_xlim(0.5, 12.5); A.set_ylim(0.5, 9); nospine(A)
A.legend(loc='upper right', frameon=False, fontsize=7); letter(A, 'b')

cols = ['#c0392b', CS, CS, CS]; B.bar(range(4), ratios, color=cols)
B.axhline(1, ls='--', c='grey', lw=0.8); B.set_yscale('log'); B.set_ylim(0.3, 25)
for i, v in enumerate(ratios): B.text(i, v * 1.07, f"{v:.2f}", ha='center', va='bottom', fontsize=7.5)
B.set_xticks(range(4)); B.set_xticklabels(["s0\n(no filter)", "s10\n(headline)", "s25", "s67"])
B.set_ylabel('Asgard : Alphaproteobacteria aELW ratio'); nospine(B); letter(B, 'a')
B.set_title('Remove the soft-core filter (s0)', fontsize=8.5, loc='left', pad=2)

for ext in ('svg', 'pdf', 'png'):
    fig.savefig(os.path.join(OUT, f"Figure1.{ext}"), dpi=300 if ext == 'png' else None, bbox_inches='tight')
print("wrote Figure1.svg/.pdf/.png to", OUT, "| panel a std", round(std[0], 2), "->", round(std[-1], 2),
      "h2h", round(np.mean(h2h), 2), "| sweep", [round(r, 2) for r in ratios])
