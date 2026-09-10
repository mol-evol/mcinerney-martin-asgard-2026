#!/usr/bin/env python3
"""Revision 2: numbers used to answer the Reply of Tobiasson et al. and Referee 2.

1. What Reply Table 1 counts: Asgard as candidate vs as best-supported sister at s0.
2. Reply's unfiltered summary statistics, reproduced from the deposited s0 file.
3. The same gene families with and without the filter (2,751 families in both runs).
4. Percentage of summed core c-ELW per group at each cut-off (Extended Data Table 1).

Usage: python reply_checks.py [DATA_DIR]      (DATA_DIR = Zenodo 15048010 deposit)
Writes: ed_table1_percent_contribution.csv
"""
import os, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(HERE)
K = ["tree_name", "euk_clade_rep"]
CORE = 0.40
ARCH = {"Euryarchaeota", "TACK group", "DPANN group", "Thermoplasmata"}
FILES = {"s0": "EPOC_data.pangenome_s0.tsv", "s10": "EPOC_data.tsv",
         "s25": "EPOC_data.pangenome_s25.tsv", "s67": "EPOC_data.pangenome_s67.tsv"}


def load(tag):
    d = pd.read_csv(os.path.join(DATA, FILES[tag]), sep="\t", low_memory=False,
                    usecols=K + ["prok_taxa", "prok_clade_size", "c-ELW"])
    d["mx"] = d.groupby(K)["c-ELW"].transform("max")
    return d


def share(d):
    s = d.groupby("prok_taxa")["c-ELW"].sum()
    return 100 * s / s.sum()


s0, s10 = load("s0"), load("s10")

# 1. Reply Table 1: a single species label cannot occur in more EPOCs than its group wins
a0 = s0[s0.prok_taxa == "Asgard"]
wins = s0.sort_values("c-ELW", ascending=False).drop_duplicates(K)
ties = s0[s0["c-ELW"] == s0["mx"]]
print("1. s0 (no filter), all tests")
print(f"   EPOCs: {s0.tree_name.nunique():,}; Asgard a candidate in {a0.tree_name.nunique():,}")
print(f"   Asgard best-supported: {wins[wins.prok_taxa == 'Asgard'].tree_name.nunique():,} EPOCs, "
      f"{(wins.prok_taxa == 'Asgard').sum():,} tests ({(ties.prok_taxa == 'Asgard').sum():,} counting ties)")
print("   Reply Table 1 lists one Asgard label in 838 EPOCs -> it counts candidate clades")

# 2. Reply's unfiltered summary statistics
per = s0.groupby("tree_name").euk_clade_rep.nunique()
best = wins
print("2. Reply statistics reproduced from s0")
print(f"   eukaryotic clades per EPOC {per.mean():.2f} (sd {per.std():.2f}); "
      f"best-sister sequences {best.prok_clade_size.sum():,.0f} "
      f"(Alphaproteobacteria {best[best.prok_taxa == 'Alphaproteobacteria'].prok_clade_size.sum():,.0f}); "
      f"median sister size {best.prok_clade_size.median():.0f}")

# 3. Same gene families with and without the filter
shared = set(s0.tree_name) & set(s10.tree_name)
m0 = share(s0[s0.tree_name.isin(shared) & (s0.mx > CORE)])
m10 = share(s10[s10.tree_name.isin(shared) & (s10.mx > CORE)])
print(f"3. {len(shared):,} gene families in both runs, share of summed core c-ELW (s0 -> s10)")
for g in ["Asgard", "Alphaproteobacteria"]:
    print(f"   {g:20s} {m0[g]:.1f}% -> {m10[g]:.1f}%")
print(f"   Asgard:Alphaproteobacteria {m0['Asgard'] / m0['Alphaproteobacteria']:.2f} -> "
      f"{m10['Asgard'] / m10['Alphaproteobacteria']:.2f}")

# 4. Extended Data Table 1
t = pd.DataFrame({tag: share((d := load(tag))[d.mx > CORE]) for tag in FILES}).fillna(0)
t = t.sort_values("s10", ascending=False)
dom = t.index.map(lambda g: "Asgard" if g == "Asgard" else ("Other archaea" if g in ARCH else "All bacteria"))
tot = t.groupby(dom).sum().loc[["All bacteria", "Other archaea"]]
out = pd.concat([t, tot]).round(1)
out.index.name = "group"
out.to_csv(os.path.join(HERE, "ed_table1_percent_contribution.csv"))
print("4. Extended Data Table 1 (% of summed core c-ELW)")
print(out.to_string())
pro = t.loc[["Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria"], "s10"].sum()
print(f"   s10: Alpha+Beta+Gamma pooled {pro:.1f}% vs Asgard {t.loc['Asgard', 's10']:.1f}%")
