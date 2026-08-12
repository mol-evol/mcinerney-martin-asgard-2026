#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functional (KEGG) enrichment of Alphaproteobacterial winners in the unfiltered (s0) set.

Tests whether the alphaproteobacterial signal that returns when the soft-core filter is
removed is function-agnostic (consistent with recent horizontal transfer) or concentrated
in the genes expected from the mitochondrial endosymbiont (oxidative phosphorylation /
electron transport, mitochondrial ribosomal proteins).

Inputs (Zenodo 15048010): EPOC_data.pangenome_s0.tsv, EPOC_annotation_KEGG.tsv,
KEGG_category_mapping.tsv. Usage: python functional_enrichment.py [DATA_DIR]
Output: s0_alpha_function_enrichment.csv (one row per KEGG pathway, n >= 10).
"""
import os, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(HERE)
CORE = 0.40

x = pd.read_csv(os.path.join(DATA, "EPOC_data.pangenome_s0.tsv"), sep="\t", low_memory=False,
                usecols=["tree_name", "euk_clade_rep", "prok_taxa", "c-ELW"])
x["celw"] = pd.to_numeric(x["c-ELW"], errors="coerce").fillna(0)
x["key"] = list(zip(x["tree_name"], x["euk_clade_rep"]))
core = set(x.groupby("key")["celw"].max().pipe(lambda s: s[s > CORE]).index)
d = x[x["key"].isin(core)]
win = d.loc[d.groupby("key")["celw"].idxmax()].copy()

a = pd.read_csv(os.path.join(DATA, "EPOC_annotation_KEGG.tsv"), sep="\t",
               usecols=["Query", "Target", "Prob"], dtype=str, low_memory=False)
a["Prob"] = pd.to_numeric(a["Prob"], errors="coerce")
# v15: stable sort. 24.5% of queries have their maximum Prob tied across >1 distinct
# target; the default (unstable) sort made this assignment irreproducible.
best = a.sort_values("Prob", kind="mergesort").drop_duplicates("Query", keep="last")
t2k = dict(zip(best["Query"], best["Target"]))
cm = pd.read_csv(os.path.join(DATA, "KEGG_category_mapping.tsv"), sep="\t", dtype=str)
cm = cm[cm["category_id"].astype(str).str.startswith("map")]
k2maps = cm.groupby("kogid")["category_id"].apply(set).to_dict()
id2name = dict(zip(cm["category_id"], cm["category_name"]))

win["maps"] = win["tree_name"].map(t2k).map(lambda k: k2maps.get(k, set()))
win["isAlpha"] = win["prok_taxa"] == "Alphaproteobacteria"
bg = win["isAlpha"].mean()

rows = []
for mp in set().union(*[m for m in win["maps"] if m]):
    sel = win["maps"].map(lambda s: mp in s)
    n = int(sel.sum())
    if n < 10:
        continue
    fa = win.loc[sel, "isAlpha"].mean()
    rows.append(dict(map_id=mp, pathway_name=id2name.get(mp, mp), n_EPOCs=n,
                     n_alpha_winner=int(win.loc[sel, "isAlpha"].sum()),
                     pct_alpha=round(100 * fa, 1), enrichment_vs_background=round(fa / bg, 2)))
t = pd.DataFrame(rows).sort_values("enrichment_vs_background", ascending=False)
out = os.path.join(HERE, "s0_alpha_function_enrichment.csv")
t.to_csv(out, index=False)
print(f"background Alphaproteobacteria-winner fraction at s0: {100*bg:.1f}%  (n_core={len(win)})")
print(f"wrote {out}  ({len(t)} pathways, n>=10)")
print(t.head(12).to_string(index=False))
