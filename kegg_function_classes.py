#!/usr/bin/env python3
"""Functional-class results for Extended Data Fig. 1e,f (revision 2).

Classes follow KEGG's top-level pathway classes (kegg_class_map in asgard_reanalysis.py):
Genetic Information Processing = informational; Metabolism = metabolic; else other.
Writes kegg_stratification.csv (panel e) and matched_family_asgard_share_by_class.csv (panel f).
Usage: python kegg_function_classes.py [DATA_DIR]
"""
import os, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(HERE)
K = ["tree_name", "euk_clade_rep"]

meta = pd.read_csv(os.path.join(DATA, "KEGG_metadata.tsv"), sep="\t", usecols=["kogid", "brite_A"], dtype=str)
cl = lambda a: "informational" if "09120" in str(a).split("|") else ("metabolic" if "09100" in str(a).split("|") else "other")
kc = dict(zip(meta.kogid, meta.brite_A.map(cl)))
ann = pd.read_csv(os.path.join(HERE, "KEGG_annotation_3col.tsv.gz"), sep="\t", dtype={"Query": str, "Target": str})
ann["Prob"] = pd.to_numeric(ann["Prob"], errors="coerce")
best = ann.sort_values("Prob", kind="mergesort").drop_duplicates("Query", keep="last")  # stable tie-break
tclass = dict(zip(best.Query, best.Target.map(kc).fillna("other")))

# panel e: Asgard-winner trees with no bacterial leaf, by class
t = pd.read_csv(os.path.join(HERE, "asgard_winner_tree_analysis.csv"))
t = t[(t.n_asg > 0) & (t.asg_mono.astype(str) != "ERR")].copy()
t["nobact"] = ~t.has_bact.astype(str).isin(["True"]); t["mono"] = t.asg_mono.astype(str).isin(["True"])
t["cls"] = t.tree.map(tclass).fillna("unannotated")
e = t.groupby("cls").apply(lambda s: pd.Series({"n": len(s), "pct_no_bacteria": 100 * s.nobact.mean(),
                                              "pct_polyphyletic": 100 * (~s.mono).mean()}))
e = e.reindex(["informational", "metabolic", "other", "unannotated"]); e.index.name = "kegg_class"
e.to_csv(os.path.join(HERE, "kegg_stratification.csv"))
print("panel e\n", e.round(1))

# panel f: same gene families with and without the filter
def load(f):
    d = pd.read_csv(os.path.join(DATA, f), sep="\t", low_memory=False, usecols=K + ["prok_taxa", "c-ELW"])
    d["mx"] = d.groupby(K)["c-ELW"].transform("max")
    d = d[d.mx > 0.4].copy(); d["cls"] = d.tree_name.map(tclass).fillna("unannotated"); return d
names = lambda f: set(pd.read_csv(os.path.join(DATA, f), sep="\t", usecols=["tree_name"]).tree_name)
shared = names("EPOC_data.pangenome_s0.tsv") & names("EPOC_data.tsv")          # 2,751 families
s0, s10 = load("EPOC_data.pangenome_s0.tsv"), load("EPOC_data.tsv")
s0, s10 = s0[s0.tree_name.isin(shared)], s10[s10.tree_name.isin(shared)]
share = lambda d: 100 * d.loc[d.prok_taxa == "Asgard", "c-ELW"].sum() / d["c-ELW"].sum()
rows = [("all", s0, s10)] + [(c, s0[s0.cls == c], s10[s10.cls == c]) for c in ["informational", "metabolic", "other"]]
f = pd.DataFrame([{"kegg_class": c, "fam_shared": b.tree_name.nunique(), "Asgard%_shared_s0": share(a),
                   "Asgard%_shared_s10": share(b)} for c, a, b in rows]).set_index("kegg_class")
f.round(2).to_csv(os.path.join(HERE, "matched_family_asgard_share_by_class.csv"))
print(f"panel f ({len(shared)} shared families)\n", f.round(1))

# survival of s10 Asgard calls when the filter is removed
w10 = s10.loc[s10.groupby(K)["c-ELW"].idxmax()]; w0 = s0.loc[s0.groupby(K)["c-ELW"].idxmax()]
m = w10.merge(w0[K + ["prok_taxa"]], on=K, suffixes=("_10", "_0")); a = m[m.prok_taxa_10 == "Asgard"]
print("Asgard s10 calls still Asgard at s0\n", a.groupby("cls").apply(lambda x: pd.Series(
    {"n": len(x), "still_Asgard_%": 100 * (x.prok_taxa_0 == "Asgard").mean()})).round(1))
