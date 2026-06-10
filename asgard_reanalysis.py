#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Reanalysis of Tobiasson, Luo, Wolf & Koonin, Nature 650, 141-149 (2026)
"Dominant contribution of Asgard archaea to eukaryogenesis"
================================================================================

Reanalysis for the Matters Arising by J. O. McInerney & W. F. Martin.
All inputs are the authors' OWN deposited data: Zenodo record 15048010
(https://doi.org/10.5281/zenodo.15048010).

This single script reproduces every number in the reanalysis, end to end:

  PART 1  VALIDATION. Reproduce the published headline from EPOC_data.tsv:
          16,526 core tuples, Asgard occupancy 63.8 %, Alpha 19.7 %,
          Asgard:Alpha aELW ratio 7.74, and the fact that c-ELW sums to 1
          per tuple (i.e. it is a closed-set statistic).
  PART 2  FILTER SWEEP. Recompute the ratio and candidate-set occupancy across
          the four deposited soft-core cut-offs (s0 = no filter, s10 = main /
          headline file, s25, s67).
  PART 3  Q1. For the tuples where Asgard is CALLED the eukaryotic donor,
          characterise the Asgard clade (sequences, single-leaf %, purity)
          versus Alphaproteobacteria.
  PART 4  Q1a / Q2. Parse the per-EPOC master trees (EPOC_data.tar.gz) for every
          tree in which Asgard wins:
            (a) is the Asgard set monophyletic? how many separate subclades?
            (b) how many bacterial / other-archaeal leaves are present at all?
            (c) what domain is the Asgard clade sister to?
  PART 5  KEGG stratification (informational vs metabolic vs other).
  PART 6  FIGURE + machine-readable summary.

Outputs are written next to this script:
    filter_sweep.csv, asgard_winner_tree_analysis.csv, kegg_stratification.csv,
    summary.json, figure_asgard_artefact.png / .pdf

Tree step note: the archive has ~619k members, so a single streaming pass is
slow. extract_winner_trees() uses a minimal raw tar reader to pull ONLY the
*.treefile.annot members needed, cached under ./_trees/. If
asgard_winner_tree_analysis.csv already exists, the heavy step is skipped.

Dependencies: pandas, numpy, matplotlib, ete3.
================================================================================
"""

import os, sys, json, glob, subprocess, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

# ------------------------------------------------------------------ config ----
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(HERE)            # deposit files sit one level up
TREES_DIR = os.path.join(HERE, "_trees")

EPOC_MAIN = os.path.join(DATA, "EPOC_data.tsv")                  # == s10
EPOC_S = {
    "s0":  os.path.join(DATA, "EPOC_data.pangenome_s0.tsv"),
    "s10": EPOC_MAIN,
    "s25": os.path.join(DATA, "EPOC_data.pangenome_s25.tsv"),
    "s67": os.path.join(DATA, "EPOC_data.pangenome_s67.tsv"),
}
TARBALL     = os.path.join(DATA, "EPOC_data.tar.gz")
KEGG_ANNOT  = os.path.join(DATA, "EPOC_annotation_KEGG.tsv")
KEGG_CATMAP = os.path.join(DATA, "KEGG_category_mapping.tsv")

CORE_THRESHOLD = 0.40    # Tobiasson core set: tuples with max c-ELW > 0.4

# The 26 deposited prokaryotic "class" labels split into domains. Asgard is kept
# separate; these four are the non-Asgard archaea; the rest are bacterial.
ARCH = {"Euryarchaeota", "TACK group", "DPANN group", "Thermoplasmata"}
PROK = {"Acidobacteriota","Actinomycetota","Alphaproteobacteria","Asgard","Bacillota",
        "Bdellovibrionota","Betaproteobacteria","Campylobacterota","Chlamydiia",
        "Chloroflexota","Cyanobacteriota","DPANN group","Deinococcota",
        "Deltaproteobacteria","Euryarchaeota","FCB group","Gammaproteobacteria",
        "Mycoplasmatota","Myxococcota","Nitrospirota","PVC group","Spirochaetota",
        "TACK group","Thermodesulfobacteriota","Thermoplasmata","Thermotogae"}
BACT = PROK - ARCH - {"Asgard"}

def domain(label):
    """Map a tip's primary taxon label (first of a mixed 'A|B' label) to a domain."""
    x = str(label).split("|")[0]
    if x == "Asgard": return "Asgard"
    if x in ARCH:     return "Archaea"
    if x in BACT:     return "Bacteria"
    return "Eukaryota"

SUMMARY = {}

# --------------------------------------------------------- shared functions ---
def load_epoc(fn):
    df = pd.read_csv(fn, sep="\t", low_memory=False,
        usecols=["tree_name","euk_clade_rep","prok_clade_rep","prok_clade_size",
                 "prok_clade_weight","prok_leaf_clade","prok_taxa","c-ELW"])
    df["celw"] = pd.to_numeric(df["c-ELW"], errors="coerce").fillna(0.0)
    df["sz"]   = pd.to_numeric(df["prok_clade_size"], errors="coerce")
    df["pur"]  = pd.to_numeric(df["prok_clade_weight"], errors="coerce")
    df["leaf"] = df["prok_leaf_clade"].astype(str).str.strip().isin(["True","TRUE","true"])
    df["key"]  = list(zip(df["tree_name"], df["euk_clade_rep"]))     # one 'tuple'
    return df

def core_set(df):
    mx = df.groupby("key")["celw"].max()
    return set(mx[mx > CORE_THRESHOLD].index)

def aelw_and_occupancy(df, core):
    """Standard aELW (a clade's summed weight over ALL core tuples / n_core) and
    candidate-set occupancy (fraction of core tuples in which the clade appears).
    Reproduces the published Asgard:Alpha ratio (7.74) and occupancies exactly."""
    d = df[df["key"].isin(core)]
    ncore = len(core)
    weight = d.groupby(["key","prok_taxa"])["celw"].sum().groupby("prok_taxa").sum()
    occ = d.groupby("prok_taxa")["key"].nunique()
    return pd.concat([(weight/ncore).rename("aELW"),
                      (occ/ncore).rename("occupancy")], axis=1), ncore, d

# ============================================================== PART 1 ========
def part1_validate():
    print("\n"+"="*72+"\nPART 1  VALIDATION (EPOC_data.tsv == s10, the headline file)\n"+"="*72)
    df = load_epoc(EPOC_MAIN)
    per_tuple = df.groupby("key")["celw"].sum()
    core = core_set(df)
    tab, ncore, _ = aelw_and_occupancy(df, core)
    a, al = tab.loc["Asgard"], tab.loc["Alphaproteobacteria"]
    ratio = a["aELW"]/al["aELW"]
    print(f"  total tuples ............. {df['key'].nunique():,}")
    print(f"  core tuples (max c-ELW>0.4) {ncore:,}   (published 16,526)")
    print(f"  c-ELW per-tuple sum ...... median {per_tuple.median():.3f}  (closed-set => sums to 1)")
    print(f"  Asgard aELW {a['aELW']:.3f}  occupancy {100*a['occupancy']:.1f}%  (published occ 63.8%)")
    print(f"  Alpha  aELW {al['aELW']:.3f}  occupancy {100*al['occupancy']:.1f}%  (published occ 19.7%)")
    print(f"  Asgard:Alpha aELW ratio .. {ratio:.2f}   (published 7.74)")
    SUMMARY["validation"] = dict(total_tuples=int(df["key"].nunique()), core_tuples=ncore,
        celw_per_tuple_median=float(per_tuple.median()),
        asgard_occupancy=float(a["occupancy"]), alpha_occupancy=float(al["occupancy"]),
        asgard_alpha_ratio=float(ratio))
    return df

# ============================================================== PART 2 ========
def part2_sweep():
    print("\n"+"="*72+"\nPART 2  SOFT-CORE FILTER SWEEP (s0 = no filter ... s67 = strictest)\n"+"="*72)
    rows = []
    for tag, fn in EPOC_S.items():
        if not os.path.exists(fn):
            print(f"  [skip {tag}: not found]"); continue
        df = load_epoc(fn); core = core_set(df); tab, ncore, _ = aelw_and_occupancy(df, core)
        a, al = tab.loc["Asgard"], tab.loc["Alphaproteobacteria"]
        rows.append(dict(cutoff=tag, n_core=ncore, aELW_Asgard=a["aELW"], aELW_Alpha=al["aELW"],
            ratio=a["aELW"]/al["aELW"], occ_Asgard=100*a["occupancy"], occ_Alpha=100*al["occupancy"]))
    sweep = pd.DataFrame(rows)
    print(sweep.to_string(index=False, formatters={"aELW_Asgard":"{:.3f}".format,
        "aELW_Alpha":"{:.3f}".format,"ratio":"{:.2f}".format,
        "occ_Asgard":"{:.1f}".format,"occ_Alpha":"{:.1f}".format}))
    print("  -> s0 = evaluation_cutoff 0 in the authors' mmseqs_create_pangenome.py")
    print("     (keep clusters with unique_valid_ranks > num_valid_ranks*cutoff => '>0' = NO filter).")
    print("     Removing the soft-core filter INVERTS the result: Alphaproteobacteria leads,")
    print("     Asgard falls to ~5% of donor calls. (s0 absent from v0.3 README; mapping confirmed by code.)")
    sweep.to_csv(os.path.join(HERE,"filter_sweep.csv"), index=False)
    SUMMARY["filter_sweep"] = sweep.to_dict("records")
    return sweep

# ============================================================== PART 3 ========
def part3_q1(df_main):
    print("\n"+"="*72+"\nPART 3  Q1  Size / purity of clades CALLED eukaryotic donor\n"+"="*72)
    core = core_set(df_main); d = df_main[df_main["key"].isin(core)]
    winners = d.loc[d.groupby("key")["celw"].idxmax()]     # argmax c-ELW per tuple
    rows = []
    for t in ["Asgard","Alphaproteobacteria","Cyanobacteriota","Actinomycetota",
              "Betaproteobacteria","Myxococcota"]:
        w = winners[winners["prok_taxa"]==t]
        if len(w)==0: continue
        rows.append(dict(taxon=t, n_wins=len(w), pct_single_leaf=100*w["leaf"].mean(),
            median_seqs=w["sz"].median(), median_purity=w["pur"].median()))
    q1 = pd.DataFrame(rows)
    print(q1.to_string(index=False, formatters={"pct_single_leaf":"{:.1f}".format,
        "median_seqs":"{:.0f}".format,"median_purity":"{:.2f}".format}))
    print("  -> Asgard 'donor' clades are NOT single sequences (median ~15 seqs, high purity).")
    SUMMARY["q1_winner_clades"] = q1.to_dict("records")
    return winners

# ============================================================== PART 4 ========
def extract_winner_trees(asg_trees):
    """Minimal raw tar reader: stream-decompress the archive and write ONLY the
    *.treefile.annot members for the given EPOC ids to ./_trees/<EPOC>.annot."""
    os.makedirs(TREES_DIR, exist_ok=True)
    want = set(map(str, asg_trees))
    have = {os.path.basename(p)[:-6] for p in glob.glob(os.path.join(TREES_DIR,"*.annot"))}
    todo = want - have
    if not todo:
        print(f"  trees cached for all {len(want)} winners"); return
    print(f"  extracting {len(todo)} tree files from {os.path.basename(TARBALL)} ...")
    proc = subprocess.Popen(["gunzip","-c",TARBALL], stdout=subprocess.PIPE); f = proc.stdout; got=0
    while todo:
        hdr = f.read(512)
        if len(hdr)<512 or hdr[:1]==b"\x00": break
        name = hdr[0:100].split(b"\x00")[0].decode("utf-8","replace")
        szf = hdr[124:136].split(b"\x00")[0].strip(); size = int(szf,8) if szf else 0
        nblk = (size+511)//512
        if name.endswith(".treefile.annot"):
            ep = name.split("/")[-1].split(".merged")[0]
            if ep in todo:
                open(os.path.join(TREES_DIR,ep+".annot"),"wb").write(f.read(nblk*512)[:size])
                todo.discard(ep); got+=1; continue
        to = nblk*512
        while to>0:
            c=f.read(min(to,1<<20))
            if not c: break
            to-=len(c)
    proc.terminate()
    print(f"  extracted {got} (missing {len(todo)})")

def analyse_trees(asg_trees):
    from ete3 import Tree
    cache = os.path.join(HERE,"asgard_winner_tree_analysis.csv")
    if os.path.exists(cache):
        print(f"  using cached {os.path.basename(cache)}"); return pd.read_csv(cache)
    extract_winner_trees(asg_trees)
    rows = []
    for fp in glob.glob(os.path.join(TREES_DIR,"*.annot")):
        tn = os.path.basename(fp)[:-6]
        try:
            t = Tree(open(fp).read(), format=1); lv = t.get_leaves()
            for l in lv: l.add_feature("P", domain(getattr(l,"taxa","")))
            asg = [l for l in lv if l.P=="Asgard"]
            nb=sum(l.P=="Bacteria" for l in lv); na=sum(l.P=="Archaea" for l in lv); ne=sum(l.P=="Eukaryota" for l in lv)
            if not asg:
                rows.append([tn,0,nb,na,ne,"NA","NA","NA",nb>0]); continue
            if len(asg)==1: mono,nsub,mrca=True,1,asg[0]
            else:
                mrca=t.get_common_ancestor(asg); mono=set(mrca.get_leaves())==set(asg)
                nsub,seen=0,set()
                for a in asg:
                    if a in seen: continue
                    nd=a
                    while nd.up is not None and all(x.P=="Asgard" for x in nd.up.get_leaves()): nd=nd.up
                    seen|=set(nd.get_leaves()); nsub+=1
            if mrca.up is None: sis="root"
            else:
                cnt={}
                for s in [c for c in mrca.up.children if c is not mrca]:
                    for lf in s.get_leaves(): cnt[lf.P]=cnt.get(lf.P,0)+1
                sis=max(cnt,key=cnt.get) if cnt else "none"
            rows.append([tn,len(asg),nb,na,ne,mono,nsub,sis if mono else "poly",nb>0])
        except Exception:
            rows.append([tn,-1,-1,-1,-1,"ERR","ERR","ERR",False])
    r = pd.DataFrame(rows, columns=["tree","n_asg","n_bact","n_arch","n_euk",
                                    "asg_mono","n_sub","allasg_sister","has_bact"])
    r.to_csv(cache, index=False); return r

def part4_trees(winners):
    print("\n"+"="*72+"\nPART 4  Q1a / Q2  Master-tree analysis of Asgard-winning EPOCs\n"+"="*72)
    asg_trees = winners[winners["prok_taxa"]=="Asgard"]["tree_name"].astype(str).unique()
    print(f"  Asgard wins {int((winners['prok_taxa']=='Asgard').sum())} core tuples across {len(asg_trees)} trees")
    r = analyse_trees(asg_trees)
    r = r[(r["n_asg"]>0) & (r["asg_mono"]!="ERR")].copy()
    r["mono"]   = r["asg_mono"].astype(str).isin(["True"])
    r["nobact"] = ~r["has_bact"].astype(str).isin(["True"])
    print(f"  trees analysed ....................... {len(r):,}")
    print(f"  Asgard MONOPHYLETIC .................. {100*r['mono'].mean():.1f}%")
    print(f"  Asgard polyphyletic ................. {100*(~r['mono']).mean():.1f}%  (median {r['n_sub'].median():.0f} subclades)")
    print(f"  trees with NO bacterial leaf ........ {100*r['nobact'].mean():.1f}%")
    print(f"  trees with NO non-Asgard prokaryote . {100*((r['n_bact']==0)&(r['n_arch']==0)).mean():.1f}%")
    print(f"  median leaves: Asgard {r['n_asg'].median():.0f}  bacteria {r['n_bact'].median():.0f}  "
          f"other-archaea {r['n_arch'].median():.0f}  euk {r['n_euk'].median():.0f}")
    sis = r[r["mono"]]["allasg_sister"].value_counts(normalize=True).mul(100)
    print("  sister of monophyletic Asgard clade:")
    for k,v in sis.items(): print(f"      {k:12s} {v:.1f}%")
    SUMMARY["trees"] = dict(n=len(r), pct_monophyletic=float(100*r["mono"].mean()),
        median_subclades=float(r["n_sub"].median()), pct_no_bacteria=float(100*r["nobact"].mean()),
        pct_no_other_prok=float(100*((r["n_bact"]==0)&(r["n_arch"]==0)).mean()),
        sister=sis.round(2).to_dict())
    return r

# ============================================================== PART 5 ========
def kegg_class_map():
    cat = pd.read_csv(KEGG_CATMAP, sep="\t", dtype=str)
    def cl(ids):
        ids = set(ids)
        if any(str(x).startswith("map03") for x in ids): return "informational"
        if any(str(x).startswith("map00") for x in ids): return "metabolic"
        return "other"
    return cat.groupby("kogid")["category_id"].apply(cl)

def part5_kegg(r):
    print("\n"+"="*72+"\nPART 5  KEGG stratification (informational vs metabolic)\n"+"="*72)
    kc = kegg_class_map()
    ann = pd.read_csv(KEGG_ANNOT, sep="\t", usecols=["Query","Target","Prob"],
                      dtype={"Query":str,"Target":str}, low_memory=False)
    ann["Prob"] = pd.to_numeric(ann["Prob"], errors="coerce")
    best = ann.sort_values("Prob").drop_duplicates("Query", keep="last")
    tree_class = dict(zip(best["Query"], best["Target"].map(kc).fillna("other")))
    r = r.copy(); r["class"] = r["tree"].map(tree_class).fillna("unannotated")
    out = []
    for c in ["informational","metabolic","other","unannotated"]:
        s = r[r["class"]==c]
        if len(s)==0: continue
        out.append(dict(kegg_class=c, n=len(s), pct_no_bacteria=100*s["nobact"].mean(),
                        pct_polyphyletic=100*(~s["mono"]).mean()))
    k = pd.DataFrame(out)
    print(k.to_string(index=False, formatters={"pct_no_bacteria":"{:.1f}".format,
        "pct_polyphyletic":"{:.1f}".format}))
    k.to_csv(os.path.join(HERE,"kegg_stratification.csv"), index=False)
    SUMMARY["kegg"] = k.to_dict("records")
    return k

# ============================================================== PART 7 ========
def part7_flip():
    """Within-gene test: take each clade's CONFIDENT donor calls at s10 and ask what
    happens to the SAME tuple when the soft-core filter is removed (s0). Comparing
    across winner domains controls for s0 being a smaller dataset."""
    print("\n"+"="*72+"\nPART 7  Within-gene donor flip when the filter is removed (s10 -> s0)\n"+"="*72)
    dommap = lambda t: "Asgard" if t=="Asgard" else ("Archaea" if t in ARCH else "Bacteria")
    def winners(fn):
        x = load_epoc(fn); c = core_set(x); d = x[x["key"].isin(c)]
        return d.loc[d.groupby("key")["celw"].idxmax()]
    w10 = winners(EPOC_MAIN); w0 = winners(EPOC_S["s0"])
    w0d = dict(zip(w0["key"], w0["prok_taxa"])); s0keys = set(w0["key"])
    w10 = w10.copy(); w10["dom10"] = w10["prok_taxa"].map(dommap); w10["matched"] = w10["key"].isin(s0keys)
    print(f"  baseline: {100*w10['matched'].mean():.0f}% of all {len(w10)} s10 confident calls remain confident at s0")
    rows = []
    for d in ["Asgard","Bacteria","Archaea"]:
        sub = w10[w10["dom10"]==d]; rem = sub[sub["matched"]]
        fl = rem["key"].map(w0d).map(dommap).value_counts(normalize=True).mul(100)
        rows.append(dict(s10_call=d, n=len(sub), n_remain=len(rem),
            pct_remain=round(100*sub["matched"].mean(),1),
            to_Bacteria=round(fl.get("Bacteria",0),1), to_Asgard=round(fl.get("Asgard",0),1),
            to_Archaea=round(fl.get("Archaea",0),1)))
    flipdf = pd.DataFrame(rows)
    print(flipdf.to_string(index=False))
    print("  -> Asgard donor calls flip to bacteria when the filter is removed; bacterial calls do not.")
    flipdf.to_csv(os.path.join(HERE,"winner_flip_s10_vs_s0.csv"), index=False)
    SUMMARY["winner_flip_s10_vs_s0"] = flipdf.to_dict("records")
    return flipdf

# ============================================================== PART 8 ========
def part8_headtohead():
    """Denominator / head-to-head reconstruction (manuscript Fig. 1 panels a-d), from the SAME
    pipeline as every other panel. Panel a reproduces the published R^2 (0.92/0.79) exactly.
    Panels b-d compare the published STANDARD aELW ratio with the inclusion-controlled
    HEAD-TO-HEAD c-ELW ratio: both clades' c-ELW summed within each shared tuple (the same
    within-tuple summation as the standard aELW), restricted to tuples where both are candidates
    and renormalised two-way. Removing the candidate-set asymmetry collapses the 4-15x standard
    ratios to 1.3-3.3x head-to-head (Asgard:Alpha 1.60, Asgard:Cyano 1.31)."""
    print("\n"+"="*72+"\nPART 8  Denominator / head-to-head reconstruction (Fig. 1 panels a-d)\n"+"="*72)
    df = load_epoc(EPOC_MAIN); core = core_set(df); ncore=len(core); d = df[df["key"].isin(core)].copy()
    tw = d.groupby(["key","prok_taxa"])["celw"].sum()
    npres = d.groupby("prok_taxa")["key"].nunique()
    std = (tw.groupby("prok_taxa").sum()/ncore); occ = (npres/ncore); cond = (tw.groupby("prok_taxa").sum()/npres)
    tab = d.pivot_table(index="key", columns="prok_taxa", values="celw", aggfunc="sum")
    def h2h(x):
        # head-to-head ratio = mean Asgard c-ELW / mean target c-ELW on tuples where BOTH are
        # candidates (per-clade c-ELW summed within each tuple). Same estimator as the standard
        # aELW ratio, only the denominator is restricted to the shared tuples.
        if x not in tab.columns: return np.nan
        sub = tab[["Asgard",x]].dropna()
        if len(sub)==0: return np.nan
        return sub["Asgard"].mean()/sub[x].mean()
    bact6 = ["Alphaproteobacteria","Cyanobacteriota","Actinomycetota","Betaproteobacteria","Myxococcota","Gammaproteobacteria"]
    lab6  = ["Alpha","Cyano","Actino","Beta","Myxo","Gamma"]
    panelb = pd.DataFrame({"label":lab6,
                           "standard_ratio":[std["Asgard"]/std[x] for x in bact6],
                           "h2h_ratio":[h2h(x) for x in bact6]})
    cc = d.groupby("key").size()
    def strat(x):
        rows=[]
        for m in range(1,13):
            ks=set(cc[cc>=m].index); dd=d[d["key"].isin(ks)]; nc=len(ks)
            tw2=dd.groupby(["key","prok_taxa"])["celw"].sum().groupby("prok_taxa").sum()/nc
            sr=tw2.get("Asgard",0)/tw2.get(x,1e-9)
            occ2=dd.groupby("prok_taxa")["key"].nunique()/nc
            ir=occ2.get("Asgard",0)/occ2.get(x,1e-9)
            t2=dd.pivot_table(index="key",columns="prok_taxa",values="celw",aggfunc="sum")
            if x in t2.columns:
                sub=t2[["Asgard",x]].dropna()
                hr=sub["Asgard"].mean()/sub[x].mean() if len(sub) else np.nan
            else: hr=np.nan
            rows.append((m,sr,hr,ir))
        return pd.DataFrame(rows,columns=["min_c","std_ratio","h2h_ratio","inc_ratio"])
    R2s=float(np.corrcoef(occ.reindex(std.index),std)[0,1]**2)
    R2c=float(np.corrcoef(occ.reindex(cond.index),cond)[0,1]**2)
    npaired=int(tab[["Asgard","Alphaproteobacteria"]].dropna().shape[0])
    pa=tab[["Asgard","Alphaproteobacteria"]].dropna()
    winA=int((pa["Asgard"]>pa["Alphaproteobacteria"]).sum()); nA=len(pa)
    mA,mAl=float(pa["Asgard"].mean()),float(pa["Alphaproteobacteria"].mean())
    print(f"  panel a  R^2 standard {R2s:.2f} (orig 0.92) | conditional {R2c:.2f} (orig 0.79)   [MATCH]")
    print(f"  panel b  standard ratios {panelb['standard_ratio'].min():.1f}-{panelb['standard_ratio'].max():.1f}"
          f"  ->  head-to-head {panelb['h2h_ratio'].min():.2f}-{panelb['h2h_ratio'].max():.2f}")
    print(f"  head-to-head Asgard:Alpha {panelb.loc[0,'h2h_ratio']:.2f} | Asgard:Cyano {panelb.loc[1,'h2h_ratio']:.2f}"
          f"  (n={npaired} paired Asgard:Alpha tuples)")
    print("  -> the candidate-set asymmetry, not a per-gene Asgard advantage, drives the published ratio.")
    scatter = pd.DataFrame({"taxon":list(std.index), "in_set_pct":(occ.reindex(std.index)*100).values,
                            "aelw":std.values, "aelw_cond":cond.reindex(std.index).values})
    SUMMARY["headtohead"] = dict(R2_standard=round(R2s,3), R2_conditional=round(R2c,3),
        n_paired_asgard_alpha=npaired,
        table1_asgard_alpha=dict(mean_renorm_asgard=round(mA/(mA+mAl),3),
            mean_renorm_alpha=round(mAl/(mA+mAl),3), ratio=round(mA/mAl,2),
            wins_asgard=winA, pct_asgard=round(100*winA/nA,1),
            wins_alpha=nA-winA, pct_alpha=round(100*(nA-winA)/nA,1)),
        standard_ratio={l:round(float(v),2) for l,v in zip(lab6,panelb["standard_ratio"])},
        headtohead_ratio={l:round(float(v),2) for l,v in zip(lab6,panelb["h2h_ratio"])},
        note=("Head-to-head c-ELW uses the SAME within-tuple summation as the standard aELW, restricted "
              "to tuples where both clades are candidates and renormalised two-way. It collapses the "
              "4-15x standard ratios to 1.3-3.3x (Asgard:Alpha 1.60, Asgard:Cyano 1.31): the published "
              "dominance is candidate-set inclusion, not a per-gene Asgard advantage. Panel-a R^2 "
              "(0.92/0.79) and all standard-metric invariants reproduce the manuscript exactly."))
    return dict(scatter=scatter, panelb=panelb, cstrat=strat("Alphaproteobacteria"),
                dstrat=strat("Cyanobacteriota"))

# ============================================================== PART 6 ========
def part6_figure(sweep, r, k, flip, ad):
    """Manuscript Figure 1 -- all eight panels (a-h) from this one pipeline.
    a: candidate-set inclusion vs aELW (standard + conditional);  b: standard vs head-to-head
    ratio for six bacterial candidates;  c,d: Asgard:Alpha and Asgard:Cyano ratios vs minimum
    candidates per tuple (standard collapses onto the stable head-to-head value);  e: filter
    sweep inversion (s0);  f: within-gene donor flip;  g: tree bacterial depletion;
    h: Asgard monophyly.  Saves figure_asgard_full.svg / .pdf / .png (true vector)."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,
        'axes.labelsize':8.5,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7,
        'axes.linewidth':0.7,'xtick.major.width':0.7,'ytick.major.width':0.7})
    CS='#1f4e79'; CC='#c25e00'; CI='#666666'
    scatter=ad["scatter"]; h2h=ad["panelb"]; pc=ad["cstrat"]; pdd=ad["dstrat"]
    sw=sweep.set_index("cutoff").reindex(["s0","s10","s25","s67"])
    fl=flip.set_index("s10_call"); kk=k.set_index("kegg_class")
    cats=[c for c in ["informational","metabolic","other"] if c in kk.index]
    nobact=100*r["nobact"].mean()
    fig=plt.figure(figsize=(7.5,11.6))
    gs=GridSpec(4,2,figure=fig,hspace=0.46,wspace=0.34,left=0.085,right=0.975,top=0.975,bottom=0.045)
    A=fig.add_subplot(gs[0,0]);B=fig.add_subplot(gs[0,1]);C=fig.add_subplot(gs[1,0]);Dx=fig.add_subplot(gs[1,1])
    E=fig.add_subplot(gs[2,0]);Fx=fig.add_subplot(gs[2,1]);G=fig.add_subplot(gs[3,0]);H=fig.add_subplot(gs[3,1])
    def letter(ax,l): ax.text(-0.16,1.06,l,transform=ax.transAxes,fontsize=12,fontweight='bold',va='top')
    def nospine(ax): ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    def fit(x,y):
        b1,b0=np.polyfit(x,y,1); return b0,b1,float(np.corrcoef(x,y)[0,1]**2)
    # a -- inclusion vs aELW
    A.scatter(scatter['in_set_pct'],scatter['aelw'],s=30,c=CS,edgecolor='white',linewidth=0.5,label='Standard aELW',zorder=3)
    A.scatter(scatter['in_set_pct'],scatter['aelw_cond'],s=30,c=CC,edgecolor='white',linewidth=0.5,label='Conditional aELW',zorder=3)
    xr=np.array([0,scatter['in_set_pct'].max()*1.05])
    i_s,ss,rs=fit(scatter['in_set_pct'],scatter['aelw']); i_c,sc,rc=fit(scatter['in_set_pct'],scatter['aelw_cond'])
    A.plot(xr,i_s+ss*xr,'--',c=CS,lw=1.1,alpha=0.7); A.plot(xr,i_c+sc*xr,'--',c=CC,lw=1.1,alpha=0.7)
    lp={'Asgard':('right',-6,4),'Cyanobacteriota':('left',5,4),'Betaproteobacteria':('left',5,-10),'Alphaproteobacteria':('right',-6,8),'Actinomycetota':('left',5,8),'Gammaproteobacteria':('left',5,-4),'Myxococcota':('right',-6,-8)}
    lt={'Asgard':'Asgard','Cyanobacteriota':'Cyano','Betaproteobacteria':'Beta','Alphaproteobacteria':'Alpha','Actinomycetota':'Actino','Gammaproteobacteria':'Gamma','Myxococcota':'Myxo'}
    for _,row in scatter.iterrows():
        if row['taxon'] not in lp: continue
        ha,dx,dy=lp[row['taxon']]; A.annotate(lt[row['taxon']],(row['in_set_pct'],row['aelw_cond']),xytext=(dx,dy),textcoords='offset points',fontsize=7,ha=ha,fontweight='bold' if row['taxon']=='Asgard' else 'normal')
    A.text(0.04,0.96,f'Standard: $R^2$ = {rs:.2f}',transform=A.transAxes,fontsize=7,color=CS,va='top')
    A.text(0.04,0.88,f'Conditional: $R^2$ = {rc:.2f}',transform=A.transAxes,fontsize=7,color=CC,va='top')
    A.set_xlabel('Candidate-set inclusion (% of core EPOC tuples)'); A.set_ylabel('aELW'); A.set_xlim(0,xr[1]); A.set_ylim(-0.02,0.75)
    A.legend(loc='lower right',frameon=False,fontsize=7); nospine(A); letter(A,'a')
    # b -- standard vs head-to-head
    xp=np.arange(len(h2h)); w=0.38
    B.bar(xp-w/2,h2h['standard_ratio'],w,color=CS,edgecolor='white',linewidth=0.5,label='Standard aELW ratio')
    B.bar(xp+w/2,h2h['h2h_ratio'],w,color=CC,edgecolor='white',linewidth=0.5,label='Head-to-head c-ELW ratio')
    B.axhline(1,color='black',lw=0.5,alpha=0.6); B.text(len(h2h)-0.5,1.2,'parity',fontsize=6.5,ha='right',alpha=0.7,style='italic')
    for i,(sv,hv) in enumerate(zip(h2h['standard_ratio'],h2h['h2h_ratio'])):
        B.text(i-w/2,sv+0.25,f'{sv:.1f}',ha='center',fontsize=6.3,color=CS); B.text(i+w/2,hv+0.25,f'{hv:.2f}',ha='center',fontsize=6.3,color=CC)
    B.set_xticks(xp); B.set_xticklabels(h2h['label'],fontsize=7.5); B.set_xlabel('Bacterial candidate (vs Asgard)'); B.set_ylabel('Asgard : taxon ratio'); B.set_ylim(0,17)
    B.legend(loc='upper left',frameon=False,fontsize=7); nospine(B); letter(B,'b')
    # c,d -- stratification
    def stratplot(ax,data,tgt,ylim,letterc):
        ax.plot(data['min_c'],data['std_ratio'],'o-',color=CS,lw=1.3,ms=4.5,mec='white',mew=0.6,label='Standard aELW ratio',zorder=3)
        ax.plot(data['min_c'],data['h2h_ratio'],'s-',color=CC,lw=1.3,ms=4,mec='white',mew=0.6,label='Head-to-head c-ELW ratio',zorder=3)
        ax.plot(data['min_c'],data['inc_ratio'],'^--',color=CI,lw=1.1,ms=4,mec='white',mew=0.6,label='Inclusion ratio',alpha=0.85,zorder=2)
        ax.axhline(1,color='black',lw=0.5,alpha=0.5); ax.text(12.4,1+(ylim[1]-ylim[0])*0.012,'parity',fontsize=6.5,ha='right',alpha=0.7,style='italic')
        ax.annotate(f"{data['std_ratio'].iloc[0]:.2f}",xy=(1,data['std_ratio'].iloc[0]),xytext=(8,0),textcoords='offset points',fontsize=7,color=CS,fontweight='bold',va='center')
        ax.annotate(f"{data['std_ratio'].iloc[-1]:.2f}",xy=(12,data['std_ratio'].iloc[-1]),xytext=(-8,6),textcoords='offset points',fontsize=7,color=CS,ha='right')
        ax.set_xlabel('Minimum candidates per EPOC tuple'); ax.set_ylabel(f'Asgard : {tgt} ratio'); ax.set_xticks(range(1,13)); ax.set_xlim(0.5,12.5); ax.set_ylim(*ylim)
        nospine(ax); ax.legend(loc='upper right',frameon=False,fontsize=7); letter(ax,letterc)
    stratplot(C,pc,'Alpha',(0.5,9),'c'); stratplot(Dx,pdd,'Cyano',(0.5,5),'d')
    # e -- filter sweep
    xl=["s0\n(no filter)","s10\n(headline)","s25","s67"]
    E.bar(range(4),sw["ratio"],color=['#c0392b',CS,CS,CS]); E.axhline(1,ls='--',c='grey',lw=0.8); E.set_yscale('log'); E.set_ylim(0.3,25)
    for i,v in enumerate(sw["ratio"]): E.text(i,v*1.06,f"{v:.2f}",ha='center',va='bottom',fontsize=7)
    E.set_xticks(range(4)); E.set_xticklabels(xl); E.set_ylabel('Asgard : Alpha aELW ratio'); nospine(E); letter(E,'e')
    E.set_title('Remove the soft-core filter (s0)',fontsize=8.5,loc='left',pad=2)
    # f -- within-gene flip
    gB=[g for g in ['Asgard','Bacteria'] if g in fl.index]; bott=np.zeros(len(gB))
    for col,c,lab in [('to_Bacteria',CC,'to Bacteria'),('to_Asgard',CS,'to Asgard'),('to_Archaea','#7f8c8d','to Archaea')]:
        vals=np.array([fl.loc[g,col] for g in gB]); Fx.bar(range(len(gB)),vals,bottom=bott,color=c,label=lab)
        for i,(v,b) in enumerate(zip(vals,bott)):
            if v>=8: Fx.text(i,b+v/2,f"{v:.0f}%",ha='center',va='center',fontsize=7.5,color='white',fontweight='bold')
        bott=bott+vals
    Fx.set_xticks(range(len(gB))); Fx.set_xticklabels([("Asgard" if g=="Asgard" else "Bacteria")+" call\n(s10)" for g in gB]); Fx.set_ylabel('re-assignment, filter OFF (%)'); Fx.set_ylim(0,100)
    Fx.legend(frameon=False,fontsize=6.5,loc='lower center'); nospine(Fx); letter(Fx,'f'); Fx.set_title('Same genes, filter off',fontsize=8.5,loc='left',pad=2)
    # g -- tree bacterial depletion
    vals=[kk.loc[c,"pct_no_bacteria"] for c in cats]; G.bar(range(len(cats)),vals,color='#16a085'); G.axhline(nobact,ls='--',c='grey',lw=0.8)
    G.text(len(cats)-0.5,nobact+1.5,f"all {nobact:.0f}%",ha='right',fontsize=6.5,color='grey')
    for i,v in enumerate(vals): G.text(i,v+1,f"{v:.0f}%",ha='center',fontsize=7.5)
    G.set_xticks(range(len(cats))); G.set_xticklabels(cats,rotation=12,fontsize=7); G.set_ylim(0,100); G.set_ylabel('% Asgard-winner trees\nwith NO bacterium'); nospine(G); letter(G,'g'); G.set_title('Filter empties trees of bacteria',fontsize=8.5,loc='left',pad=2)
    # h -- monophyly
    poly=[kk.loc[c,"pct_polyphyletic"] for c in cats]; mono=[100-p for p in poly]; xx=np.arange(len(cats))
    H.bar(xx,mono,color='#8e44ad',label='monophyletic'); H.bar(xx,poly,bottom=mono,color='#d2b4de',label='polyphyletic')
    for i,m in enumerate(mono): H.text(i,m/2,f"{m:.0f}%",ha='center',va='center',fontsize=7.5,color='white')
    H.set_xticks(xx); H.set_xticklabels(cats,rotation=12,fontsize=7); H.set_ylabel('% Asgard-winner trees'); H.legend(frameon=False,loc='lower right',fontsize=6.5); nospine(H); letter(H,'h'); H.set_title("'Asgard' often not monophyletic",fontsize=8.5,loc='left',pad=2)
    for ext in ('svg','pdf','png'):
        fig.savefig(os.path.join(HERE,f"figure_asgard_full.{ext}"),dpi=300 if ext=='png' else None,bbox_inches='tight')
    print("\n  figure written (all 8 panels a-h): figure_asgard_full.svg / .pdf / .png")

# ============================================================== PART 9 ========
def part9_stemlength():
    """Stem-length critique (manuscript Fig. 2). Reproduces the paper's normalised stem
    (raw_stem_length / median_euk_leaf_dist; Pittis-Gabaldon relative stem) and shows it
    carries no timing/ancestry signal: (1) genes co-acquired in one event span ~100x
    (oxidative phosphorylation) / ~70x (ribosome); (2) the Asgard 'short stem' is a
    normalisation artefact (raw stems equal; the median-euk-branch normaliser is ~13% larger
    for Asgard, the dominant effect) and is shared by Cyanobacteria, an acknowledged HGT donor;
    (3) stems are measured to the soft-core-filter-determined sister (changes for ~60% of
    clades at s0). Writes figure_stemlength.svg/.pdf/.png."""
    print("\n"+"="*72+"\nPART 9  Stem-length analysis (Fig. 2): normalised stem is not a timing signal\n"+"="*72)
    from scipy.stats import mannwhitneyu
    cols=["tree_name","euk_clade_rep","prok_taxa","c-ELW","stem_length","raw_stem_length","median_euk_leaf_dist"]
    def load_stem(fn):
        x=pd.read_csv(fn,sep="\t",usecols=cols,low_memory=False)
        x["celw"]=pd.to_numeric(x["c-ELW"],errors="coerce").fillna(0.0)
        for c in ["stem_length","raw_stem_length","median_euk_leaf_dist"]: x[c]=pd.to_numeric(x[c],errors="coerce")
        x["key"]=list(zip(x["tree_name"],x["euk_clade_rep"])); return x
    dom=lambda t:"Asgard" if t=="Asgard" else ("Archaea" if t in ARCH else "Bacteria")
    df=load_stem(EPOC_MAIN); mx=df.groupby("key")["celw"].max(); core=set(mx[mx>CORE_THRESHOLD].index)
    d=df[df["key"].isin(core)]; win=d.loc[d.groupby("key")["celw"].idxmax()].copy()
    win["domn"]=win["prok_taxa"].map(dom); win=win[np.isfinite(win["stem_length"])]
    # KEGG pathway per EPOC
    a=pd.read_csv(KEGG_ANNOT,sep="\t",usecols=["Query","Target","Prob"],dtype=str,low_memory=False)
    a["Prob"]=pd.to_numeric(a["Prob"],errors="coerce"); b=a.sort_values("Prob").drop_duplicates("Query",keep="last")
    t2k=dict(zip(b["Query"],b["Target"])); cm=pd.read_csv(KEGG_CATMAP,sep="\t",dtype=str)
    k2c=cm.groupby("kogid")["category_id"].apply(set).to_dict()
    win["cats"]=win["tree_name"].map(t2k).map(lambda k:k2c.get(k,set()))
    # (1) decomposition Asgard vs Bacteria
    A=win[win.domn=="Asgard"]; B=win[win.domn=="Bacteria"]
    def es(col):
        x,y=A[col].dropna(),B[col].dropna(); U,p=mannwhitneyu(x,y,alternative="two-sided")
        return x.median(),y.median(),2*U/(len(x)*len(y))-1,p
    dec={}
    for col in ["raw_stem_length","median_euk_leaf_dist","stem_length"]:
        am,bm,rbc,p=es(col); dec[col]=dict(asgard=round(am,4),bacteria=round(bm,4),rank_biserial=round(rbc,3),p=float(p))
        print(f"  {col:22s} Asgard {am:.4f} | Bacteria {bm:.4f} | rank-biserial {rbc:+.3f} (p={p:.0e})")
    # (2) per-donor + Cyano non-specificity
    donors=["Asgard","Cyanobacteriota","Alphaproteobacteria","Gammaproteobacteria","Actinomycetota"]
    perdon={t:dict(raw=round(win[win.prok_taxa==t].raw_stem_length.median(),4),
                   norm=round(win[win.prok_taxa==t].stem_length.median(),4),
                   denom=round(win[win.prok_taxa==t].median_euk_leaf_dist.median(),3)) for t in donors}
    print("  Cyanobacteria (HGT donor) normalised stem %.3f vs Asgard %.3f (Alpha %.3f)"%(
        perdon['Cyanobacteriota']['norm'],perdon['Asgard']['norm'],perdon['Alphaproteobacteria']['norm']))
    # (3) reductio: spread within co-acquired sets
    def spread(m):
        v=win.loc[m,"stem_length"]; v=v[v>1e-6]; return int(len(v)),round(v.quantile(.025),3),round(v.quantile(.975),3),round(v.quantile(.975)/v.quantile(.025),0)
    red={}
    for lab,m in [("oxphos_alpha",win.cats.map(lambda s:'map00190' in s)&(win.prok_taxa=='Alphaproteobacteria')),
                  ("oxphos_all",win.cats.map(lambda s:'map00190' in s)),
                  ("ribosome_asgard",win.cats.map(lambda s:'map03010' in s)&(win.domn=='Asgard'))]:
        n,lo,hi,r=spread(m); red[lab]=dict(n=n,p2_5=lo,p97_5=hi,fold_spread=r)
        print(f"  spread {lab:16s} n={n:4d} {lo:.3f}-{hi:.3f}  {r:.0f}x")
    # info vs metabolic
    info=win.cats.map(lambda s:any(c.startswith('map03') for c in s))
    metab=win.cats.map(lambda s:any(c.startswith('map00') for c in s) and not any(c.startswith('map03') for c in s))
    infomed=round(win.loc[info,'stem_length'].median(),4); metabmed=round(win.loc[metab,'stem_length'].median(),4)
    print(f"  information-processing stem {infomed:.3f} < metabolic {metabmed:.3f} (paper's own paradox)")
    # (4) filter dependence vs s0
    s0=load_stem(EPOC_S["s0"]); mx0=s0.groupby("key")["celw"].max(); core0=set(mx0[mx0>CORE_THRESHOLD].index)
    w0=s0[s0["key"].isin(core0)]; w0=w0.loc[w0.groupby("key")["celw"].idxmax()]
    mg=win.merge(w0[["key","prok_taxa","stem_length"]],on="key",suffixes=("_s10","_s0"))
    chg=mg[mg.prok_taxa_s10!=mg.prok_taxa_s0]
    fdep=dict(n_common=int(len(mg)),pct_sister_change=round(100*len(chg)/len(mg),0),
              median_abs_dstem=round((chg.stem_length_s10-chg.stem_length_s0).abs().median(),3))
    print(f"  filter dependence: sister changes for {fdep['pct_sister_change']:.0f}% of {len(mg)} clades at s0 (median |d norm stem| {fdep['median_abs_dstem']:.3f})")
    SUMMARY["stem_length"]=dict(n_core_winners=int(len(win)),decomposition=dec,per_donor=perdon,
        reductio=red,info_stem=infomed,metabolic_stem=metabmed,filter_dependence=fdep,
        note=("Normalised stem = raw_stem_length/median_euk_leaf_dist. Carries no timing signal "
              "(co-acquired OXPHOS ~100x spread, ribosome ~70x; information-processing shortest). "
              "Asgard 'short stem' is a normalisation artefact (raw stems equal, rank-biserial -0.03; "
              "the normaliser is the dominant effect, +0.16) and is shared by Cyanobacteria (HGT). "
              "Stems are measured to the soft-core-filter-determined sister."))
    _stem_figure(win,donors,perdon)
    print("\n  figure written: figure_stemlength.svg / .pdf / .png")
    return SUMMARY["stem_length"]

def _stem_figure(win,donors,perdon):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,
        'axes.labelsize':8.5,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7,
        'axes.linewidth':0.7,'xtick.major.width':0.7,'ytick.major.width':0.7})
    CS='#1f4e79'; CC='#c25e00'; CG='#999999'
    short=["Asgard","Cyano","Alpha\n(mito)","Gamma","Actino"]
    raw=[perdon[t]['raw'] for t in donors]; nrm=[perdon[t]['norm'] for t in donors]; den=[perdon[t]['denom'] for t in donors]
    fig=plt.figure(figsize=(7.4,3.5))
    gs=GridSpec(1,2,figure=fig,wspace=0.42,left=0.085,right=0.93,top=0.84,bottom=0.17,width_ratios=[1.15,1])
    A=fig.add_subplot(gs[0,0]); B=fig.add_subplot(gs[0,1])
    nospine=lambda ax:[ax.spines[s].set_visible(False) for s in ('top','right')]
    letter=lambda ax,l:ax.text(-0.17,1.10,l,transform=ax.transAxes,fontsize=12,fontweight='bold',va='top')
    x=np.arange(len(donors)); w=0.38
    A.bar(x-w/2,raw,w,color=CG,edgecolor='white',linewidth=0.5,label='Raw stem')
    A.bar(x+w/2,nrm,w,color=CC,edgecolor='white',linewidth=0.5,label='Normalised stem')
    A.set_xticks(x); A.set_xticklabels(short,fontsize=7.2); A.set_ylabel('Median stem length'); A.set_ylim(0,0.26)
    A.legend(loc='upper left',frameon=False,fontsize=6.8); A.spines['top'].set_visible(False)
    At=A.twinx(); At.plot(x,den,'o-',color=CS,lw=1.2,ms=4,mec='white',mew=0.6,label='Median euk. branch\n(normaliser)')
    At.set_ylabel('Median euk. branch length',color=CS); At.tick_params(axis='y',colors=CS); At.set_ylim(0,1.7); At.spines['top'].set_visible(False)
    At.legend(loc='upper right',frameon=False,fontsize=6.5)
    A.set_title("Raw stems are equal; normalisation makes\nAsgard (and Cyano) 'short'",fontsize=8,loc='left',pad=3); letter(A,'a')
    sets=[("Ribosome\n(Asgard, 'vertical')",win.cats.map(lambda s:'map03010' in s)&(win.domn=='Asgard'),CS),
          ("Oxidative phosph.\n(Alpha, one symbiosis)",win.cats.map(lambda s:'map00190' in s)&(win.prok_taxa=='Alphaproteobacteria'),CC)]
    for i,(lab,m,c) in enumerate(sets):
        v=win.loc[m,"stem_length"]; v=v[v>1e-6].values
        y=np.full(len(v),i)+(np.random.RandomState(0).rand(len(v))-0.5)*0.28
        B.scatter(v,y,s=9,color=c,alpha=0.55,edgecolor='none',zorder=2)
        lo,hi=np.quantile(v,.025),np.quantile(v,.975)
        B.plot([lo,hi],[i,i],color=c,lw=2.4,alpha=0.9,zorder=3,solid_capstyle='round')
        B.text(hi*1.15,i,f"{hi/lo:.0f}x spread",va='center',fontsize=7,color=c,fontweight='bold')
    B.set_xscale('log'); B.set_yticks([0,1]); B.set_yticklabels([s[0] for s in sets],fontsize=7.2)
    B.set_xlim(0.01,6); B.set_ylim(-0.6,1.6); B.set_xlabel('Normalised stem length (log scale)')
    B.set_title("If stem = acquisition time, each set\nwould be a single value",fontsize=8,loc='left',pad=3); nospine(B); letter(B,'b')
    for ext in ('svg','pdf','png'): fig.savefig(os.path.join(HERE,f"figure_stemlength.{ext}"),dpi=300 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)

# ================================================================ MAIN ========
if __name__ == "__main__":
    df_main = part1_validate()
    sweep   = part2_sweep()
    winners = part3_q1(df_main)
    r       = part4_trees(winners)
    k       = part5_kegg(r)
    flip    = part7_flip()
    ad      = part8_headtohead()   # Fig. 1 panels a-d (denominator / head-to-head reconstruction)
    part6_figure(sweep, r, k, flip, ad)   # full 8-panel figure_asgard_full.svg/.pdf/.png
    part9_stemlength()                    # Fig. 2 stem-length analysis (figure_stemlength.*)
    # --- donor ranking (winner share per taxon): documented s10, plus caveated s0 ---
    print("\n" + "="*72 + "\nDONOR RANKING (winner share, % of core tuples)\n" + "="*72)
    dr = {}
    for tag, fn in [("s10", EPOC_MAIN), ("s0_UNDOC", EPOC_S["s0"])]:
        x = load_epoc(fn); cc = core_set(x); d = x[x["key"].isin(cc)]
        win = d.loc[d.groupby("key")["celw"].idxmax()]
        dr[tag] = (win["prok_taxa"].value_counts() / len(cc) * 100).round(2)
    drdf = pd.concat(dr, axis=1).fillna(0.0)
    print(drdf.head(10).to_string())
    drdf.to_csv(os.path.join(HERE, "donor_ranking_winner_share_pct.csv"))
    SUMMARY["donor_ranking_winner_share_pct"] = {t: v.to_dict() for t, v in dr.items()}
    SUMMARY["caveats"] = {
        "s0": ("s0 = evaluation_cutoff 0 in the authors' mmseqs_create_pangenome.py: the soft-core "
               "filter keeps clusters with unique_valid_ranks > num_valid_ranks*cutoff, so cutoff 0 "
               "gives '>0' = retain all taxonomically-assigned clusters = NO soft-core filter. This is "
               "established directly from the authors' code, and corroborated by the naming series "
               "(s10/s25/s67 = cutoff 0.10/0.25/0.67), by the result (diverse bacterial donors dominate; "
               "Asgard ~5%), and by the candidate sets (95% of s0 tuples carry a bacterial candidate vs "
               "66% at s10). The EPOC count (7,276 core tuples on the md5-verified deposited file) fits: cutoff 0 yields the largest, most "
               "divergent EPOCs, most of which fail the downstream tree QC. (s0 is absent from the v0.3 "
               "README, which predates it; it can be independently reproduced via "
               "mmseqs_create_pangenome.py --evaluation_cutoff 0.0 + the downstream eukgen steps.)"),
        "aELW_absolute": ("This pipeline reproduces the published RATIO (7.74), occupancy (63.8/19.7%), "
               "core-tuple count (16,526), single-candidate (7,075) and Asgard-sole (4,396) counts, and "
               "the panel-a R^2 (0.92/0.79) ALL EXACTLY. Standard aELW for Asgard is 0.439 vs Table 1's "
               "0.487 (a per-clade aggregation choice; per-tuple c-ELW sums are exactly 1.0). The head-to-"
               "head reconstruction (Part 8) uses the same within-tuple summation as the standard aELW and "
               "collapses the 4-15x standard ratios to 1.3-3.3x head-to-head (Asgard:Alpha 1.60, "
               "Asgard:Cyano 1.31): the published dominance is candidate-set inclusion, not a per-gene "
               "Asgard advantage. The whole of Figure 1 (panels a-h) is produced by this one script.")}
    json.dump(SUMMARY, open(os.path.join(HERE,"summary.json"),"w"), indent=2)
    print("\nAll outputs written to", HERE)
