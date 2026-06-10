#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recompute the aELW table for the Matters Arising on ONE explicit core set.

Purpose
-------
Produce an internally self-consistent Extended Data Table 1 directly from the
deposited EPOC TSVs, so the manuscript table and the deposited code cannot drift.
Every quantity is computed on the SAME core set (max c-ELW > 0.4 -> 16,526 tuples).

Definition (stated explicitly, computed only from deposited data)
-----------------------------------------------------------------
c-ELW partitions each tuple's total weight (= 1.0) among that tuple's candidate
clades. Within a tuple, a CLASS's weight is the SUM of the c-ELW of its candidate
clades (its share of that tuple). aELW then averages this per-tuple class weight:
    standard aELW     = sum over all core tuples / (number of core tuples)
    conditional aELW  = sum over all core tuples / (tuples where the class is a candidate)
    occupancy         = (tuples where the class is a candidate) / (number of core tuples)
By construction:  standard aELW = conditional aELW x occupancy,
and the standard aELW summed over all classes = mean per-tuple total = 1.0.

Head-to-head c-ELW (inclusion-controlled): on tuples where BOTH classes are
candidates, take each class's per-tuple summed c-ELW, average over those tuples,
and report the ratio (and the two-way renormalised means that sum to 1).

NOTE ON PROVENANCE: this recomputation does NOT reproduce the previously tabulated
absolute aELW (Asgard 0.487 / Alpha 0.063). Those values correspond to a smaller
(~14,880-tuple) denominator and are almost certainly the authors' published aELW on
their own core; the script prints the implied denominator so this can be confirmed.
The RATIOS, occupancy, head-to-head, sole-candidate counts and the cut-off sweep are
unaffected and reproduce the manuscript exactly.

Usage
-----
    python recompute_aelw.py [DATA_DIR]
DATA_DIR defaults to the parent of this script (the deposit layout). Expects
EPOC_data.tsv (== s10) and, for the sweep, EPOC_data.pangenome_s{0,25,67}.tsv.
Outputs: table1_recomputed.csv, aelw_per_taxon.csv, table1_recomputed.json
"""
import os, sys, json
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(HERE)
CORE = 0.40
FOCAL = "Asgard"
BACT6 = ["Alphaproteobacteria", "Cyanobacteriota", "Actinomycetota",
         "Betaproteobacteria", "Myxococcota", "Gammaproteobacteria"]
SFILES = {"s0":  "EPOC_data.pangenome_s0.tsv", "s10": "EPOC_data.tsv",
          "s25": "EPOC_data.pangenome_s25.tsv", "s67": "EPOC_data.pangenome_s67.tsv"}


def load(fn):
    df = pd.read_csv(fn, sep="\t", low_memory=False,
                     usecols=["tree_name", "euk_clade_rep", "prok_taxa", "c-ELW"])
    df["celw"] = pd.to_numeric(df["c-ELW"], errors="coerce").fillna(0.0)
    df["key"] = list(zip(df["tree_name"], df["euk_clade_rep"]))
    return df


def core_rows(df):
    mx = df.groupby("key")["celw"].max()
    core = set(mx[mx > CORE].index)
    return df[df["key"].isin(core)], len(core)


def aelw_table(d, ncore):
    """Per-class occupancy, standard and conditional aELW (one core set)."""
    tw = d.groupby(["key", "prok_taxa"])["celw"].sum()            # per-tuple class weight
    num = tw.groupby("prok_taxa").sum()                           # numerator (sum over tuples)
    npres = tw.reset_index().groupby("prok_taxa")["key"].nunique()  # tuples where present
    out = pd.DataFrame({"occupancy": npres / ncore,
                        "standard_aELW": num / ncore,
                        "conditional_aELW": num / npres,
                        "n_present": npres, "weight_sum": num})
    return out.sort_values("standard_aELW", ascending=False)


def head_to_head(d, focal, others):
    wide = d.pivot_table(index="key", columns="prok_taxa", values="celw", aggfunc="sum")
    rows = []
    for x in others:
        if x not in wide.columns:
            continue
        sub = wide[[focal, x]].dropna()
        mf, mx = sub[focal].mean(), sub[x].mean()
        wins = int((sub[focal] > sub[x]).sum())
        rows.append(dict(vs=x, n_paired=len(sub),
                         mean_focal=round(mf, 3), mean_other=round(mx, 3),
                         renorm_focal=round(mf / (mf + mx), 3),
                         renorm_other=round(mx / (mf + mx), 3),
                         h2h_ratio=round(mf / mx, 2),
                         wins_focal=wins, pct_focal=round(100 * wins / len(sub), 1),
                         wins_other=len(sub) - wins))
    return pd.DataFrame(rows)


def sole_candidate(d):
    cc = d.groupby("key").size()
    single = set(cc[cc == 1].index)
    sd = d[d["key"].isin(single)]
    sole = sd.groupby("key")["prok_taxa"].first()
    return len(single), sole.value_counts()


def main():
    main_fn = os.path.join(DATA, SFILES["s10"])
    df = load(main_fn)
    d, ncore = core_rows(df)
    tab = aelw_table(d, ncore)
    h2h = head_to_head(d, FOCAL, BACT6)
    n_single, sole = sole_candidate(d)

    a, al = tab.loc["Asgard"], tab.loc["Alphaproteobacteria"]
    print("=" * 70)
    print(f"Recomputed aELW on a single core set (max c-ELW > {CORE}); n_core = {ncore:,}")
    print("=" * 70)
    print("\nAsgard vs Alphaproteobacteria")
    print(f"  occupancy           {100*a.occupancy:5.1f}% / {100*al.occupancy:5.1f}%"
          f"   ratio {a.occupancy/al.occupancy:.2f}")
    print(f"  standard aELW       {a.standard_aELW:.3f} / {al.standard_aELW:.3f}"
          f"   ratio {a.standard_aELW/al.standard_aELW:.2f}")
    print(f"  conditional aELW    {a.conditional_aELW:.3f} / {al.conditional_aELW:.3f}"
          f"   ratio {a.conditional_aELW/al.conditional_aELW:.2f}")
    hh = h2h[h2h.vs == "Alphaproteobacteria"].iloc[0]
    print(f"  head-to-head c-ELW  {hh.renorm_focal:.3f} / {hh.renorm_other:.3f}"
          f"   ratio {hh.h2h_ratio:.2f}  (n={hh.n_paired}, wins {hh.wins_focal}/{hh.wins_other}={hh.pct_focal}%)")
    print(f"  sole-candidate      Asgard {int(sole.get('Asgard',0)):,} / "
          f"Alpha {int(sole.get('Alphaproteobacteria',0)):,}  of {n_single:,} single-candidate tuples")

    # internal-consistency diagnostics (informational, not assertions)
    chk = abs(a.standard_aELW - a.conditional_aELW * a.occupancy) < 1e-9
    total = tab["standard_aELW"].sum()
    print(f"\n  consistency: standard = conditional x occupancy ? {chk};  "
          f"sum of standard aELW over all classes = {total:.4f} (~1.0 if every tuple's c-ELW sums to 1)")

    # provenance diagnostic for the previously published absolute aELW
    for name, pub in [("Asgard", 0.487), ("Alphaproteobacteria", 0.063)]:
        implied = tab.loc[name, "weight_sum"] / pub
        print(f"  to obtain published {name} aELW {pub}: implied denominator "
              f"{implied:,.0f}  (this core = {ncore:,})")

    # head-to-head table
    print("\nHead-to-head vs six bacterial candidates (standard ratio -> head-to-head ratio):")
    for _, r in h2h.iterrows():
        std = tab.loc["Asgard", "standard_aELW"] / tab.loc[r.vs, "standard_aELW"]
        print(f"  {r.vs:22s} {std:5.2f}  ->  {r.h2h_ratio:.2f}   (n={r.n_paired})")

    # cut-off sweep
    print("\nSensitivity to soft-core cut-off (standard Asgard:Alpha ratio, occupancy):")
    sweep = []
    for tag, fn in SFILES.items():
        p = os.path.join(DATA, fn)
        if not os.path.exists(p):
            print(f"  {tag}: [not found]"); continue
        dd, nc = core_rows(load(p))
        t = aelw_table(dd, nc)
        r = t.loc["Asgard", "standard_aELW"] / t.loc["Alphaproteobacteria", "standard_aELW"]
        sweep.append(dict(cutoff=tag, n_core=nc,
                          aELW_Asgard=round(t.loc["Asgard", "standard_aELW"], 3),
                          aELW_Alpha=round(t.loc["Alphaproteobacteria", "standard_aELW"], 3),
                          ratio=round(r, 2),
                          occ_Asgard=round(100 * t.loc["Asgard", "occupancy"], 1),
                          occ_Alpha=round(100 * t.loc["Alphaproteobacteria", "occupancy"], 1)))
        print(f"  {tag:4s} n_core {nc:>6,}  aELW {sweep[-1]['aELW_Asgard']:.3f}/"
              f"{sweep[-1]['aELW_Alpha']:.3f}  ratio {r:5.2f}  occ {sweep[-1]['occ_Asgard']:.1f}/{sweep[-1]['occ_Alpha']:.1f}")

    # write outputs
    tab.round(4).to_csv(os.path.join(HERE, "aelw_per_taxon.csv"))
    h2h.to_csv(os.path.join(HERE, "table1_recomputed_headtohead.csv"), index=False)
    pd.DataFrame(sweep).to_csv(os.path.join(HERE, "table1_recomputed_sweep.csv"), index=False)
    json.dump(dict(n_core=ncore,
                   asgard=tab.loc["Asgard"].round(4).to_dict(),
                   alpha=tab.loc["Alphaproteobacteria"].round(4).to_dict(),
                   head_to_head=h2h.to_dict("records"),
                   sole_candidate=dict(n_single=n_single,
                                       asgard=int(sole.get("Asgard", 0)),
                                       alpha=int(sole.get("Alphaproteobacteria", 0))),
                   sweep=sweep),
              open(os.path.join(HERE, "table1_recomputed.json"), "w"), indent=2)
    print("\nWrote aelw_per_taxon.csv, table1_recomputed_headtohead.csv, "
          "table1_recomputed_sweep.csv, table1_recomputed.json to", HERE)


if __name__ == "__main__":
    main()
