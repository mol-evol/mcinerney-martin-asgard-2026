# Re-analysis of Tobiasson et al. (2026) — Asgard contribution to eukaryogenesis

Code and derived data for the Matters Arising by **J. O. McInerney & W. F. Martin** on
Tobiasson, Luo, Wolf & Koonin, "Dominant contribution of Asgard archaea to eukaryogenesis",
*Nature* **650**, 141–149 (2026).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20628634.svg)](https://doi.org/10.5281/zenodo.20628634)

This release (**v2.0.0**) accompanies the revised Comment. Every analysis uses **only the
authors' own deposited data** (Zenodo record
[10.5281/zenodo.15048010](https://doi.org/10.5281/zenodo.15048010)) and their analytical
pipeline; no external data are introduced.

## What changed at revision (v1.0.0 → v2.0.0)

| Change | Where |
|---|---|
| Main Fig. 1 panels swapped: filter effect = **a**, occupancy diagnostic = **b** | `make_main_figure.py`, `rebuild_v15_figures.py` |
| "Head-to-head" relabelled **"Conditional"** throughout (it is a diagnostic, not a summary statistic) | both figure scripts |
| Extended Data Fig. 1 reduced from 8 panels to 5 (a–e) | `asgard_reanalysis.py: part6_figure`, `rebuild_v15_figures.py` |
| — deleted: Asgard:Alpha stratification and filter sweep (now main Fig. 1a,b) | |
| — deleted: Asgard-monophyly panel (Asgard paraphyly is expected if eukaryotes branch within Asgard) | |
| — relettered d→c, f→d, g→e; full clade names with leader lines in panel a | |
| Extended Data Fig. 2b retitled: the spread argument now targets **post-acquisition adaptation**, not acquisition timing, which Tobiasson et al. themselves reject | `asgard_reanalysis.py: _stem_figure` |
| **KEGG best-hit tie-break made deterministic** (`kind='mergesort'`). 24.5% of queries have their maximum `Prob` tied across >1 distinct target; the previous default (unstable quicksort) made those assignments irreproducible between runs | `asgard_reanalysis.py` (Parts 5 and 9), `functional_enrichment.py`, `rebuild_ed_fig2.py` |
| Numbers that moved as a result: ED Fig. 2b spreads 72×/47× → **70×/50×** (n = 65/50 → **58/54**); SI functional enrichment recomputed | manuscript, Supplementary Information |

## Data

Download from Zenodo 15048010 and place in a single directory:

| File | Used by |
|------|---------|
| `EPOC_data.tsv` (= s10) | all |
| `EPOC_data.pangenome_s0.tsv`, `…_s25.tsv`, `…_s67.tsv` | filter sweep, s0 tests |
| `EPOC_data.tar.gz` (master trees, ~1.2 GB) | tree analysis (Parts 4–5) |
| `EPOC_annotation_KEGG.tsv`, `KEGG_category_mapping.tsv` | KEGG stratification, functional test, stem sets |

By default the scripts expect these one directory above the script (the deposit layout);
otherwise pass the data directory as the first argument.

### Reduced inputs included here

To allow the figures to be rebuilt without the full deposit, two column subsets are provided:

- `EPOC_slim_columns.tsv` — `tree_name, euk_clade_rep, prok_taxa, raw_stem_length,
  median_euk_leaf_dist, stem_length, c-ELW` from `EPOC_data.tsv`. Note that
  `EPOC_data.tsv` carries an unnamed leading index column, so these are fields 2,3,14,17,18,19,28.
- `KEGG_annotation_3col.tsv.gz` — `Query, Target, Prob` from `EPOC_annotation_KEGG.tsv`.
  The full three columns are kept rather than a pre-reduced best-hit table, because the
  best-hit selection depends on the tie-break rule and must be reproducible from the input.
- `EPOC_s0_slim_columns.tsv` — `tree_name, euk_clade_rep, prok_taxa, c-ELW` from
  `EPOC_data.pangenome_s0.tsv` (fields 2,3,14,28), for the s0 functional enrichment.

### A note on reproducibility

Nearly half of all KEGG queries score `Prob = 100`, and 24.5% have that maximum tied
across more than one distinct target. Selecting a single best hit therefore requires a
tie-break. All scripts here now use a stable sort, so the selection is deterministic.
The scientific conclusion does not depend on the choice: the within-event normalised-stem
spread is 45- to 71-fold under every rule we tested (first-max, last-max, alphabetical),
against inter-donor differences roughly an order of magnitude smaller.

## Scripts

| Script | Output |
|--------|--------|
| `asgard_reanalysis.py` | full pipeline (Parts 1–9): validation, filter sweep, within-gene flip, tree analysis, KEGG stratification, conditional c-ELW, stem length; writes Extended Data Fig. 1 (5 panels) and the stem figure |
| `recompute_aelw.py [DATA_DIR]` | recomputed aELW table on the 16,526-tuple core |
| `make_main_figure.py [DATA_DIR] [OUT_DIR]` | main-text Figure 1 (two panels, v15 order) |
| `functional_enrichment.py [DATA_DIR]` | KEGG functional enrichment of s0 alphaproteobacterial winners |
| `rebuild_v15_figures.py` | rebuilds main Fig. 1 and Extended Data Fig. 1 from the reduced inputs above, without the 1.2 GB tree archive |
| `rebuild_ed_fig2.py` | rebuilds Extended Data Fig. 2 from the reduced inputs, without the 1.2 GB tree archive |

## Requirements

Python 3.10+. `pip install -r requirements.txt` (pandas, numpy, scipy, matplotlib, ete3).
`asgard_reanalysis.py` also calls `gunzip` to stream the tree archive.

## Reproducing the headline numbers

```bash
python recompute_aelw.py /path/to/zenodo_15048010      # 7.74 = 3.24 x 2.39; conditional 1.60
python make_main_figure.py /path/to/zenodo_15048010 .  # Figure 1
python functional_enrichment.py /path/to/zenodo_15048010
```

All of the following reproduce exactly on the current deposit: 16,526 core tuples
(max c-ELW > 0.4); aELW ratio 7.74 = occupancy 3.24 × conditional 2.39; occupancy
63.8% / 19.7%; 1,509 paired tuples with conditional ratio 1.601; 7,075 single-candidate
tuples of which 4,396 (62.1%) are Asgard-only; stratification 7.74 → 1.62; panel-a
R² 0.92 / 0.79; filter sweep 0.44 / 7.74 / 6.13 / 14.02; 5,381 Asgard-winner trees,
63.3% with no bacterial leaf. See `../RECONCILIATION.md`.

## Citation

Please cite the Matters Arising (McInerney & Martin), this code archive, and the original
data deposit:

- **This code**, version of record for the published analysis:
  [10.5281/zenodo.21905181](https://doi.org/10.5281/zenodo.21905181) (v2.0.0).
  To cite all versions and always resolve to the latest, use the concept DOI
  [10.5281/zenodo.20628634](https://doi.org/10.5281/zenodo.20628634).
- **Original data** (Tobiasson et al.):
  [10.5281/zenodo.15048010](https://doi.org/10.5281/zenodo.15048010).

## License

MIT — see `LICENSE`.
