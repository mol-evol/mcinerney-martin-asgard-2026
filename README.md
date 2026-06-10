# Re-analysis of Tobiasson et al. (2026) — Asgard contribution to eukaryogenesis

Code for the Matters Arising by **J. O. McInerney & W. F. Martin** on Tobiasson, Luo,
Wolf & Koonin, "Dominant contribution of Asgard archaea to eukaryogenesis", *Nature*
**650**, 141–149 (2026).

Every analysis uses **only the authors' own deposited data** (Zenodo record
[10.5281/zenodo.15048010](https://doi.org/10.5281/zenodo.15048010)) and their analytical
pipeline; no external data are introduced.

## Data (not included here — download from Zenodo)

Place the following files from Zenodo 15048010 in a single directory:

| File | Used by |
|------|---------|
| `EPOC_data.tsv` (= s10) | all |
| `EPOC_data.pangenome_s0.tsv`, `…_s25.tsv`, `…_s67.tsv` | filter sweep, s0 tests |
| `EPOC_data.tar.gz` (master trees, ~1.2 GB) | tree analysis (Parts 4–5) |
| `EPOC_annotation_KEGG.tsv`, `KEGG_category_mapping.tsv` | KEGG stratification, functional test |

By default the scripts expect these one directory above the script (the deposit layout);
otherwise pass the data directory as the first argument.

## Scripts

| Script | Output |
|--------|--------|
| `asgard_reanalysis.py` | full pipeline (Parts 1–9): validation, filter sweep, within-gene flip, tree analysis, KEGG stratification, head-to-head, stem length; writes the 8-panel and stem figures |
| `recompute_aelw.py [DATA_DIR]` | recomputed aELW table on the 16,526-tuple core |
| `make_main_figure.py [DATA_DIR] [OUT_DIR]` | main-text Figure 1 (two panels) |
| `functional_enrichment.py [DATA_DIR]` | KEGG functional enrichment of s0 alphaproteobacterial winners |

## Requirements

Python 3.10+. `pip install -r requirements.txt` (pandas, numpy, scipy, matplotlib, ete3).
`asgard_reanalysis.py` also calls `gunzip` to stream the tree archive.

## Reproducing the headline numbers

```bash
python recompute_aelw.py /path/to/zenodo_15048010      # 7.74 = 3.24 x 2.39; head-to-head 1.60
python make_main_figure.py /path/to/zenodo_15048010 .  # Figure 1
python functional_enrichment.py /path/to/zenodo_15048010
```

## Notes

The in-code comment in `asgard_reanalysis.py` that the s0 set has "5,846" EPOCs is
superseded: the md5-verified deposited s0 file yields 7,276 core tuples.

## Citation

Please cite the Matters Arising (McInerney & Martin) and the original data deposit
(Tobiasson et al., Zenodo 10.5281/zenodo.15048010).

## License

MIT — see `LICENSE`.
