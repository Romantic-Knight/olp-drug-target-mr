# OLP drug-target Mendelian randomization — analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24052532.svg)](https://doi.org/10.5281/zenodo.24052532)

Archived as `v1.0.0` at DOI **10.5281/zenodo.24052533**; the badge above resolves to the latest archived version.

Code accompanying the manuscript *"From germline risk to lesional tissue:
druggable-genome Mendelian randomization and single-cell profiling nominate
GDF11 inhibition for oral lichen planus"*.

This repository is the analysis code only. No data are redistributed here:
all inputs are public summary statistics or public single-cell data, listed
below with their accessions.

- **[`FIGURE_PROVENANCE.md`](FIGURE_PROVENANCE.md)** — which script draws each
  figure, panel by panel, and which files it reads. Start here.
  Supplementary figures are numbered in order of first citation in the
  manuscript: Supplementary Figure 1 is the regional-association plate and
  Supplementary Figure 2 the dataset-QC plate. The stems written by the
  builders therefore differ from the delivered filenames in two places, and
  the renaming is documented under each figure.
- `analysis/` — the pipeline and the figure builders, in their original layout.
- `env/` — the pixi environment specifications used for each analysis domain.
- `provenance/` — machine-readable versions of the above
  (`figure_provenance.csv`, `script_manifest.csv`).

## What the study does

A druggable-genome MR screen in oral lichen planus (OLP), then colocalisation,
directionality and pleiotropy checks, then protein-level and docking evidence,
then a single-cell atlas analysis of lesional tissue.

| stage | what it produces |
|---|---|
| `01_druggable_set` | the druggable-gene union and its cis-eQTL intersection |
| `02_mr` | instruments, Wald-ratio MR, replication, IV strength, Steiger test |
| `03_coloc` | coloc.abf and coloc.susie over 22 hits × 2 case definitions |
| `04_pleiotropy` | PheWAS and pleiotropy validation plots |
| `05_pqtl` | pQTL assessment (deCODE / INTERVAL / Knight ADRC) |
| `06_docking` | receptor and ligand preparation, Vina docking, decoys, re-docking |
| `07_molecular_dynamics` | GROMACS build/run/analysis, MM-GBSA, fpocket |
| `08_single_cell` | GSE211630 preprocessing, clustering, annotation, composition, erosion, LIANA, MR back-mapping |
| `09_figures` | the figure plates and their panel scripts |
| `10_manuscript` | supplementary tables and the Word deliverables |

Stages are a reading order, not a directory layout: files keep their original
paths (`analysis/scripts/...`, `analysis/md/scripts/...`) so that every relative
path inside them still resolves. `provenance/script_manifest.csv` maps each file
to its stage.

## Order of operations

```bash
# 1. environment (pick the domain you need)
cd env/scrna  && pixi install     # single-cell analysis + single-cell figure plates
cd env/md     && pixi install     # GROMACS / AmberTools / fpocket
cd env/docking && pixi install    # AutoDock Vina toolchain

# 2. data preparation, then each analysis stage
Rscript analysis/scripts/prepare_coloc_data.R
Rscript analysis/scripts/run_coloc.R
Rscript analysis/scripts/run_susie_coloc.R
...

# 3. figures
Rscript analysis/scripts/build_restructured_figures.R          # Figures 1-3
cd analysis/scrna && pixi run python ../../analysis/scripts/\
build_restructured_sc_figures.py all                           # Figures 5-7, Supp 1/3/4
python analysis/md/scripts/fig_md_integrated.py                # Figure 4
```

The single-cell scripts must run inside `env/scrna`; the R figures were run
under R 4.4.2 with the package set pinned in `env/md/pixi.toml`. Exact tool
versions are recorded in `env/`.

## Data

Nothing in this repository downloads or contains participant-level data.

| input | accession / source |
|---|---|
| Whole-blood cis-eQTL instruments | eQTLGen consortium (Võsa et al. 2021), <https://www.eqtlgen.org> |
| OLP outcome summary statistics | FinnGen release R13, endpoints `K11_LICHEN_PLANUS_WIDE` and `K11_ORAL_LICHEN_PLANUS_INCLAVO` |
| Linkage-disequilibrium reference | 1000 Genomes Phase 3, European (GRCh37) |
| Plasma pQTL | deCODE (Ferkingstad et al. 2021), INTERVAL (Sun et al. 2018), Knight ADRC |
| Druggable-gene resources | Pharos, DGIdb, Finan et al. 2017, Open Targets |
| Single-cell atlas | GEO **GSE211630** (buccal mucosa). Files used: GSM6481637 (normal), GSM6481638–GSM6481640 (non-erosive lesions), GSM6481641–GSM6481642 (erosive lesions) |
| Docking structures | PDB **2A9A** (SUOX, chicken) and **5E4G** (GDF11) |
| Docking ligands | DrugSig / PubChem, prepared as described in Methods 2.8 |

Derived tables under `results/` are not included, because several of them are
third-party summary statistics at the gene level. `FIGURE_PROVENANCE.md` lists
every table each figure needs, and which script produces it — including the
handful whose generating script was **not** retained (see the last section
there; the gap is stated rather than papered over).

## Reproducing on your own checkout

The scripts were developed in a Windows working copy and were copied here
**byte for byte**, so they are the exact code that produced the manuscript.
Two consequences, both of which need a one-line edit before a run:

1. **A few files hard-code the project root.** Where you see
   `ROOT = r"D:\AI-Wisp-Science\Drug Target MR"` (or a `DEST_DIR`, `FG_DIR`,
   `BASE` variant of it) near the top of a script, point it at your checkout.
   `provenance/script_manifest.csv` lists every file, and searching for
   `AI-Wisp-Science` finds all of them:

   ```bash
   grep -rln 'AI-Wisp-Science' analysis/
   ```
   and edit the constant in place. A blind `sed` is risky here: the path sits
   inside Python raw strings and R string literals with differing slash
   conventions.

2. **Seven figure scripts load a plotting sidecar from an absolute path.** The
   single-cell and docking plates call helper modules (`figure-style/runtime.py`
   and, in two cases, `figure-style/kernel.py`) that provide
   `apply_figure_style()`, `set_frame()` and `panel_letter()`. Those modules are
   part of the authors' local tooling and are **not** distributed here. The
   scripts will fail at import until you supply modules exposing those three
   functions, or delete the calls and accept default matplotlib styling. The
   affected files are:

   ```
   analysis/md/scripts/fig_md_integrated.py                           runtime
   analysis/scripts/plot_scrna_candidates_gse211630.py                runtime + kernel
   analysis/scripts/plot_scrna_cellqc_gse211630.py                    kernel
   analysis/scripts/scrna_figD1_build_gse211630.py                    runtime
   analysis/scripts/scrna_suppfig5_dataset_qc_gse211630.py            runtime
   analysis/scripts/scrna_suppfig7_marker_dotplot_gse211630.py        runtime
   analysis/scripts/scrna_suppfig8_candidate_expression_gse211630.py  runtime
   ```

Neither of these affects the analysis itself, with one exception worth naming:
`analysis/scripts/plot_scrna_candidates_gse211630.py` writes
`results/tables/geo_gse211630_candidate_loc.csv` as a side product while also
calling the sidecar, so that one table needs the script to run.

## Not included, deliberately

| left out | why |
|---|---|
| `_`-prefixed scripts, `qa_*`, `check_*`, `test_*`, `debug_*`, `probe_*` | one-off inspection, audit and fix-up drivers used while building the paper; they assert on intermediate states that no longer exist |
| `analysis/scripts/MR/drugGene*.R` | third-party tutorial templates (Ke Biowolf), superseded by the retained MR scripts |
| `build_options.R` | renders alternative layout proposals that were not used |
| the server transfer scripts | infrastructure for a rented GPU host, not analysis; two of them contained a credential and are excluded for that reason |
| `analysis/scripts/_retired/`, `__pycache__`, pixi environments | build artefacts |

## Software versions

The pinned specifications are in `env/`:

| environment | platform pin | provides | lock file |
|---|---|---|---|
| `env/scrna` | `win-64` | Python 3.12, scanpy, anndata, leidenalg, igraph, LIANA | yes |
| `env/md` | `linux-64` | GROMACS, AmberTools, gmx_MMPBSA, MDAnalysis | not retained |
| `env/docking` | `win-64` | AutoDock Vina toolchain, RDKit | yes |

`env/md` ships its `pixi.toml` and its stage README only. No solved lock file
survived for that environment, so re-solving it on a new host may resolve to
different package builds; the GROMACS and AmberTools versions actually used are
recorded below.

The two `win-64` pins reflect the machine the analysis was run on; change the
`platforms` line for Linux. `env/md` is `linux-64` because AmberTools is
Linux-only.

Measured tool versions, recorded on the compute host at the time the results
were produced: GROMACS 2026.3 (`conda_forge` build), AmberTools 26.0,
gmx_MMPBSA 1.5.0.3, Ubuntu 22.04.5, NVIDIA driver 560.35.03. R 4.4.2 provided
`data.table`, `ggplot2`, `ggrepel`, `patchwork`, `VennDiagram`, `coloc` (with
`coloc.susie`), `susieR`, `cowplot`, `scales`, `GenomicRanges` and `Rsamtools`.
The Python analysis used `scanpy`, `anndata`, `scipy`, `pandas`, `numpy`,
`rdkit`, `liana`, `scrublet` and `python-docx`.

## Contact

All authors: Yahui Lu, Xiteng Yin, Jianchuan Ran, Wenguang Xu, Yufeng Wang,
Wei Han.

Corresponding author: Wei Han — Nanjing Stomatological Hospital, Affiliated
Hospital of Medical School, Institute of Stomatology, Nanjing University,
Nanjing 210008, China. doctorhanwei@nju.edu.cn

## Licence

MIT — see [`LICENSE`](LICENSE). The licence covers this analysis code. It does
not extend to the third-party summary statistics or single-cell data analysed
with it, which remain under their own terms.
