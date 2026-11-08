# Figure provenance

Which script draws each manuscript figure, and which files it reads.

All seven main figures and all four supplementary figures in `figures/` are listed.
Every entry below was established by reading the code: the plate builders for the
panel-to-object mapping, the panel scripts for their data reads, and each input
table's producer by locating the actual write call. Nothing here is inferred from
filenames.

## How the figures were assembled

The plates are **composed from re-rendered panels**, never cut from earlier images.
Three programs do all the work:

| program | figures | mechanism |
|---|---|---|
| `analysis/scripts/build_restructured_figures.R` | 1, 2, 3 | sources each panel script into an isolated R environment with every figure-writing call neutralised, then composes the returned ggplot objects with `patchwork` |
| `analysis/scripts/build_restructured_sc_figures.py` | 5, 6, 7, Supp 1, 3, 4 | calls each panel module's panel function on a freshly created axes of the new plate, so the drawn content is the original panel at its new position |
| `analysis/md/scripts/fig_md_integrated.py` | 4 | builds the 9-panel docking/MD plate directly from the result tables |

Main figures land in `figures/main_figures_v2/`, supplementary figures in
`figures/supplementary_v2/`, in both `.png` (300 dpi) and `.pdf`. All three
programs expect to be run from the project root, because every data path in the
panel scripts is written relative to it (`results/tables/...`,
`data/processed/...`).

Supplementary Figures 2 and 4 are the two exceptions to the naming above: their
delivered filenames were assigned after generation, from sources named differently.
See the note under each.

Reproduce in this order (each plate writes only into `figures/`):

```bash
# 1-3, R plate
Rscript analysis/scripts/build_restructured_figures.R
# 5-7 and Supp 1/3/4, single-cell plate
cd analysis/scrna && pixi run python \
    ../../analysis/scripts/build_restructured_sc_figures.py all
# 4, docking/MD plate
python analysis/md/scripts/fig_md_integrated.py
```

## Figure 1 — `fig01_design_screen`

**Output** `figures/main_figures_v2/fig01_design_screen.png` (300 dpi) and `.pdf`, 180 x 210 mm

**Built by** `analysis/scripts/build_restructured_figures.R`

```
save_v2("fig01_design_screen", fig1_v2, 180, 210)
```

Rendered by sourcing the panel scripts into isolated environments with every figure-writing call neutralised, so the plates are composed from re-rendered panels and never from cropped images.

| Panel | Content | Drawing code | Input file(s) |
|---|---|---|---|
| **A** | Study-design schematic: druggable union to MR screen to tiering | `fig1_parts.R -> render_design()` | `data/processed/targets/druggable_genes_final.csv` |
| **B** | Four-source druggable-gene Venn (Pharos, DGIdb, Finan et al., Open Targets) | `fig1_parts.R -> render_venn()` | `data/processed/targets/druggable_genes_final.csv` |
| **C** | Druggable genes carrying a cis-eQTL instrument (5,801 of 9,261) | `fig1_parts.R -> render_funnel()` | `data/processed/targets/druggable_genes_final.csv`<br>`data/processed/targets/eqtlgen_genes.csv`<br>`data/processed/targets/druggable_x_eqtlgen.csv` |
| **D** | Manhattan plot, WIDE case definition (5,040 genes) | `journal_fig2.R -> png_a  (re-scaled x axis in the plate script)` | `results/tables/mr_olp_wide_results.csv` |
| **E** | Manhattan plot, INCLAVO case definition | `journal_fig2.R -> png_b` | `results/tables/mr_olp_inclavo_results.csv` |
| **F** | Volcano plot of the WIDE screen | `journal_fig2.R -> vol` | `results/tables/mr_olp_wide_results.csv` |

## Figure 2 — `fig02_replication_colocalisation`

**Output** `figures/main_figures_v2/fig02_replication_colocalisation.png` (300 dpi) and `.pdf`, 180 x 280 mm

**Built by** `analysis/scripts/build_restructured_figures.R`

```
save_v2("fig02_replication_colocalisation", fig2_v2, 180, 280)
```

Panels A, E and F reuse objects that were originally standalone supplementary figures; they are re-rendered here, not copied as images.

| Panel | Content | Drawing code | Input file(s) |
|---|---|---|---|
| **A** | Dual-cohort forest plot of the 22 WIDE Bonferroni-significant genes | `build_supplementary_figures.R -> p1` | `results/tables/supplementary/suppl_table_02_wide_hits_mr_replication_coloc.csv` |
| **B** | Forest plot of the five tiered candidates | `journal_fig3.R -> p` | `results/tables/final_target_evidence.csv` |
| **C** | Multi-signal (SuSiE) vs single-signal (coloc.abf) PP.H4 | `journal_fig4.R -> pb` | `results/tables/coloc_olp_results.csv`<br>`results/tables/susie_coloc_results.csv` |
| **D** | coloc.abf PP.H4 for the 22 WIDE hits, both case definitions | `journal_fig4.R -> pa` | `results/tables/coloc_olp_results.csv`<br>`results/tables/mr_olp_wide_results.csv` |
| **E** | Full posterior decomposition across the 44 gene-definition runs | `build_supplementary_figures.R -> p2` | `results/tables/supplementary/_retired/coloc_abf_susie_44runs.csv` |
| **F** | Bonferroni vs BH-FDR multiplicity sensitivity | `build_supplementary_figures.R -> p3` | `results/tables/supplementary/suppl_table_01_full_mr_screen_5040.csv` |

## Figure 3 — `fig03_pleiotropy_docking`

**Output** `figures/main_figures_v2/fig03_pleiotropy_docking.png` (300 dpi) and `.pdf`, 180 x 207 mm

**Built by** `analysis/scripts/build_restructured_figures.R`

```
save_v2("fig03_pleiotropy_docking", fig3_v2, 180, 207)
```

All four panels are objects of build_final_mr_figures.R, which itself sources journal_fig5.R for the PheWAS panel.

| Panel | Content | Drawing code | Input file(s) |
|---|---|---|---|
| **A** | PheWAS of each candidate's top cis-eQTL across FinnGen R13 phenotypes | `build_final_mr_figures.R -> phewas  (= journal_fig5.R -> pa)` | `results/tables/phewas_olp_all.csv` |
| **B** | Steiger directionality test | `build_final_mr_figures.R -> steiger` | `results/tables/steiger_direction_test.csv` |
| **C** | SuSiE pleiotropy for the chr12 candidates | `build_final_mr_figures.R -> susie` | `results/tables/susie_pleiotropy_results.csv` |
| **D** | PP.H4 matrix across the pleiotropy colocalisation grid | `build_final_mr_figures.R -> heatmap` | `results/tables/pleiotropy_coloc_results.csv` |

## Figure 4 — `fig04_docking_md`

**Output** `figures/main_figures_v2/fig04_docking_md.png` (300 dpi) and `.pdf`, single canvas, 9 panels

**Built by** `analysis/md/scripts/fig_md_integrated.py`

```
python analysis/md/scripts/fig_md_integrated.py
```

The only plate written directly to figures/main_figures_v2/ by a standalone script (the others are composed by the two restructure builders).

| Panel | Content | Drawing code | Input file(s) |
|---|---|---|---|
| **A** | Redocking all 22 pairs: previous vs re-scored affinity | `fig_md_integrated.py (panel A block)` | `results/tables/docking_redo_comparison.csv` |
| **B** | Three-seed reproducibility and SO4 self-docking control | `fig_md_integrated.py (panel B block)` | `analysis/docking/redo/docking_redo_summary.json`<br>`analysis/docking/redo/validation_so4.json` |
| **C** | Drift, backbone RMSD, RMSF and SASA triad over 200 ns | `fig_md_integrated.py (panel C block)` | `results/tables/md_curves.npy`<br>`results/tables/md_rmsf.csv`<br>`results/tables/md_sasa.csv` |
| **D** | Cofactor effect on Spiperone drift (apo vs holo) | `fig_md_integrated.py (panel D block)` | `results/tables/md_curves.npy` |
| **E** | Residue-contact hotspots | `fig_md_integrated.py (panel E block)` | `results/tables/md_residue_contacts.csv` |
| **F** | Drift versus conformational RMSD | `fig_md_integrated.py (panel F block)` | `results/tables/md_curves.npy` |
| **G** | Hydrogen bonds per frame | `fig_md_integrated.py (panel G block)` | `results/tables/md_hbonds.csv` |
| **H** | MM-GBSA binding free energy per replicate | `fig_md_integrated.py (panel H block)` | `results/tables/md_mmgbsa.csv` |
| **I** | fpocket druggability of the two receptors | `fig_md_integrated.py (panel I block)` | `results/tables/fpocket_summary.json` |

## Figure 5 — `fig05_atlas_composition`

**Output** `figures/main_figures_v2/fig05_atlas_composition.png` (300 dpi) and `.pdf`, 180 x 250 mm

**Built by** `analysis/scripts/build_restructured_sc_figures.py`

```
cd analysis/scrna && pixi run python ../../analysis/scripts/build_restructured_sc_figures.py main4
```

The single-cell plates are matplotlib; each panel is a function of a panel module called on an axes created on the new plate, so nothing is cropped.

| Panel | Content | Drawing code | Input file(s) |
|---|---|---|---|
| **A** | UMAP of the 46,478-cell atlas | `scrna_erosion_gse211630.panel_A` | `data/processed/geo_GSE211630/scanpy/gse211630_olp_lognorm.h5ad`<br>`data/processed/geo_GSE211630/scanpy/gse211630_olp_obs_annotation.csv` |
| **B** | Same embedding coloured by condition | `scrna_erosion_gse211630.panel_B` | `same as A` |
| **C** | Raw Leiden clusters c0-c20 | `scrna_suppfig5_umap_gse211630.panel_a` | `data/processed/geo_GSE211630/scanpy/gse211630_olp_lognorm.h5ad` |
| **D** | Per-sample coverage on the shared embedding | `scrna_suppfig5_umap_gse211630.panel_b` | `same as C` |
| **E** | Per-sample cell-type composition at family level | `scrna_plot_composition_tnk_gse211630.panel_composition_bars` | `results/tables/geo_gse211630_composition_per_sample.csv`<br>`results/tables/geo_gse211630_composition_shift.csv` |
| **F** | Lineage-marker dot plot (51 markers) | `scrna_suppfig7_marker_dotplot_gse211630.panel_dotplot + panel_coverage` | `data/processed/geo_GSE211630/ann/GSM648163{8..42}_filtered.h5ad`<br>`data/processed/geo_GSE211630/scanpy/gse211630_olp_lognorm.h5ad` |

## Figure 6 — `fig06_candidate_expression`

**Output** `figures/main_figures_v2/fig06_candidate_expression.png` (300 dpi) and `.pdf`, 180 x 250 mm

**Built by** `analysis/scripts/build_restructured_sc_figures.py`

```
cd analysis/scrna && pixi run python ../../analysis/scripts/build_restructured_sc_figures.py main5
```

Panel objects come from plot_scrna_candidates_gse211630.py (executed with writes redirected) and from the erosion module.

| Panel | Content | Drawing code | Input file(s) |
|---|---|---|---|
| **A-E** | Expression of the five candidates on the UMAP embedding | `plot_scrna_candidates_gse211630.py -> panel_candidate_umap(ax, gene)` | `data/processed/geo_GSE211630/scanpy/gse211630_olp_lognorm.h5ad`<br>`results/tables/geo_gse211630_celltype_annotation.csv` |
| **F-J** | Per-cell-type dotplots for the same five genes | `scrna_erosion_gse211630.erosion_grid` | `results/tables/geo_gse211630_erosion_candidate_celltype.csv`<br>`results/tables/geo_gse211630_erosion_per_sample.csv` |

## Figure 7 — `fig07_target_mapping`

**Output** `figures/main_figures_v2/fig07_target_mapping.png` (300 dpi) and `.pdf`, 240 x 270 mm

**Built by** `analysis/scripts/build_restructured_sc_figures.py`

```
cd analysis/scrna && pixi run python ../../analysis/scripts/build_restructured_sc_figures.py main6
```

240 mm wide, above the usual 180 mm single-page limit: it would be scaled to about 75% for submission.

| Panel | Content | Drawing code | Input file(s) |
|---|---|---|---|
| **A** | Target-message direction per cell subtype | `build_final_sc_fig8.panel_A` | `results/tables/geo_gse211630_mr_backmap_gene.csv`<br>`results/tables/geo_gse211630_mr_backmap_subtype.csv` |
| **B** | Oaxaca decomposition of the tissue delta | `build_final_sc_fig8.panel_B` | `same as A` |
| **C** | Seven T/NK functional module scores | `scrna_plot_composition_tnk_gse211630.panel_tnk_signature` | `results/tables/geo_gse211630_tnk_module_per_sample.csv`<br>`results/tables/geo_gse211630_tnk_module_shift.csv` |
| **D** | Incoming focus signalling heatmap | `build_final_sc_fig8.panel_C` | `results/tables/geo_gse211630_liana_incoming_focus.csv` |
| **E** | Sustained T/NK signalling shift | `build_final_sc_fig8.panel_D` | `results/tables/geo_gse211630_liana_signal_shift.csv` |
| **F** | Ligand-level extremes of the sustained losses | `build_final_sc_fig8.panel_E` | `results/tables/geo_gse211630_liana_signal_shift.csv` |

## Supplementary figures

### Supplementary Figure 1 — `suppfig01_dataset_qc`

**Output** `figures/supplementary_v2/suppfig01_dataset_qc.png` / `.pdf`, 180 x 240 mm

**Panels** A dataset provenance and per-cell QC; B QC metrics distribution; C clustering overview

**Built by** `analysis/scripts/build_restructured_sc_figures.py  (builder key: supp1)`

```
cd analysis/scrna && pixi run python ../../analysis/scripts/build_restructured_sc_figures.py supp1
```

**Panel objects** `scrna_suppfig5_dataset_qc_gse211630.py`

| Input file |
|---|
| `results/tables/geo_gse211630_qc_summary.csv` |
| `data/processed/geo_GSE211630/qc_cell_metrics.csv.gz` |
| `data/processed/geo_GSE211630/preprocess_summary.json` |
| `data/processed/geo_GSE211630/scanpy/gse211630_all_analyzed.h5ad` |

### Supplementary Figure 2 — `suppfig02_regional_plots`

**Output** `figures/supplementary_v2/suppfig02_regional_plots.png` / `.pdf`, 180 x 165 mm

**Panels** 5 candidate genes x 2 case definitions of regional association plots

**Built by** `analysis/scripts/build_supplementary_figures.R  (object p4)`

```
save_supp("supp_fig04_regional_plots", p4, height_mm = 165)
```

**Panel objects** `build_supplementary_figures.R -> p4`

| Input file |
|---|
| `data/processed/coloc/coloc_eqtl_windows.csv` |
| `data/processed/coloc/coloc_fg_wide.csv` |
| `data/processed/coloc/coloc_fg_inclavo.csv` |
| `data/processed/coloc/ld_1000g/<gene>_ld.csv.gz` |

This plate is the R object `p4`, written as `supp_fig04_regional_plots` and
delivered under the name `suppfig02_regional_plots`. The delivered files are
byte-identical to `figures/supplementary/supp_fig04_regional_plots.{png,pdf}`
(MD5 verified), so the mapping is confirmed by hash and not merely asserted.

### Supplementary Figure 3 — `suppfig03_detection_composition`

**Output** `figures/supplementary_v2/suppfig03_detection_composition.png` / `.pdf`, 180 x 200 mm

**Panels** A per-gene detection summary; B erosive composition shift

**Built by** `analysis/scripts/build_restructured_sc_figures.py  (builder key: supp3)`

```
cd analysis/scrna && pixi run python ../../analysis/scripts/build_restructured_sc_figures.py supp3
```

**Panel objects** `scrna_suppfig8_candidate_expression_gse211630.panel_C + scrna_plot_composition_tnk_gse211630.panel_subtype_dots`

| Input file |
|---|
| `results/tables/geo_gse211630_candidate_celltype_enrichment.csv` |
| `results/tables/geo_gse211630_composition_shift.csv` |

### Supplementary Figure 4 — `suppfig04_candidate_landscape`

**Output** `figures/supplementary_v2/suppfig04_candidate_landscape.png` / `.pdf`, 180 x 215 mm

**Panels** A candidate expression across 19 cell types; B per-cell-type UMAP facets

**Built by** `analysis/scripts/build_restructured_sc_figures.py  (builder key: supp2)`

```
cd analysis/scrna && pixi run python ../../analysis/scripts/build_restructured_sc_figures.py supp2
```

**Panel objects** `scrna_suppfig8_candidate_expression_gse211630.panel_A / panel_B`

| Input file |
|---|
| `data/processed/geo_GSE211630/scanpy/gse211630_olp_lognorm.h5ad` |
| `results/tables/geo_gse211630_celltype_annotation.csv` |
| `results/tables/geo_gse211630_candidate_celltype_enrichment.csv` |
| `results/tables/geo_gse211630_candidate_loc.csv` |

The builder writes this plate under the stem `suppfig02_candidate_landscape`;
it is delivered as `suppfig04_candidate_landscape`. Unlike Supplementary
Figure 2, no byte-identical copy under the original stem survives in the
repository, so this rename is recorded from the figure-renumbering step
rather than re-verified by hash. Its companion `suppfig03_...` was written by
the same builder run at the same timestamp, which is consistent with the two
having been renamed together.

## Input tables and the scripts that produce them

`NOT RETAINED` means no generating script survives in this repository: the table
was produced during an earlier analysis stage whose driver was not kept. Those
tables are inputs here, and their provenance is recorded in the manuscript
Methods rather than in code.

| Input file | Produced by | Contents |
|---|---|---|
| `results/tables/mr_olp_wide_results.csv` | **NOT RETAINED** | WIDE-definition MR screen, 5,040 genes |
| `results/tables/mr_olp_inclavo_results.csv` | **NOT RETAINED** | INCLAVO-definition MR screen, 5,040 genes |
| `results/tables/mr_olp_hits_crosscohort.csv` | `analysis/scripts/mr_wald_ratio.R` | cross-cohort hit table |
| `results/tables/final_target_evidence.csv` | **NOT RETAINED** | single authoritative table for the five tiered candidates |
| `results/tables/phewas_olp_all.csv` | **NOT RETAINED** | PheWAS of the five candidate instruments |
| `results/tables/coloc_olp_results.csv` | `analysis/scripts/run_coloc.R` | coloc.abf, 44 runs |
| `results/tables/susie_coloc_results.csv` | `analysis/scripts/run_susie_coloc.R` | coloc.susie, 44 runs |
| `results/tables/steiger_direction_test.csv` | `analysis/scripts/run_steiger_test.py` | Steiger test |
| `results/tables/pleiotropy_coloc_results.csv` | `analysis/scripts/run_pleiotropy_coloc.R` | pleiotropy colocalisation grid |
| `results/tables/susie_pleiotropy_results.csv` | `analysis/scripts/run_susie_pleiotropy.R` | SuSiE pleiotropy |
| `results/tables/supplementary/suppl_table_01_full_mr_screen_5040.csv` | `analysis/scripts/build_supplementary_tables.py` | full 5,040-gene screen |
| `results/tables/supplementary/suppl_table_02_wide_hits_mr_replication_coloc.csv` | `analysis/scripts/build_supplementary_tables.py` | 22 WIDE hits with replication and coloc |
| `results/tables/supplementary/_retired/coloc_abf_susie_44runs.csv` | `analysis/scripts/build_supplementary_tables.py` | 44-run posterior table |
| `data/processed/targets/druggable_genes_final.csv` | **NOT RETAINED** | 9,261-gene druggable union |
| `data/processed/targets/eqtlgen_genes.csv` | `analysis/scripts/dump_eqtlgen_symbols.R` | eQTLGen gene list |
| `data/processed/targets/druggable_x_eqtlgen.csv` | **NOT RETAINED** | union x eQTLGen intersection |
| `data/processed/coloc/coloc_eqtl_windows.csv` | `analysis/scripts/prepare_coloc_data.R` | eQTL regional windows |
| `data/processed/coloc/coloc_fg_wide.csv` | `analysis/scripts/prepare_coloc_data.R` | FinnGen WIDE regional windows |
| `data/processed/coloc/coloc_fg_inclavo.csv` | `analysis/scripts/prepare_coloc_data.R` | FinnGen INCLAVO regional windows |
| `analysis/docking/redo/docking_redo_summary.json` | `analysis/scripts/redock.py` | 3-seed redocking summary |
| `analysis/docking/redo/validation_so4.json` | `analysis/scripts/validate_docking.py` | crystallographic SO4 self-docking control |
| `results/tables/docking_redo_comparison.csv` | `analysis/scripts/compare_docking_redo.py` | previous vs re-scored affinity, 22 pairs |
| `results/tables/md_curves.npy` | `analysis/md/scripts/analyse_md_local.py` | raw per-frame curves |
| `results/tables/md_stability_all.csv` | `analysis/md/scripts/analyse_md_local.py` | one row per run |
| `results/tables/md_rmsf.csv` | `analysis/md/scripts/rmsf_sasa_md.py` | per-residue RMSF |
| `results/tables/md_sasa.csv` | **NOT RETAINED** | per-run SASA from gmx sasa output |
| `results/tables/md_residue_contacts.csv` | `analysis/md/scripts/interactions_md.py` | residue contacts |
| `results/tables/md_hbonds.csv` | `analysis/md/scripts/interactions_md.py` | hydrogen bonds |
| `results/tables/md_mmgbsa.csv` | **NOT RETAINED** | parsed from the 12 Amber MMPBSA.py FINAL_RESULTS files |
| `results/tables/fpocket_summary.json` | **NOT RETAINED** | summarised from the fpocket runs |
| `results/tables/geo_gse211630_composition_per_sample.csv` | `analysis/scripts/scrna_composition_tnk_gse211630.py` | per-sample composition |
| `results/tables/geo_gse211630_composition_shift.csv` | `analysis/scripts/scrna_composition_tnk_gse211630.py` | composition shift |
| `results/tables/geo_gse211630_tnk_module_per_sample.csv` | `analysis/scripts/scrna_composition_tnk_gse211630.py` | per-sample T/NK module scores |
| `results/tables/geo_gse211630_tnk_module_shift.csv` | `analysis/scripts/scrna_composition_tnk_gse211630.py` | T/NK module shift |
| `results/tables/geo_gse211630_celltype_annotation.csv` | `analysis/scripts/scrna_annotate_gse211630.py` | 19-subtype annotation |
| `results/tables/geo_gse211630_candidate_celltype_enrichment.csv` | `analysis/scripts/scrna_annotate_gse211630.py` | candidate enrichment |
| `results/tables/geo_gse211630_candidate_loc.csv` | `analysis/scripts/plot_scrna_candidates_gse211630.py` | candidate localisation (side product) |
| `results/tables/geo_gse211630_erosion_candidate_celltype.csv` | `analysis/scripts/scrna_erosion_gse211630.py` | 95 subtype + 40 family rows |
| `results/tables/geo_gse211630_erosion_per_sample.csv` | `analysis/scripts/scrna_erosion_gse211630.py` | 470 per-sample rows |
| `results/tables/geo_gse211630_mr_backmap_gene.csv` | `analysis/scripts/scrna_mr_backmap_gse211630.py` | gene-level back-mapping |
| `results/tables/geo_gse211630_mr_backmap_subtype.csv` | `analysis/scripts/scrna_mr_backmap_gse211630.py` | subtype-level back-mapping |
| `results/tables/geo_gse211630_liana_signal_shift.csv` | `analysis/scripts/scrna_liana_integrate_gse211630.py` | LIANA interaction shift |
| `results/tables/geo_gse211630_liana_incoming_focus.csv` | `analysis/scripts/scrna_liana_integrate_gse211630.py` | incoming focus signalling |
| `results/tables/geo_gse211630_qc_summary.csv` | `analysis/scripts/plot_gse211630_qc.py` | dataset QC summary |
| `data/processed/geo_GSE211630/scanpy/gse211630_olp_lognorm.h5ad` | `analysis/scripts/build_lognorm_allgenes_gse211630.py` | 46,478 x 17,669 log-normalised |
| `data/processed/geo_GSE211630/scanpy/gse211630_olp_obs_annotation.csv` | `analysis/scripts/scrna_annotate_gse211630.py` | per-cell annotation |
| `data/processed/geo_GSE211630/scanpy/gse211630_all_analyzed.h5ad` | `analysis/scripts/scanpy_cluster_gse211630.py` | all-sample analysed object |

### Tables without a retained producer (9 of 47)

- `results/tables/mr_olp_wide_results.csv`
- `results/tables/mr_olp_inclavo_results.csv`
- `results/tables/final_target_evidence.csv`
- `results/tables/phewas_olp_all.csv`
- `data/processed/targets/druggable_genes_final.csv`
- `data/processed/targets/druggable_x_eqtlgen.csv`
- `results/tables/md_sasa.csv`
- `results/tables/md_mmgbsa.csv`
- `results/tables/fpocket_summary.json`

These fall into three groups, and the gap is worth stating plainly:

1. **The MR screen itself** (`mr_olp_wide_results.csv`, `mr_olp_inclavo_results.csv`)
   and the druggable-set assembly (`druggable_genes_final.csv`,
   `druggable_x_eqtlgen.csv`) were produced before the scripts in this folder were
   formalised. The MR model is fully described in Methods 2.4, and the complete
   5,040-gene screen ships as Supplementary Table 1, so the numbers can be
   checked against it even though the driver script is missing.
2. **Two summary tables** (`final_target_evidence.csv`, `phewas_olp_all.csv`) were
   assembled by hand from the per-analysis outputs and cross-checked against them.
3. **Three post-processing summaries** (`md_mmgbsa.csv`, `md_sasa.csv`,
   `fpocket_summary.json`) were parsed from raw tool output (`MMPBSA.py`
   `FINAL_RESULTS` files, `gmx sasa`, `fpocket -o`) with throwaway commands. The
   underlying raw outputs are archived, so they can be rebuilt from those.

Everything downstream of these tables is scripted end to end.
