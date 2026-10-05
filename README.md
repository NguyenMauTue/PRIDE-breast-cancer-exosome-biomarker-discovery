# Transparent Weighting of Heterogeneous Evidence for Auditable Candidate Prioritization

Code for **AHP-CDS**, the method introduced in the preprint of the same title.

<!-- [?] One sentence: what this repository is. Author's own words. -->

- **Preprint:** https://doi.org/10.64898/2026.05.14.725271 (this README describes the code as of **v7**)
- **Earlier version:** the README and code of preprint v1 are preserved in archived, preprint-v1 folder. They describe a different analysis (AGRN-centred, weight-perturbation robustness labels) and do not match the current preprint.

---

## What this repository does

<!-- [?] Author's section. The contribution as he states it (explicit, auditable multi-criteria integration rule), and the three research questions in his words. -->

---

## Reproducing the analysis

**Environment.** R 4.5.2 / Bioconductor 3.22 (`renv.lock`, 205 packages pinned).

```r
# open AHP-CDS.Rproj, then
renv::restore()
source("Master_run_all.R")   # runs the preamp, the main pipeline and the ablation pipeline
```

**Order.**
0. `R/Preamp/fetch_annotation_table.R`, the *preamp* (short for preamplifier: the step that runs first). It has to run before the main pipeline, because the main pipeline reads the annotation table it builds (biomaRt query). `Master run_all.R` calls it first.
1. `PXD056161/scripts/run_all.R`, scripts 01–12: data → QC → imputation → limma → GSEA → STRING network → AHP-CDS scoring.
2. `Ablation_Table/scripts/run_all.R`: baselines, RWR comparator, RQ1, RQ2, and the tier/G* analysis. It reads the outputs of the main pipeline in `PXD056161/results/`.

**Input data.** PRIDE PXD056161, MDA-MB-231 vs MCF-10A, nanoLC-MS arm, n = 3 per condition. Script 01 downloads the file from PRIDE and applies the MaxQuant-based filtering, so nothing has to be placed by hand.

**External resources.** STRING, Reactome, Ensembl BioMart and CTD were queried in August 2026. Live services change, so a re-run later may not reproduce the same pool. The CTD export in `Ablation_Table/data/CTD_curated_genes_diseases.csv` is the frozen snapshot used for the paper.

**Random seeds.** `30032026` for the main pipeline (imputation, GSEA); `21082026` for the tier-3 permutation test and the null model in `RQ5_tier_audit_CTD.R`.

---

## Repository structure

```
├── Master run_all.R           # entry point
├── AHP-CDS.Rproj
├── renv.lock
├── PXD056161/                 # main pipeline
│   ├── data/
│   ├── scripts/               # 01–12 + run_all.R
│   └── results/{figures,tables}
├── Ablation_Table/            # RQ1, RQ2, RQ3 analyses
│   ├── data/                  # candidate_pool_merged.csv, CTD_curated_genes_diseases.csv
│   ├── scripts/               # 01, 02, RQ1, RQ2, RQ5 + run_all.R
│   └── results/{figures,tables}
├── R/
│   ├── Helper/                # ahp_weights, network_helper, string_api_helper, ranking_comparison_utils, ...
│   └── Preamp/                # fetch_annotation_table.R, runs before the main pipeline
├── Visualization/             # Fig 1, Fig 2, Fig S1 scripts (figures exported by hand)
└── Archived/                  # superseded work, see below
```

`R/Helper/cross_dataset_helper.R` is only used by scripts in `Archived/`. Some other files under `data/` and `results/` (for example `PXD056161/results/tables/CDS_candidates_themed.csv`) are not produced or read by the current scripts; they are left in place on purpose.

**Archived/** holds work that is not part of the current preprint's results: earlier ablation and null-model experiments (GOSemSim / ECM gene-set nulls, network randomization, a 2×2 double null, a G–R circularity check), the v1 theme-classification scripts, the v1 cross-dataset validation on PXD012162, and early drafts of the figure scripts. Kept for the record, not maintained.

---

## From script to paper

Mapping below is **proposed, unconfirmed**; S-numbers refer to the Supplementary Methods.

| Script | Does | Paper |
|---|---|---|
| `R/Preamp/fetch_annotation_table.R` | biomaRt annotation table (GO terms, gene names) | S1 |
| `PXD056161/scripts/01_load_and_filter_data.R` | downloads the PRIDE file, filters by MaxQuant output | S1 |
| `02_metadata_and_matrix.R` | sample metadata, LFQ matrix | S1 |
| `03_quality_control.R` | correlation matrix, PCA of the variance structure | S1 |
| `04_missingness_analysis.R` | missingness patterns | S2 |
| `05_imputation.R` | left-shifted Gaussian imputation, imputation QC | S2.5, Fig S1 |
| `06_Differential_expression_analysis.R` | limma, Tumor − Normal | S3 |
| `07_Robustness_analysis.R` | complete-case vs imputed comparison. "Robustness" here means robustness to imputation; it has nothing to do with the weight-perturbation labels of v1 | S3 |
| `08_Pathway_enrichment_analysis.R` | Reactome GSEA | S4 |
| `09_Extract_genes_for_STRING.R` | candidate pool, contaminant filter | S5 |
| `10_Protein-protein_interaction_network.R`, `11_Network-expression_integration.R` | STRING network; 11 joins the limma table to the network table before scoring | S6 |
| `12_CDS_via_AHP.R` | AHP weights, CDS | S7, Methods |
| `Ablation_Table/scripts/01_ablation_conditions.R` | baseline rankings (FC-only, Centrality-only, Equal-weight), shared by all three RQs | RQ1–RQ3 |
| `02_ren2019_gr.R` | RWR comparator | S8 |
| `RQ1_ranking_difference.R` | rank-level comparison | RQ1, §5.1 |
| `RQ2_mechanism_analysis.R` | contribution shares, LOCO | RQ2, §5.2 |
| `RQ5_tier_audit_CTD.R` | CTD tiers, G*, null model | RQ3, §5.3, S9 |
| `Visualization/Mega_Figure_1.R` | | Fig 1 |
| `Visualization/Mega_Figure_2.R` | | Fig 2 |
| `Visualization/Supplementary_Figure.R` | | Fig S1 |

> **Naming note:** file names `RQ5_*` correspond to **RQ3** in the paper. The files kept their original names.

---

## Main output files

`PXD056161/results/tables/`
- `annotated_gene_pool.csv` — candidate pool after annotation and contaminant filtering
- `network_summary.csv` — candidates retained in the STRING network
- `CDS_candidates_final_result.csv` — AHP-CDS scores and ranks
- `differential_expression_imputed.csv`, `complete_case_DEA.csv`, `reactome_gsea_full.csv`, `reactome_gsea_filtered.csv`, `limma_network_table.csv`, and QC tables

`Ablation_Table/results/tables/`
- `ablation_ranks_partial.csv` — ranks under all baselines
- `RQ1_*.csv`, `RQ2A_*.csv`, `RQ2B_*.csv` — RQ1 and RQ2 tables
- `candidates_tiered_CTD.csv`, `RQ5_Gstar_sweep_results.csv`, `RQ5_null_distribution.csv` — tiers, G* sweep, null
- `rwr_ranks_full.csv`, `ren2019_seed_genes.csv` — RWR comparator

Figures are drawn by the scripts in `Visualization/` and exported by hand, so that resolution (dpi) and size can be set per figure. For that reason the `results/figures/` folders are mostly empty; they are kept only because the scripts still point to those paths.

---

## Scope and limitations

**Why there is no weight-sensitivity analysis.** Preprint v1 perturbed the AHP weights and labelled candidates by how stable their rank stayed. That analysis was removed: shaking the weights tests the researcher's original judgment, and the premise of AHP-CDS is that this judgment is made a priori and stated openly. The leave-one-criterion-out (LOCO) analysis that remains is also a form of sensitivity check, but it asks a different question: how much each criterion contributes to the rank of each protein, not whether the researcher's judgment was right.

<!-- [?] Author's section, rest. Facts available to build from: n = 3 per condition; exploratory framing; no wet-lab validation; the Case-Study Selection Audit exists in Archived/ but is exploratory and not used in the paper. -->

---

## Citation, disclosure, license

- **Cite:** the preprint (DOI at the top of this file)
- **LLM use:** Large language models (Claude, GPT, and Gemini) were used during manuscript preparation in an adversarial-review capacity: to stress-test analytical claims, check internal consistency between reported statistics and their interpretation, and assist with language editing. All scientific content, analyses, and conclusions are the author's own; no manuscript text was generated from a prompt without direct authorship and verification. Further detail on this workflow is provided in the Supplementary Material.
  This README was drafted with the assistance of Claude (Anthropic), which proposed its structure and the script-to-paper table from the repository's file listing. The author supplied the content and decisions and is responsible for it.
- **License:** MIT, see `LICENSE`
- **Author:** Nguyễn Mậu Tuệ, Faculty of Biology, VNU University of Science, Hanoi
