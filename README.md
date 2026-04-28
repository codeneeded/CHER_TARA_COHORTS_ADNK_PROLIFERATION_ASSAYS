# CHER / TARA Cohorts — ADNK Assay QC and Normalization

End-to-end quality control and normalization pipeline for **Antibody-Dependent NK cell (ADNK)** activation assays measured by flow cytometry across two HIV cohorts:

- **CHER** — longitudinal pediatric cohort with a **treatment-interruption** design.
- **TARA** — cross-sectional cohort comparing **virally suppressed vs unsuppressed** participants.

The input is a single Excel workbook of FlowJo-exported frequencies (one sheet per plate / batch × antigen). The pipeline produces QC reports, batch-corrected and dual-normalized data, per-cohort analysis-ready tables, cleaned metadata, and a self-contained HTML summary report. Sample identifiers are auto-classified as TARA (`PID` starts with `CP`) or CHER otherwise.

---

## Assay design

| | |
|---|---|
| **Modality** | Flow cytometry (FlowJo exports) |
| **Cell populations measured** | NK CD56-bright, NK CD56-dim, non-NK |
| **Functional markers** | CD107a (degranulation), FasL, IFN-γ, MIP-1β |
| **Antigens** | p24 (Gag), gp120 (Env) — inferred from sheet name |
| **Plate / batch structure** | "Batch 1" vs "Batch 2" — inferred from sheet name |
| **Controls** | PBS (no-antigen background), HIV-negative plasma, **N6 bNAb** (positive control), NK stained, Unstained |

The **primary readouts** the pipeline focuses on are NK CD56-dim CD107a and IFN-γ — the canonical ADCC degranulation and cytokine readouts.

### QC thresholds

Hard-coded thresholds applied during QC:

| Threshold | Value |
|---|---|
| Minimum N6 CD107a response | 0.05 |
| Maximum PBS background | 0.15 |
| Minimum event count | 5,000 |
| Minimum lymphocyte frequency | 0.50 |
| Maximum batch fold-difference | 3.0 |

---

## Repository structure

```
CHER_TARA_COHORTS_ADNK_PROLIFERATION_ASSAYS/
├── ADNK_Data_CHER_and_TARA_with_metadata_040826.xlsx   # Input: FlowJo exports + metadata
├── R_Scripts/
│   └── ADNK_QC.R                                       # Single end-to-end pipeline (~1,500 lines)
└── ADNK_QC_Output/                                     # All pipeline outputs
    ├── 01_QC_Reports/
    │   ├── Control_Statistics_Full.csv                 # Per-plate control means / CVs
    │   ├── QC1_Control_Hierarchy.{csv,txt}             # PBS < HIV-neg < N6 ordering check
    │   ├── QC2_N6_Performance.csv                      # Positive-control sensitivity
    │   ├── QC3_Batch_Effects.csv                       # Batch 1 vs Batch 2 fold-difference
    │   └── QC4_Sample_Summary.csv                      # Sample-level pass/fail
    ├── 02_QC_Plots/
    │   ├── 01_Control_Performance.png
    │   ├── 02_Batch_Distribution_RAW.png
    │   ├── 03_Event_Counts.png
    │   └── 04_Lymphocyte_Freq.png
    ├── 03_Batch_Correction/
    │   ├── 01_Batch_Correction_Comparison.png
    │   ├── 02_Batch_Effect_Statistics.csv
    │   └── 03_Normalized_Density_Overlay.png
    ├── 04_Normalized_Data/
    │   ├── ADNK_All_Normalized.xlsx                    # Combined cohorts, all normalization columns
    │   ├── ADNK_Analysis_Ready.xlsx                    # Combined, analysis-ready
    │   ├── CHER_Analysis_Ready.xlsx                    # CHER-only analysis-ready table
    │   └── TARA_Analysis_Ready.xlsx                    # TARA-only analysis-ready table
    ├── 05_Cleaned_Metadata/
    │   ├── CHER_Metadata.csv                           # Per-PID cleaned metadata
    │   ├── CHER_Timepoint_Mapping.csv                  # Raw timepoint → numeric mapping
    │   ├── TARA_Metadata.csv
    │   └── Missing_Data_Summary.csv
    └── 06_Summary/
        ├── ADNK_QC_Report.html                         # Self-contained report (embedded images)
        └── QC_Summary_Table.csv
```

---

## Pipeline overview

The full workflow is implemented in [`R_Scripts/ADNK_QC.R`](R_Scripts/ADNK_QC.R).

### 1 – Data ingest
All sheets of `ADNK_Data_CHER_and_TARA_with_metadata_040826.xlsx` are read with `readxl`, type-coerced, and concatenated. Each row is annotated with:
- `Sheet` / `Plate` (the original FlowJo sheet name)
- `Batch` (`Batch 1` / `Batch 2`, inferred from sheet name)
- `Antigen` (`p24` / `gp120`, inferred from sheet name)
- `Cohort` (`TARA` if PID begins with `CP`, else `CHER`)
- `IsControl` / `ControlType` (PBS, HIV_neg, N6_bnAb, Staining)

Cohort-specific timepoint handling: TARA's "Timepoint relative to phase" column is treated as **age in months**; CHER timepoints (e.g. `0 B`, `12 WK`, `14 VST`, `0WK`) are parsed into a numeric scale.

### 2 – Control statistics
Per-plate (`Sheet × Batch × Antigen × ControlType`) means and dispersion are computed across all 12 functional readouts and exported as `Control_Statistics_Full.csv`.

### 3 – Four QC checks

1. **QC1 — Control hierarchy.** Verifies the expected biological ordering on each plate: `PBS < HIV-neg plasma < N6 bNAb`. Plate-level pass/fail with reason codes → `QC1_Control_Hierarchy.{csv,txt}`.
2. **QC2 — N6 performance.** Confirms the N6 broadly neutralizing antibody positive control reaches the minimum CD107a response (0.05). Plates where N6 fails to activate are flagged. → `QC2_N6_Performance.csv`.
3. **QC3 — Batch effects.** Per-marker fold difference between Batch 1 and Batch 2 against the 3.0× threshold. → `QC3_Batch_Effects.csv`.
4. **QC4 — Sample-level.** Event-count and lymphocyte-frequency filters per sample, with cohort-stratified failure summaries. → `QC4_Sample_Summary.csv`.

### 4 – Normalization (two stages)

1. **Background subtraction** (`*_BGsub`). For every sample, the **plate-matched PBS control mean** is subtracted from each functional marker; negatives are clamped to zero. This removes plate-specific autofluorescence/spontaneous background.
2. **N6 normalization** (`*_NormN6`). Each background-subtracted marker is expressed as a **percentage of the N6 maximum response on the same plate**:

   ```
   *_NormN6 = ((sample - PBS) / (N6 - PBS)) * 100
   ```

   Values are clipped to [0, 500]. This yields a within-plate, positive-control–scaled metric that is comparable across batches and antigens — the headline column for downstream analysis.

### 5 – Batch correction confirmation
Diagnostic plots and statistics show the residual batch effect after background subtraction and N6 normalization (`03_Batch_Correction/`), so it can be verified before any downstream modeling.

### 6 – QC plots
Plate-level control performance, raw batch distribution, event counts, and lymphocyte frequencies are saved to `02_QC_Plots/` as ggplot/patchwork figures.

### 7 – Cleaned metadata export
Per-cohort metadata tables (`CHER_Metadata.csv`, `TARA_Metadata.csv`) are written alongside the CHER timepoint mapping and a `Missing_Data_Summary.csv` flagging participants with incomplete records.

### 8 – Analysis-ready tables
The cleaned, normalized data are exported in three forms:
- `ADNK_All_Normalized.xlsx` — every column (raw, BGsub, N6 means, NormN6) for full audit.
- `ADNK_Analysis_Ready.xlsx` — combined cohorts, downstream-ready.
- `CHER_Analysis_Ready.xlsx` / `TARA_Analysis_Ready.xlsx` — cohort-split versions for cohort-specific modeling (longitudinal mixed effects for CHER, cross-sectional comparisons for TARA).

### 9 – HTML summary
A single self-contained `ADNK_QC_Report.html` (with embedded images, no external file dependencies) summarizes the run for sharing.

---

## Dependencies

Base R plus:

- [`readxl`](https://CRAN.R-project.org/package=readxl), [`writexl`](https://CRAN.R-project.org/package=writexl) — Excel I/O
- [`dplyr`](https://CRAN.R-project.org/package=dplyr), [`tidyr`](https://CRAN.R-project.org/package=tidyr), [`stringr`](https://CRAN.R-project.org/package=stringr) — wrangling
- [`ggplot2`](https://ggplot2.tidyverse.org/), [`patchwork`](https://patchwork.data-imaginist.com/), [`scales`](https://CRAN.R-project.org/package=scales) — visualization

Missing packages are auto-installed at the top of the script.

---

## Reproducing the analysis

1. Clone the repo and open `CHER_TARA_COHORTS_ADNK_PROLIFERATION_ASSAYS.Rproj` in RStudio.
2. At the top of `R_Scripts/ADNK_QC.R`, update `BASE_DIR` to point at your local clone (the script currently has `/home/akshay-iyer/Documents/CHER_TARA_COHORTS_ADNK_PROLIFERATION_ASSAYS`).
3. Source `R_Scripts/ADNK_QC.R`. The full `ADNK_QC_Output/` tree will be regenerated, including the HTML summary in `06_Summary/`.

The pipeline is single-script and idempotent — re-running overwrites previous outputs without manual intervention.
