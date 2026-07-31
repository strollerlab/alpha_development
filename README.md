# Alpha Burst & Rhythm Development — Analysis Code

Reproducibility code for *"The emergence and maturation of the infant alpha peak
reflect a transition from transient bursts to sustained oscillations."*

Python, R, and MATLAB code that extracts the EEG metrics, builds the analysis
datasets, fits the statistical models, and renders every figure and table in the
manuscript and supplement. All the data to run the analysis can be found in the 
companion OSF data repository (https://osf.io/zhqke/).

This README covers (1) getting it running, (2) the output taxonomy, (3) the
supplementary-statistics pipeline, and (4) the naming convention. 
For the identifier ↔ manuscript-term lookup tables, see `NOMENCLATURE_CROSSWALK.md`.

---

## 1. Quick start

### Prerequisites
- **Python** 3.10 with `specparam`, MNE, and `bycycle`, and complementary Python codes. 
- **R** ≥ 4.2 with the packages installed by `00_Setup_PackageInstallation.R`
  (mgcv, gratia, lme4, lmerTest, emmeans, performance, flextable, openxlsx,
  tidyverse, etc.).
- **MATLAB** with EEGLAB and FieldTrip — topomaps only.
- The anonymized data from the OSF repository. 

### Configuration for the analysis

**Edit `config_paths.R`.** Every R script sources it as its first
action. No path appears anywhere else in the R codebase.

| Variable | Points at |
|---|---|
| `path2code` | The folder where the code lies (i.e., this folder)
| `path2data` | `Data/` — the merged cross-visit CSVs the analyses read. This is provided in the OSF repository |
| `path2sets` | `Data/` — per-visit CSVs written by `DatasetCreation_*` and the Python stage |
| `path2root` | where `Results/` should be created (need not exist yet) |

`config_paths.R` expands `~`, tolerates a missing trailing slash, creates
`path2root` if absent, and **stops with a named error** if any input directory is
missing. Any mistyped path fails immediately rather than halfway through a model
fit.

Each script locates `config_paths.R` itself: it checks the working directory, the
parent folder, and a `Code/` subfolder. If none of those work, it crashes. If you want
to store the code in different folders (i.e., from an unrelated directory), set `CODE_FOLDER` 
on the second line of the script, or `setwd()` to the folder where `config_paths.R` is. 
The script stops with a message naming both options and printing the directory it searched from.

`AlphaBurstRhythm.Rproj` is a convenience to avoid the error of not finding `config_paths.R`: double-clicking it m
akes RStudio set the working directory to the same folder. If all the code is in the same folder, it automatically 
finds `config_paths.R` even though the path is not specified. It is not required and has no effect outside RStudio.

Python and MATLAB have their own small config blocks at the top of the entry-point
files (`CODE_DIR` / `SET_DIR` / `DATA_DIR`, and `path2*` respectively). Set
`DATA_DIR` to the same folder as `path2sets` so the R stage finds the extractor
output.

### Run order

Scripts are numbered for execution order. For data privacy, clean EEG sets after MADE preprocessing are not available in the repository, so the Python stage is not runnable. This also limits running DatasetCreation_1a/1b/2a/2b/3, which require intermediate or Python outputs, and SupplementaryMethods_1. This doesn't mean these codes are not functional. Instead, our Python codes and DatasetCreation scripts are easily adaptable. With very few modifications, our Python codes can process clean-epoched EEG data from EEGlab, generating single-child outputs with the metrics of interest (and more). The rest of the pipeline can be run with the merged CSVs found in our OSF repository.

```
# --- Python feature extraction (§4a; run once per subject, before R) ---
EEG_metrics_rest_NN.py                          # imports jrpc + lagged_autocoherence
EEG_metrics_rest_NN_basedonbursts.py            # imports jrpc

00_Setup_PackageInstallation.R                  # install/verify packages, global constants
01_Utils_ProcFunctions.R                        # helpers (sourced by analysis scripts)

# --- Dataset creation (writes Data/ and Data/Merged/) ---
DatasetCreation_1a_ParametrizedPSD.R            # specparam params + PSDs -> long, per visit
DatasetCreation_1b_ParametrizedPSD_by_Burst.R   # same, conditioned on burst presence. It also runs alpha_lifespan by visit and burst.
DatasetCreation_2a_Burst_and_Lifespan.R         # cycle-by-cycle bursts + LAcH lifespan
DatasetCreation_2b_LAcH_Figures.R               # LAcH summary CSVs for figures, including merged and by burst.
DatasetCreation_3_Merge_Across_Visits.R         # merge per-visit files -> unified long sets
DatasetCreation_4_Topomaps.R                    # build topological_heatmaps_data.csv

# --- Analysis + main outputs ---
DataAnalysis_0_Sample_and_Metrics_Descriptives.R
DataAnalysis_1_WholeBraind_and_ROI_GAMM_Development.R
DataAnalysis_2_BurstImpact_MLM_Models.R
DataAnalysis_3_AlphaPeakPrediction_and_MLM_models.R
DataAnalysis_4_AlphaPeakPrediction_and_MLM_models_voltamp.R   # Extended Data Fig. 5 robustness check

# --- Figures ---
Figure_1_Code.R
Figures_2_3_4_Code.R
Figures_S1_S2_S3_topomaps_fieldtrip.m           # MATLAB; needs DatasetCreation_4 output
                                                # (topomaps now numbered Fig S1a/b, S2, S3 in the SI)
                         

# --- Supplementary methods/results ---
SupplementaryMethods_1_ParametrizedPSD_RangeSelection.R
SupplementaryMethods_2_ParametrizedPSD_FinalRange_Descriptives.R
SupplementaryMethods_3_Lifespan_RangeSelection.R
SupplementaryResuls_2_BurstDevelopment_adjEpoch.R

# --- Assemble the supplementary statistics file (LAST; §3) ---
SupplementaryData_BuildWorkbook.R
```

Every R script analysis/figure sources `config_paths.R`, then `00_Setup` and (where
needed) `01_Utils`, so any one can be run individually once the merged datasets
(i.e., the ones provided) exist. `SupplementaryData_BuildWorkbook.R` reads a cache of
other scripts, so it must run after them to generate the _Supplementary_Statistical_Tables.xlsx.

---

## 2. Output Folder Organization
Everything is written under a single `Results/` root, organised into three
top-level tiers that mirror the contents of the manuscript — **`MainText/`**, **`ExtendedData/`**,
and **`SupplementaryInformation/`** — each split into `Figures/` and `Tables/`.
`SupplementaryInformation/` is further divided into `GeneralInformation/`,
`Methods/`, `Results/` and `StatisticalTables/` to match the structure of the Supplementary Information
document. Sub-folders are derived from `path2root` and created on first run.

```
Results/
├── MainText/
│   ├── Figures/                    Fig 1, 2, 3, 4, 5
│   └── Tables/
│       ├── (Table 2a, 2b)          main-text sociodemographic tables
│       ├── Development/            Whole-brain GAMM development result tables
│       ├── BurstImpact/            Burst-impact MLM result tables
│       └── AlphaPrediction/        Peak-prediction + MLM-contribution tables
├── ExtendedData/
│   ├── Figures/                    Extended Data Fig. 1, 3, 4, 5
│   └── Tables/                     Extended Data Table 1, 2 (+ Fig. 5 robustness HTML)
└── SupplementaryInformation/
    ├── GeneralInformation/
    │   ├── Figures/                (topomap panels: Fig. S1a/b, S2, S3 — from MATLAB)
    │   └── Tables/                 Table S1, S2, S3, S4, S5, S6 + complementary 
    ├── Methods/
    │   ├── Figures/                Fig SM1, SM2, SM4   (+ SM3 quality topomaps but these need to be moved manually)
    │   └── Tables/                 Table SM1, range-comparison tables
    └── Results/
        ├── Figures/                Fig SR1, SR2
        └── Tables/
            ├── (ROI-stratified GAMM tables, adj-epoch tables)
            └── SupplementaryTables/   ★ Table_SR1.csv … Table_SM3.csv
    └── StatisticalTables/
            └── Supplementary_Statistical_Tables.xlsx   ★ the submitted file
```

### Organization

- **The three tiers match the manuscript's own division.** Items in
  Extended Data (Fig. 1, 3, 4, 5; Table 1, 2) are filed under `ExtendedData/`,
  separately from the Supplementary Information items.
- **Descriptive S-tables live in `SupplementaryInformation/GeneralInformation/Tables/`**
  carrying the Supplementary Information document's own numbering (Table S1–S6),
  independent of the main-text and Extended Data numbering.
- **Main-text result tables** split by analysis into `Development`, `BurstImpact`
  and `AlphaPrediction`, matching the three Results sections. The two main-text
  sociodemographic tables (Table 2a, 2b) sit directly under `MainText/Tables/` and need
  to be merged.
- **`SupplementaryInformation/Results/Tables/SupplementaryTables/`** holds the
  numbered statistical tables that accompany the submission (the `.xlsx`). The HTML
  tables elsewhere are unchanged and remain the at-a-glance view.
- A script emitting more than one category writes to more than one sink and declares
  the tiered path variables it needs (`path2main_tab`, `path2ed_tab`, `path2ed_fig`,
  `path2si_gi_tab`, and the existing `path2tabs` / `path2figs` / `path2suppres`),
  all derived from `path2root`.

### Where each script writes

| Script | Figures | Tables |
|---|---|---|
| `Figure_1_Code.R` | `MainText/Figures` (Fig 1) | — |
| `Figures_2_3_4_Code.R` | `MainText/Figures` (Fig 2–4), `ExtendedData/Figures` (Extended Data Fig. 3, 4) | — |
| `Figures_S3_S4_S5_topomaps_fieldtrip.m` † | `SupplementaryInformation/GeneralInformation/Figures` (Fig S1a/b, S2, S3; SM3 panels) | — |
| `DataAnalysis_0_…Descriptives` | `ExtendedData/Figures` (Extended Data Fig. 1) | `MainText/Tables` (Table 2a, 2b), `ExtendedData/Tables` (Extended Data Table 1a and 1b), `SupplementaryInformation/GeneralInformation/Tables` (Table S1, S2, S3, S5, S6) |
| `DataAnalysis_1_…GAMM_Development` | `SupplementaryInformation/Results/Figures` (Fig SR1) | `MainText/Tables/Development`, `SupplementaryInformation/Results/Tables`, **SR1, SR2, SR7, SR8** |
| `DataAnalysis_2_…BurstImpact_MLM` | `MainText/Figures` | `MainText/Tables/BurstImpact`, `ExtendedData/Tables` (Extended Data Table 2), **SR3, SR4** |
| `DataAnalysis_3_…AlphaPeakPrediction` | `MainText/Figures` (Fig 5) | `MainText/Tables/AlphaPrediction`, **SR5, SR6** |
| `DataAnalysis_4_…voltamp` (Extended Data Fig. 5) | `ExtendedData/Figures` (Extended Data Fig. 5) | `ExtendedData/Tables`, **SR11** |
| `SupplementaryMethods_1_…RangeSelection` | `SupplementaryInformation/Methods/Figures` (Fig SM1) | `SupplementaryInformation/Methods/Tables` |
| `SupplementaryMethods_2_…FinalRange` | `SupplementaryInformation/Methods/Figures` (Fig SM2) | `SupplementaryInformation/GeneralInformation/Tables` (Table S4) |
| `SupplementaryMethods_3_…Lifespan` | `SupplementaryInformation/Methods/Figures` (Fig SM4) | `SupplementaryInformation/Methods/Tables` (Table SM1), **SM2, SM3** |
| `SupplementaryResuls_2_…adjEpoch` | `SupplementaryInformation/Results/Figures` (Fig SR2) | `SupplementaryInformation/Results/Tables`, **SR9, SR10** |
| `SupplementaryData_BuildWorkbook.R` | — | `SupplementaryInformation/Results/Tables/SupplementaryTables/` (assembles all of the above) |


> **Note on the topomaps.** The MATLAB topomap script saves one JPEG per metric into
> `SupplementaryInformation/GeneralInformation/Figures/`. The aperiodic/oscillatory,
> burst and lifespan panels are Fig. S1a/b, S2 and S3; the R²/MAE quality panels
> belong to Fig. SM3.
---

## 3. Auxiliary Code to Generate Supplementary Table Statistics

The model-output tables are too large to typeset in the manuscript, so they ship as
a data file.

### `SupplementaryTables_Helper.R` — library, sourced (never run separatedly; it only contains functions).

| Function | Purpose |
|---|---|
| `nn_supp_table(data, id, cols, note)` | Registers **one numbered table**. Writes `Table_<id>.csv` and caches it for the workbook. Returns its input unchanged, so it drops into an existing pipe without disturbing the flextable/HTML branch. |
| `nn_gam_df(gam_obj)` | Residual df of a `gamm` fit, with three fallbacks — the denominator for parametric *t* tests. |
| `nn_lmer_df(model, term)` | Satterthwaite df for one `lmerTest` term; warns loudly if the model was fitted with `lme4::lmer` instead. |
| `nn_save_table(ft, dir, stem)` | Saves a flextable as HTML and optionally Word in one call. |

`NN_TABLE_INDEX` at the top of the file is the **single place** where table numbers,
captions, and column groupings are declared. Renumbering is one edit there.

### `SupplementaryData_BuildWorkbook.R` — run last

Reads the cache and typesets `Supplementary_Statistical_Tables.xlsx`: one sheet per
table with a caption row, spanning sub-headers where a table reports more than one
thing, bold rules, vertically merged key columns, and a Note. Applied automatically
from the column headers:

- `95% CI lower` + `95% CI upper` → a single `95% CI` column as `[l, u]`
- P columns → `<0.0001` rather than `0.0000`
- a significance column (`*/**/***`) after each FDR-corrected P
- rows ordered as the manuscript presents them — **power spectrum → bursts → alpha
  lifespan**, aperiodic before oscillatory (`NN_METRIC_ORDER`); Fig. 5d predictors
  keep their own panel order (`NN_PREDICTOR_ORDER`)

### The tables

| ID | Contents | Produced by |
|---|---|---|
| SR1 | Whole-brain GAMM developmental trajectories | `DataAnalysis_1` |
| SR2 | Burst × age difference smooths | `DataAnalysis_1` |
| SR3 | Burst-presence effects (Fig. 4) | `DataAnalysis_2` |
| SR4 | Burst vs non-burst contrasts by visit | `DataAnalysis_2` |
| SR5 | Alpha-peak classifier, 1 month (Fig. 5a–c) | `DataAnalysis_3` |
| SR6 | Contributions to oscillatory alpha (Fig. 5d) | `DataAnalysis_3` |
| SR7 | Regional (ROI) differences | `DataAnalysis_1` |
| SR8 | Regional (ROI) age smooths | `DataAnalysis_1` |
| SR9 / SR10 | Adjusted-epoch robustness | `SupplementaryResuls_2` |
| SR11 | Voltage-amplitude robustness | `DataAnalysis_4` |
| SM2 / SM3 | Alpha-lifespan frequency-range selection | `SupplementaryMethods_3` |

Every model reports its test statistic, degrees of freedom, and exact *P*.
Mixed-model df are Satterthwaite-approximated; GAMM smooths are tested with *F* on
the reference df, and their parametric terms based on the residual df, which is the same as EDF.

---

## 4.  Python Code - Files Generated

Two scripts to compute most metrics of interest based on the processed EEG extracting per-subject CSVs that `DatasetCreation_` 
and `SupplementaryMethods_1_RangeSelection.R` use. These codes are functional, and with some adjustments
can be run on other datasets.

| Module | Role |
|---|---|
| `convenience_functions_jrp.py` (imported as `jrpc`) | Shared toolkit: PSD estimation, lagged coherence, specparam trial selection, CSV writer. No `__main__`. |
| `lagged_autocoherence.py` (imported as `la`) | LAcH reference implementation, adapted from Zhang et al. (2025), *Imaging Neuroscience*. Imported by `EEG_metrics_rest_NN.py`. |
| `EEG_metrics_rest_NN.py` | **Data generated.** `aperosc_parameters*.csv`, `psds*.csv`, `burst_properties_bycycle.csv`, `lagged_coh_py_hilb*.csv` → feeds `DatasetCreation_1a`, `2a`. |
| `EEG_metrics_rest_NN_basedonbursts.py` | **Data generated.** Burst-conditioned `aperosc_parameters_burst*.csv`, `psds_burst*.csv`, `lagged_coh_py_hilb*_burst.csv` → feeds `DatasetCreation_1b`. |

---

## 5. Naming convention

`NOMENCLATURE_CROSSWALK.md` provides the correspondence among the column names and the manuscript terms. Below you can find the general rules:

| Domain | Canonical | Rationale |
|---|---|---|
| Proportions | `prop_*` | Every one is a ratio in **[0, 1]** |
| Unit of segmentation | `epoch` |
| Epoch count | `epochs` | Epoched resting-state EEG - clean epochs after processing. |
| ROI | `region` | The manuscript says *region of interest* / *region*. |
| Fit error | `mae` | The stored value is specparam's `fm.error_`, and `compute_error()` defaults to `error_metric='mae'` |
| Rhythmicity metric | `lifespan`  |
| Burst-corrected ratios | `<metric>_corrected` | Suffix, so the metric sorts next to its absolute counterpart (`volt_amp`, `volt_amp_corrected`). |

Display labels are the **exact Table 1 term**, so a reader can map a panel to a row
of Table 1.

### Note about the naming of the epoch-related variables

**`epochs` vs `epochsprop` vs `prop_epochs`.** All epoch-related, all different.
`epochs` is a *count* per electrode (model covariate, filtered against
`EPOCHS_THRESHOLD`). `epochsprop` is the *proportion* kept by Specparam's interim R²
selection, which was dropped because it did not apply, and `prop_epochs` is the *proportion* of
clean epochs per session if, from `CrossVisit_EEG_CleaningDescriptives.csv`, z-scored.

**`epoch` vs `epoch_block`.** `epoch` is the ordinal position within the **whole
recording**; `epoch_block` is the position **within its own condition block** and
restarts at 0 for every block, so it equals `epoch` only when the data are unblocked
(`blocks == 0`). Whichever is used must match how the PSD array was sliced: use
`epoch_block` *only* after filtering on `['block']`. Mixing them produces no error,
just a silently mislabelled burst/no-burst mask.


### Acronyms

LAcH = Lagged Autocoherence Hilbert · ROI = region of interest · GAMM = generalized
additive mixed model · MLM = multilevel mixed model · PSD = power spectral density ·
specparam = power-spectrum parametrisation · MAE = mean absolute error · EDF =
effective degrees of freedom · FDR = false discovery rate.
