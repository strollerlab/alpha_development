# Nomenclature Crosswalk

Equivalence tables only: **code identifier ↔ display label ↔ manuscript term**.
For the reasoning behind the convention, see `README.md` §5.

**Source of truth:** Rico-Picó et al., **Table 1** ("Aperiodic and oscillatory alpha
metrics used in this study") + Online Methods.

---

## A. Parametrised power spectrum

| Code identifier | Display label | Manuscript (Table 1) |
|---|---|---|---|
| `offset` | "Offset" | **Offset** |
| `slope` | "Slope" | **Slope** (aperiodic exponent) |
| `alpha_freq` to `peak_freq` | "Peak Freq." | **Peak Frequency** |
| `alpha_ampl` to `peak_ampl` | "Peak Amp." | **Peak Amplitude** |
| `alpha_peak` to `peak_prop` | "Prop. of Peaks" | **Prop. of Peaks** |
| `alpha_osc` to `osc_ampl` | "Band Power" | **Band Power** |
| `alpha_widt` | Gaussian peak width |
| `prop_nopeak` | Prop. of electrodes with no peak |
| `prop_thetapeak` | Prop. of electrodes with theta peak |
| `prop_alphapeak` | Prop. of electrodes with alpha peak |
| `prop_bothpeaks` | Prop. of electrodes with both peaks |
| `r2value` | "Model Fit (R²)" | specparam fit R² |
| `mae` | "MAE" | **MAE** (mean absolute error) |

> `alpha_*` is the per-electrode raw column; the right-hand name is the aggregated
> summary used in models, tables and figures.

## B. Cycle-by-cycle burst analysis

| Code identifier | Kind | Display label | Manuscript (Table 1) |
|---|---|---|---|
| `volt_amp` | "Volt. Amp." | **Voltage Amplitude** (Absolute) |
| `volt_amp_corrected` | "Corrected Volt. Amp." | **Voltage Amplitude** (Corrected) |
| `band_amp` | "Band Amp." | **Band Amplitude** (Absolute) |
| `band_amp_corrected` | "Corrected Band Amp." | **Band Amplitude** (Corrected) |
| `prop_bursty_epochs` | "Prop. of Epochs w/ Burst" | **Prop. of Epochs w/ Burst** |
| `prop_bursty_cycles_burst` | "Prop. of Cycles w/ Burst" | **Prop. of Cycles w/ Burst** |
| `prop_bursty_cycles` | Prop. of cycles, averaged over *all* epochs |
| `avg_burst_duration` | "Burst Duration" | **Burst Duration** |
| `is_burst` / `burst` / `burst_type` | "Cycle Type" ("Burst" / "NoBurst") | burst vs non-burst cycle |
| `bursty_epochs` | Per-epoch binary burst flag (intermediate) |
| `rise_decay_asym` / `peak_trough_asym` | Cycle waveform asymmetry |

## C. Rhythmicity (Lagged Autocoherence Hilbert)

| Code identifier | Display label | Manuscript (Table 1) |
|---|---|---|---|
| `alpha_LAcH` | "Lifespan" | **Lifespan** |
| `LAcH` | "LAcH" | **LAcH** (Lagged Autocoherence Hilbert) |
| `cLAcH` | LAcH normalised per electrode × frequency |
| `LAcH_lifespan` | Lifespan container (in-memory) |
| `lagged_hilbert_autocoherence()` | **LAcH** (Methods; Fig. 1E) |

## D. Epoch counts and proportions

Four related quantities. They are **not** interchangeable.

| Code identifier | Source | Meaning |
|---|---|---|---|
| `epochs` | Python Code | **Count** of retained artefact-free epochs, per electrode. Model covariate; filtered against `EPOCHS_THRESHOLD`. |
| `epochsprop` | Python Code | **Proportion** kept by the interim R² selection (`epochs_included / num_epochs`). Not used in this study|
| `totalepochs` | Python Code | Epoch count **before** selection. Not used in this study. |
| `prop_epochs` | `CrossVisit_EEG_CleaningDescriptives.csv` | **Proportion** of clean epochs, per session. Used as z-scored covariate. |
| `clean_epochs_rest` | `CrossVisit_EEG_CleaningDescriptives.csv` | Count of clean epochs, per visit. |
| `nepochs_group` | Python Code | Epochs per fitting group (Specparam setting) if two-stage fitting was indicated. Not used in this study. |

## E. Cross-cutting

| Code identifier | Display label | Meaning |
|---|---|---|---|
| `sujid` | Participant identifier |
| `session_age` | "Visit" | Visit age in months (1, 6, 12, 15→18, 30, 36, 42, 48) |
| `age_months` | "Age (months)" | Exact chronological age at visit. Used in the GAMM smooth term |
| `ch` | Electrode label |
| `region` | "Region" | ROI: Central, Frontal, Occipital, Parietal, Temporal |
| `hemis` | Hemisphere. Not used in this study|
| `chinclu` | Electrode-inclusion flag. Filter to 1 retain the 60 electrodes of the study |
| `goodch` | Count of good electrodes per participant. Used mostly in Specparam |
| `epoch` | Epoch index **within the whole recording**. Used in Python mostly to select burst and non-burst epochs |
| `epoch_block` | Epoch index **within its own condition block** |
| `block` | Condition block label |
| `freq` | "Frequency (Hz)" | Frequency bin |
| `inclusion_lmm` | Participant retained for longitudinal models |

## F. Statistical output columns

Emitted by `nn_supp_table()` into the numbered supplementary tables (`README.md` §3).
Display labels are what appears in `Supplementary_Statistical_Tables.xlsx`.

| Code identifier | Display label | Meaning |
|---|---|---|
| `age_edf`, `diff_edf`, `k_edf` | "EDF" | Effective degrees of freedom of a smooth (>1 = non-linear) |
| `age_refdf`, `diff_refdf` | "Ref. df" | Reference df — the numerator df of the smooth's *F* test |
| `age_F`, `diff_F` | "F" | *F* statistic for the smooth |
| `epochs_estimate`, `burst_estimate`, `Burst_estimate`, `estimate` | "Beta …" | Fixed-effect coefficient |
| `epochs_SE`, `Burst_SE`, `std.error` | "s.e." | Standard error |
| `epochs_t`, `Burst_t`, `Visit_t`, `statistic` | "t" | *t* statistic (lmerTest / mgcv parametric term) |
| `epochs_df`, `Burst_df`, `Visit_df`, `df` | "df" | Satterthwaite df (MLM) or residual df (GAMM parametric term) |
| `z_val` | "z (mean)" | Wald *z* — logistic classifier only; **no df applies** |
| `Boot_CI_Lower` / `Boot_CI_Upper` | "95% CI" | Merged to `[l, u]` by the workbook builder |
| `p.value`, `*_pval` | "P" | Uncorrected *P* |
| `p_fdr`, `*_pval_fdr` | "P (FDR)" | Benjamini–Hochberg corrected *P* |
| `partial_r2`, `effect_size` | "Partial R2" | (RSS_reduced − RSS_full)/RSS_reduced for the age smooth |
| `R2_Marginal` / `R2_Conditional` | "Marginal R2" / "Conditional R2" | Nakagawa *R*² (fixed / fixed + random) |
| `n`, `total_obs` | "n (children)" / "n (observations)" | Effective sample size for that model |

## G. Supplementary table IDs

Declared in `NN_TABLE_INDEX` (`SupplementaryTables_Helper.R`) — the single place to
renumber. `SR*` continues the Supplementary Results series (Fig. SR1–SR2); `SM*`
continues Supplementary Methods (Table SM1).

| ID | Contents | Script |
|---|---|---|
| SR1 | Whole-brain GAMM developmental trajectories | `DataAnalysis_1` |
| SR2 | Burst × age difference smooths | `DataAnalysis_1` |
| SR3 | Burst-presence effects (Fig. 4) | `DataAnalysis_2` |
| SR4 | Burst vs non-burst contrasts by visit | `DataAnalysis_2` |
| SR5 | Alpha-peak classifier at 1 month (Fig. 5a–c) | `DataAnalysis_3` |
| SR6 | Contributions to oscillatory alpha by visit (Fig. 5d) | `DataAnalysis_3` |
| SR7 | Regional (ROI) differences | `DataAnalysis_1` |
| SR8 | Regional (ROI) age smooths | `DataAnalysis_1` |
| SR9 | Adjusted-epoch robustness: trajectories | `SupplementaryResuls_2` |
| SR10 | Adjusted-epoch robustness: difference smooths | `SupplementaryResuls_2` |
| SR11 | Voltage-amplitude robustness: contributions | `DataAnalysis_4` |
| SM2 | Alpha-lifespan range selection (whole sample) | `SupplementaryMethods_3` |
| SM3 | Alpha-lifespan range selection by visit | `SupplementaryMethods_3` |

## H. Path variables

Set once in `config_paths.R`; every other path is derived from these four.

| Variable | Points at | Aliases used inside scripts |
|---|---|---|
| `path2code` | the folder where the code lies |
| `path2data` | `Data/Merged/` |
| `path2sets` | `Data/` | `path2save`, `path2freqs`, `path2psd`, `path2subsets` |
| `path2root` | `Results/` | `path2figs`, `path2tabs`, `path2desc`, `path2suppfig`, `path2suppres`, `path2save` (results context) |

## I. Constants (`00_Setup_PackageInstallation.R`)

| Constant | Value | Applied to |
|---|---|---|
| `EPOCHS_THRESHOLD` | 5 | `epochs >= EPOCHS_THRESHOLD` |
| `R2_THRESH` | 0.900 | `r2value > R2_THRESH` |
| `MAE_THRESH` | 0.10 | `mae < MAE_THRESH` |
| `CH_THRESHOLD` | 12 | `goodch >= CH_THRESHOLD` |

## J. Convention

| Domain |
|---|---|
| Proportions | `prop_*` (a ratio in [0, 1]) |
| Unit of segmentation | `epoch` |
| Epoch count | `epochs` |
| ROI | `region` |
| Fit error in Specparam | `mae` |
| Rhythmicity metric | `lifespan` |
| Rhythmicity acronym during computation | `LAcH` |
| Burst-corrected ratios | `<metric>_corrected` |
| Identifier style | `snake_case` (except the `LAcH` acronym) |
| Supplementary table IDs | `SR*` / `SM*`, declared in `NN_TABLE_INDEX` |
| Helper functions | `nn_*` prefix, to keep them distinct from package functions |

## K. Exceptions

| Scope | Note |
|---|---|
| `lagged_autocoherence.py` | Vendored from Zhang et al. (2025), *Imaging Neuroscience*; byte-identical to upstream apart from comments. Internal names follow the original publication, not this convention. |
| `lcoh` / `lcoh_avg` | Created by the Python code only. Harmonised to `LAcH` on read in `DatasetCreation_2a` |
