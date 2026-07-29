# ---------------------------------------------------------------------------
# PATHS: edit config_paths.R once; nothing in this file needs changing.
# ---------------------------------------------------------------------------
# config_paths.R is looked for in the working directory. If R was started
# somewhere else, set CODE_FOLDER on the next line to this script's folder.
CODE_FOLDER = ""          # e.g. "~/AlphaBurstRhythm/Code"  If you have open the code from the project, you don't need to modify this line. Otherwise, select where the code folder that contains the config_paths.R is

local({
  cand = c(if (nzchar(CODE_FOLDER)) file.path(path.expand(CODE_FOLDER), "config_paths.R"),
           "config_paths.R", "../config_paths.R", "Code/config_paths.R")
  hit  = cand[file.exists(cand)]
  if (!length(hit)) # If the code cannot find the config_path.R it will stop the execution avoiding crashing.
    stop("config_paths.R not found.\n",
         "Fix either way:\n",
         "  1. setwd(\"/path/to/AlphaBurstRhythm/Code\")   then re-run, or\n", # Solution proposed 1: just add the directory of the code
         "  2. set CODE_FOLDER at the top of this script to that same path.\n", #Solution proposed 2: set the CODE_FOLDER variable to the directory of the code
         "Currently looking from: ", getwd(), call. = FALSE)
  source(hit[1], local = FALSE)
})

source(file.path(path2code, "00_Setup_PackageInstallation.R"))

# =============================================================================
# SECTION 1: ANALYSIS PARAMETERS
# =============================================================================
# Local copies of global thresholds from 00_Setup_PackageInstallation.R.
# These allow the script to be run standalone without sourcing 00_Setup first.

r2_thresh        = R2_THRESH         # Minimum Specparam R-squared (0.900)
mae_thresh       = MAE_THRESH        # Maximum Specparam MAE (0.10)
epochs_threshold = EPOCHS_THRESHOLD  # Minimum clean epochs per electrode (5)


# =============================================================================
# SECTION 2: PATH DEFINITIONS
# =============================================================================
# IMPORTANT: Update only `path2root` to match your local directory/
# All output sub-folders below derive from it and mirror the manuscript's own
# organisation (see README / Output folder taxonomy).
#
# path2data      : merged/analysis-ready input data (read-only here)
# path2main_tab  : MainText/Tables/                         .../ Table 2a, 2b
# path2ed_tab    : ExtendedData/Tables/                     .../ Extended Data Table 1
# path2ed_fig    : ExtendedData/Figures/                    .../ Extended Data Fig. 1
# path2si_gi_tab : SupplementaryInformation/GeneralInformation/Tables/
#                                                           .../ Table S1, S2, S3, S5, S6 + complementary

path2data = path2sets   # EDIT   # = Data/; this script prepends "" to each file name
# path2root comes from config_paths.R

# Manuscript-congruent output destinations (see README, Output folder taxonomy).
path2main_tab  = file.path(path2root, "MainText", "Tables")
path2ed_tab    = file.path(path2root, "ExtendedData", "Tables")
path2ed_fig    = file.path(path2root, "ExtendedData", "Figures")
path2si_gi_tab = file.path(path2root, "SupplementaryInformation", "GeneralInformation", "Tables")

# Backward-compatible aliases retained so downstream calls that were not repointed
# still resolve; new calls should use the tiered names above.
path2desc    = path2main_tab   # main-text Table 2 sink
path2supptab = path2si_gi_tab  # SI general-information table sink
path2suppfig = path2ed_fig     # this script's only figure (Extended Data Fig. 1)
path2save = path2main_tab
path2figs = path2ed_fig

for (p in c(path2main_tab, path2ed_tab, path2ed_fig, path2si_gi_tab))
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)

# =============================================================================
# SECTION 3: DATA LOADING & PREPARATION
# =============================================================================

# --- Electrode map ---
# Expands short region codes to full names; retains only the 60 analysis channels.
electrodes = read_csv(file.path(path2data, "electrodes.csv")) |> # EDIT the path
  dplyr::select(label, region, hemis, chinclu) |>
  # Recode raw 2-letter region codes to full anatomical region names
  mutate(region = case_when(
    region == "Fr" ~ "Frontal",
    region == "P"  ~ "Parietal",
    region == "T"  ~ "Temporal",
    region == "O"  ~ "Occipital",
    TRUE         ~ "Central"
  )) |>
  filter(chinclu == 1) |>
  # Rename 'region' to 'region' for consistency across downstream operations
  rename(region = region)

# Variables of interest from the Specparam output file.
aper_voi = c("sujid", "session_age", "ch", "region", "chinclu", "epochs",
             "r2value", "alpha_peak", "inclusion_final_dummy",
             "alpha_ampl", "alpha_osc", "alpha_freq", "slope", "offset", "mae", "goodch")

# --- Aperiodic / oscillatory data ---
# Filters: R² > threshold, MAE < threshold, included channel, minimum clean epochs,
# at least 12 good channels per subject, and excluded subjects removed.
aperiodic_data = read_csv(file.path(path2data, "Aperiodic_Oscillatory_ByCycle_Long.csv")) |>
  dplyr::select(all_of(aper_voi)) |>
  filter(
    goodch >= CH_THRESHOLD,
    r2value > r2_thresh,
    chinclu == 1,
    epochs >= epochs_threshold,
    mae < mae_thresh,
    !sujid %in% EXCLUDED_SUBJECTS
  ) |>
  # Remap session labels: 15 to 18 mo., 40 to 42 mo.
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  ) |>
  # Rename 'region' to 'region' for consistency across datasets
  rename(region = region) |>
  group_by(sujid, session_age, region) |>
  # Collapse electrode-level values to subject×visit×region level by averaging
  summarise(
    slope      = mean(slope,                       na.rm = TRUE),
    offset     = mean(offset,                      na.rm = TRUE),
    peak_ampl  = mean(alpha_ampl,                  na.rm = TRUE),
    # peak_freq computed only from cycles with alpha peak present (alpha_peak == 1)
    peak_freq  = mean(alpha_freq[alpha_peak == 1], na.rm = TRUE),
    osc_ampl   = mean(alpha_osc,                   na.rm = TRUE),
    # peak_prop: proportion of cycles with alpha peak
    peak_prop   = mean(alpha_peak,                  na.rm = TRUE),
    epochs     = mean(epochs,                      na.rm = TRUE),
    r2value    = mean(r2value,                     na.rm = TRUE),
    goodch     = n(),
    goodch_prop = n() / length(unique(ch)),
    .groups    = "drop"
  ) |>
  group_by(sujid, session_age) |>
  mutate(goodch_sum = sum(goodch)) |>
  ungroup()

suj_aper = aperiodic_data |> distinct(sujid) |> pull(sujid)

# --- Lagged alpha LAcH data ---
# alpha_LAcH: Hilbert-based Lagged alpha LAcH measuring oscillatory temporal stability.
hlcoh_data = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long.csv")) |>
  # Merge electrode region info (from mapping above) by channel label
  merge(electrodes |> dplyr::select(label, chinclu, region), by.x = "ch", by.y = "label") |>
  filter(chinclu == 1, epochs >= epochs_threshold, !sujid %in% EXCLUDED_SUBJECTS) |>
  # Harmonize session_age labels: 15 to 18 mo., 40 to 42 mo.
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  ) |>
  group_by(sujid, session_age, region) |>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE), .groups = "drop")

suj_lcoh = hlcoh_data |> distinct(sujid) |> pull(sujid)

# --- Burst properties data ---
# is_burst column identifies Burst vs. NoBurst cycle rows.
# Corrected amplitude ratios (Burst / NoBurst) remove non-specific power changes.
burst_data = read_csv(file.path(path2data, "BurstProperties_ByCycle_Long.csv")) |>
  # Merge full electrode metadata
  merge(electrodes, by.x = "ch", by.y = "label") |>
  filter(epochs >= epochs_threshold, chinclu == 1) |>
  # Harmonize session_age labels: 15 to 18 mo., 40 to 42 mo.
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  ) |>
  group_by(session_age, sujid, region, is_burst) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Compute corrected amplitude ratios in wide format: ratio of Burst to NoBurst to isolate burst-specific effects
burst_data_wide = burst_data |>
  dplyr::select(sujid, session_age, region, is_burst, volt_amp, band_amp) |>
  pivot_wider(names_from = is_burst, values_from = c(volt_amp, band_amp)) |>
  # Corrected amplitudes = Burst-condition / Non-Burst-condition (removes non-specific power changes)
  mutate(
    corrected_volt_amp = volt_amp_Burst / volt_amp_NoBurst,
    corrected_band_amp = band_amp_Burst / band_amp_NoBurst
  )

# Keep NoBurst rows as the amplitude baseline and attach corrected ratios via left join
burst_data = burst_data |>
  filter(is_burst == "NoBurst") |>
  dplyr::select(-is_burst) |>
  left_join(burst_data_wide)

suj_burst = burst_data |> distinct(sujid) |> pull(sujid)

# --- EEG preprocessing descriptives (defines the analysis sample) ---
# inclusion_lmm = 1 marks participants retained for the longitudinal models.
eeg_desc = read_csv(file.path(path2data, "CrossVisit_EEG_CleaningDescriptives.csv")) |>
  # Filter to visits with minimum clean epochs; inclusion_lmm flag for downstream LMM models
  filter(clean_epochs_rest >= epochs_threshold, !sujid %in% EXCLUDED_SUBJECTS) |>
  mutate(inclusion_lmm = 1) |>
  dplyr::select(sujid, session_age, prop_epochs, nsessions_withdata,
                inclusion_lmm, total_bch, total_ics, clean_epochs_rest) |>
  # Harmonize session_age labels: 15 to 18 mo., 40 to 42 mo.
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

suj_desc_eeg = eeg_desc |> distinct(sujid) |> pull(sujid)

# --- Sociodemographic / age data ---
desc_and_ages = read_csv(file.path(path2data, "Sociodemographic_Descriptives_Long_Updated.csv"))

# --- Sociodemographic variables ----------------------------------------------
# This release ships two files, and they cannot be joined to each other:
#
#   Sociodemographic_Descriptives_Long_Updated.csv   keyed by sujid
#       visit structure and the analysis covariates - session_age, age_months,
#       dev_filter, Cohort, GestationalAge_weeks, child_sex_dc.
#
#   Sociodemographic_SES_delinked.csv                keyed by demo_id
#       one row per participant: Cohort plus the six sociodemographic variables
#       (income-to-needs, family size, maternal education, race, ethnicity,
#       sex). Row order is shuffled and no sujid is present, so these values
#       cannot be attributed to an individual EEG record.
#
SES_VARS = c("ITN_mean", "FamilySize_mean", "mother_education_mean",
             "child_race", "child_ethnicity", "child_sex_dc")

ses_file = file.path(path2data, "Sociodemographic_SES_delinked.csv")
if (!file.exists(ses_file))
  stop("Missing ", basename(ses_file), " in:\n  ", path2data,
       "\nIt ships alongside the other CSVs in the data download.", call. = FALSE)

# One row per participant already; nothing to collapse.
descriptives_sample_SES = readr::read_csv(ses_file, show_col_types = FALSE)

missing_ses = setdiff(SES_VARS, names(descriptives_sample_SES))
if (length(missing_ses))
  stop(basename(ses_file), " is missing: ", paste(missing_ses, collapse = ", "),
       call. = FALSE)
if (!"demo_id" %in% names(descriptives_sample_SES))
  stop(basename(ses_file), " must have a demo_id column.", call. = FALSE)

message(sprintf("[demographics] %d participants from %s (de-linked; not joinable to EEG).",
                nrow(descriptives_sample_SES), basename(ses_file)))

cat(paste0("Participants with developmental disorders excluded: ",
           length(unique(desc_and_ages |> filter(dev_filter == 0) |> pull(sujid))), "\n"))

desc_and_ages = desc_and_ages |>
  # Inclusion criteria: dev_filter==1 (no developmental disorders), exclude pilot children (BUDDY001/006) and other children with no usable EEG data,
  filter(dev_filter == 1, !sujid %in% c("BUDDY001", "BUDDY006"), sujid %in% suj_desc_eeg) |>
  # Harmonize session_age labels; impute missing Cohort as 2 (Cohort 3 in final factor labels)
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age),
    Cohort      = if_else(is.na(Cohort), 2, Cohort) # Participants w missing Cohort belonged to cohort 2
  ) |>
  # Mean-impute missing SES covariates so that subjects with one missing covariate
  # are not excluded from the balance model.
  mutate(across(any_of(c("ITN_mean", "FamilySize_mean", "mother_education_mean",
                         "GestationalAge_weeks")),
                ~ if_else(is.na(.x), mean(.x, na.rm = TRUE), .x)))


# Recode raw numeric codes to labeled factor levels for demographic presentation in Tables
# child_sex_dc: 0=Male, 1=Female
# child_ethnicity: 1=Hispanic/Latino, 2=Non-Hispanic/Latino
# child_race: 1=White, 2=Black/AA, 3=Asian, 4=AI/AN, 5=NH/PI, 6=Other, 7=Not Answer
# Cohort: 0=Cohort1, 1=Cohort2, 2=Cohort3
nn_recode = function(d, var, levels, labels) {
  if (!var %in% names(d)) return(d)
  if (!is.numeric(d[[var]])) return(d)          # already labelled
  d[[var]] = factor(d[[var]], levels = levels, labels = labels)
  d
}
RACE_LEVELS = c(1,2,3,4,5,6,7)
RACE_LABELS = c("White", "Black/African American", "Asian",
                "American Indian/Alaska Native",
                "Native Hawaiian/Pacific Islander", "Other", "Not Answer")

# Analysis tier: sex (kept for Table S5 by visit) and Cohort.
desc_and_ages = desc_and_ages |>
  nn_recode("child_sex_dc", c(0,1), c("Male", "Female")) |>
  nn_recode("Cohort",       c(0,1,2), c("Cohort 1", "Cohort 2", "Cohort 3"))

# De-linked tier: all three categorical sociodemographic variables, plus Cohort.
descriptives_sample_SES = descriptives_sample_SES |>
  nn_recode("child_sex_dc",    c(0,1), c("Male", "Female")) |>
  nn_recode("child_ethnicity", c(1,2), c("Hispanic/Latino", "Non-Hispanic/Latino")) |>
  nn_recode("child_race",      RACE_LEVELS, RACE_LABELS) |>
  nn_recode("Cohort",          c(0,1,2), c("Cohort 1", "Cohort 2", "Cohort 3"))

# Subject-level wide SES data (one row per participant) for table display
desc_and_ages_wide = desc_and_ages |>
  dplyr::select(sujid, contains("mean"), GestationalAge_weeks) |>
  distinct() |>
  filter(sujid %in% eeg_desc$sujid) |>
  mutate(across(any_of("ITN_mean"), ~ if_else(is.na(.x), mean(.x, na.rm = TRUE), .x)))

# Combined descriptives frame: EEG preprocessing + sociodemographic data for descriptive Table
# Left-anchored on eeg_desc (preprocessing metrics and inclusion flags)
descriptives = left_join(
  eeg_desc,
  desc_and_ages |>
    filter(sujid %in% eeg_desc$sujid) |>
    dplyr::select(any_of(c("sujid", "age_months", "session_age", "dev_filter",
                           "Cohort", "child_race", "child_sex_dc",
                           "child_ethnicity")))
) |>
  left_join(desc_and_ages_wide) |>
  # Final QC: drop rows missing preprocessing info, keep only non-excluded subjects (dev_filter==1)
  filter(!is.na(prop_epochs), dev_filter == 1) |>
  # Coerce Cohort factor to character for final label harmonization
  mutate(Cohort = if_else(is.na(Cohort), "Cohort 3", as.character(Cohort)))

included_participants = unique(descriptives$sujid)


# =============================================================================
# SECTION 4: SOCIODEMOGRAPHIC SUMMARY Tables
# =============================================================================

reset_gtsummary_theme()

# --- Table 1: Whole-sample sociodemographics ---
ses_labels = list(
  ITN_mean              ~ "Income to Needs Ratio",
  FamilySize_mean       ~ "Family Size (# people in household)",
  mother_education_mean ~ "Maternal Education (years)",
  child_race            ~ "Child Race",
  child_ethnicity       ~ "Child Ethnicity",
  child_sex_dc          ~ "Child Sex"
)

descriptives_sample_SES |>
  dplyr::select(-any_of(c("demo_id", "Cohort"))) |>
  tbl_summary(
    statistic = list(all_continuous()   ~ "{mean} ({sd})\n[{p2.5}, {p97.5}]",
                     all_categorical()  ~ "{n} ({p}%)"),
    label = ses_labels,
    missing = "ifany", 
    digits  = all_continuous() ~ 2
  ) |>
  modify_header(label ~ "**Socio-Demographic Variable**") |>
  bold_labels() |>
  modify_caption(paste0("**Table 2a.** Sociodemographic characteristics of ",
                        "participants included in the study. Gestational age is ",
                        "reported in Table 2b.")) |>
  as_flex_table() |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(j = 2, align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  flextable::save_as_docx(path = file.path(path2desc, "Table_2a_Sociodemographics_Characteristics.docx"))

# --- Table 2b: gestational age ------------------------------------------------
ga_person = desc_and_ages |>
  dplyr::distinct(sujid, GestationalAge_weeks) |>
  dplyr::filter(!is.na(GestationalAge_weeks))

ga_summary = tibble::tibble(
  Variable = "Gestational Age (weeks)",
  N        = nrow(ga_person),
  Value    = sprintf("%.2f (%.2f)\n[%.2f, %.2f]",
                     mean(ga_person$GestationalAge_weeks),
                     sd(ga_person$GestationalAge_weeks),
                     quantile(ga_person$GestationalAge_weeks, 0.025),
                     quantile(ga_person$GestationalAge_weeks, 0.975))
)

ga_summary |>
  flextable() |>
  set_header_labels(Variable = "Socio-Demographic Variable",
                    N        = "N",
                    Value    = paste0("Overall\nN = ", nrow(ga_person))) |>
  set_caption(paste0("Table 2b. Gestational age of participants included in the study",
                     "Values are mean (s.d.) [2.5th, 97.5th percentile].")) |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(j = 2:3, align = "center", part = "all") |>
  bold(part = "header") |>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  flextable::save_as_docx(path = file.path(path2desc, "Table_2b_GestationalAge.docx"))

message(sprintf(
  "[Table 2b] Gestational age (weeks): M = %.2f (SD = %.2f) [%.2f, %.2f], n = %d",
  mean(ga_person$GestationalAge_weeks), sd(ga_person$GestationalAge_weeks),
  quantile(ga_person$GestationalAge_weeks, 0.025),
  quantile(ga_person$GestationalAge_weeks, 0.975), nrow(ga_person)))

# --- Table S1a: Sociodemographics by cohort (balance check) ---
descriptives_sample_SES |>
  dplyr::select(-any_of(c("demo_id", "sujid", "session_age", "age_months",
                          "dev_filter"))) |>
  tbl_summary(
    by = Cohort,
    statistic = list(all_continuous()  ~ "{mean} ({sd})\n[{p2.5}, {p97.5}]",
                     all_categorical() ~ "{n} ({p}%)"),
    label   = ses_labels,
    missing = "ifany", 
    digits  = all_continuous() ~ 2
  ) |>
  add_p() |>
  add_overall() |>
  add_significance_stars() |>
  bold_p(t = 0.05) |>
  modify_header(label ~ "**Socio-Demographic Variable**") |>
  bold_labels() |>
  modify_caption(paste0("**Table S1a.** Sociodemographic characteristics by ",
                        "cohort. Gestational age is reported in Table S1b.")) |>
  modify_footnote(all_stat_cols() ~ "ANOVA for continuous variables; chi-square for categorical variables.") |>
  as_flex_table() |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(j = 2:6, align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  flextable::save_as_docx(path = file.path(path2ed_tab, "Extended_Data_Table_1a_Sociodemographics_Characteristics_by_Cohort.docx"))


# --- Table S1b: gestational age by cohort -------------------------------------
ga_long = desc_and_ages |>
  dplyr::distinct(sujid, Cohort, GestationalAge_weeks) |>
  dplyr::filter(!is.na(GestationalAge_weeks)) |>
  dplyr::group_by(Cohort) |>
  dplyr::summarise(
    N     = dplyr::n(),
    Value = sprintf("%.2f (%.2f)\n[%.2f, %.2f]",
                    mean(GestationalAge_weeks), sd(GestationalAge_weeks),
                    quantile(GestationalAge_weeks, 0.025),
                    quantile(GestationalAge_weeks, 0.975)),
    .groups = "drop")

# One column per cohort with N in the header, matching Table S1a's layout so the
# row appends cleanly rather than arriving with N_Cohort/Value_Cohort columns.
ga_by_cohort = ga_long |>
  dplyr::select(Cohort, Value) |>
  tidyr::pivot_wider(names_from = Cohort, values_from = Value) |>
  dplyr::mutate(Variable = "Gestational Age (weeks)", .before = 1)

ga_cohort_headers = stats::setNames(
  as.list(sprintf("%s\nN = %d", ga_long$Cohort, ga_long$N)),
  as.character(ga_long$Cohort))

ga_by_cohort |>
  flextable() |>
  set_header_labels(values = c(list(Variable = "Socio-Demographic Variable"),
                               ga_cohort_headers)) |>
  set_caption(paste0("Table S1b. Gestational age by cohort. It has to be merged with Extended_Table_1",
                     "[2.5th, 97.5th percentile].")) |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  bold(part = "header") |>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  flextable::save_as_docx(path = file.path(path2ed_tab,
                                           "Extended_Data_Table_1b_GestationalAge_by_Cohort.docx"))

message("[Table S1b] gestational age by cohort written.")


# =============================================================================
# SECTION 5: COHORT BALANCE MODEL (MULTINOMIAL LOGISTIC REGRESSION)
# =============================================================================
# Tests whether sociodemographic variables differ between cohorts (balance check).

ses_variables = descriptives_sample_SES |> drop_na()

# Full model: Cohort ~ all SES predictors.
# Six predictors. The published model additionally included gestational age,
# which is not in the de-linked file; its coefficients therefore differ slightly
# from Table S1 in the paper. The published model is the seven-predictor one. This 
# does not affect the results: Public x2(20) = 87.53, p < .0001, Manuscript χ2(22)= 93.720, P < 0.0001
balance_model = VGAM::vglm(
  Cohort ~ ITN_mean + FamilySize_mean + mother_education_mean +
    child_race + child_ethnicity + child_sex_dc,
  data   = ses_variables,
  family = "multinomial"
)
summary(balance_model)

# Likelihood Ratio Test against null model
null_model = VGAM::vglm(Cohort ~ 1, family = "multinomial", data = ses_variables)
lrt_anova  = VGAM::lrtest(null_model, balance_model)

rtf_file = rtf::RTF(file.path(path2supptab, "Complementary_Analysis_CohortDifferencies.doc"))
rtf::addParagraph(rtf_file, "Likelihood Ratio Test for the Cohort Balance Model", bold = TRUE, underline = TRUE)
rtf::addTable(rtf_file, lrt_anova@Body)
rtf::done(rtf_file)

# =============================================================================
# SECTION 6: EEG PREPROCESSING DESCRIPTIVES
# =============================================================================

# --- Table S5: EEG preprocessing descriptives by visit ---
descriptives |>
  dplyr::select(session_age, child_sex_dc, age_months, total_bch, total_ics,
                prop_epochs, clean_epochs_rest) |>
  mutate(session_age = factor(paste(session_age, 'mo.'), levels = c('1 mo.', '6 mo.', '12 mo.', '18 mo.', '30 mo.', '36 mo.', '42 mo.', '48 mo.'))) |>
  tbl_summary(
    by = session_age,
    statistic = list(all_continuous()  ~ "{mean} ({sd})\n[{p2.5}, {p97.5}]",
                     all_categorical() ~ "{n} ({p}%)"),
    label = list(
      age_months        ~ "Age (months)",
      total_bch         ~ "Bad Channels",
      total_ics         ~ "ICA Components Removed",
      clean_epochs_rest ~ "Clean Epochs (N)",
      prop_epochs        ~ "Clean Epochs (%)",
      child_sex_dc      ~ "Child Sex"
    ),
    missing = "no",
    digits  = all_continuous() ~ 3
  ) |>
  modify_header(label ~ "") |>
  bold_labels() |>
  modify_caption(paste0("**EEG preprocessing descriptives by visit (N = ",
                        length(unique(descriptives$sujid)), ").**")) |>
  modify_footnote(all_stat_cols() ~ "Mean (SD) [2.5th, 97.5th centile] for continuous; N (%) for categorical.") |>
  as_flex_table() |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(j = 2:8, align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  flextable::save_as_docx(path = file.path(path2supptab, "Table_S3_EEG_Preprocessing_Metrics.docx"))

# --- Complementary Table: Distribution of number of visits with usable data ---
descriptives |>
  dplyr::select(sujid, nsessions_withdata) |>
  mutate(nsessions_withdata = if_else(nsessions_withdata == 0, 1, nsessions_withdata)) |>
  distinct() |>
  group_by(nsessions_withdata) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(percent = round(n / sum(n) * 100, 1)) |>
  flextable() |>
  set_header_labels(nsessions_withdata = "Visits with Usable Data", n = "N Participants",
                    percent = "%") |>
  set_caption("**Distribution of Participants by Number of Usable Visits**") |>
  theme_booktabs() |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(j = 2, align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  save_as_html(path = file.path(path2supptab, "Complementary_Table_Sample_by_NumberOfVisitsProvided.html"))



# =============================================================================
# SECTION 7: SAMPLE PLOT (LONGITUDINAL COVERAGE BY COHORT)
# =============================================================================

descriptives_plot_long = descriptives |> dplyr::select(age_months, sujid, Cohort)

plot_sample = ggplot(
  descriptives_plot_long,
  aes(x     = age_months,
      y     = sujid,
      color = factor(Cohort, levels = c("Cohort 1", "Cohort 2", "Cohort 3")))
) +
  facet_grid(Cohort ~ ., scales = "free_y") +
  geom_line(aes(group = sujid), alpha = 0.75, linewidth = 0.75) +
  geom_point(size = 2) +
  scale_color_manual(values = COLORS_COHORT) +  # From 00_Setup_PackageInstallation.R
  labs(x = "Age (months)", y = "") +
  theme(
    axis.text.y    = element_blank(),
    axis.ticks.y   = element_blank(),
    legend.title   = element_blank(),
    legend.text    = element_text(size = 12),
    legend.position = "top",
    panel.background = element_blank()
  )

# SES distribution plots per metric (categorical and continuous)
# Configuration defines which metrics to plot and their display labels
metrics_config = list(
  list(metric = "child_race",            labels = "Child Race"),
  list(metric = "child_sex_dc",          labels = "Child Sex"),
  list(metric = "ITN_mean",              labels = "Income to Needs Ratio"),
  list(metric = "mother_education_mean", labels = "Maternal Education (years)")
)

plot_ses_list = list()
for (item in metrics_config) {
  v     = item$metric
  label = item$labels
  
  # These four panels are participant-level sociodemographic distributions, so
  # they come from descriptives_sample_SES. `descriptives` is the EEG frame and
  # holds only child_sex_dc of the four: race and the SES means are in the
  # de-linked file, which cannot be joined to it.
  plot_data = descriptives_sample_SES |>
    dplyr::select(any_of(c("demo_id", "sujid")), all_of(v), Cohort) |>
    distinct() |>
    filter(!is.na(!!sym(v)))
  
  if (v == "child_race") {
    # Abbreviate race labels for plot x-axis to avoid crowding (Table_S retain full labels)
    plot_data <- plot_data|>
      mutate(child_race = factor(child_race, levels= c("White", "Black/African American", "Asian",
                                                       "American Indian/Alaska Native",
                                                       "Native Hawaiian/Pacific Islander", "Other", "Not Answer"),
                                 labels = c("White", "Black/AA", "Asian",
                                            "AI/AN", "NH/PI", "Other", "NA")))}
  
  if (v %in% c("child_race", "child_sex_dc")) {
    # Categorical variable: bar chart by cohort
    p_ses = ggplot(plot_data, aes(factor(!!sym(v)), fill = factor(Cohort))) +
      geom_bar(color = "black", alpha = 0.5, linewidth = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(values = COLORS_COHORT) +
      THEME_BASE + THEME_TEXT_COMPACT +
      labs(x = label, y = "Count", fill = "Cohort")
  } else {
    # Continuous variable: histogram by cohort (bins=30)
    p_ses = ggplot(plot_data, aes(x = !!sym(v), fill = factor(Cohort))) +
      geom_histogram(color = "black", alpha = 0.5, linewidth = 0.5, bins = 30) +
      scale_fill_manual(values = COLORS_COHORT) +
      THEME_BASE + THEME_TEXT +
      labs(x = label, y = "Count", fill = "Cohort")
  }
  plot_ses_list[[v]] = p_ses
}

plot_topright = ggpubr::ggarrange(
  plotlist     = plot_ses_list[1:4],
  ncol         = 2, nrow = 2,
  common.legend = TRUE, legend = "bottom"
)
plot_final = ggpubr::ggarrange(plot_sample, plot_topright, ncol = 2, widths = c(0.4, 0.6))
ggsave(file.path(path2suppfig, "Extended_Data_Fig_1_LongitudinalCohort_Characteristics.jpeg"), plot_final, width = 10, height = 7)


# =============================================================================
# SECTION 8: MISSING DATA SUMMARIES
# =============================================================================

# Variable label maps for display in Table_S: maps code names to manuscript display labels
# Sociodemographic variable labels
socio_var_labels = c(
  ITN_mean              = "Income to Needs Ratio",
  FamilySize_mean       = "Family Size",
  mother_education_mean = "Maternal Education (Years)",
  child_race            = "Child Race",
  child_ethnicity       = "Child Ethnicity",
  child_sex_dc          = "Child Sex",
  GestationalAge_weeks  = "Gestational Age (Weeks)"
)

# EEG metric variable labels (used in recode() for pivot operations)
# Maps internal column names to display labels for descriptive Table_S
eeg_var_labels = c(
  slope                   = "Slope",
  offset                  = "Offset",
  peak_prop                = "Prop. of Peaks",
  peak_ampl               = "Peak Amp.",
  osc_ampl                = "Band Power",
  peak_freq               = "Peak Freq.",
  alpha_LAcH              = "Lifespan",
  prop_bursty_epochs     = "Prop. of Epochs w/ Bursts",
  prop_bursty_cycles_burst = "Prop. of Cycles w/ Bursts",
  avg_burst_duration      = "Burst Duration",
  volt_amp_Burst          = "Volt. Amp. (Burst Cycles)",
  volt_amp_NoBurst        = "Volt. Amp. (Non-Burst Cycles)",
  band_amp_Burst          = "Band Amp. (Burst Cycles)",
  band_amp_NoBurst        = "Band Amp. (Non-Burst Cycles)"
)


# --- Complementary Table: Missing sociodemographic variables ---
# Load raw sociodemographic data (before imputation) to show true missingness patterns
# The seven variables are split across the two files, so missingness is counted
# in each and the results stacked. Counting only the analysis tier would silently
# report on two variables instead of seven.
raw_socio_analysis = read_csv(file.path(path2data, "Sociodemographic_Descriptives_Long_Updated.csv")) |>
  filter(sujid %in% suj_desc_eeg, dev_filter == 1) |>
  # Keep only one row per subject (sociodemographics are stable within-subject)
  distinct(sujid, .keep_all = TRUE) |>
  dplyr::select(any_of(setdiff(names(socio_var_labels), names(descriptives_sample_SES))))

raw_socio_ses = descriptives_sample_SES |>
  dplyr::select(any_of(names(socio_var_labels)))

raw_socio = dplyr::bind_cols(
  raw_socio_ses,
  if (ncol(raw_socio_analysis) && nrow(raw_socio_analysis) == nrow(raw_socio_ses))
    raw_socio_analysis else NULL
)
if (nrow(raw_socio_analysis) != nrow(raw_socio_ses))
  message("[missingness] the two sources have different participant counts (",
          nrow(raw_socio_ses), " vs ", nrow(raw_socio_analysis),
          "); reporting the de-linked variables only.")

raw_socio |>
  # Count missingness per variable
  summarise(across(everything(), ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "Variable", values_to = "Missing_N") |>
  # Compute percentages and remap variable names using socio_var_labels
  mutate(
    Total        = nrow(raw_socio),
    Missing_Pct  = round(Missing_N / Total * 100, 1),
    Variable     = recode(Variable, !!!socio_var_labels)
  ) |>
  flextable() |>
  set_header_labels(Variable = "Metric", Missing_N = "N Missing",
                    Total = "N Total", Missing_Pct = "% Missing") |>
  set_caption("Missing Sociodemographic Data") |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(j = 2:4, align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  save_as_html(path = file.path(path2supptab, "Complementary_Table_MissinSociodemographic_Information.html"))


# =============================================================================
# SECTION 9: EEG METRIC DESCRIPTIVES BY VISIT
# =============================================================================

# Subset all EEG datasets to only the final analysis sample (included_participants)
# ensures consistency across all descriptive summaries
aperiodic_data = filter(aperiodic_data, sujid %in% included_participants)
hlcoh_data     = filter(hlcoh_data,     sujid %in% included_participants)
burst_data     = filter(burst_data,     sujid %in% included_participants)

# Sanity check: we have to have 207 participants in all-epochs metrics even though some datasets got discarded due to processing problems or data quality.
cat(sprintf("Participants in aperiodic data:  %d\n", n_distinct(aperiodic_data$sujid)))
cat(sprintf("Participants in alpha LAcH data:  %d\n", n_distinct(hlcoh_data$sujid)))
cat(sprintf("Participants in burst data:      %d\n", n_distinct(burst_data$sujid)))

# Create a common subject × session grid (left anchor for merging)
base_grid = eeg_desc |> filter(sujid %in% included_participants) |>
  dplyr::select(sujid, session_age)

# Summarise each dataset to one row per subject × visit
aper_summ = aperiodic_data |>
  group_by(sujid, session_age) |>
  summarise(across(c(slope, offset, peak_prop, peak_ampl, osc_ampl, peak_freq),
                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")

coh_summ = hlcoh_data |>
  group_by(sujid, session_age) |>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE), .groups = "drop")

burst_summ = burst_data |>
  group_by(sujid, session_age) |>
  summarise(across(c(prop_bursty_epochs, prop_bursty_cycles_burst, avg_burst_duration,
                     volt_amp_Burst, volt_amp_NoBurst, band_amp_Burst, band_amp_NoBurst),
                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Define factor levels for consistent ordering in output Tables
# visit_levels: canonical visit ages (months) in longitudinal sequence
visit_levels = paste0(c(1, 6, 12, 18, 30, 36, 42, 48), " mo.")
# var_levels: variable display order (power spectrum, burst properties, then lifespan)
var_levels  = c("Slope", "Offset", "Prop. of Peaks", "Peak Amp.", "Band Power",
                "Peak Freq.", "Prop. of Epochs w/ Bursts", "Prop. of Cycles w/ Bursts", "Burst Duration",
                "Volt. Amp. (Burst Cycles)", "Volt. Amp. (Non-Burst Cycles)",
                "Band Amp. (Burst Cycles)", "Band Amp. (Non-Burst Cycles)", "Lifespan")

# Combine all EEG metric summaries into long format with consistent variable labels
# Left-anchored on base_grid (subject×visit skeleton)
eeg_long = base_grid |>
  merge(coh_summ,   by = c("sujid", "session_age"), all.x = TRUE) |>
  merge(burst_summ, by = c("sujid", "session_age"), all.x = TRUE) |>
  merge(aper_summ,  by = c("sujid", "session_age"), all.x = TRUE) |>
  pivot_longer(cols = -c(sujid, session_age), names_to = "Variable", values_to = "Value") |>
  # Remap variable names to display labels; append " mo." to session_age for display
  mutate(
    Variable    = recode(Variable, !!!eeg_var_labels),
    session_age = paste0(session_age, " mo.")
  ) |>
  # Coerce to ordered factors for consistent table/plot ordering
  mutate(session_age = factor(session_age, levels = visit_levels),
         Variable    = factor(Variable,    levels = var_levels))

# --- Table S1: EEG descriptives by visit (Mean (SD)) ---
eeg_long |>
  group_by(session_age, Variable) |>
  summarise(
    Mean     = mean(Value, na.rm = TRUE),
    SD       = sd(Value,   na.rm = TRUE),
    n_visit  = n(),
    .groups  = "drop"
  ) |>
  mutate(
    Display_Stat = sprintf("%.3f\n(%.3f)", Mean, SD),
    visit_label  = paste0(session_age)
  ) |>
  dplyr::select(Variable, visit_label, Display_Stat) |>
  pivot_wider(names_from = visit_label, values_from = Display_Stat) |>
  rename(Metric = Variable)|>
  flextable() |>
  set_caption("**EEG Metric Descriptives by Visit (Mean (SD))**") |>
  add_footer_lines("Values are Mean (SD) across participants at each visit.") |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(j = 2:9, align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  save_as_docx(path = file.path(path2si_gi_tab, "Table_S1_WholeBrain_Descriptives_by_visit.docx"))

# --- Table: EEG descriptives by brain region (Mean (SD)) ---
aper_summ_reg  = aperiodic_data |>
  group_by(sujid, session_age, region) |>
  summarise(across(c(slope, offset, peak_prop, peak_ampl, osc_ampl, peak_freq),
                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")

coh_summ_reg   = hlcoh_data |>
  group_by(sujid, session_age, region) |>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE), .groups = "drop")

burst_summ_reg = burst_data |>
  group_by(sujid, session_age, region) |>
  summarise(across(c(prop_bursty_epochs, prop_bursty_cycles_burst, avg_burst_duration,
                     volt_amp_Burst, volt_amp_NoBurst, band_amp_Burst, band_amp_NoBurst),
                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Combine all EEG metric summaries by brain region into long format
# Full outer joins to retain all region observations
eeg_long_reg = aper_summ_reg |>
  full_join(coh_summ_reg,   by = c("sujid", "session_age", "region")) |>
  full_join(burst_summ_reg, by = c("sujid", "session_age", "region")) |>
  filter(!is.na(region)) |>
  pivot_longer(cols = -c(sujid, session_age, region), names_to = "Variable", values_to = "Value") |>
  # Remap variable names to display labels; apply ordered factor for consistent presentation
  mutate(Variable = recode(Variable, !!!eeg_var_labels),
         Variable = factor(Variable, levels = var_levels))


# Loop over metric classes (power spectrum, burst properties, lifespan) to generate separate Table_S S3XX
for (t in c('powerspectrum', 'burst', 'lifespan')) {
  
  # Variable list for this metric class (subset of var_levels)
  vois = switch(t,
                powerspectrum = c("Slope", "Offset", "Prop. of Peaks", "Peak Freq.", " Peak Amp.", "Band Power"),
                burst         = c("Prop. of Epochs w/ Bursts", "Prop. of Cycles w/ Bursts", "Burst Duration",
                                  "Volt. Amp. (Burst Cycles)", "Volt. Amp. (Non-Burst Cycles)",
                                  "Band Amp. (Burst Cycles)", "Band Amp. (Non-Burst Cycles)"),
                lifespan       = c("Lifespan"))
  eeg_long_reg |>
    filter(Variable %in% vois)|>
    group_by(region, Variable) |>
    summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop") |>
    mutate(Display_Stat = sprintf("%.3f\n(%.3f)", Mean, SD)) |>
    dplyr::select(Variable, region, Display_Stat) |>
    pivot_wider(names_from = region, values_from = Display_Stat) |>
    flextable() |>
    set_caption("**EEG Metric Descriptives by Brain Region (Mean (SD))**") |>
    add_footer_lines("Values are Mean (SD) across all participants and visits.") |>
    theme_booktabs() |> fix_border_issues() |> autofit() |>
    set_table_properties(layout = "fixed") |>
    align(align = "left", part = "all") |>
    align(j = 2:6, align = "center", part = "all") |>
    bold(part = "header")|>
    padding(padding = 1, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
    set_header_labels(Variable = "Metric") |>
    save_as_docx(path = file.path(path2si_gi_tab, paste0("Complementary_EEG_ROIs_Descriptives_", t, ".docx")))
  
  eeg_long_reg |>
    filter(Variable %in% vois)|>
    mutate(session_age = paste0(session_age, " mo."),
           session_age = factor(session_age, levels = visit_levels),
           Variable    = factor(Variable,    levels = var_levels)) |>
    group_by(region, session_age, Variable) |>
    summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop") |>
    mutate(Display_Stat = sprintf("%.3f\n(%.3f)", Mean, SD)) |>
    dplyr::select(Variable, region, session_age, Display_Stat) |>
    pivot_wider(names_from = session_age, values_from = Display_Stat) |>
    rename(ROI = region) |>
    rename(Metric = Variable) |>
    arrange(Metric, ROI) |>
    flextable() |>
    merge_v(j = "Metric") |>
    valign(j = "Metric", valign = "top") |>
    set_caption("**EEG Metric Descriptives by Region and Visit (Mean (SD))**") |>
    add_footer_lines("Values are Mean (SD). Rows grouped by metric; columns are visit ages.") |>
    theme_booktabs() |> fix_border_issues() |> autofit() |>
    set_table_properties(layout = "fixed") |>
    align(align = "left", part = "all") |>
    align(j = 2:9, align = "center", part = "all") |>
    bold(part = "header")|>
    padding(padding = 1, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
    set_header_labels(Variable = "Metric") |>
    save_as_docx(path = file.path(path2si_gi_tab, paste0("Table_S2_ROIs_Descriptives_by_visit_", t, ".docx")))
}

# =============================================================================
# SECTION 10: BURST VS. NON-BURST DESCRIPTIVES (MODEL FIT AND METRICS)
# =============================================================================
# Reloads the burst-conditioned aperiodic and alpha LAcH files.
# These contain separate rows for burst and non-burst cycle estimates.

# Variable list for burst-conditioned aperiodic data (includes burst label for stratification)
aper_voi_burst = c("sujid", "session_age", "ch", "region", "chinclu", "r2value", "goodch",
                   "offset", "slope", "alpha_freq", "alpha_ampl", "alpha_osc",
                   "alpha_peak", "inclusion_final_dummy", "epochs", "burst", "mae")

aperiodic_data_burst = read_csv(file.path(path2data, "Aperiodic_Oscillatory_ByCycle_Long_Burst.csv")) |>
  dplyr::select(all_of(aper_voi_burst)) |>
  # Rename epochs to epochs for consistency with other datasets
  rename(epochs = epochs) |>
  # Apply same quality filters as aperiodic_data
  filter(
    goodch >= CH_THRESHOLD,
    r2value > r2_thresh,
    chinclu == 1,
    epochs >= epochs_threshold,
    mae < mae_thresh,
    !sujid %in% EXCLUDED_SUBJECTS,
    sujid %in% included_participants
  ) |>
  # Harmonize session_age labels: 15 to 18 mo., 40 to 42 mo.
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  ) |>
  # Aggregate electrode-level to subject×visit×region×burst-state level
  group_by(session_age, sujid, burst, region) |>
  summarise(
    slope     = mean(slope,                       na.rm = TRUE),
    offset    = mean(offset,                      na.rm = TRUE),
    peak_ampl = mean(alpha_ampl,                  na.rm = TRUE),
    # peak_freq computed only from cycles with alpha peak present
    peak_freq = mean(alpha_freq[alpha_peak == 1], na.rm = TRUE),
    osc_ampl  = mean(alpha_osc,                   na.rm = TRUE),
    peak_prop  = mean(alpha_peak,                  na.rm = TRUE),
    epochs    = mean(epochs,                      na.rm = TRUE),
    r2value   = mean(r2value,                     na.rm = TRUE),
    # Rename mae to mae for clarity (Mean Absolute Error from Specparam)
    mae       = mean(mae,                         na.rm = TRUE),
    goodch    = n(),
    .groups   = "drop"
  ) |>
  rename(region = region) |>
  # Compute total good channels per subject×visit×burst-state combination
  mutate(goodch_sum = sum(goodch), .by = c(session_age, sujid, burst))

hlcoh_burst = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long_Burst.csv"))

hlcoh_burst = hlcoh_burst|>
  # Merge electrode inclusion info
  merge(electrodes |> dplyr::select(label, chinclu), by.x = "ch", by.y = "label") |>
  filter(chinclu == 1, epochs >= epochs_threshold) |>
  # Count good channels per subject×visit×burst and apply inclusion threshold (>75% of 60 channels)
  mutate(
    goodch = n(), .by = c(session_age, sujid, burst),
    inclusion_final_dummy = as.integer(goodch > CH_THRESHOLD)
  ) |>
  filter(inclusion_final_dummy == 1, sujid %in% included_participants) |>
  # Harmonize session_age labels: 15 to 18 mo., 40 to 42 mo.
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  ) 


# --- Table S8 ---

hlcoh_pvals = hlcoh_burst|>
  # Aggregate to subject×visit×burst level for paired tests
  group_by(session_age, sujid, burst) |>
  summarise(
    epochs     = mean(epochs,     na.rm = TRUE),
    goodch     = n(),
    .groups    = "drop"
  ) |>
  mutate(sesion_age = paste0(session_age, " mo."))|>
  # Reshape to wide format (separate columns for Burst vs. NoBurst)
  pivot_wider(names_from = burst, values_from = c(epochs, goodch))|>
  # Paired t-tests per visit: Burst vs. NoBurst for epochs and channel count
  group_by(session_age) |>
  summarise(
    p_epochs     = t.test(epochs_burst,     epochs_noburst,     paired = TRUE)$p.value,
    p_goodch     = t.test(goodch_burst,     goodch_noburst,     paired = TRUE)$p.value,
    .groups      = "drop"
  ) |>
  # Reshape to long format and recode variable names
  pivot_longer(starts_with("p_"), names_to = "Variable", values_to = "p_value") |>
  mutate(
    Variable = case_when(
      Variable == "p_epochs"     ~ "# of Trials",
      Variable == "p_goodch"     ~ "# Electrodes"
    )
  ) |>
  # Apply FDR correction within each variable, then format p-values
  mutate(p_value = p.adjust(p_value, method = "fdr"), .by = Variable) |>
  mutate(p_value = style_pvalue(p_value))

hlcoh_table = hlcoh_burst|>
  # Apply final inclusion filters (re-applied for clarity in context of burst file)
  filter(chinclu == 1, epochs >= epochs_threshold) |>
  mutate(
    goodch = n(), .by = c(session_age, sujid, burst),
    inclusion_final_dummy = as.integer(goodch > CH_THRESHOLD)
  ) |>
  filter(inclusion_final_dummy == 1, sujid %in% included_participants) |>
  # Harmonize session_age labels: 15 to 18 mo., 40 to 42 mo.
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  ) |>
  # Aggregate to subject×visit×burst level
  group_by(session_age, sujid, burst) |>
  summarise(
    epochs     = mean(epochs,     na.rm = TRUE),
    goodch     = n(),
    .groups    = "drop"
  )|>
  # Reshape to long format and recode variable names
  pivot_longer(c(epochs, goodch), names_to = "Variable", values_to = "Value") |>
  mutate(Variable = case_when(
    Variable == "epochs"  ~ "# of Trials",
    Variable == "goodch"  ~ "# Electrodes"
  )) |>
  # Compute descriptive statistics (Mean, SD) per metric×visit×burst
  group_by(Variable, session_age, burst) |>
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop") |>
  mutate(Display_Stat = sprintf("%.3f (%.3f)", Mean, SD)) |>
  dplyr::select(Variable, session_age, burst, Display_Stat) |>
  # Reshape to wide format (separate columns for Burst vs. NoBurst)
  pivot_wider(names_from = burst, values_from = Display_Stat) |>
  # Attach p-values from paired t-tests
  left_join(hlcoh_pvals, by = c("Variable", "session_age")) |>
  # Format visit ages and apply ordered factor
  mutate(session_age = paste0(session_age, " mo.")) |>
  mutate(session_age = factor(session_age, levels = visit_levels)) |>
  arrange(Variable, session_age) |>
  rename(Metric = Variable) |>
  flextable() |>
  merge_v(j = "Metric") |> valign(j = "Metric", valign = "top") |>
  set_header_labels(Variable = "Metric", session_age = "Visit",
                    burst = "Burst Mean\n(SD)", noburst = "Non-Burst\nMean (SD)",
                    p_value = "p (FDR)") |>
  add_footer_lines("Paired t-tests per visit; FDR correction (Benjamini-Hochberg) within each metric") |>
  set_caption("# Trials and Ch in burst- and non-burst containing sets") |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(j = 2:4, align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  save_as_docx(path = file.path(path2supptab, "Table_S6_BurstvsNonBurst_Datasets_Lifespan_PreprocessingComparison.docx"))


# --- Table S5: Parametrization Model Fit by Burst State and Visit ---
# Compares Specparam R² and MAE (model fit), epoch counts, and channel availability between Burst and Non-Burst conditions
aper_model_fit = aperiodic_data_burst |>
  # Aggregate to subject×visit×burst level
  group_by(sujid, session_age, burst) |>
  summarise(
    r2value = mean(r2value, na.rm = TRUE),
    epochs  = mean(epochs,  na.rm = TRUE),
    mae     = mean(mae,     na.rm = TRUE),
    goodch  = mean(goodch_sum,  na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(session_age = paste0(session_age, " mo."))

# Reshape to wide format for paired tests
aper_wide_test = aper_model_fit |>
  pivot_wider(names_from = burst, values_from = c(r2value, epochs, mae, goodch))

# Paired t-tests per visit; FDR-corrected within each metric
p_vals_aper = aper_wide_test |>
  group_by(session_age) |>
  summarise(
    p_r2value = t.test(r2value_burst, r2value_noburst, paired = TRUE)$p.value,
    p_epochs  = t.test(epochs_burst,  epochs_noburst,  paired = TRUE)$p.value,
    p_mae     = t.test(mae_burst,     mae_noburst,     paired = TRUE)$p.value,
    p_goodch  = t.test(goodch_burst,  goodch_noburst,  paired = TRUE)$p.value,
    .groups   = "drop"
  ) |>
  pivot_longer(starts_with("p_"), names_to = "Variable", values_to = "p_value") |>
  mutate(
    Variable = case_when(
      Variable == "p_r2value" ~ "R\u00b2",
      Variable == "p_epochs"  ~ "# of Trials",
      Variable == "p_mae"     ~ "MAE",
      Variable == "p_goodch"  ~ "# Electrodes"
    )
  ) |>
  # Apply FDR correction within each metric, then format p-values
  mutate(p_value = p.adjust(p_value, method = "fdr"), .by = Variable) |>
  mutate(p_value = style_pvalue(p_value))

# Compute descriptive statistics (Mean, SD) for model fit metrics by visit×burst state
aper_model_fit |>
  # Reshape to long format
  pivot_longer(c(r2value, epochs, mae, goodch), names_to = "Variable", values_to = "Value") |>
  # Recode variable names to display labels
  mutate(Variable = case_when(
    Variable == "r2value" ~ "R\u00b2",
    Variable == "epochs"  ~ "# of Trials",
    Variable == "mae"     ~ "MAE",
    Variable == "goodch"  ~ "# Electrodes"
  )) |>
  # Compute Mean and SD per metric, visit, and burst state
  group_by(Variable, session_age, burst) |>
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop") |>
  mutate(Display_Stat = sprintf("%.3f (%.3f)", Mean, SD)) |>
  dplyr::select(Variable, session_age, burst, Display_Stat) |>
  # Reshape to wide format (separate columns for Burst vs. NoBurst)
  pivot_wider(names_from = burst, values_from = Display_Stat) |>
  # Attach p-values from paired t-tests
  left_join(p_vals_aper, by = c("Variable", "session_age")) |>
  # Apply ordered factor for consistent table ordering
  mutate(session_age = factor(session_age, levels = visit_levels)) |>
  arrange(Variable, session_age) |>
  rename(Metric = Variable) |>
  flextable() |>
  merge_v(j = "Metric") |> valign(j = "Metric", valign = "top") |>
  set_header_labels(Variable = "Metric", session_age = "Visit",
                    burst = "Burst Mean (SD)", noburst = "Non-Burst Mean (SD)",
                    p_value = "p (FDR)") |>
  add_footer_lines("Paired t-tests per visit; FDR correction (Benjamini-Hochberg) within each metric. MAE = Mean Absolute Error from Specparam.") |>
  set_caption("Parametrization Model Fit (R\u00b2 and MAE) and Trial Count by Burst State and Visit") |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  theme_booktabs() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(j = 2:4, align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  save_as_docx(path = file.path(path2supptab, "Table_S5_BurstvsNonBurst_Datasets_ParametrizedPSD_PreprocessingComparison.docx"))



