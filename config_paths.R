# =============================================================================
# config_paths.R — THE ONLY FILE YOU NEED TO EDIT
# -----------------------------------------------------------------------------
# Reproducibility package for Rico-Picó et al., "The emergence and maturation of
# the infant alpha peak reflect a transition from transient bursts to sustained
# oscillations."
#
# Set the three paths below once and the whole pipeline runs.
#
# HOW SCRIPTS FIND THIS FILE
#   Each script looks for config_paths.R in the working directory, the parent
#   folder, and a Code/ subfolder. If none of those work - running with Rscript
#   from an unrelated directory, say - set CODE_FOLDER on the second line of the
#   script you are running, or setwd() to this folder. The script stops with a
#   message naming both options and printing where it searched.
#
#   AlphaBurstRhythm.Rproj is a convenience only: opening it makes RStudio set
#   the working directory here. It has no effect outside RStudio.
# =============================================================================

# --- 1. Where the code lives (this folder) -----------------------------------
path2code = ""

# --- 2. The data folder ------------------------------------------------------
# ONE FLAT FOLDER holding every CSV from the data download. There is no Merged/
# subdirectory: that name belonged to the dataset-creation stage, which a
# released dataset does not have.
#
# It should contain, among others:
#   Aperiodic_Oscillatory_ByCycle_Long.csv
#   BurstProperties_ByCycle_Long.csv
#   LaggedCoh_Hilb_ByCycle_Long.csv
#   CrossVisit_EEG_CleaningDescriptives.csv
#   electrodes.csv
#   Sociodemographic_Descriptives_Long_Updated.csv   (keyed by sujid)
#   Sociodemographic_SES_delinked.csv                (keyed by demo_id)
path2data = ""

# --- 3. Where results are written --------------------------------------------
# Results/Figures/... and Results/Tables/... are created automatically beneath
# this root on first run. It does not need to exist beforehand.
path2root = ""

# =============================================================================
# Nothing below this line needs editing.
# =============================================================================

# Normalise so trailing slashes are optional above and file.path() never
# produces a doubled separator.
for (.p in c("path2code", "path2data", "path2root")) {
  assign(.p, normalizePath(path.expand(get(.p)), winslash = "/", mustWork = FALSE))
}

# Fail loudly and early rather than halfway through a model fit.
for (.p in c("path2code", "path2data")) {
  if (!dir.exists(get(.p)))
    stop("config_paths.R: ", .p, " does not exist:\n  ", get(.p),
         "\nEdit config_paths.R before running any script.", call. = FALSE)
}
if (!dir.exists(path2root)) dir.create(path2root, recursive = TRUE, showWarnings = FALSE)
rm(.p)

# `path2sets` is retained as an alias of path2data so that the dataset-creation
# scripts, which predate the flat layout, keep working unchanged. Analysis
# scripts should use path2data.
path2sets = path2data

message("[config] code   : ", path2code)
message("[config] data   : ", path2data)
message("[config] results: ", path2root)
