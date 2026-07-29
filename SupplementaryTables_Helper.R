# =============================================================================
# SupplementaryTables_Helper.R
# -----------------------------------------------------------------------------
# Produces SUBMISSION-READY supplementary CSV/XLSX alongside the existing HTML.
#
# THE PROBLEM THIS SOLVES
#   Transform into model-output tables, which were too large to typeset inside the manuscript, and
#   they currently exist only as HTML previews with no table number.
#
# HOW IT WORKS
#   Each analysis script builds ONE display-ready data frame (already filtered,
#   selected and renamed to publication headers), then calls:
#
#       nn_supp_table(tbl_data, id = "SR1",
#                     note = "Smooth terms tested with ...")
#
#   That single call
#     (a) writes  Table_SR1.csv  into Tables/SupplementaryTables/
#     (b) caches the table so SupplementaryData_BuildWorkbook.R can assemble
#         every table into one workbook, one sheet per table.
#   The same frame is then piped into flextable() for the HTML preview, so the
#   HTML and the submitted table can never drift apart.
#
# TABLE NUMBERING
#   Numbers live in NN_TABLE_INDEX below. Edit that block (and only that block)
#   to renumber; the scripts pass an id and the index supplies the caption.
#
# =============================================================================

# This file is a library: it is source()d by the analysis scripts and never run
# on its own. It expects path2root to already exist, which config_paths.R does.
if (!exists("path2root"))
  stop("SupplementaryTables_Helper.R: path2root not found. Source config_paths.R first.",
       call. = FALSE)

suppressPackageStartupMessages({ library(dplyr) })

# -----------------------------------------------------------------------------
# MASTER TABLE INDEX - the single place to change numbering or captions.
# The SR series continues the Supplementary Results convention already used for
# Fig. SR1-SR2; the SM series continues Table SM1.
# -----------------------------------------------------------------------------
# NOTE: the power-spectrum range-comparison tables (Supplementary Methods) are
# produced directly by gtsummary with their statistics embedded in the rendered
# table, so they are not part of this series; they remain HTML/typeset SI tables
# and now carry the Friedman df in their source note.
NN_TABLE_INDEX = list(
  SR1  = list(order =  1,
    title = "Whole-brain GAMM developmental trajectories of aperiodic and alpha metrics.",
    groups = list(
      "Age smooth" = c("EDF", "Ref. df", "F", "P (age)", "P (age, FDR)"),
      "Proportion of retained epochs (covariate)" =
        c("Beta (prop. epochs)", "s.e.", "t", "df", "P (epochs)", "P (epochs, FDR)"))),

  SR2  = list(order =  2,
    title = "Whole-brain GAMM burst x age difference smooths for cycle amplitude.",
    groups = list(
      "Difference smooth (burst x age)" =
        c("EDF (difference smooth)", "Ref. df", "F", "P", "P (FDR)"),
      "Cycle-type contrast" =
        c("Beta (burst vs non-burst)", "s.e.", "t", "df", "P (burst)"))),

  SR3  = list(order =  3,
    title = "Effect of alpha-burst presence on power-spectrum and rhythmicity metrics.",
    groups = list(
      "Burst main effect" = c("Beta (burst vs non-burst)", "s.e.", "95% CI", "t", "df", "P", "P (FDR)"),
      "Burst x visit interaction" = c("Beta (burst x visit)", "s.e. (interaction)",
        "95% CI (interaction)", "t (interaction)", "df (interaction)",
        "P (interaction)", "P (interaction, FDR)"),
      "Model fit" = c("Marginal R2", "Conditional R2"))),

  SR4  = list(order =  4,
    title = "Burst versus non-burst contrasts stratified by visit.",
    groups = list(
      "Descriptives" = c("Burst, mean (s.d.)", "Non-burst, mean (s.d.)"),
      "Contrast"     = c("Beta", "s.e.", "95% CI", "t", "df", "P", "P (FDR)"),
      "Model fit"    = c("Marginal R2", "Conditional R2"))),

  SR5  = list(order =  5,
    title = "Logistic classification of alpha-peak presence at the 1-month visit."),

  SR6  = list(order =  6,
    title = "Contributions of burst properties and rhythmicity to oscillatory alpha, by visit.",
    groups = list(
      "Fixed effect" = c("Beta", "s.e.", "t", "df", "95% CI", "P", "P (FDR)"),
      "Model fit"    = c("Marginal R2", "Conditional R2"))),

  SR7  = list(order =  7,
    title = "Regional (ROI) differences in aperiodic, alpha, burst and rhythmicity metrics."),

  SR8  = list(order =  8,
    title = "Regional (ROI) age smooths for aperiodic, alpha, burst and rhythmicity metrics."),

  SR9  = list(order =  9,
    title = "Robustness check: burst metric trajectories with visit-adjusted epoch length.",
    groups = list(
      "Age smooth" = c("EDF", "Ref. df", "F", "P (age)", "P (age, FDR)"),
      "Proportion of retained epochs (covariate)" =
        c("Beta (prop. epochs)", "s.e.", "t", "df", "P (epochs)", "P (epochs, FDR)"))),

  SR10 = list(order = 10,
    title = "Robustness check: burst x age difference smooths with visit-adjusted epoch length.",
    groups = list(
      "Difference smooth (burst x age)" =
        c("EDF (difference smooth)", "Ref. df", "F", "P", "P (FDR)"),
      "Cycle-type contrast" =
        c("Beta (burst vs non-burst)", "s.e.", "t", "df", "P (burst)"))),

  SR11 = list(order = 11,
    title = "Robustness check: contributions to oscillatory alpha using corrected voltage amplitude.",
    groups = list(
      "Fixed effect" = c("Beta", "s.e.", "t", "df", "95% CI", "P", "P (FDR)"),
      "Model fit"    = c("Marginal R2", "Conditional R2"))),

  SM2  = list(order = 13,
    title = "Alpha-lifespan frequency-range selection: oscillation-based versus all-cycles range.",
    groups = list(
      "Range definition (OB vs AC)" = c("Beta (OB vs AC)", "s.e.", "95% CI", "t", "df", "P"),
      "Visit"                       = c("Beta (visit)", "s.e. (visit)", "t (visit)", "df (visit)", "P (visit)"),
      "Interaction"                 = c("Beta (interaction)", "s.e. (interaction)",
                                        "t (interaction)", "df (interaction)", "P (interaction)"),
      "Model fit"                   = c("Marginal R2", "Conditional R2"))),

  SR12 = list(order = 12,
    title = "Basis-dimension and residual diagnostics for the developmental GAMMs.",
    groups = list(
      "Basis dimension check (mgcv::k.check)" = c("k'", "EDF", "k-index", "P (k-check)"),
      "Residual diagnostics"                  = c("n (residuals)", "Skewness", "Excess kurtosis"))),

  SM3  = list(order = 14,
    title = "Alpha-lifespan frequency-range selection stratified by visit.",
    groups = list(
      "Descriptives" = c("OB, mean (s.d.)", "AC, mean (s.d.)"),
      "Contrast"     = c("Beta (OB vs AC)", "s.e.", "95% CI", "t", "df", "P", "P (FDR)"),
      "Model fit"    = c("Marginal R2", "Conditional R2")))
)

# -----------------------------------------------------------------------------
# Output locations
# -----------------------------------------------------------------------------
nn_supp_dir = function() {
  root = if (exists("path2root", envir = globalenv())) get("path2root", envir = globalenv()) else getwd()
  d = file.path(root, "SupplementaryInformation", "Results", "Tables", "SupplementaryTables")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
nn_cache_dir = function() {
  d = file.path(nn_supp_dir(), "_cache")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# -----------------------------------------------------------------------------
# Degrees-of-freedom helpers
# -----------------------------------------------------------------------------
# Smooth terms are F tests on Ref.df (already in summary()$s.table). Parametric
# terms are t tests on the model residual df, which is what this returns.
nn_gam_df = function(gam_obj) {
  val = suppressWarnings(tryCatch(as.numeric(stats::df.residual(gam_obj)), error = function(e) NA_real_))
  if (length(val) != 1 || !is.finite(val))
    val = suppressWarnings(tryCatch(as.numeric(gam_obj$df.residual), error = function(e) NA_real_))
  if (length(val) != 1 || !is.finite(val))
    val = suppressWarnings(tryCatch(as.numeric(gam_obj$df.null) - as.numeric(gam_obj$rank),
                                    error = function(e) NA_real_))
  if (length(val) != 1 || !is.finite(val)) {
    warning("nn_gam_df(): residual df could not be recovered; returning NA.")
    return(NA_real_)
  }
  val
}

# Satterthwaite df for one lmerTest term. Warns loudly if the model was fitted
# with lme4::lmer (no df column) rather than silently returning an empty column.
nn_lmer_df = function(model, term) {
  ct = as.data.frame(coef(summary(model)))
  if (!"df" %in% names(ct)) {
    warning("nn_lmer_df(): no df column - was the model fitted with lmerTest::lmer()?")
    return(NA_real_)
  }
  if (!term %in% rownames(ct)) return(NA_real_)
  as.numeric(ct[term, "df"])
}

# -----------------------------------------------------------------------------
# nn_supp_table(): register + write ONE numbered supplementary table.
#
#   data  display-ready data frame: already filtered, already reduced to the
#         columns a reader needs, already renamed to publication headers.
#   id    entry in NN_TABLE_INDEX, e.g. "SR1"
#   note  the table's Note line (test, df method, correction, abbreviations)
#   csv   write the per-table .csv immediately (default TRUE)
#
# Returns `data` invisibly so it can be piped straight into flextable().
# -----------------------------------------------------------------------------
nn_supp_table = function(data, id, cols = NULL, title = NULL, note = "", csv = TRUE) {
  stopifnot(is.data.frame(data))

  passthrough = data

  out = as.data.frame(data)
  if (!is.null(cols)) {
    have = intersect(names(cols), names(out))
    miss = setdiff(names(cols), names(out))
    if (length(miss)) {
      warning(sprintf("nn_supp_table('%s'): column(s) not found and skipped: %s",
                      id, paste(miss, collapse = ", ")))
    }
    out = out[, have, drop = FALSE]
    names(out) = unname(cols[have])
  }
  if (is.null(NN_TABLE_INDEX[[id]])) {
    warning(sprintf("nn_supp_table(): id '%s' is not in NN_TABLE_INDEX; add it there.", id))
    idx = list(order = 999, title = "")
  } else {
    idx = NN_TABLE_INDEX[[id]]
  }
  if (is.null(title)) title = idx$title

  # Never ship a column that is empty for every row.
  keep = vapply(out, function(cc) !all(is.na(cc)), logical(1))
  if (any(!keep)) {
    message(sprintf("[Table %s] dropping all-NA column(s): %s",
                    id, paste(names(out)[!keep], collapse = ", ")))
    out = out[, keep, drop = FALSE]
  }

  saveRDS(list(id = id, order = idx$order, title = title, note = note,
               groups = idx$groups,
               data = out, stamp = Sys.time()),
          file.path(nn_cache_dir(), paste0("Table_", id, ".rds")))

  if (isTRUE(csv)) {
    f   = file.path(nn_supp_dir(), paste0("Table_", id, ".csv"))
    con = file(f, open = "w", encoding = "UTF-8")
    # Caption and Note travel with the CSV so the file is self-describing.
    writeLines(paste0("Table ", id, ". ", title), con)
    utils::write.csv(out, con, row.names = FALSE, na = "")
    if (nzchar(note)) writeLines(paste0("Note. ", note), con)
    close(con)
  }

  message(sprintf("[Table %s] %d rows x %d cols%s",
                  id, nrow(out), ncol(out),
                  if (csv) paste0(" -> Table_", id, ".csv") else " (cached)"))
  # Return the ORIGINAL frame so the pipe continues unchanged.
  invisible(passthrough)
}

# -----------------------------------------------------------------------------
# nn_save_table(): save a flextable as HTML and (optionally) Word in one call.
# -----------------------------------------------------------------------------
nn_save_table = function(ft, dir, stem, docx = FALSE) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  flextable::save_as_html(ft, path = file.path(dir, paste0(stem, ".html")))
  if (isTRUE(docx)) {
    tryCatch(flextable::save_as_docx(ft, path = file.path(dir, paste0(stem, ".docx"))),
             error = function(e) warning(sprintf("nn_save_table(): Word export failed for '%s': %s",
                                                 stem, conditionMessage(e))))
  }
  invisible(ft)
}
