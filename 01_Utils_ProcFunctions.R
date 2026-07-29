# ---------------------------------------------------------------------------
# PATHS: edit config_paths.R once; nothing in this file needs changing.
# ---------------------------------------------------------------------------
# config_paths.R is looked for in the working directory. If R was started
# somewhere else, set CODE_FOLDER on the next line to this script's folder.
CODE_FOLDER = ""          # e.g. "~/AlphaBurstRhythm/Code"  (leave "" if unsure)

local({
  cand = c(if (nzchar(CODE_FOLDER)) file.path(path.expand(CODE_FOLDER), "config_paths.R"),
           "config_paths.R", "../config_paths.R", "Code/config_paths.R")
  hit  = cand[file.exists(cand)]
  if (!length(hit))
    stop("config_paths.R not found.\n",
         "Fix either way:\n",
         "  1. setwd(\"/path/to/AlphaBurstRhythm/Code\")   then re-run, or\n",
         "  2. set CODE_FOLDER at the top of this script to that same path.\n",
         "Currently looking from: ", getwd(), call. = FALSE)
  source(hit[1], local = FALSE)
})

# -----------------------------------------------------------------------------
# Script: 01_Utils_ProcFunctions.R
# Purpose: Utility functions shared across analysis scripts.
#          - has_converged()             : checks lmer convergence
#          - proc_iter()                : iterative GLMER classification with ROC 
#          - median_se()                 : median ± SE stat summary for ggplot2
#          - proc_graph_ggplot2()        : plots mean ROC curve from proc_iter output
#          - add_stat_pairwise_wilcox() : paired Wilcoxon test for gtsummary
#         - add_stat_overall_friedman() : Friedman test for gtsummary
# =============================================================================
# Inputs:  None (defines functions only — source this script, do not run standalone)
# Outputs: None
# Dependencies: 00_Setup_PackageInstallation.R
# =============================================================================

# =============================================================================
# FUNCTION 1: has_converged
# =============================================================================
# Checks whether a fitted lmer model has converged without singularity.
# Returns 1 (converged) or 0 (singular / failed).

has_converged = function(model) {
  if (!inherits(model, "merMod")) {
    stop("Input must be a fitted lmer model object.")
  }
  # A NULL value for conv$lme4 indicates successful convergence
  if (is.null(unlist(model@optinfo$conv$lme4))) {
    return(1)  # Converged cleanly
  } else {
    # Treat singular fit as non-converged; all other warnings are accepted
    return(ifelse(isSingular(model), 0, 1))
  }
}

# =============================================================================
# FUNCTION 2: proc_iter
# =============================================================================
# Runs n_iter iterations of a balanced train/test split GLMER classification.
# At each iteration:
#   - Creates an 80/20 train/test split stratified by dep_var
#   - Balances both halves by down-sampling the majority class so classification is not affected by an uneven distribution of the events
#   - Fits a GLMER logistic model
#   - Evaluates AUC, accuracy, sensitivity, specificity, F1, Brier score, MCC
# 
# Returns a list:
#   $auc, $accuracy, $sensitivity, $specificity, $precision, $f1, $brier, $mcc
#   $results : data frame of per-iteration coefficient estimates
#   $roc     : list of pROC roc objects (one per iteration)

proc_iter = function(dataset, n_iter, dep_var, fixed_eff, rand_eff, covariates = NULL) {

  # Pre-allocate metric storage vectors
  auc_values         = numeric(n_iter)
  accuracy_values    = numeric(n_iter)
  sensitivity_values = numeric(n_iter)
  specificity_values = numeric(n_iter)
  precision_values   = numeric(n_iter)
  f1_values          = numeric(n_iter)
  brier_scores       = numeric(n_iter)
  mcc_values         = numeric(n_iter)

  results_list = vector("list", n_iter)
  roc_list     = vector("list", n_iter)

  # Build model formula from supplied fixed effects (rename dep_var + fixed_eff strings to formula terms)
  fixed_eff_string = paste(fixed_eff, collapse = " + ")
  model_formula    = if (is.null(covariates)) {
    as.formula(paste(dep_var, "~", fixed_eff_string, "+", rand_eff))
  } else {
    cov_string = paste(covariates, collapse = " + ")
    as.formula(paste(dep_var, "~", fixed_eff_string, "+", cov_string, "+", rand_eff))
  }

  for (i in seq_len(n_iter)) {

    # --- Create 80/20 stratified split (stratified by dep_var to preserve class balance in train/test) ---
    train_idx  = createDataPartition(dataset[[dep_var]], p = 0.8, list = FALSE)
    train_data = dataset[ train_idx, ]
    test_data  = dataset[-train_idx, ]

    # --- Balance training set: down-sample majority class to match minority class size ---
    train_pos = train_data |> filter(.data[[dep_var]] == 1)
    train_neg = train_data |> filter(.data[[dep_var]] == 0)
    if (nrow(train_pos) < nrow(train_neg)) {
      train_data_bal = bind_rows(train_pos, sample_n(train_neg, nrow(train_pos)))
    } else {
      train_data_bal = bind_rows(sample_n(train_pos, nrow(train_neg)), train_neg)
    }

    # --- Balance test set by down-sampling majority class (same approach as training) ---
    test_pos = test_data |> filter(.data[[dep_var]] == 1)
    test_neg = test_data |> filter(.data[[dep_var]] == 0)
    if (nrow(test_pos) < nrow(test_neg)) {
      test_data_bal = bind_rows(test_pos, sample_n(test_neg, nrow(test_pos)))
    } else {
      test_data_bal = bind_rows(sample_n(test_pos, nrow(test_neg)), test_neg)
    }

    common_ids = intersect(train_data_bal$sujid, test_data_bal$sujid)
    if (!setequal(train_data_bal$sujid, test_data_bal$sujid)) {
      train_data_bal = train_data_bal[train_data_bal$sujid %in% common_ids, ]
      test_data_bal  = test_data_bal[test_data_bal$sujid  %in% common_ids, ]
    }

    # --- Fit GLMER logistic regression and extract coefficients/fit statistics ---
    model_glmer   = glmer(model_formula, data = train_data_bal, family = binomial(link = "logit"))
    summary_glmer = summary(model_glmer)
    r2_values     = r.squaredGLMM(model_glmer)

    # Extract fixed-effect coefficients for specified terms (beta, SE, z-value, p-value) and rename
    iter_results_df = as.data.frame(summary_glmer$coefficients) |>
      rownames_to_column("term") |>
      filter(term %in% fixed_eff) |>
      rename(beta  = "Estimate",
             SE    = "Std. Error",
             z_val = "z value",
             p_val = "Pr(>|z|)") |>
      mutate(iteration      = i,
             marginal_r2    = r2_values[[1]],
             conditional_r2 = r2_values[[2]])

    results_list[[i]] = iter_results_df

    # --- Predict on balanced test set and compute classification metrics (AUC, accuracy, etc.) ---
    test_data_bal$predicted_prob = predict(
      model_glmer, newdata = test_data_bal, type = "response", allow.new.levels = TRUE
    )

    roc_obj       = roc(test_data_bal[[dep_var]], test_data_bal$predicted_prob)
    auc_values[i] = auc(roc_obj)
    roc_list[[i]] = roc_obj

    # Classify predictions: threshold at 0.5
    predicted_class = ifelse(test_data_bal$predicted_prob > 0.5, 1, 0)

    # Compute confusion matrix cells: TP (True Positive), TN, FP, FN
    TP = sum(predicted_class == 1 & test_data_bal[[dep_var]] == 1, na.rm = TRUE)
    TN = sum(predicted_class == 0 & test_data_bal[[dep_var]] == 0, na.rm = TRUE)
    FP = sum(predicted_class == 1 & test_data_bal[[dep_var]] == 0, na.rm = TRUE)
    FN = sum(predicted_class == 0 & test_data_bal[[dep_var]] == 1, na.rm = TRUE)

    # Matthews Correlation Coefficient
    mcc_den = sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) *
                   as.numeric(TN + FP) * as.numeric(TN + FN))
    mcc_values[i] = if (mcc_den == 0) 0 else
      (as.numeric(TP) * as.numeric(TN) - as.numeric(FP) * as.numeric(FN)) / mcc_den

    # Compute other classification metrics
    accuracy_values[i]    = (TP + TN) / (TP + TN + FP + FN)
    sensitivity_values[i] = if ((TP + FN) > 0) TP / (TP + FN) else NA  # Recall: proportion of true positives identified
    specificity_values[i] = if ((TN + FP) > 0) TN / (TN + FP) else NA  # Proportion of true negatives identified
    precision_values[i]   = if ((TP + FP) > 0) TP / (TP + FP) else NA  # Proportion of predicted positives that are true
    f1_values[i] = {
      pr = precision_values[i]; re = sensitivity_values[i]
      if (!is.na(pr) && !is.na(re) && (pr + re) > 0) 2 * pr * re / (pr + re) else NA  # Harmonic mean of precision and recall
    }

    # Brier score: mean squared error between predicted probability and true binary class (lower is better)
    valid_idx     = !is.na(test_data_bal$predicted_prob) & !is.na(test_data_bal[[dep_var]])
    brier_scores[i] = mean(
      (test_data_bal$predicted_prob[valid_idx] - test_data_bal[[dep_var]][valid_idx])^2
    )
  }

  return(list(
    auc         = auc_values,
    roc         = roc_list,
    accuracy    = accuracy_values,
    sensitivity = sensitivity_values,
    specificity = specificity_values,
    precision   = precision_values,
    f1          = f1_values,
    brier       = brier_scores,
    mcc         = mcc_values,
    results     = bind_rows(results_list)
  ))
}


# =============================================================================
# FUNCTION 3: median_se
# =============================================================================
# Computes median ± SE of the median for use as a ggplot2 stat_summary function.
# The factor 1.2533 is the asymptotic correction for the SE of the sample median.
# Usage: stat_summary(fun.data = median_se, geom = "errorbar")

median_se = function(x) {
  x      = na.omit(x)
  med    = median(x)
  se     = sd(x) / sqrt(length(x))
  # Constant 1.2533 ≈ sqrt(pi/2): asymptotic correction for SE of the median relative to SE of the mean
  se_med = 1.2533 * se
  data.frame(y = med, ymin = med - se_med, ymax = med + se_med)
}


# =============================================================================
# FUNCTION 4: proc_graph_ggplot2
# =============================================================================
# Creates a ggplot2 ROC curve from a list of pROC roc() objects (from proc_iter).
# Individual iteration curves are plotted in grey; the average curve is in black.

proc_graph_ggplot2 = function(roc_list) {

  # Extract ROC coordinates (specificity/sensitivity) from each pROC roc object into a single data frame
  df_plot = bind_rows(lapply(seq_along(roc_list), function(i) {
    data.frame(
      specificity = roc_list[[i]]$specificities,
      sensitivity = roc_list[[i]]$sensitivities,
      iter        = i  # Tag each iteration to allow per-iteration line coloring
    )
  }))

  # Interpolate each ROC curve onto a common specificity grid (0 to 1, 100 points); compute mean sensitivity
  grid     = seq(0, 1, length.out = 100)
  avg_sens = rowMeans(do.call(cbind, lapply(roc_list, function(r) {
    approx(x = r$specificities, y = r$sensitivities, xout = grid)$y
  })), na.rm = TRUE)

  # Create averaged ROC curve data frame
  df_avg = data.frame(specificity = grid, sensitivity = avg_sens)

  ggplot(df_plot, aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(aes(group = iter), alpha = 0.05, color = "grey70") +
    geom_line(data = df_avg, aes(x = 1 - specificity, y = sensitivity),
              color = "black", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    labs(x = "1 - Specificity", y = "Sensitivity")
}


# =============================================================================
# FUNCTION 5: add_stat_pairwise_wilcox
# =============================================================================
# Computes all pairwise paired Wilcoxon tests for a given continuous variable
# across two or more groups, with FDR correction across comparisons.
# Designed for use inside gtsummary::add_stat().
#
# Returns a one-row data frame formatted as "V=X; p=Y" for each comparison.
# The column names are formatted as "**Group1 vs. Group2**" so they render
# as bold in the gtsummary table.
#
# Parameters:
#   data    : data frame (passed automatically by gtsummary)
#   variable: name of the dependent variable (string)
#   by      : name of the grouping variable (string)
#   id_var  : name of the subject identifier column (string)
#   method  : p-value correction method (default "fdr")


add_stat_pairwise_wilcox = function(data, variable, by, id_var, method = "fdr", ...) {
  # Extract all unique levels from grouping variable (by) and generate all pairwise combinations
  levels      = sort(unique(data[[by]]))
  comparisons = utils::combn(levels, 2, simplify = FALSE)

  # Perform paired Wilcoxon test for each pairwise comparison
  raw_results = lapply(comparisons, function(pair) {
    g1   = pair[1];  g2   = pair[2]
    d1   = data[data[[by]] == g1, ];  d2   = data[data[[by]] == g2, ]
    # Sort both vectors by subject ID to ensure correct pairing
    vec1 = d1[[variable]][order(d1[[id_var]])]
    vec2 = d2[[variable]][order(d2[[id_var]])]

    test_res = tryCatch(
      wilcox.test(vec1, vec2, paired = TRUE),
      error = function(e) NULL
    )

    if (is.null(test_res)) {
      list(V = NA, p = NA, label = glue::glue("**{g1} vs. {g2}**"))
    } else {
      list(V = test_res$statistic, p = test_res$p.value,
           label = glue::glue("**{g1} vs. {g2}**"))
    }
  })

  # Apply p-value adjustment (FDR default) across all comparisons
  p_raw_vec = sapply(raw_results, function(x) x$p)
  p_adj_vec = p.adjust(p_raw_vec, method = method)

  # Format results as "V=X; p=Y" strings with adjusted p-values; use comparison labels as column names
  formatted_list = lapply(seq_along(raw_results), function(i) {
    res = raw_results[[i]]
    if (is.na(res$V)) return(NA)
    glue::glue("V={round(res$V, 2)}\np={gtsummary::style_pvalue(p_adj_vec[i])}") |>
      setNames(res$label)
  })

  formatted_list |> unlist() |> t() |> as.data.frame()
}


# =============================================================================
# FUNCTION 6: add_stat_overall_friedman
# =============================================================================
# Runs a Friedman repeated-measures test for a given variable across two or
# more ordered groups. Designed for use inside gtsummary::add_stat().
#
# Only subjects present in ALL groups (complete cases) are included in the
# Friedman test to satisfy the paired/balanced design requirement.
#
# Returns a one-row data frame with formatted test result "Chi²=X; p=Y".
#
# Parameters:
#   data    : data frame (passed automatically by gtsummary)
#   variable: name of the dependent variable (string)
#   by      : name of the grouping variable (string)
#   id_var  : name of the subject identifier column (string)
#

add_stat_overall_friedman = function(data, variable, by, id_var, ...) {

  # Identify subjects present in every group level (complete cases for balanced design)
  complete_ids = data |>
    group_by(.data[[id_var]]) |>
    summarise(n_groups = n_distinct(.data[[by]]), .groups = "drop") |>
    filter(n_groups == length(unique(data[[by]]))) |>  # Only subjects with all group levels
    pull(.data[[id_var]])

  # Filter to complete cases: Friedman test requires balanced repeated-measures data
  clean_data = data |> filter(.data[[id_var]] %in% complete_ids)

  # Build Friedman test formula: DV ~ grouping_var | subject_id
  f_formula = as.formula(paste(variable, "~", by, "|", id_var))

  test_res = tryCatch(
    friedman.test(f_formula, data = clean_data),
    error = function(e) NULL
  )

  if (is.null(test_res)) return(data.frame(Overall = "Test error"))

  stat_str = round(test_res$statistic, 2)
  p_str    = gtsummary::style_pvalue(test_res$p.value)

  # Format output with chi-squared statistic and p-value; Unicode χ² = \u03c7\u00b2
  data.frame(Overall = glue::glue("\u03c7\u00b2={stat_str}\np={p_str}"))
}
