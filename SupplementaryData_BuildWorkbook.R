# =============================================================================
# SupplementaryData_BuildWorkbook.R
# -----------------------------------------------------------------------------
# Typesets every numbered supplementary table into ONE workbook:
#
#   Supplementary_Tables.xlsx
#     - "Contents" sheet
#     - one sheet per table, formatted as a journal table:
#         caption row
#         spanning sub-header row (when a table reports more than one thing,
#           e.g. the age smooth AND the epochs covariate in the whole-brain
#           models) with a rule under each span
#         bold header row with rules above and below
#         data, with repeated key values merged vertically
#         Note: Tables went through manual typesetting after the document creation. Format may slightly differ. 
#
# Formatting applied automatically from the publication headers, so no call
# site needs to change:
#   * "95% CI lower" + "95% CI upper"  ->  single "95% CI" column as [l, u]
#   * P columns                         ->  "<0.0001" instead of 0.0000
#   * a significance column (***/**/*)  ->  added after the FDR-corrected P
#   * embedded newlines in labels       ->  replaced with a space
#
# RUN ORDER: last, after all analysis scripts.
# =============================================================================

# ---- PATHS --------------------------------------------
# path2root comes from config_paths.R Do not need modify
# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# PATHS: edit config_paths.R once; nothing in this file needs changing.
# ---------------------------------------------------------------------------
# config_paths.R is looked for in the working directory. If R was started
# somewhere else, set CODE_FOLDER on the next line to this script's folder.
CODE_FOLDER = ""          # e.g. "~/AlphaBurstRhythm/Code"  If you have opened the code from the project, you don't need to modify this line. Otherwise, select where the code folder that contains the config_paths.R is

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

source(file.path(path2code, "SupplementaryTables_Helper.R")) # We call the helper function

suppressPackageStartupMessages({
  library(openxlsx); library(dplyr)
})

DIGITS   = 3      # estimates, s.e., CI, R2. This follows journal standards and only affects how the rounding works. 
DIGITS_P = 4      # P values
P_FLOOR  = 1e-4   # below this, print "<0.0001"

# -----------------------------------------------------------------------------
#  Helpers
# -----------------------------------------------------------------------------
is_p_col   = function(h) grepl("^P\\b|^P \\(|^p-value|^P value", h) #Determine if a column is a P value statistic
is_fdr_col = function(h) grepl("FDR", h, fixed = TRUE) # Determine if a column is a P value FDR corrected 
# t, df, F and all estimates print at DIGITS (3) decimals; only counts and the
# visit label are whole numbers. "^t" is deliberately NOT here.
is_int_col = function(h) grepl("^n \\(|^N \\(|^Iterations$|^Visit", h)

fmt_p = function(p) {
  ifelse(is.na(p), "",
    ifelse(p < P_FLOOR, paste0("<", format(P_FLOOR, scientific = FALSE)),
           formatC(p, format = "f", digits = DIGITS_P)))
}
stars = function(p) { # P to stars marking
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}
clean_txt = function(x) {
  x = gsub("[\r\n]+", " ", as.character(x))   # labels built with "\n" in flextable
  gsub("\\s+", " ", trimws(x))
}

# -----------------------------------------------------------------------------
# ROW ORDER
# Every table is sorted to mirror the order in which the manuscript presents the
# results and the metrics presented in Table 1. 
#
#     1) Power spectrum, 2) Bursts, 3) Alpha lifespan
#
# Within the power spectrum, aperiodic before oscillatory.
# -----------------------------------------------------------------------------
NN_FAMILY_ORDER = c(
  "powerspec|power spec|parametr|paramet|psd",   # 1. parameterized power spectrum
  "burst",                                       # 2. burst properties
  "lifespan|lach|rhythm"                         # 3. alpha lifespan / rhythmicity
)

NN_METRIC_ORDER = c(
  # --- power spectrum: aperiodic, then oscillatory ---
  "offset",
  "slope",
  "prop.*peak",
  "peak freq",
  "^peak amp|^corrected peak amp",
  "band power",
  "peak-to-peak|volt\\.? ?amp",
  "band amp",
  "prop.*epoch",
  "prop.*cycle",
  "duration",
  "lifespan|lach"
)

# The Fig. 5d / S8d order
NN_PREDICTOR_ORDER = c(
  "band amp|volt\\.? ?amp",
  "duration",
  "prop.*epoch",
  "prop.*cycle",
  "lifespan|lach"
)

nn_rank = function(x, patterns) {
  x = tolower(trimws(as.character(x)))
  vapply(x, function(v) {
    if (is.na(v) || !nzchar(v)) return(as.numeric(length(patterns) + 1))
    hit = which(vapply(patterns, function(p) grepl(p, v), logical(1)))
    if (length(hit)) as.numeric(hit[1]) else as.numeric(length(patterns) + 1)
  }, numeric(1), USE.NAMES = FALSE)
}

# Visits print as "1 mo.", "6 mo." ... so a plain sort would give 1, 12, 18, 6.
nn_visit_num = function(x) {
  v = suppressWarnings(as.numeric(sub("^[^0-9-]*(-?[0-9.]+).*$", "\\1", as.character(x))))
  ifelse(is.na(v), Inf, v)
}

# The model-family column arrives as internal codes ("powerspectrum", "burst",
# "alphalifespan"). Relabel for print AFTER ordering, so the patterns above are
# matched against the raw codes.
NN_FAMILY_LABELS = c(
  powerspectrum = "Power spectrum",
  parametrized  = "Power spectrum",
  parameterized = "Power spectrum",
  burst         = "Burst properties",
  alphalifespan = "Alpha lifespan"
)
                
nn_relabel_family = function(d) {
  for (col in intersect(c("Model family", "Component"), names(d))) {
    if (col != "Model family") next
    key = tolower(trimws(as.character(d[[col]])))
    hit = key %in% names(NN_FAMILY_LABELS)
    d[[col]][hit] = unname(NN_FAMILY_LABELS[key[hit]])
  }
  d
}

nn_order_rows = function(d) {
  if (!nrow(d)) return(d)
  keys = list()
  if ("Model family"       %in% names(d)) keys[[length(keys)+1]] = nn_rank(d[["Model family"]], NN_FAMILY_ORDER)
  if ("Dependent variable" %in% names(d)) keys[[length(keys)+1]] = nn_rank(d[["Dependent variable"]], NN_METRIC_ORDER)
  if ("Metric"             %in% names(d)) keys[[length(keys)+1]] = nn_rank(d[["Metric"]], NN_METRIC_ORDER)
  if ("Visit (months)"     %in% names(d)) keys[[length(keys)+1]] = nn_visit_num(d[["Visit (months)"]])
  if ("Predictor"          %in% names(d)) keys[[length(keys)+1]] = nn_rank(d[["Predictor"]], NN_PREDICTOR_ORDER)
  if (!length(keys)) return(d)
  keys[[length(keys)+1]] = seq_len(nrow(d))   # stable: preserve original order within ties
  d[do.call(order, keys), , drop = FALSE]
}

# Insert `newcol` immediately after column `pos`, by name assignment only.
# cbind(df, df, NULL) can collapse to a zero-row frame inside cbind.data.frame,
# which is what produced "arguments imply differing number of rows: n, 0".
insert_after = function(d, pos, newcol, name) {
  out = d[, seq_len(pos), drop = FALSE]
  out[[name]] = newcol
  if (pos < ncol(d)) {
    rest = d[, (pos + 1):ncol(d), drop = FALSE]
    for (nm in names(rest)) out[[nm]] = rest[[nm]]
  }
  out
}

# Combine "<x> lower" / "<x> upper" pairs into one bracketed column, in place.
combine_ci = function(d) {
  if (!nrow(d) || !ncol(d)) return(d)
  # Matches "95% CI lower" and suffixed variants such as
  # "95% CI lower (interaction)", pairing each with its "upper" twin.
  lows = grep(" lower", names(d), fixed = TRUE)
  for (lo in rev(lows)) {
    hi = which(names(d) == sub(" lower", " upper", names(d)[lo], fixed = TRUE))
    if (!length(hi)) next
    vals = ifelse(is.na(d[[lo]]) | is.na(d[[hi]]), "",
                  sprintf("[%.3f, %.3f]", d[[lo]], d[[hi]]))
    d[[lo]] = vals
    names(d)[lo] = sub(" lower", "", names(d)[lo], fixed = TRUE)
    d = d[, -hi, drop = FALSE]
  }
  d
}

# -----------------------------------------------------------------------------
# Load cache files with the table information
# -----------------------------------------------------------------------------
cache = list.files(nn_cache_dir(), pattern = "\\.rds$", full.names = TRUE)
if (!length(cache)) stop("No tables in ", nn_cache_dir(), " - run the analysis scripts first.")
tabs = lapply(cache, readRDS)
tabs = tabs[order(vapply(tabs, function(x) x$order, numeric(1)))]
missing = setdiff(names(NN_TABLE_INDEX), vapply(tabs, function(x) x$id, character(1)))
if (length(missing)) warning("Declared but not produced: ", paste(missing, collapse = ", "))
message(sprintf("[nn] typesetting %d supplementary tables", length(tabs)))

# -----------------------------------------------------------------------------
# Styles
# -----------------------------------------------------------------------------
FONT = "Aptos"
s_caption = createStyle(fontName = FONT, fontSize = 11, textDecoration = "bold",
                        halign = "left", valign = "center", wrapText = FALSE)
s_group   = createStyle(fontName = FONT, fontSize = 10, textDecoration = "bold",
                        halign = "center", valign = "center",
                        border = "bottom", borderStyle = "thin")
s_head    = createStyle(fontName = FONT, fontSize = 10, textDecoration = "bold",
                        halign = "center", valign = "center", wrapText = TRUE,
                        border = "TopBottom", borderStyle = c("medium", "medium"))
s_key1    = createStyle(fontName = FONT, fontSize = 10, textDecoration = "bold",
                        halign = "center", valign = "center", wrapText = TRUE)
s_key2    = createStyle(fontName = FONT, fontSize = 10, textDecoration = "italic",
                        halign = "right",  valign = "center", wrapText = TRUE)
s_txt     = createStyle(fontName = FONT, fontSize = 10, halign = "left",   valign = "center")
s_ctr     = createStyle(fontName = FONT, fontSize = 10, halign = "center", valign = "center")
s_num     = createStyle(fontName = FONT, fontSize = 10, halign = "center", valign = "center",
                        numFmt = paste0("0.", strrep("0", DIGITS)))
s_int     = createStyle(fontName = FONT, fontSize = 10, halign = "center", valign = "center",
                        numFmt = "0")
s_star    = createStyle(fontName = FONT, fontSize = 10, halign = "left", valign = "center")
s_rule    = createStyle(border = "top", borderStyle = "thin")
s_bottom  = createStyle(border = "bottom", borderStyle = "medium")
s_note    = createStyle(fontName = FONT, fontSize = 9, fontColour = "#444444",
                        wrapText = TRUE, valign = "top")

wb = createWorkbook()

# -----------------------------------------------------------------------------
# Contents: This is only created in the public data so readers can have the contents of the book
# -----------------------------------------------------------------------------
addWorksheet(wb, "Contents")
writeData(wb, "Contents", "Supplementary Tables", startRow = 1, startCol = 1)
addStyle(wb, "Contents", createStyle(fontName = FONT, fontSize = 13, textDecoration = "bold"),
         rows = 1, cols = 1)
writeData(wb, "Contents", paste(
  "The emergence and maturation of the infant alpha peak reflect a transition from",
  "transient bursts to sustained oscillations. These tables support main text analysis",
  "Every comparison reported in the manuscript appears here with its",
  "test statistic, degrees of freedom and exact P value."),
  startRow = 3, startCol = 1)
addStyle(wb, "Contents", s_note, rows = 3, cols = 1)
cont = data.frame(Table = paste("Table", vapply(tabs, function(x) x$id, character(1))),
                  Title = vapply(tabs, function(x) x$title, character(1)),
                  Rows  = vapply(tabs, function(x) nrow(x$data), integer(1)),
                  stringsAsFactors = FALSE)
writeData(wb, "Contents", cont, startRow = 6, startCol = 1, headerStyle = s_head)
setColWidths(wb, "Contents", cols = 1:3, widths = c(14, 100, 8))
freezePane(wb, "Contents", firstActiveRow = 7)

# -----------------------------------------------------------------------------
# One typeset sheet per table
# -----------------------------------------------------------------------------
for (tb in tabs) {
  sh  = paste("Table", tb$id)
  message("[nn] typesetting ", sh, " (", nrow(tb$data), " rows)")
  dat = tb$data
  names(dat) = clean_txt(names(dat))
  for (j in seq_along(dat)) if (is.character(dat[[j]]) || is.factor(dat[[j]]))
    dat[[j]] = clean_txt(dat[[j]])
  dat = combine_ci(dat)
  dat = nn_order_rows(dat)      # manuscript order: power spectrum -> burst -> lifespan
  dat = nn_relabel_family(dat)  # internal codes -> printable family names

  # A significance column after EVERY FDR-corrected P (or after the plain P when
  # a table reports no FDR column), so tables with two effect blocks - e.g. the
  # burst main effect and the burst x visit interaction - get a star for each.
  if (nrow(dat) > 0) {
    p_cols  = which(vapply(names(dat), is_p_col, logical(1)))
    anchors = p_cols[vapply(names(dat)[p_cols], is_fdr_col, logical(1))]
    if (!length(anchors) && length(p_cols)) anchors = max(p_cols)
    if (length(anchors))
      anchors = anchors[vapply(anchors, function(a) is.numeric(dat[[a]]), logical(1))]
    anchors = sort(as.integer(anchors), decreasing = TRUE)   # insert right-to-left
    for (k in seq_along(anchors)) {
      a = anchors[k]
      # Star columns carry unique blank names (" ", "  ", ...) because Excel
      # cannot hold two columns with the same header.
      dat = insert_after(dat, a, stars(dat[[a]]), strrep(" ", k))
    }
  }
  is_star_col = function(h) grepl("^ +$", h)
  # P columns as formatted text so 0.0000 never appears
  for (j in which(vapply(names(dat), is_p_col, logical(1)))) {
    if (is.numeric(dat[[j]])) dat[[j]] = fmt_p(dat[[j]])
  }
  for (j in seq_along(dat)) if (is.character(dat[[j]]) || is.factor(dat[[j]]))
    dat[[j]] = clean_txt(dat[[j]])

  ncol_t = ncol(dat)
  groups = tb$groups
  has_g  = !is.null(groups) && length(groups)

  cap_row = 1
  grp_row = if (has_g) 2 else NA
  hdr_row = if (has_g) 3 else 2
  first   = hdr_row + 1
  last    = hdr_row + nrow(dat)

  addWorksheet(wb, sh)

  # caption across the full width
  writeData(wb, sh, paste0("Table ", tb$id, ". ", tb$title), startRow = cap_row, startCol = 1)
  mergeCells(wb, sh, cols = 1:ncol_t, rows = cap_row)
  addStyle(wb, sh, s_caption, rows = cap_row, cols = 1:ncol_t, gridExpand = TRUE)

  # spanning sub-headers
  if (has_g) {
    for (gname in names(groups)) {
      cols_in = which(names(dat) %in% groups[[gname]])
      if (!length(cols_in)) next
      rng = min(cols_in):max(cols_in)
      # absorb a trailing significance column so the span stays contiguous
      nxt = max(rng) + 1
      if (nxt <= ncol_t && is_star_col(names(dat)[nxt])) rng = min(rng):nxt
      writeData(wb, sh, gname, startRow = grp_row, startCol = min(rng))
      if (length(rng) > 1) mergeCells(wb, sh, cols = rng, rows = grp_row)
      addStyle(wb, sh, s_group, rows = grp_row, cols = rng, gridExpand = TRUE)
    }
  }

  # header + data
  writeData(wb, sh, dat, startRow = hdr_row, startCol = 1, headerStyle = s_head)

  key_n = 0
  for (j in seq_len(ncol_t)) {
    h = names(dat)[j]
    col_is_key = j <= 3 && !is.numeric(dat[[j]]) && !is_star_col(h)
    if (col_is_key) key_n = max(key_n, j)
    style =
      if (is_star_col(h))                s_star
      else if (col_is_key && j == 1)     s_key1
      else if (col_is_key)               s_key2
      else if (!is.numeric(dat[[j]]))    s_ctr
      else if (is_int_col(h))            s_int
      else                               s_num
    addStyle(wb, sh, style, rows = first:last, cols = j, gridExpand = TRUE)
  }

  # merge repeated key values vertically and rule between blocks
  if (key_n >= 1 && nrow(dat) > 1) {
    for (j in seq_len(key_n)) {
      key = do.call(paste, c(lapply(seq_len(j), function(k) as.character(dat[[k]])), sep = "\r"))
      run_start = 1
      for (i in 2:(nrow(dat) + 1)) {
        if (i > nrow(dat) || key[i] != key[run_start]) {
          if (i - 1 > run_start)
            mergeCells(wb, sh, cols = j, rows = (first + run_start - 1):(first + i - 2))
          if (j == 1 && i <= nrow(dat))
            addStyle(wb, sh, s_rule, rows = first + i - 1, cols = 1:ncol_t,
                     gridExpand = TRUE, stack = TRUE)
          run_start = i
        }
      }
    }
  }

  addStyle(wb, sh, s_bottom, rows = last, cols = 1:ncol_t, gridExpand = TRUE, stack = TRUE)

  if (nzchar(tb$note)) {
    nrow_note = last + 1
    writeData(wb, sh, paste0("Note. ", tb$note,
                             " *P < 0.05; **P < 0.01; ***P < 0.001."),
              startRow = nrow_note, startCol = 1)
    mergeCells(wb, sh, cols = 1:ncol_t, rows = nrow_note)
    addStyle(wb, sh, s_note, rows = nrow_note, cols = 1:ncol_t, gridExpand = TRUE)
    setRowHeights(wb, sh, rows = nrow_note, heights = 46)
  }

  w = pmin(pmax(nchar(names(dat)) + 3, 9), 30)
  if (key_n >= 1) w[1:key_n] = pmin(pmax(nchar(names(dat)[1:key_n]) + 6, 16), 30)
  w[vapply(names(dat), is_star_col, logical(1))] = 5
  setColWidths(wb, sh, cols = seq_len(ncol_t), widths = w)
  setRowHeights(wb, sh, rows = hdr_row, heights = 30)
  freezePane(wb, sh, firstActiveRow = first, firstActiveCol = min(key_n + 1, ncol_t))
}

out = file.path(nn_supp_dir(), "Supplementary_Statistical_Tables.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)
message("[nn] wrote ", out)
message("[nn] tables: ", paste(vapply(tabs, function(x) x$id, character(1)), collapse = ", "))
