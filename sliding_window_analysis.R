# -----------------------------------------------------------------------------
# Repository location and shared utilities
# -----------------------------------------------------------------------------
repo <- local({
  candidates <- vapply(sys.frames(), function(fr) {
    if (is.null(fr$ofile)) "" else as.character(fr$ofile)[1L]
  }, character(1))
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  candidates <- c(rev(candidates), sub("^--file=", "", args))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    candidates <- c(candidates, rstudioapi::getSourceEditorContext()$path)
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  roots <- c(dirname(normalizePath(candidates, winslash = "/", mustWork = TRUE)), getwd())
  roots <- roots[file.exists(file.path(roots, "analysis_helpers.R"))]
  if (!length(roots)) stop("Open a saved repository script or set the working directory to the repository root.")
  normalizePath(roots[1L], winslash = "/", mustWork = TRUE)
})
source(file.path(repo, "analysis_helpers.R"), local = TRUE)

# -----------------------------------------------------------------------------
# Inputs and analysis settings
# -----------------------------------------------------------------------------
require_packages(c("lme4", "lmerTest", "openxlsx"))
window_sizes <- c(50L, 60L, 70L, 80L)
step_size <- 1L
# Preserve the original two fitted models. Published main windows use Main.
model_names <- c("Main", "Pathogen_adjusted")
formulas <- list(
  Main = Multifunctionality ~ Multidiversity + (1 | Location) + (1 | Vegetation),
  Pathogen_adjusted = Multifunctionality ~ Multidiversity + Log10PlantPathogenRelativeAbundance +
    (1 | Location) + (1 | Vegetation)
)
# -----------------------------------------------------------------------------
# Parallel settings
# -----------------------------------------------------------------------------
workers <- 2L
# Deterministic fits; no window bootstrap or location-omission analysis.

# -----------------------------------------------------------------------------
# Data preparation: already calculated sample-level indices
# -----------------------------------------------------------------------------
required <- c("SampleID", "Location", "Vegetation", "Multidiversity",
              "Multifunctionality", "PlantPathogenRelativeAbundance")
files <- c(RS = "sliding_samples_RS.csv", BS = "sliding_samples_BS.csv")
datasets <- lapply(files, function(f) read_input(repo, f, required))
scaling <- list()
for (label in names(datasets)) {
  d <- datasets[[label]]
  unique_ids(d$SampleID, paste(label, "SampleID"))
  numeric_columns(d, required[4:6])
  if (any(d$PlantPathogenRelativeAbundance <= 0 | d$PlantPathogenRelativeAbundance > 100))
    stop("Pathogen abundance must be percent in (0, 100]; no pseudocount is added.")
  for (v in c("Location", "Vegetation")) {
    d[[v]] <- trimws(as.character(d[[v]]))
    if (anyNA(d[[v]]) || any(!nzchar(d[[v]]))) stop("Missing grouping identifiers.")
  }
  if (nrow(d) < max(window_sizes)) stop("Too few samples for the specified windows.")
  sx <- sd(d$Multidiversity); sy <- sd(d$Multifunctionality)
  if (!is.finite(sx / sy) || min(sx, sy) <= 0) stop("Invalid full-compartment scaling.")
  scaling[[label]] <- data.frame(Dataset = label, N = nrow(d), SD_X = sx, SD_Y = sy, ScaleFactor = sx / sy)
  d$OriginalRow <- seq_len(nrow(d))
  d$Log10PlantPathogenRelativeAbundance <- log10(d$PlantPathogenRelativeAbundance)
  d <- d[order(d$PlantPathogenRelativeAbundance, d$SampleID), , drop = FALSE]
  datasets[[label]] <- d
}
tasks <- list()
for (label in names(datasets)) for (model in model_names) for (size in window_sizes) {
  d <- datasets[[label]]; starts <- seq.int(1L, nrow(d) - size + 1L, by = step_size)
  for (k in seq_along(starts)) tasks[[length(tasks) + 1L]] <- list(
    Dataset = label, Model = model, WindowSize = size, WindowID = k,
    data = d[seq.int(starts[k], length.out = size), , drop = FALSE],
    scale = scaling[[label]]$ScaleFactor)
}

# -----------------------------------------------------------------------------
# Statistical analysis: REML and Satterthwaite slope tests
# -----------------------------------------------------------------------------
fit_window <- function(task, formulas) {
  d <- task$data; size <- task$WindowSize
  row <- data.frame(Dataset = task$Dataset, Model = task$Model, WindowSize = size,
    WindowID = task$WindowID, N = nrow(d), StartSampleID = d$SampleID[1],
    EndSampleID = tail(d$SampleID, 1), PathogenMean = mean(d$PlantPathogenRelativeAbundance),
    MeanLog10Pathogen = mean(d$Log10PlantPathogenRelativeAbundance),
    Beta = NA_real_, SE = NA_real_, DF = NA_real_, P = NA_real_,
    StandardizedBeta = NA_real_, StandardizedSE = NA_real_,
    StandardizedCI_low = NA_real_, StandardizedCI_high = NA_real_,
    MarginalR2 = NA_real_, Singular = NA, Converged = FALSE, Status = "pending", Messages = "")
  messages <- character()
  answer <- tryCatch(withCallingHandlers({
    for (v in c("Location", "Vegetation")) d[[v]] <- factor(d[[v]])
    ng <- vapply(d[c("Location", "Vegetation")], nlevels, integer(1))
    if (min(ng) < 2L || max(ng) >= size) stop("Random-effect groups not estimable in this window.")
    if (min(sd(d$Multidiversity), sd(d$Multifunctionality)) <= 0) stop("Zero within-window variance.")
    form <- formulas[[task$Model]]
    design <- model.matrix(lme4::nobars(form), d)
    if (qr(design)$rank < ncol(design)) stop("Rank-deficient fixed-effect design.")
    fit <- lmerTest::lmer(form, data = d, REML = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa", check.rankX = "stop.deficient",
                                 optCtrl = list(maxfun = 200000)))
    co <- coef(summary(fit))["Multidiversity", ]
    row$Beta <- co["Estimate"]; row$SE <- co["Std. Error"]; row$DF <- co["df"]; row$P <- co["Pr(>|t|)"]
    row$StandardizedBeta <- row$Beta * task$scale
    row$StandardizedSE <- row$SE * task$scale
    row$StandardizedCI_low <- (row$Beta - qt(.975, row$DF) * row$SE) * task$scale
    row$StandardizedCI_high <- (row$Beta + qt(.975, row$DF) * row$SE) * task$scale
    vc <- as.data.frame(lme4::VarCorr(fit))
    vf <- var(as.numeric(lme4::getME(fit, "X") %*% lme4::fixef(fit)))
    row$MarginalR2 <- vf / (vf + sum(vc$vcov[vc$grp != "Residual"]) + sigma(fit)^2)
    row$Singular <- lme4::isSingular(fit, tol = 1e-4)
    cm <- fit@optinfo$conv$lme4$messages
    row$Converged <- all(fit@optinfo$conv$opt == 0) && !length(cm[!grepl("singular", cm, ignore.case = TRUE)])
    row$Status <- if (!row$Converged) "not_converged" else if (row$Singular) "singular" else "ok"
    if (any(!is.finite(c(row$Beta, row$SE, row$DF, row$P, row$StandardizedCI_low, row$StandardizedCI_high))))
      row$Status <- "invalid_inference"
    row
  }, warning = function(e) {messages <<- c(messages, conditionMessage(e)); invokeRestart("muffleWarning")},
     message = function(e) {messages <<- c(messages, conditionMessage(e)); invokeRestart("muffleMessage")}),
  error = function(e) {row$Status <- "failed"; messages <<- c(messages, conditionMessage(e)); row})
  answer$Messages <- paste(unique(messages), collapse = " | ")
  answer
}
output <- new_output(repo, "sliding_windows")
execution <- map_fits(tasks, fit_window, workers, formulas = formulas)
windows <- do.call(rbind, execution$values)
trends <- lapply(split(windows, interaction(windows$Dataset, windows$Model, windows$WindowSize, drop = TRUE)), function(z) {
  good <- z[z$Status %in% c("ok", "singular") & is.finite(z$StandardizedBeta), , drop = FALSE]
  row <- data.frame(Dataset = z$Dataset[1], Model = z$Model[1], WindowSize = z$WindowSize[1],
    NWindows = nrow(good), Intercept = NA_real_, Slope = NA_real_, SE = NA_real_,
    NominalCI_low = NA_real_, NominalCI_high = NA_real_, NominalP = NA_real_, R2 = NA_real_, Status = "insufficient_data")
  if (nrow(good) < 3L || length(unique(good$MeanLog10Pathogen)) < 2L) return(row)
  fit <- lm(StandardizedBeta ~ MeanLog10Pathogen, data = good)
  sm <- summary(fit); cf <- coef(sm)["MeanLog10Pathogen", ]; ci <- confint(fit, "MeanLog10Pathogen")
  row$Intercept <- unname(coef(fit)[1]); row$Slope <- cf["Estimate"]; row$SE <- cf["Std. Error"]
  row$NominalCI_low <- ci[1]; row$NominalCI_high <- ci[2]
  row$NominalP <- cf["Pr(>|t|)"]; row$R2 <- sm$r.squared; row$Status <- "ok"
  row
})

# -----------------------------------------------------------------------------
# Result export; all estimates retain their fit status
# -----------------------------------------------------------------------------
write_table(windows, output, "Fig2_S3_sliding_window_results.csv")
for (label in names(datasets)) write_table(
  windows[windows$Dataset == label & windows$Model == "Main" & windows$WindowSize == 60L, ],
  output, paste0("Fig2ab_Sliding_window_60_", label, ".csv"))
write_table(do.call(rbind, trends), output, "Fig2_S3_OLS_trends.csv")
status <- if (all(windows$Status %in% c("ok", "singular"))) "COMPLETE" else "NEEDS_REVIEW"
save_audit(output, list(Standardization = do.call(rbind, scaling),
  Fit_diagnostics = windows[c("Dataset", "Model", "WindowSize", "WindowID", "Status", "Singular", "Messages")]),
  file.path(repo, "data", files), list(WindowSizes = window_sizes, Step = step_size,
  Backend = execution$backend, Workers = execution$workers),
  c("Slopes use full-compartment SD(X)/SD(Y), not within-window SDs.",
    "PathogenMean is the arithmetic mean of percentages; MeanLog10Pathogen is mean(log10(percent)).",
    "Singular converged models are retained and flagged. Failed fits are not silently replaced.",
    "OLS trends and their inference are unadjusted for overlapping windows.",
    "Pathogen_adjusted adds log10 pathogen abundance; it is NOT the nested sensitivity model of Figure S4."), status)
if (status != "COMPLETE") stop("Window fitting requires review; inspect ", output)
message("Sliding-window results: ", output)
