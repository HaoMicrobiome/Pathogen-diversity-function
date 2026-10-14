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
require_packages(c("mediation", "openxlsx"))
bootstrap_replicates <- 5000L
seed <- 2025L
expected_windows <- 245L  # Set NULL only for an explicitly different input dataset.
# root_effect_strength is the window mean of paired rhizosphere-bulk soil
# community dissimilarity: the equal-weight mean of bacterial and fungal
# Bray-Curtis dissimilarities, calculated from ASV relative abundances.
# root_effect_strength_without_pathogens uses the same definition after removing
# plant-pathogenic fungal ASVs and renormalizing the fungal abundances.
# X = pathogen abundance; M = root-effect strength; Y = diversity effect.
required <- c("WindowID", "WindowSize", "PlantPathogenRelativeAbundance", "root_effect_strength", "DiversityEffect")
w <- read_input(repo, "mediation_windows.csv", required)
unique_ids(w$WindowID, "WindowID")
numeric_columns(w, c("WindowID", "WindowSize"))
if (!identical(as.numeric(w$WindowID), as.numeric(seq_len(nrow(w)))))
  stop("WindowID must retain the original consecutive order; do not drop failed windows silently.")
if (anyNA(w$WindowSize) || any(w$WindowSize != 60L)) stop("Expected 60-sample windows.")
if (!is.null(expected_windows) && nrow(w) != expected_windows) stop("Expected 245 paired-sample windows, not the 247 Figure 2 RS windows.")
if (is.unsorted(w$WindowID)) stop("Preserve ascending WindowID order.")
cases <- c(AllFungi = "root_effect_strength")
if ("root_effect_strength_without_pathogens" %in% names(w)) cases <- c(cases, PathogenExcluded = "root_effect_strength_without_pathogens")
numeric_columns(w, c("PlantPathogenRelativeAbundance", unname(cases), "DiversityEffect"), allow_missing = TRUE)
if (any(w$PlantPathogenRelativeAbundance <= 0 | w$PlantPathogenRelativeAbundance > 100, na.rm = TRUE))
  stop("Exposure must be arithmetic-mean pathogen percentage in (0, 100].")
for (col in cases) if (any(w[[col]] < 0 | w[[col]] > 1, na.rm = TRUE)) stop("Bray-Curtis means must lie in [0, 1].")

# -----------------------------------------------------------------------------
# Parallel settings
# -----------------------------------------------------------------------------
# Keep the original serial bootstrap and its random-number stream unchanged.

# -----------------------------------------------------------------------------
# Statistical analysis: mediation of window-level associations
# -----------------------------------------------------------------------------
output <- new_output(repo, "mediation")
tables <- list(); audit <- list(); all_paths <- list(); all_effects <- list(); proportions <- list()
for (case in names(cases)) {
  raw <- data.frame(X = w$PlantPathogenRelativeAbundance, M = w[[cases[[case]]]], Y = w$DiversityEffect)
  ok <- complete.cases(raw); raw <- raw[ok, , drop = FALSE]
  if (nrow(raw) < 10L) stop("Too few complete windows.")
  mu <- vapply(raw, mean, numeric(1)); sig <- vapply(raw, sd, numeric(1))
  if (any(!is.finite(sig) | sig <= 0)) stop("Invalid window-level standardization.")
  d <- as.data.frame(Map(function(x, m, s) (x - m) / s, raw, mu, sig))
  med_fit <- lm(M ~ X, data = d)
  out_fit <- lm(Y ~ X + M, data = d)
  set.seed(seed)
  med <- mediation::mediate(model.m = med_fit, model.y = out_fit, treat = "X", mediator = "M",
    boot = TRUE, sims = bootstrap_replicates, control.value = 0, treat.value = 1,
    conf.level = .95, boot.ci.type = "perc")
  effects <- data.frame(AnalysisCase = case, Effect = c("Indirect_ACME", "Direct_ADE", "Total"),
    Estimate = c(med$d.avg, med$z.avg, med$tau.coef),
    CI_low = c(med$d.avg.ci[1], med$z.avg.ci[1], med$tau.ci[1]),
    CI_high = c(med$d.avg.ci[2], med$z.avg.ci[2], med$tau.ci[2]),
    P = c(med$d.avg.p, med$z.avg.p, med$tau.p))
  effects$Significance <- stars(effects$P); all_effects[[case]] <- effects
  path_row <- function(fit, term, from, to, type) {
    cf <- coef(summary(fit)); ci <- confint(fit, term)
    data.frame(AnalysisCase = case, From = from, To = to, PathType = type,
      StandardizedEstimate = cf[term, "Estimate"], SE = cf[term, "Std. Error"],
      DF = df.residual(fit), CI_low = ci[1], CI_high = ci[2], P = cf[term, "Pr(>|t|)"])
  }
  paths <- rbind(path_row(med_fit, "X", "PlantPathogen", case, "a"),
    path_row(out_fit, "M", case, "DiversityEffect", "b"),
    path_row(out_fit, "X", "PlantPathogen", "DiversityEffect", "c_prime"))
  paths$Significance <- stars(paths$P); all_paths[[case]] <- paths
  ratio <- if (is.finite(med$tau.coef) && abs(med$tau.coef) > .Machine$double.eps) med$d.avg / med$tau.coef else NA_real_
  proportions[[case]] <- data.frame(AnalysisCase = case, Measure = "ACME / Total", Proportion = ratio, Percent = 100 * ratio)
  audit[[paste0(case, "_Scaling")]] <- data.frame(Variable = names(mu), Mean = mu, SD = sig)
  audit[[paste0(case, "_Excluded")]] <- data.frame(WindowID = w$WindowID[!ok])
  audit[[paste0(case, "_R2")]] <- data.frame(Model = c("M ~ X", "Y ~ X + M"), R2 = c(summary(med_fit)$r.squared, summary(out_fit)$r.squared))
  sens <- tryCatch(mediation::medsens(med, rho.by = .1), error = function(e) e)
  audit[[paste0(case, "_Sensitivity")]] <- data.frame(Text = if (inherits(sens, "error")) conditionMessage(sens) else capture.output(summary(sens)))
  if (!inherits(sens, "error") && length(sens$rho) == length(sens$d0))
    write_table(data.frame(Rho = sens$rho, ACME = sens$d0), output, paste0("Fig4a_", case, "_sensitivity.csv"))
  write_table(cbind(WindowID = w$WindowID[ok], d), output, paste0("Fig4a_", case, "_standardized_windows.csv"))
}

# -----------------------------------------------------------------------------
# Result export
# -----------------------------------------------------------------------------
write_table(do.call(rbind, all_paths), output, "Fig4a_mediation_paths.csv")
write_table(do.call(rbind, all_effects), output, "Fig4a_mediation_effects.csv")
write_table(do.call(rbind, proportions), output, "Fig4a_mediated_proportion.csv")
save_audit(output, audit, file.path(repo, "data", "mediation_windows.csv"),
  list(Seed = seed, BootstrapReplicates = bootstrap_replicates, Cases = names(cases), NWindows = nrow(w)),
  c("X is arithmetic-mean pathogen percentage, without log transformation. X, M, Y are z-scored across complete windows.",
    "Path P values are OLS t-tests. Effect CIs are 95% percentile bootstrap intervals.",
    "The reported mediated proportion is ACME / Total; no alternative proportion definitions are exported.",
    "The bootstrap and OLS inference do not correct dependence among overlapping windows.",
    "These are observational statistical associations, not proof of causal mediation."))
message("Mediation results: ", output)
