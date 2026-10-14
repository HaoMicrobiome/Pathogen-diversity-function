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
# Inputs and model specification
# -----------------------------------------------------------------------------
require_packages(c("lavaan", "openxlsx"))
variables <- c("MAT", "MAP", "pH", "TC", "EctomycorrhizaeRelativeAbundance",
               "PlantPathogenRelativeAbundance", "Multidiversity", "Multifunctionality")
d <- read_input(repo, "sem_groups.csv", c("SampleID", "Scheme", "PathogenGroup", variables))
numeric_columns(d, variables, allow_missing = TRUE)
schemes <- c("Threshold_main", "TopBottom102_sensitivity")
if (anyNA(d[c("SampleID", "Scheme", "PathogenGroup")]) ||
    any(!d$Scheme %in% schemes) || any(!d$PathogenGroup %in% c("Low", "High"))) stop("Invalid group identifiers.")
if (!"Threshold_main" %in% d$Scheme) stop("The primary grouping is required.")
for (scheme in unique(d$Scheme)) {
  part <- d[d$Scheme == scheme, , drop = FALSE]
  unique_ids(part$SampleID, paste(scheme, "SampleID"))
  if (!setequal(unique(part$PathogenGroup), c("Low", "High"))) stop("Supply both allocated groups for each scheme.")
}
model <- '
  pH ~ MAT + MAP
  TC ~ MAT + MAP
  EctomycorrhizaeRelativeAbundance ~ pH + TC + MAT + MAP
  PlantPathogenRelativeAbundance ~ pH + TC + MAT + EctomycorrhizaeRelativeAbundance + MAP
  Multidiversity ~ PlantPathogenRelativeAbundance + pH + TC + MAT + EctomycorrhizaeRelativeAbundance + MAP
  Multifunctionality ~ PlantPathogenRelativeAbundance + Multidiversity + pH + TC + MAT + EctomycorrhizaeRelativeAbundance + MAP
'

# -----------------------------------------------------------------------------
# Parallel settings
# -----------------------------------------------------------------------------
workers <- 2L  # Independent deterministic group models; PSOCK on Windows.

# -----------------------------------------------------------------------------
# Group-specific scaling before complete-case filtering
# -----------------------------------------------------------------------------
tasks <- list(); counts <- list(); scaling <- list(); exclusions <- list()
for (scheme in schemes[schemes %in% d$Scheme]) for (group in c("Low", "High")) {
  a <- d[d$Scheme == scheme & d$PathogenGroup == group, , drop = FALSE]
  x <- a[variables]
  mu <- vapply(x, mean, numeric(1), na.rm = TRUE)
  sig <- vapply(x, sd, numeric(1), na.rm = TRUE)
  if (any(!is.finite(sig) | sig <= 0)) stop("Cannot scale an empty or constant group variable.")
  ok <- complete.cases(x)
  if (sum(ok) < 20L) stop("Fewer than 20 complete cases in ", scheme, "/", group)
  z <- as.data.frame(Map(function(v, m, s) (v - m) / s, x, mu, sig))
  k <- length(tasks) + 1L
  tasks[[k]] <- list(Scheme = scheme, Group = group, data = z[ok, , drop = FALSE], SampleID = a$SampleID[ok])
  counts[[k]] <- data.frame(Scheme = scheme, Group = group, NAllocated = nrow(a), NAnalyzed = sum(ok), NExcluded = sum(!ok))
  scaling[[k]] <- data.frame(Scheme = scheme, Group = group, Variable = variables, Mean = mu, SD = sig)
  if (any(!ok)) exclusions[[k]] <- data.frame(Scheme = scheme, Group = group, SampleID = a$SampleID[!ok],
    MissingVariables = apply(is.na(x[!ok, , drop = FALSE]), 1, function(v) paste(variables[v], collapse = "; ")))
}

# -----------------------------------------------------------------------------
# Statistical analysis: original six-equation lavaan SEM
# -----------------------------------------------------------------------------
fit_group <- function(task, model) {
  messages <- character()
  result <- tryCatch(withCallingHandlers({
    fit <- lavaan::sem(model, data = task$data, estimator = "ML", meanstructure = TRUE, missing = "listwise")
    if (!isTRUE(lavaan::lavInspect(fit, "converged"))) stop("SEM did not converge.")
    admissible <- isTRUE(lavaan::lavInspect(fit, "post.check"))
    ss <- lavaan::standardizedSolution(fit, type = "std.all", ci = TRUE)
    ss <- ss[ss$op == "~", , drop = FALSE]
    if (nrow(ss) != 26L || anyDuplicated(paste(ss$lhs, ss$rhs))) stop("Expected 26 unique directed paths.")
    if (any(!is.finite(as.matrix(ss[c("est.std", "se", "pvalue", "ci.lower", "ci.upper")]))))
      stop("Non-finite standardized path inference; inspect the original model and input.")
    list(Status = if (admissible) "COMPLETE" else "INADMISSIBLE", paths = ss,
      r2 = lavaan::lavInspect(fit, "r2"),
      fit = lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", "rmsea.ci.lower",
        "rmsea.ci.upper", "rmsea.pvalue", "cfi", "tli", "srmr", "aic", "bic")))
  }, warning = function(e) {messages <<- c(messages, conditionMessage(e)); invokeRestart("muffleWarning")}),
  error = function(e) {messages <<- c(messages, conditionMessage(e)); list(Status = "FAILED")})
  result$Messages <- paste(unique(messages), collapse = " | "); result
}
output <- new_output(repo, "SEM_Fig3")
execution <- map_fits(tasks, fit_group, workers, model = model)
paths <- list(); responses <- list(); global <- list(); diagnostics <- list()
for (i in seq_along(tasks)) {
  task <- tasks[[i]]; result <- execution$values[[i]]
  diagnostics[[i]] <- data.frame(Scheme = task$Scheme, Group = task$Group,
    Status = result$Status, Messages = result$Messages)
  if (is.null(result$paths)) next
  ss <- result$paths
  paths[[i]] <- data.frame(Scheme = task$Scheme, Group = task$Group, Response = ss$lhs, Predictor = ss$rhs,
    StandardizedBeta = ss$est.std, StandardizedSE = ss$se, StandardizedP = ss$pvalue,
    StandardizedCI_low = ss$ci.lower, StandardizedCI_high = ss$ci.upper,
    Significance = stars(ss$pvalue), FitStatus = result$Status)
  responses[[i]] <- data.frame(Scheme = task$Scheme, Group = task$Group,
    Response = names(result$r2), R2 = unname(result$r2), FitStatus = result$Status)
  global[[i]] <- cbind(data.frame(Scheme = task$Scheme, Group = task$Group,
    N = nrow(task$data), FitStatus = result$Status), as.data.frame(as.list(result$fit)))
  write_table(cbind(SampleID = task$SampleID, task$data), output,
    paste0("Fig3_S5_", task$Scheme, "_", task$Group, "_standardized_inputs.csv"))
}

# -----------------------------------------------------------------------------
# Result export: std.all is the only path coefficient scale
# -----------------------------------------------------------------------------
if (length(paths)) {
  write_table(do.call(rbind, paths), output, "Fig3_S5_SEM_path_coefficients.csv")
  write_table(do.call(rbind, responses), output, "Fig3_S5_SEM_response_R2.csv")
  write_table(do.call(rbind, global), output, "Fig3_S5_SEM_global_fit.csv")
}
writeLines(model, file.path(output, "SEM_model.txt"))
status <- if (all(vapply(execution$values, function(x) x$Status == "COMPLETE", logical(1)))) "COMPLETE" else "NEEDS_REVIEW"
audit <- list(Sample_counts = do.call(rbind, counts), Standardization = do.call(rbind, scaling),
  Diagnostics = do.call(rbind, diagnostics))
if (length(exclusions)) audit$Excluded_samples <- do.call(rbind, exclusions)
save_audit(output, audit, file.path(repo, "data", "sem_groups.csv"),
  list(Estimator = "ML", MeanStructure = TRUE, Backend = execution$backend, Workers = execution$workers),
  c("Groups are supplied, not redefined here. No median split or fixed 6-percent regrouping is performed.",
    "Scaling uses all allocated observations per group, na.rm=TRUE, before complete-case filtering.",
    "All path coefficients, SEs, P values and CIs use standardizedSolution(type='std.all').",
    "Random effects, automatic path deletion and group-difference tests are not part of this SEM."), status)
if (status != "COMPLETE") stop("SEM requires review; inspect ", output)
message("SEM results: ", output)
