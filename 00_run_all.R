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
# Preflight and ordered execution
# -----------------------------------------------------------------------------
required <- c("sliding_samples_RS.csv", "sliding_samples_BS.csv", "mediation_windows.csv", "sem_groups.csv")
missing <- required[!file.exists(file.path(repo, "data", required))]
if (length(missing)) stop("Provide the documented local inputs first: ", paste(missing, collapse = ", "),
                          ". No study data are included in the repository.")
require_packages(c("lme4", "lmerTest", "mediation", "lavaan", "openxlsx"))
# These analyses use separate prepared inputs; mediation does not consume Figure 2 windows.
for (script in c("sliding_window_analysis.R", "mediation.R", "lavaan_SEM_Fig3.R")) {
  message("Running ", script)
  source(file.path(repo, script), local = new.env(parent = globalenv()))
}
