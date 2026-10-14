# -----------------------------------------------------------------------------
# Input validation and package requirements
# -----------------------------------------------------------------------------
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages first: ", paste(missing, collapse = ", "))
}
read_input <- function(repo, name, required) {
  path <- file.path(repo, "data", name)
  if (!file.exists(path)) stop("Input not distributed with this repository: data/", name,
                             ". See data/README.md for the schema.")
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
  if (!nrow(d) || anyDuplicated(names(d)) || !all(required %in% names(d)))
    stop("Empty table, duplicate columns or missing required columns in ", name)
  d
}
numeric_columns <- function(d, columns, allow_missing = FALSE) {
  for (nm in columns) {
    x <- d[[nm]]
    if (!is.numeric(x) || any(is.infinite(x)) || (!allow_missing && anyNA(x)))
      stop("Invalid numeric values in ", nm)
  }
}
unique_ids <- function(x, label) {
  if (anyNA(x) || any(!nzchar(trimws(as.character(x)))) || anyDuplicated(x))
    stop("Missing, empty or duplicate ", label)
}
stars <- function(p) ifelse(!is.finite(p), "NA", ifelse(p < .001, "***",
                         ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))

# -----------------------------------------------------------------------------
# Parallel execution of independent deterministic fits
# -----------------------------------------------------------------------------
map_fits <- function(tasks, fun, workers = 2L, ...) {
  if (length(workers) != 1L || !is.finite(workers) || workers < 1L || workers != floor(workers))
    stop("workers must be a positive integer.")
  workers <- min(as.integer(workers), length(tasks))
  cl <- NULL
  if (workers > 1L) cl <- tryCatch(parallel::makePSOCKcluster(workers), error = function(e) {
    warning("PSOCK unavailable; using serial execution: ", conditionMessage(e)); NULL
  })
  if (is.null(cl)) return(list(values = lapply(tasks, fun, ...), workers = 1L, backend = "serial"))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    NULL
  })
  list(values = parallel::parLapply(cl, tasks, fun, ...), workers = workers, backend = "PSOCK")
}

# -----------------------------------------------------------------------------
# Output and audit utilities
# -----------------------------------------------------------------------------
new_output <- function(repo, analysis) {
  parent <- file.path(repo, "outputs")
  dir.create(parent, showWarnings = FALSE, recursive = TRUE)
  stem <- paste0(analysis, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  for (i in 0:9999) {
    target <- file.path(parent, if (i == 0L) stem else paste0(stem, "_", i))
    if (dir.create(target, showWarnings = FALSE)) {
      writeLines("RUNNING", file.path(target, "STATUS.txt")); return(target)
    }
    if (!dir.exists(target)) stop("Cannot create output directory.")
  }
  stop("Cannot allocate output directory.")
}
write_table <- function(d, folder, name) {
  previous <- options(digits = 17); on.exit(options(previous), add = TRUE)
  write.csv(d, file.path(folder, name), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
save_audit <- function(folder, tables, inputs, settings, notes, status = "COMPLETE") {
  input_info <- data.frame(File = basename(inputs), MD5 = unname(tools::md5sum(inputs)))
  tables$Input_checksums <- input_info
  tables$Settings <- data.frame(Setting = names(settings), Value = vapply(settings, toString, character(1)))
  tables$Notes <- data.frame(Note = notes)
  tables$Session <- data.frame(Text = capture.output(sessionInfo()))
  openxlsx::write.xlsx(tables, file.path(folder, "Analysis_audit.xlsx"), overwrite = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(folder, "sessionInfo.txt"))
  writeLines(status, file.path(folder, "STATUS.txt"))
}
