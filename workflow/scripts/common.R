project_root <- function() {
  root <- Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd())
  normalizePath(root, mustWork = TRUE)
}

project_path <- function(path) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, mustWork = FALSE))
  }
  file.path(project_root(), path)
}

read_traits <- function(path) {
  path <- project_path(path)
  if (!file.exists(path)) stop("Trait configuration does not exist: ", path)

  traits <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required <- c("trait", "file", "N", "sample_prev", "population_prev")
  missing <- setdiff(required, names(traits))
  if (length(missing)) stop("Missing trait columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(traits$trait)) stop("Trait names must be unique")
  if (any(!nzchar(traits$trait))) stop("Trait names cannot be empty")

  traits$file <- vapply(traits$file, project_path, character(1))
  traits
}

read_workflow_settings <- function(path = "config/workflow.tsv") {
  path <- project_path(path)
  settings <- read.delim(
    path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!identical(names(settings), c("key", "value"))) {
    stop("Workflow settings must have exactly two columns: key and value")
  }
  if (anyDuplicated(settings$key)) stop("Workflow setting keys must be unique")
  stats::setNames(as.list(settings$value), settings$key)
}

read_models <- function(path) {
  path <- project_path(path)
  models <- read.delim(
    path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  required <- c("label", "model_file", "role", "factor_name", "snp_regression")
  missing <- setdiff(required, names(models))
  if (length(missing)) stop("Missing model columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(models$label)) stop("Model labels must be unique")
  if (sum(models$role == "primary") != 1L) stop("Exactly one model must have role=primary")
  models$model_file <- vapply(models$model_file, project_path, character(1))
  models
}

logical_column <- function(values, name) {
  normalized <- toupper(trimws(as.character(values)))
  valid <- normalized %in% c("TRUE", "FALSE")
  if (any(!valid)) stop("Column ", name, " must contain TRUE or FALSE")
  normalized == "TRUE"
}

numeric_or_na <- function(values) {
  values[values == ""] <- NA
  suppressWarnings(as.numeric(values))
}

requested_cores <- function() {
  candidates <- c(
    Sys.getenv("GENOMICSEM_CORES", unset = ""),
    Sys.getenv("NSLOTS", unset = ""),
    Sys.getenv("SLURM_CPUS_PER_TASK", unset = ""),
    "1"
  )
  value <- candidates[nzchar(candidates)][1]
  cores <- suppressWarnings(as.integer(value))
  if (is.na(cores) || cores < 1L) 1L else cores
}

ensure_detectable_cores <- function() {
  detected <- suppressWarnings(parallel::detectCores())
  if (length(detected) == 1L && !is.na(detected) && detected >= 1L) return(invisible(detected))

  # lavaan builds its default options from detectCores(). Some restricted local
  # environments return NA, which prevents even a single-core model from fitting.
  parallel_namespace <- asNamespace("parallel")
  unlockBinding("detectCores", parallel_namespace)
  assign(
    "detectCores",
    function(all.tests = FALSE, logical = TRUE) 1L,
    envir = parallel_namespace
  )
  lockBinding("detectCores", parallel_namespace)
  invisible(1L)
}

ensure_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

# Use one explicit date across stages, including jobs that finish after midnight.
results_directory <- function() {
  run_date <- Sys.getenv("GENOMICSEM_RUN_DATE", unset = format(Sys.time(), "%Y-%m-%d", tz = "America/Los_Angeles"))
  parsed <- suppressWarnings(as.Date(run_date, format = "%Y-%m-%d"))
  if (is.na(parsed) || !identical(format(parsed, "%Y-%m-%d"), run_date)) stop("GENOMICSEM_RUN_DATE must be YYYY-MM-DD")
  ensure_directory(file.path(project_root(), "results", run_date))
}
