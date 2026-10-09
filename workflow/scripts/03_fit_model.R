args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) {
  stop("Usage: 03_fit_model.R <ldsc-output.rds> <model.txt> [output-label]")
}

source(file.path(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), "workflow", "scripts", "common.R"))
ensure_detectable_cores()
library(GenomicSEM)

ldsc_path <- project_path(args[[1]])
model_path <- project_path(args[[2]])
output_label <- if (length(args) == 3L) args[[3]] else ""
if (nzchar(output_label) && !grepl("^[A-Za-z0-9_.-]+$", output_label)) {
  stop("Output label may contain only letters, numbers, dots, underscores, and hyphens")
}
if (!file.exists(ldsc_path)) stop("LDSC output does not exist: ", ldsc_path)
if (!file.exists(model_path)) stop("Model file does not exist: ", model_path)

covariance_structure <- readRDS(ldsc_path)
model <- paste(readLines(model_path, warn = FALSE), collapse = "\n")
if (!nzchar(trimws(model))) stop("Model file is empty")

fit_warnings <- character()
fit <- withCallingHandlers(
  usermodel(
    covstruc = covariance_structure,
    estimation = "DWLS",
    model = model,
    CFIcalc = TRUE,
    std.lv = TRUE,
    imp_cov = TRUE
  ),
  warning = function(condition) {
    fit_warnings <<- c(fit_warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }
)

output_dir <- if (nzchar(output_label)) {
  ensure_directory(file.path(project_root(), ".work", "models", output_label))
} else {
  ensure_directory(file.path(project_root(), ".work", "models", "model"))
}
saveRDS(fit, file.path(output_dir, "model_fit.rds"))
if (!is.null(fit$modelfit)) {
  write.table(fit$modelfit, file.path(output_dir, "model_fit_statistics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
if (!is.null(fit$results)) {
  write.table(fit$results, file.path(output_dir, "model_parameters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
writeLines(capture.output(fit), file.path(output_dir, "model_fit.txt"))
writeLines(unique(fit_warnings), file.path(output_dir, "model_warnings.txt"))
