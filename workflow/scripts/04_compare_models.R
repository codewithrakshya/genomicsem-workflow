args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 2L)) {
  stop("Usage: 04_compare_models.R <models.tsv> [output-dir]")
}

source(file.path(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), "workflow", "scripts", "common.R"))
models <- read_models(args[[1]])
output_dir <- ensure_directory(project_path(if (length(args) == 2L) args[[2]] else intermediate_results_directory()))

read_result <- function(label, filename) {
  path <- project_path(file.path(".work", "models", label, filename))
  if (!file.exists(path)) stop("Missing model result: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

parameter_columns <- c(
  "lhs", "op", "rhs", "Unstand_Est", "Unstand_SE",
  "STD_Genotype", "STD_Genotype_SE", "STD_All", "p_value"
)

parameter_tables <- lapply(models$label, function(label) {
  values <- read_result(label, "model_parameters.tsv")
  missing <- setdiff(parameter_columns, names(values))
  if (length(missing)) stop("Missing parameter columns for ", label, ": ", paste(missing, collapse = ", "))
  values <- values[, parameter_columns]
  names(values)[-(1:3)] <- paste0(names(values)[-(1:3)], "_", label)
  values
})

parameter_comparison <- Reduce(
  function(left, right) merge(left, right, by = c("lhs", "op", "rhs"), all = TRUE, sort = FALSE),
  parameter_tables
)

fit_tables <- lapply(models$label, function(label) {
  values <- read_result(label, "model_fit_statistics.tsv")
  values$model <- label
  values[, c("model", setdiff(names(values), "model"))]
})
fit_comparison <- do.call(rbind, fit_tables)

write.table(
  parameter_comparison, file.path(output_dir, "model_parameter_comparison.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  fit_comparison, file.path(output_dir, "model_fit_comparison.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

summary_lines <- c("GenomicSEM model comparison")
for (index in seq_len(nrow(models))) {
  label <- models$label[index]
  fit <- fit_tables[[index]]
  parameters <- read_result(label, "model_parameters.tsv")
  residuals <- parameters[parameters$op == "~~" & parameters$lhs == parameters$rhs, ]
  negative <- residuals[!is.na(residuals$Unstand_Est) & residuals$Unstand_Est < 0, ]
  summary_lines <- c(
    summary_lines,
    paste0("Model: ", label, " (", models$role[index], ")"),
    paste0("  df: ", fit$df),
    paste0("  chi-square: ", fit$chisq),
    paste0("  p-value: ", fit$p_chisq),
    paste0("  CFI: ", fit$CFI),
    paste0("  SRMR: ", fit$SRMR),
    if (nrow(negative)) paste0("  negative residuals: ", paste(paste0(negative$lhs, "=", negative$Unstand_Est), collapse = ", ")) else "  negative residuals: none"
  )
}
writeLines(summary_lines, file.path(output_dir, "model_comparison_summary.txt"))
