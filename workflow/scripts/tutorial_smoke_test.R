root <- normalizePath(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
library_dir <- Sys.getenv("R_LIBS_USER", unset = file.path(root, "workflow", "R", "library"))
.libPaths(c(library_dir, .libPaths()))

if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  stop("GenomicSEM is not installed. Run workflow/scripts/install_genomicsem.R first.")
}

traits <- c("trait_1", "trait_2", "trait_3")
S <- matrix(
  c(
    0.40, 0.18, 0.14,
    0.18, 0.35, 0.12,
    0.14, 0.12, 0.30
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(traits, traits)
)

# Three traits yield six unique covariance elements. A small diagonal sampling
# covariance matrix is sufficient for this installation smoke test.
V <- diag(c(0.0004, 0.0002, 0.0002, 0.0004, 0.0002, 0.0004))
covariance_structure <- list(V = V, S = S)

model <- "CommonFactor =~ trait_1 + trait_2 + trait_3"
fit <- GenomicSEM::usermodel(
  covstruc = covariance_structure,
  estimation = "DWLS",
  model = model,
  CFIcalc = TRUE,
  std.lv = TRUE,
  imp_cov = TRUE
)

output_dir <- file.path(root, ".work", "setup", "tutorial_smoke_test")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, file.path(output_dir, "smoke_test_fit.rds"))
if (!is.null(fit$modelfit)) {
  write.table(fit$modelfit, file.path(output_dir, "fit_statistics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
if (!is.null(fit$results)) {
  write.table(fit$results, file.path(output_dir, "parameters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
writeLines(capture.output(fit), file.path(output_dir, "smoke_test_output.txt"))

cat("GenomicSEM smoke test completed. Results:", output_dir, "\n")
