root <- normalizePath(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
library_dir <- Sys.getenv("R_LIBS_USER", unset = file.path(root, "workflow", "R", "library"))
.libPaths(c(library_dir, .libPaths()))

required <- c("GenomicSEM", "lavaan", "data.table", "Matrix")
status <- vapply(required, requireNamespace, logical(1), quietly = TRUE)

cat("Project:", root, "\n")
cat("R:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("R library:", library_dir, "\n")
for (package in required) {
  version <- if (status[[package]]) as.character(packageVersion(package)) else "MISSING"
  cat(package, version, sep = ": ", "\n")
}

if (!all(status)) {
  stop("The environment is incomplete. Run workflow/scripts/install_genomicsem.R first.")
}

output_dir <- file.path(root, ".work", "setup", "environment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("Environment check passed.\n")
