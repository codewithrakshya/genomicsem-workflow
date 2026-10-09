root <- normalizePath(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
library_dir <- Sys.getenv("R_LIBS_USER", unset = file.path(root, "workflow", "R", "library"))
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))

`%||%` <- function(x, y) if (is.null(x) || !nzchar(x)) y else x

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = library_dir)
}

current_r <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
if (dir.exists(file.path(library_dir, "lavaan")) && requireNamespace("lavaan", quietly = TRUE)) {
  lavaan_description <- packageDescription("lavaan")
  lavaan_built_r <- sub("^R ([0-9]+\\.[0-9]+).*$", "\\1", lavaan_description$Built)
  if (nzchar(lavaan_built_r) && !identical(lavaan_built_r, current_r)) {
    message("Reinstalling lavaan for R ", current_r, " (installed package was built for R ", lavaan_built_r, ")")
    remove.packages("lavaan", lib = library_dir)
    install.packages("lavaan", lib = library_dir)
  }
}

pin_file <- file.path(root, "config", "genomicsem_ref.txt")
pinned_ref <- if (file.exists(pin_file)) trimws(readLines(pin_file, warn = FALSE)[1]) else ""
ref <- Sys.getenv("GENOMICSEM_REF", unset = pinned_ref)
if (!nzchar(ref)) stop("Set GENOMICSEM_REF or config/genomicsem_ref.txt to a commit SHA")
repository <- "GenomicSEM/GenomicSEM"
if (nzchar(ref)) repository <- paste0(repository, "@", ref)

reinstall <- identical(tolower(Sys.getenv("GENOMICSEM_REINSTALL", unset = "false")), "true")
installed <- requireNamespace("GenomicSEM", quietly = TRUE)
description <- if (installed) packageDescription("GenomicSEM") else NULL
built_r <- if (installed) sub("^R ([0-9]+\\.[0-9]+).*$", "\\1", description$Built) else ""
reinstall_for_r_version <- installed && nzchar(built_r) && !identical(built_r, current_r)
if (reinstall_for_r_version) {
  message("Rebuilding GenomicSEM for R ", current_r, " (installed package was built for R ", built_r, ")")
}
reinstall <- reinstall || reinstall_for_r_version
if (installed && !reinstall && !identical(description$RemoteSha, ref)) {
  stop("Installed GenomicSEM differs from the requested commit. Use GENOMICSEM_REINSTALL=true to install the pin.")
}
if (!installed || reinstall) {
  remotes::install_github(
    repository,
    lib = library_dir,
    dependencies = TRUE,
    upgrade = "never",
    build_vignettes = FALSE,
    force = TRUE
  )
}


required_packages <- c("GenomicSEM", "lavaan")
unavailable <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
output_dir <- file.path(root, ".work", "setup", "environment")
if (length(unavailable)) {
  stop("Required package(s) could not be loaded after setup: ", paste(unavailable, collapse = ", "))
}
lavaan_built_r <- sub("^R ([0-9]+\\.[0-9]+).*$", "\\1", packageDescription("lavaan")$Built)
if (nzchar(lavaan_built_r) && !identical(lavaan_built_r, current_r)) {
  stop("lavaan was built for R ", lavaan_built_r, " but this workflow is using R ", current_r)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

description <- packageDescription("GenomicSEM")
built_r <- sub("^R ([0-9]+\\.[0-9]+).*$", "\\1", description$Built)
if (nzchar(built_r) && !identical(built_r, current_r)) {
  stop("GenomicSEM was built for R ", built_r, " but this workflow is using R ", current_r)
}
details <- c(
  paste("Installed:", format(Sys.time(), tz = "UTC")),
  paste("R:", R.version.string),
  paste("GenomicSEM version:", as.character(packageVersion("GenomicSEM"))),
  paste("Remote ref:", description$RemoteRef %||% "not recorded"),
  paste("Remote SHA:", description$RemoteSha %||% "not recorded"),
  paste("Library:", find.package("GenomicSEM"))
)
writeLines(details, file.path(output_dir, "installation.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat(paste(details, collapse = "\n"), "\n")
