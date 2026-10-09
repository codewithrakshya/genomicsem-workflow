root <- normalizePath(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
library_dir <- Sys.getenv("R_LIBS_USER", unset = file.path(root, "workflow", "R", "library"))
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))

`%||%` <- function(x, y) if (is.null(x) || !nzchar(x)) y else x

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = library_dir)
}

pin_file <- file.path(root, "config", "genomicsem_ref.txt")
pinned_ref <- if (file.exists(pin_file)) trimws(readLines(pin_file, warn = FALSE)[1]) else ""
ref <- Sys.getenv("GENOMICSEM_REF", unset = pinned_ref)
if (!nzchar(ref)) stop("Set GENOMICSEM_REF or config/genomicsem_ref.txt to a commit SHA")
repository <- "GenomicSEM/GenomicSEM"
if (nzchar(ref)) repository <- paste0(repository, "@", ref)

reinstall <- identical(tolower(Sys.getenv("GENOMICSEM_REINSTALL", unset = "false")), "true")
installed <- requireNamespace("GenomicSEM", quietly = TRUE)
if (installed && !reinstall && !identical(packageDescription("GenomicSEM")$RemoteSha, ref)) {
  stop("Installed GenomicSEM differs from the requested commit. Use GENOMICSEM_REINSTALL=true to install the pin.")
}
if (!installed || reinstall) {
  remotes::install_github(
    repository,
    lib = library_dir,
    dependencies = TRUE,
    upgrade = "never",
    build_vignettes = FALSE
  )
}

output_dir <- file.path(root, ".work", "setup", "environment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

description <- packageDescription("GenomicSEM")
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
