args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 01_munge.R <traits.tsv> <hapmap3-reference>")
}

source(file.path(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), "workflow", "scripts", "common.R"))
library(GenomicSEM)

traits <- read_traits(args[[1]])
missing_files <- traits$file[!file.exists(traits$file)]
if (length(missing_files)) stop("Missing GWAS files:\n", paste(missing_files, collapse = "\n"))

hm3 <- project_path(args[[2]])
if (!file.exists(hm3)) stop("HapMap3 reference does not exist: ", hm3)

sample_n <- numeric_or_na(as.character(traits$N))
cores <- requested_cores()
output_dir <- ensure_directory(file.path(project_root(), ".work", "munged"))
old_dir <- setwd(output_dir)
on.exit(setwd(old_dir), add = TRUE)

munge(
  files = as.vector(traits$file),
  hm3 = hm3,
  trait.names = as.vector(traits$trait),
  N = sample_n,
  info.filter = as.numeric(Sys.getenv("GENOMICSEM_INFO_FILTER", unset = "0.9")),
  maf.filter = as.numeric(Sys.getenv("GENOMICSEM_MAF_FILTER", unset = "0.01")),
  parallel = cores > 1L,
  cores = cores,
  overwrite = TRUE
)
