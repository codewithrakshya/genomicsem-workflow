args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 02_ldsc.R <traits.tsv> <ld-directory> <weight-directory>")
}

source(file.path(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), "workflow", "scripts", "common.R"))
library(GenomicSEM)

traits <- read_traits(args[[1]])
trait_files <- file.path(project_root(), ".work", "munged", paste0(traits$trait, ".sumstats.gz"))
missing_files <- trait_files[!file.exists(trait_files)]
if (length(missing_files)) stop("Missing munged files:\n", paste(missing_files, collapse = "\n"))

ld_dir <- project_path(args[[2]])
weight_dir <- project_path(args[[3]])
if (!dir.exists(ld_dir)) stop("LD-score directory does not exist: ", ld_dir)
if (!dir.exists(weight_dir)) stop("LD-weight directory does not exist: ", weight_dir)

sample_prev <- numeric_or_na(as.character(traits$sample_prev))
population_prev <- numeric_or_na(as.character(traits$population_prev))

work_dir <- ensure_directory(file.path(project_root(), ".work", "ldsc"))
old_dir <- setwd(work_dir)
on.exit(setwd(old_dir), add = TRUE)

ldsc_output <- ldsc(
  traits = trait_files,
  sample.prev = sample_prev,
  population.prev = population_prev,
  ld = paste0(normalizePath(ld_dir), "/"),
  wld = paste0(normalizePath(weight_dir), "/"),
  trait.names = as.vector(traits$trait)
)

saveRDS(ldsc_output, file.path(work_dir, "ldsc_output.rds"))
write.table(ldsc_output[[1]], file.path(work_dir, "sampling_covariance_V.tsv"), sep = "\t", quote = FALSE)
write.table(ldsc_output[[2]], file.path(work_dir, "genetic_covariance_S.tsv"), sep = "\t", quote = FALSE)
genetic_correlation <- cov2cor(as.matrix(ldsc_output[[2]]))
rownames(genetic_correlation) <- colnames(genetic_correlation)
genetic_correlation_table <- data.frame(
  trait = rownames(genetic_correlation),
  genetic_correlation,
  check.names = FALSE
)
write.table(
  genetic_correlation_table,
  file.path(results_directory(), "genetic_correlations.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
genetic_covariance <- as.matrix(ldsc_output[[2]])
rownames(genetic_covariance) <- colnames(genetic_covariance)
genetic_covariance_table <- data.frame(
  trait = rownames(genetic_covariance),
  genetic_covariance,
  check.names = FALSE
)
write.table(
  genetic_covariance_table,
  file.path(results_directory(), "genetic_covariances.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(capture.output(str(ldsc_output)), file.path(work_dir, "ldsc_output_structure.txt"))
