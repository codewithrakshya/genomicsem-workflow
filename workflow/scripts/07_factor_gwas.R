args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else "pilot"
settings_path <- if (length(args) >= 2L) args[[2]] else "config/workflow.tsv"
if (!mode %in% c("pilot", "full")) stop("Usage: 07_factor_gwas.R [pilot|full] [workflow.tsv]")

source(file.path(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), "workflow", "scripts", "common.R"))
ensure_detectable_cores()
suppressPackageStartupMessages(library(GenomicSEM))

settings <- read_workflow_settings(settings_path)
traits <- read_traits(settings$traits_config)
models <- read_models(settings$models_config)
primary <- models[models$role == "primary", , drop = FALSE]
root <- project_root()
covariance_structure <- readRDS(project_path(".work/ldsc/ldsc_output.rds"))
model <- paste(c(readLines(primary$model_file, warn = FALSE), primary$snp_regression), collapse = "\n")

if (mode == "pilot") {
  input_dir <- project_path(".work/factor_gwas/pilot")
  files <- file.path(input_dir, paste0(traits$trait, ".tsv.gz"))
  reference <- file.path(input_dir, "reference.tsv")
  output_dir <- ensure_directory(input_dir)
  output_stem <- paste0(settings$analysis_id, "_factor_gwas_pilot")
} else {
  files <- vapply(traits$full_file, project_path, character(1))
  reference <- project_path(settings$gwas_reference)
  output_dir <- ensure_directory(file.path(root, ".work", "factor_gwas", "full"))
  output_stem <- paste0(settings$analysis_id, "_factor_gwas")
}

missing <- c(files[!file.exists(files)], reference[!file.exists(reference)])
if (length(missing)) stop("Missing factor-GWAS input(s): ", paste(missing, collapse = ", "))

old_working_directory <- getwd()
on.exit(setwd(old_working_directory), add = TRUE)
setwd(output_dir)

message("Harmonizing ", mode, " SNP inputs with the configured allele/MAF reference")
snps <- sumstats(
  files = files,
  ref = reference,
  trait.names = traits$trait,
  se.logit = logical_column(traits$se_logit, "se_logit"),
  OLS = logical_column(traits$OLS, "OLS"),
  linprob = rep(FALSE, nrow(traits)),
  N = numeric_or_na(as.character(traits$N)),
  betas = rep(TRUE, nrow(traits)),
  info.filter = as.numeric(settings$info_filter),
  maf.filter = as.numeric(settings$maf_filter),
  keep.indel = FALSE,
  parallel = FALSE,
  direct.filter = TRUE
)
if (!nrow(snps)) stop("No SNPs remained after multivariate harmonization")
saveRDS(snps, file.path(output_dir, paste0(output_stem, "_harmonized.rds")))

message("Running ", primary$factor_name, " ~ SNP for ", nrow(snps), " SNPs using model ", primary$label)
if (mode == "full") {
  measurement_path <- project_path(file.path(".work", "models", primary$label, "model_fit.rds"))
  measurement_fit <- readRDS(measurement_path)$results
  fit <- userGWAS(
    covstruc = covariance_structure,
    SNPs = snps,
    model = model,
    analytic = TRUE,
    usermod = measurement_fit,
    batch_size = 100000
  )
  result <- fit
} else {
  fit <- userGWAS(
    covstruc = covariance_structure,
    SNPs = snps,
    estimation = "DWLS",
    model = model,
    printwarn = TRUE,
    sub = gsub(" ", "", primary$snp_regression, fixed = TRUE),
    cores = requested_cores(),
    parallel = FALSE,
    GC = "standard",
    std.lv = TRUE,
    fix_measurement = TRUE,
    Q_SNP = TRUE
  )
  result <- if (is.list(fit) && length(fit) == 1L && is.data.frame(fit[[1]])) fit[[1]] else fit
}
if (!is.data.frame(result)) stop("Unexpected userGWAS output structure")

if (mode == "full") {
  result_path <- file.path(output_dir, paste0(output_stem, ".tsv.gz"))
  connection <- gzfile(result_path, open = "wt")
  write.table(result, connection, sep = "\t", quote = FALSE, row.names = FALSE)
  close(connection)
  report_path <- project_path(file.path(results_directory(), paste0(settings$analysis_id, "_factor_gwas.tsv.gz")))
} else {
  result_path <- file.path(output_dir, paste0(output_stem, ".tsv"))
  write.table(result, result_path, sep = "\t", quote = FALSE, row.names = FALSE)
  report_path <- project_path(file.path(results_directory(), paste0(settings$analysis_id, "_factor_gwas_pilot.tsv")))
}
saveRDS(fit, file.path(output_dir, paste0(output_stem, ".rds")))
file.copy(result_path, report_path, overwrite = TRUE)

if (mode == "full") {
  p_column <- paste0("p_val_", primary$factor_name)
  beta_column <- paste0("beta_", primary$factor_name)
  se_column <- paste0("SE_", primary$factor_name)
  required <- c(p_column, beta_column, se_column, "Q_omnibus_pval")
  if (!all(required %in% names(result))) stop("Missing expected analytic output columns")
  top <- which.min(result[[p_column]])
  summary <- data.frame(
    metric = c(
      "variants", "genome_wide_significant_rows", "q_significant_rows", "top_snp",
      "top_chr_bp", "top_beta", "top_se", "top_p", "top_q_p"
    ),
    value = c(
      nrow(result), sum(result[[p_column]] < 5e-8, na.rm = TRUE),
      sum(result$Q_omnibus_pval < 5e-8, na.rm = TRUE), result$SNP[top],
      paste0(result$CHR[top], ":", result$BP[top]), result[[beta_column]][top],
      result[[se_column]][top], result[[p_column]][top], result$Q_omnibus_pval[top]
    )
  )
  write.table(
    summary, project_path(file.path(results_directory(), paste0(settings$analysis_id, "_factor_gwas_summary.tsv"))),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
message("Factor-GWAS results written to ", report_path)
