args <- commandArgs(trailingOnly = TRUE)
settings_path <- if (length(args)) args[[1]] else ".work/config/workflow.tsv"
source(file.path(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), "workflow", "scripts", "common.R"))

settings <- read_workflow_settings(settings_path)
required_settings <- c(
  "analysis_id", "traits_config", "models_config", "hm3_reference",
  "ld_reference", "weight_reference", "gwas_reference", "pilot_snps",
  "info_filter", "maf_filter"
)
missing_settings <- setdiff(required_settings, names(settings))
if (length(missing_settings)) stop("Missing workflow settings: ", paste(missing_settings, collapse = ", "))

traits <- read_traits(settings$traits_config)
required_traits <- c(
  "display_name", "raw_file", "full_file", "chr_col", "bp_col", "snp_col",
  "a1_col", "a2_col", "af_col", "beta_col", "se_col", "p_col", "lp_col",
  "n_col", "effect_multiplier", "se_logit", "OLS"
)
missing_traits <- setdiff(required_traits, names(traits))
if (length(missing_traits)) stop("Missing trait columns: ", paste(missing_traits, collapse = ", "))
invisible(logical_column(traits$se_logit, "se_logit"))
invisible(logical_column(traits$OLS, "OLS"))
if (any(!nzchar(traits$p_col) & !nzchar(traits$lp_col))) {
  stop("Each trait needs either p_col or lp_col")
}
if (length(unique(traits$build)) != 1L) stop("All traits must use the same genome build")

models <- read_models(settings$models_config)
missing_models <- models$model_file[!file.exists(models$model_file)]
if (length(missing_models)) stop("Missing model files: ", paste(missing_models, collapse = ", "))

references <- c(settings$hm3_reference, settings$ld_reference, settings$weight_reference, settings$gwas_reference)
missing_references <- references[!file.exists(vapply(references, project_path, character(1))) &
                                   !dir.exists(vapply(references, project_path, character(1)))]
if (length(missing_references)) stop("Missing references: ", paste(missing_references, collapse = ", "))

cat("Configuration valid\n")
cat("Analysis:", settings$analysis_id, "\n")
cat("Traits:", paste(traits$trait, collapse = ", "), "\n")
cat("Models:", paste(models$label, collapse = ", "), "\n")
cat("Primary model:", models$label[models$role == "primary"], "\n")
