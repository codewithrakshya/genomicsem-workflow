# ============================================================
# Calculate N-hat for GenomicSEM Sleep Duration Factor GWAS
# ============================================================

# Install data.table if needed
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

library(data.table)


# ------------------------------------------------------------
# 1. FILE PATH
# ------------------------------------------------------------

# CHANGE THIS to the location of your factor GWAS file.
# Example:
# gwas_file <- "/path/to/sleep_duration_factor_gwas.tsv.gz"

gwas_file <- "sleep_duration_factor_gwas.tsv.gz"


# ------------------------------------------------------------
# 2. READ THE GENOMICSEM FACTOR GWAS
# ------------------------------------------------------------

cat("\nReading factor GWAS...\n")

gwas <- fread(
  gwas_file,
  showProgress = TRUE
)

cat("GWAS loaded successfully.\n")
cat("Total SNP rows:", format(nrow(gwas), big.mark = ","), "\n\n")


# ------------------------------------------------------------
# 3. CHECK REQUIRED COLUMNS
# ------------------------------------------------------------

required_columns <- c(
  "MAF",
  "SE_SleepDuration"
)

missing_columns <- setdiff(required_columns, names(gwas))

if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    "\n\nAvailable columns are:\n",
    paste(names(gwas), collapse = ", ")
  )
}

cat("Required columns found:\n")
cat("  MAF\n")
cat("  SE_SleepDuration\n\n")


# ------------------------------------------------------------
# 4. CONVERT TO NUMERIC
# ------------------------------------------------------------

gwas[, MAF := as.numeric(MAF)]
gwas[, SE_SleepDuration := as.numeric(SE_SleepDuration)]


# ------------------------------------------------------------
# 5. APPLY MAF FILTER
#
# Keep SNPs with:
#
#       0.10 < MAF < 0.40
#
# and valid positive standard errors.
# ------------------------------------------------------------

n_total <- nrow(gwas)

gwas_nhat <- gwas[
  !is.na(MAF) &
  !is.na(SE_SleepDuration) &
  MAF > 0.10 &
  MAF < 0.40 &
  SE_SleepDuration > 0
]

n_filtered <- nrow(gwas_nhat)

cat("MAF/SE filtering\n")
cat("----------------\n")
cat(
  "Original SNPs: ",
  format(n_total, big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "SNPs retained (0.10 < MAF < 0.40): ",
  format(n_filtered, big.mark = ","),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. CALCULATE SNP-SPECIFIC N-HAT
#
# Formula:
#
#               1
# N_hat = -----------------
#         2*p*(1-p)*SE^2
#
# where:
# p  = MAF
# SE = factor GWAS standard error
# ------------------------------------------------------------

gwas_nhat[, N_hat :=
  1 / (
    2 *
    MAF *
    (1 - MAF) *
    SE_SleepDuration^2
  )
]


# ------------------------------------------------------------
# 7. REMOVE INVALID N-HAT VALUES
# ------------------------------------------------------------

gwas_nhat <- gwas_nhat[
  is.finite(N_hat) &
  N_hat > 0
]

cat(
  "SNPs with valid N-hat: ",
  format(nrow(gwas_nhat), big.mark = ","),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. CALCULATE FACTOR-LEVEL N-HAT
#
# IMPORTANT:
# Use the MEAN of SNP-specific N-hat values.
# ------------------------------------------------------------

factor_N_hat <- mean(gwas_nhat$N_hat)


# ------------------------------------------------------------
# 9. QC / DISTRIBUTION
# ------------------------------------------------------------

median_N_hat <- median(gwas_nhat$N_hat)

sd_N_hat <- sd(gwas_nhat$N_hat)

N_hat_quantiles <- quantile(
  gwas_nhat$N_hat,
  probs = c(
    0.01,
    0.05,
    0.25,
    0.50,
    0.75,
    0.95,
    0.99
  ),
  na.rm = TRUE
)


# ------------------------------------------------------------
# 10. PRINT RESULTS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("        GenomicSEM Sleep Duration Factor N-hat\n")
cat("============================================================\n\n")

cat(
  "Number of SNPs used: ",
  format(nrow(gwas_nhat), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "Mean N-hat:           ",
  format(round(factor_N_hat), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "Median N-hat:         ",
  format(round(median_N_hat), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "SD of N-hat:          ",
  format(round(sd_N_hat), big.mark = ","),
  "\n\n",
  sep = ""
)

cat("N-hat quantiles:\n")

print(
  round(N_hat_quantiles)
)

cat("\n")

cat(
  "FINAL FACTOR N-HAT = ",
  format(round(factor_N_hat), big.mark = ","),
  "\n",
  sep = ""
)

cat("\n============================================================\n")


# ------------------------------------------------------------
# 11. SAVE SNP-LEVEL N-HAT RESULTS
# ------------------------------------------------------------

fwrite(
  gwas_nhat,
  file = "SleepDuration_factor_Nhat_SNPs.tsv.gz",
  sep = "\t"
)


# ------------------------------------------------------------
# 12. SAVE SUMMARY
# ------------------------------------------------------------

summary_table <- data.table(
  metric = c(
    "total_GWAS_SNPs",
    "SNPs_used_for_Nhat",
    "MAF_lower_bound",
    "MAF_upper_bound",
    "mean_Nhat",
    "median_Nhat",
    "SD_Nhat",
    "Nhat_1_percentile",
    "Nhat_5_percentile",
    "Nhat_25_percentile",
    "Nhat_50_percentile",
    "Nhat_75_percentile",
    "Nhat_95_percentile",
    "Nhat_99_percentile"
  ),

  value = c(
    n_total,
    nrow(gwas_nhat),
    0.10,
    0.40,
    factor_N_hat,
    median_N_hat,
    sd_N_hat,
    N_hat_quantiles[1],
    N_hat_quantiles[2],
    N_hat_quantiles[3],
    N_hat_quantiles[4],
    N_hat_quantiles[5],
    N_hat_quantiles[6],
    N_hat_quantiles[7]
  )
)

fwrite(
  summary_table,
  file = "SleepDuration_factor_Nhat_summary.tsv",
  sep = "\t"
)


# ------------------------------------------------------------
# 13. FINISHED
# ------------------------------------------------------------

cat("\nFiles saved:\n")
cat("  SleepDuration_factor_Nhat_summary.tsv\n")
cat("  SleepDuration_factor_Nhat_SNPs.tsv.gz\n\n")

cat("Done.\n")