args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: compare_factor_scaling.R <standardized-pilot.tsv> <unit-loading-pilot.tsv>")
}

library(data.table)

# ============================================================
# OLD: std.lv = TRUE pilot
# ============================================================

old <- fread(args[[1]])

# See columns
cat("\nOLD columns:\n")
print(names(old))


# ============================================================
# NEW: Dashti-anchored std.lv = FALSE pilot
# ============================================================

new <- fread(args[[2]])

cat("\nNEW columns:\n")
print(names(new))


# ============================================================
# Standardize OLD column names
# ============================================================

# Original pilot userGWAS output uses:
# est, SE, Z_Estimate, Pval_Estimate

old_factor <- old[
  lhs == "SleepDuration" &
  op == "~" &
  rhs == "SNP"
]

old_factor <- old_factor[, .(
  SNP,
  beta_old = est,
  SE_old = SE,
  Z_old = Z_Estimate,
  P_old = Pval_Estimate
)]


# ============================================================
# Standardize NEW column names
# ============================================================

new_factor <- new[, .(
  SNP,
  beta_new = beta_SleepDuration,
  SE_new = SE_SleepDuration,
  Z_new = Z_SleepDuration,
  P_new = p_val_SleepDuration
)]


# ============================================================
# MERGE SAME SNPs
# ============================================================

comparison <- merge(
  old_factor,
  new_factor,
  by = "SNP"
)

cat(
  "\nSNPs compared:",
  nrow(comparison),
  "\n\n"
)


# ============================================================
# Calculate differences
# ============================================================

comparison[, beta_ratio := beta_new / beta_old]

comparison[, SE_ratio := SE_new / SE_old]

comparison[, Z_difference := Z_new - Z_old]

comparison[, P_difference := P_new - P_old]


# ============================================================
# SUMMARY
# ============================================================

cat("============================================================\n")
cat("       OLD vs NEW FACTOR SCALING COMPARISON\n")
cat("============================================================\n\n")

cat("OLD:\n")
cat("  std.lv = TRUE\n")
cat("  factor variance fixed to 1\n\n")

cat("NEW:\n")
cat("  std.lv = FALSE\n")
cat("  Dashti loading fixed to 1\n\n")


cat("Number of SNPs compared:\n")
print(nrow(comparison))


cat("\nBeta correlation:\n")
print(
  cor(
    comparison$beta_old,
    comparison$beta_new,
    use = "complete.obs"
  )
)


cat("\nSE correlation:\n")
print(
  cor(
    comparison$SE_old,
    comparison$SE_new,
    use = "complete.obs"
  )
)


cat("\nZ correlation:\n")
print(
  cor(
    comparison$Z_old,
    comparison$Z_new,
    use = "complete.obs"
  )
)


cat("\n-log10(P) correlation:\n")

valid_p <- (
  comparison$P_old > 0 &
  comparison$P_new > 0
)

print(
  cor(
    -log10(comparison$P_old[valid_p]),
    -log10(comparison$P_new[valid_p]),
    use = "complete.obs"
  )
)


cat("\nMedian beta scaling ratio (NEW / OLD):\n")
print(
  median(
    comparison$beta_ratio[
      is.finite(comparison$beta_ratio)
    ],
    na.rm = TRUE
  )
)


cat("\nMedian SE scaling ratio (NEW / OLD):\n")
print(
  median(
    comparison$SE_ratio[
      is.finite(comparison$SE_ratio)
    ],
    na.rm = TRUE
  )
)


cat("\nMaximum absolute Z difference:\n")
print(
  max(
    abs(comparison$Z_difference),
    na.rm = TRUE
  )
)


cat("\nMaximum absolute P difference:\n")
print(
  max(
    abs(comparison$P_difference),
    na.rm = TRUE
  )
)


# ============================================================
# Show first 10 SNPs
# ============================================================

cat("\n============================================================\n")
cat("First 10 SNP comparisons\n")
cat("============================================================\n\n")

print(
  comparison[
    1:min(10, .N),
    .(
      SNP,
      beta_old,
      beta_new,
      SE_old,
      SE_new,
      Z_old,
      Z_new,
      P_old,
      P_new
    )
  ]
)


# ============================================================
# SAVE
# ============================================================

fwrite(
  comparison,
  "SleepDuration_scaling_comparison.tsv",
  sep = "\t"
)

cat(
  "\nComparison saved as:\n",
  "SleepDuration_scaling_comparison.tsv\n"
)