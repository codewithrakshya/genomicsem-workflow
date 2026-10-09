args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 2L, 3L)) {
  stop("Usage: 05_visualize_results.R <project-root> [output-directory] [workflow.tsv]")
}

root <- normalizePath(args[[1]], mustWork = TRUE)
Sys.setenv(GENOMICSEM_PROJECT_ROOT = root)
source(file.path(root, "workflow", "scripts", "common.R"))
output_dir <- if (length(args) >= 2L) args[[2]] else results_directory()
output_dir <- ensure_directory(project_path(output_dir))
settings <- read_workflow_settings(if (length(args) == 3L) args[[3]] else "config/workflow.tsv")
traits <- read_traits(settings$traits_config)
models <- read_models(settings$models_config)

baseline <- models$label[models$role == "baseline"]
primary <- models$label[models$role == "primary"]
if (length(baseline) != 1L) baseline <- models$label[1]

read_result <- function(filename) {
  path <- file.path(output_dir, filename)
  if (!file.exists(path)) stop("Missing result: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

correlations <- read_result("genetic_correlations.tsv")
parameters <- read_result("model_parameter_comparison.tsv")
fit <- read_result("model_fit_comparison.tsv")
trait_order <- traits$trait
trait_labels <- stats::setNames(traits$display_name, traits$trait)
n_traits <- length(trait_order)

correlation_matrix <- as.matrix(correlations[, trait_order, drop = FALSE])
storage.mode(correlation_matrix) <- "numeric"
rownames(correlation_matrix) <- correlations$trait
correlation_matrix <- correlation_matrix[trait_order, trait_order, drop = FALSE]

loading_rows <- parameters[parameters$op == "=~", ]
loading_rows <- loading_rows[match(trait_order, loading_rows$rhs), ]
residual_rows <- parameters[
  parameters$op == "~~" & parameters$lhs == parameters$rhs & parameters$lhs %in% trait_order,
]
residual_rows <- residual_rows[match(trait_order, residual_rows$lhs), ]
fit_primary <- fit[fit$model == primary, ]
if (nrow(fit_primary) != 1L) stop("Expected one primary-model fit row")

column <- function(prefix, label) {
  name <- paste0(prefix, "_", label)
  if (!name %in% names(parameters)) stop("Missing comparison column: ", name)
  name
}

blue <- "#0072B2"
orange <- "#D55E00"
neutral <- "#767676"
grid <- "#D9D9D9"

draw_errorbar <- function(x, y, se, color) {
  if (is.na(se)) return(invisible(NULL))
  lower <- x - 1.96 * se
  upper <- x + 1.96 * se
  segments(lower, y, upper, y, col = color, lwd = 1.5)
  segments(lower, y - 0.04, lower, y + 0.04, col = color, lwd = 1.5)
  segments(upper, y - 0.04, upper, y + 0.04, col = color, lwd = 1.5)
}

draw_figure <- function() {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1.25, 1.25))
  par(family = "sans", las = 1, xaxs = "i", yaxs = "i")

  par(mar = c(6.5, 6.5, 3.2, 1.2))
  plot.new()
  plot.window(xlim = c(0.5, n_traits + 0.5), ylim = c(n_traits + 0.5, 0.5), asp = 1)
  palette <- hcl.colors(201, "Blue-Red 3")
  for (i in seq_len(n_traits)) {
    for (j in seq_len(n_traits)) {
      value <- correlation_matrix[i, j]
      color_index <- max(1L, min(201L, round((value + 1) * 100) + 1L))
      rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = palette[color_index], border = "white")
      text(j, i, sprintf("%.2f", value), cex = min(1.05, 3 / sqrt(n_traits)), font = 2)
    }
  }
  axis(1, at = seq_len(n_traits), labels = unname(trait_labels[trait_order]), tick = FALSE, las = 2, cex.axis = 0.75)
  axis(2, at = seq_len(n_traits), labels = unname(trait_labels[trait_order]), tick = FALSE, cex.axis = 0.75)
  box(col = neutral)
  title("A  Genetic correlations", adj = 0, font.main = 2, cex.main = 1.05)

  y <- rev(seq_len(n_traits))
  baseline_loading <- loading_rows[[column("STD_Genotype", baseline)]]
  primary_loading <- loading_rows[[column("STD_Genotype", primary)]]
  baseline_se <- loading_rows[[column("STD_Genotype_SE", baseline)]]
  primary_se <- loading_rows[[column("STD_Genotype_SE", primary)]]
  loading_limit <- range(c(0, baseline_loading - 1.96 * baseline_se, baseline_loading + 1.96 * baseline_se,
                           primary_loading - 1.96 * primary_se, primary_loading + 1.96 * primary_se), na.rm = TRUE)
  par(mar = c(5.0, 7.0, 3.2, 1.2))
  plot(NA, xlim = loading_limit, ylim = c(0.5, n_traits + 0.5), xlab = "Standardized factor loading", ylab = "", axes = FALSE)
  abline(v = c(0, 1), col = c(neutral, grid), lty = c(1, 2))
  axis(1); axis(2, at = y, labels = unname(trait_labels[trait_order]), tick = FALSE); box(col = neutral)
  for (i in seq_along(y)) {
    segments(baseline_loading[i], y[i] + 0.10, primary_loading[i], y[i] - 0.10, col = grid, lwd = 2)
    draw_errorbar(baseline_loading[i], y[i] + 0.10, baseline_se[i], blue)
    draw_errorbar(primary_loading[i], y[i] - 0.10, primary_se[i], orange)
  }
  points(baseline_loading, y + 0.10, pch = 16, col = blue, cex = 1.1)
  points(primary_loading, y - 0.10, pch = 17, col = orange, cex = 1.1)
  legend("bottomright", legend = c(baseline, primary), col = c(blue, orange), pch = c(16, 17), bty = "n", cex = 0.82)
  title("B  Factor loadings", adj = 0, font.main = 2, cex.main = 1.05)

  baseline_residual <- residual_rows[[column("Unstand_Est", baseline)]]
  primary_residual <- residual_rows[[column("Unstand_Est", primary)]]
  x_limits <- range(c(baseline_residual, primary_residual, 0), finite = TRUE)
  padding <- max(diff(x_limits) * 0.15, 0.01)
  par(mar = c(5.0, 7.0, 3.2, 1.2))
  plot(NA, xlim = x_limits + c(-padding, padding), ylim = c(0.5, n_traits + 0.5), xlab = "Residual genetic variance", ylab = "", axes = FALSE)
  abline(v = 0, col = neutral, lty = 2)
  axis(1); axis(2, at = y, labels = unname(trait_labels[trait_order]), tick = FALSE); box(col = neutral)
  segments(baseline_residual, y + 0.10, primary_residual, y - 0.10, col = grid, lwd = 2)
  points(baseline_residual, y + 0.10, pch = 16, col = blue, cex = 1.1)
  points(primary_residual, y - 0.10, pch = 17, col = orange, cex = 1.1)
  title("C  Residual variances", adj = 0, font.main = 2, cex.main = 1.05)

  mtext(
    sprintf("Primary model (%s): chi-square(%s) = %s, p = %s; CFI = %s; SRMR = %s",
      primary, fit_primary$df, signif(fit_primary$chisq, 3), signif(fit_primary$p_chisq, 3),
      signif(fit_primary$CFI, 3), signif(fit_primary$SRMR, 3)),
    side = 1, outer = TRUE, line = -1.1, cex = 0.82, col = neutral
  )
}

height <- max(920, 720 + n_traits * 60)
png(file.path(output_dir, "genomicsem_summary.png"), width = 2400, height = height, res = 200, bg = "white")
par(oma = c(2.2, 0, 1.4, 0)); draw_figure(); dev.off()
pdf(file.path(output_dir, "genomicsem_summary.pdf"), width = 12, height = max(4.8, 3.6 + n_traits * 0.35), useDingbats = FALSE)
par(oma = c(2.2, 0, 1.4, 0)); draw_figure(); dev.off()
cat("Figures written to", output_dir, "\n")
