args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 10_plot_manhattan.R <gwas.tsv.gz> <output.png> <factor-name>")
}

source(file.path(Sys.getenv("GENOMICSEM_PROJECT_ROOT", unset = getwd()), "workflow", "scripts", "common.R"))
suppressPackageStartupMessages(library(data.table))

gwas_path <- project_path(args[[1]])
output_path <- project_path(args[[2]])
factor_name <- args[[3]]
p_column <- paste0("p_val_", factor_name)
if (!file.exists(gwas_path)) stop("Full GWAS result does not exist: ", gwas_path)

gwas <- data.table::fread(
  gwas_path,
  select = c("CHR", "BP", p_column),
  showProgress = FALSE,
  nThread = 1L
)
chromosome <- suppressWarnings(as.integer(sub("^chr", "", gwas$CHR, ignore.case = TRUE)))
position <- suppressWarnings(as.numeric(gwas$BP))
p_value <- suppressWarnings(as.numeric(gwas[[p_column]]))
valid <- chromosome %in% 1:22 & is.finite(position) & position > 0 &
  is.finite(p_value) & p_value >= 0 & p_value <= 1
if (!any(valid)) stop("No valid autosomal SNPs found for the Manhattan plot")

chromosome <- chromosome[valid]
position <- position[valid]
p_value <- p_value[valid]
minus_log10_p <- -log10(pmax(p_value, .Machine$double.xmin))
chromosomes <- sort(unique(chromosome))
chromosome_lengths <- vapply(chromosomes, function(value) max(position[chromosome == value]), numeric(1))
gap <- max(chromosome_lengths) * 0.01
offsets <- c(0, cumsum(head(chromosome_lengths + gap, -1L)))
x <- position + offsets[match(chromosome, chromosomes)]
x_ticks <- offsets + chromosome_lengths / 2
x_limit <- max(offsets + chromosome_lengths)
y_limit <- max(8, min(50, ceiling(max(minus_log10_p))))
colors <- c("#0072B2", "#D55E00")[(match(chromosome, chromosomes) - 1L) %% 2L + 1L]

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
png(output_path, width = 2600, height = 1100, res = 200)
par(mar = c(5, 5, 3.5, 1.5), family = "sans", las = 1)
plot(
  x, pmin(minus_log10_p, y_limit),
  pch = ".", col = colors, xlim = c(0, x_limit), ylim = c(0, y_limit),
  xaxt = "n", xlab = "Chromosome", ylab = expression(-log[10](P)),
  main = paste(factor_name, "factor GWAS"), bty = "l"
)
axis(1, at = x_ticks, labels = chromosomes, tick = FALSE)
abline(h = -log10(5e-8), col = "#222222", lty = 2, lwd = 1.5)
legend(
  "topright", legend = "P = 5 x 10^-8", lty = 2, lwd = 1.5,
  col = "#222222", bty = "n", cex = 0.9
)
if (max(minus_log10_p) > y_limit) {
  mtext("Values above 50 are capped", side = 3, adj = 1, line = 0.3, cex = 0.8)
}
dev.off()
cat("Manhattan plot written to", output_path, "\n")