# Configurable GenomicSEM workflow

A summary-statistic workflow for multivariable LD score regression, genetic
factor modeling, model comparison, and optional factor GWAS. Results include
a self-contained HTML scientific report with figures, methods, interpretation,
and references.

**This repository contains fictional configuration examples only. No GWAS data,
reference panels, fitted study models, results, installed R libraries, or run
history are distributed.** Example sample sizes are placeholders, not estimates.

## Requirements

- Bash, Python 3.9+ with IANA timezone data, and R/Rscript on PATH.
- An R build and compiler/system libraries compatible with GenomicSEM dependencies.
- Internet access for initial R package installation and reference downloads.
- Optional Snakemake: use `workflow/snakemake/environment.yml`.
- Sufficient memory/storage for your GWAS. Full-genome requirements depend on
  variant count and model size; no universal minimum is claimed.

## Configure your analysis

1. Obtain GWAS summary statistics with appropriate permission and place them in
   `data/raw/`. All traits must have compatible genome builds and ancestries.
2. Edit `config/traits.tsv`: source column names, effect direction, sample sizes,
   trait types, and binary-trait prevalence assumptions. Examples are fictional.
3. Download HapMap3, ancestry-matched LD scores/weights, and the SNP allele/MAF
   reference described in the [GenomicSEM tutorials](https://github.com/GenomicSEM/GenomicSEM/wiki).
   Set their paths in `config/workflow.tsv`. The supplied paths target European
   resources and must be changed for other appropriate resources.
4. Define your model in `config/model_common_factor.txt` and configure
   `config/models.tsv`. Exactly one model is designated primary for the SNP stage.
5. Fill in cohort metadata and overlap notes. For private configuration, use
   ignored `config/local/` and set `GENOMICSEM_SETTINGS=config/local/workflow.tsv`.

## Run

```bash
./run_pipeline.sh setup
./run_pipeline.sh check
./run_pipeline.sh validate
./run_pipeline.sh all
# Review the model before running SNP associations:
./run_pipeline.sh gwas-pilot
./run_pipeline.sh gwas
```

The example cannot run an analysis until actual inputs/references are supplied.
`all` stops after model fitting and reporting. It does not install packages or run
the full GWAS. Normal analysis runs do not reinstall GenomicSEM.

Outputs: `results/YYYY-MM-DD/<analysis_id>_report.html`, QC tables, genetic
covariances/correlations, fit comparisons, figures, and optional GWAS outputs.
Rebuild HTML with `./run_pipeline.sh html-report`. Set `GENOMICSEM_RUN_DATE` to
keep multiple invocations in one dated folder. Same-date runs overwrite outputs;
intermediates in `.work/` are shared. Use separate copies for independent analyses.

Snakemake provides the same stages:

```bash
snakemake setup --cores 1
snakemake check --cores 1
snakemake --cores 1
```

## Scientific limits

The example three-indicator one-factor model has zero degrees of freedom and
cannot provide an overall model-fit test. It is a syntax example, not a validated
measurement model. No residual variance is automatically fixed to zero. Any
boundary constraint must be justified for the specific analysis. With one model,
the comparison figure uses that model as both baseline and primary.

The full SNP stage uses analytic userGWAS for first-order measurement models;
more complex models may need a different estimation strategy. Significant SNP
rows are not independent loci. Heterogeneity and model assumptions require review.

The preparation stage does not retain INFO scores, so its downstream INFO
threshold does not establish filtering on original imputation quality. Effective
sample sizes and prevalence assumptions require study-specific justification.

## Documentation

- [Configuration and reuse](docs/REUSE.md)
- [Snakemake](docs/SNAKEMAKE.md)
- [Dated results](docs/RESULTS_LAYOUT.md)
- [Reproducibility](docs/REPRODUCIBILITY.md)
- [Publication checklist](PUBLISHING.md)

## References

Grotzinger et al. (2019), *Genomic structural equation modelling provides insights
into the multivariate genetic architecture of complex traits*.
[doi:10.1038/s41562-019-0566-x](https://doi.org/10.1038/s41562-019-0566-x).
GenomicSEM remains a separate dependency with its own license and citation requirements.

## License

This workflow is available under the [MIT License](LICENSE). Dependencies retain
their own licenses and citation requirements.
