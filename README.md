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
- Python dependency: PyYAML, installed with `python -m pip install -r requirements.txt`.
- Sufficient memory/storage for your GWAS. Full-genome requirements depend on
  variant count and model size; no universal minimum is claimed.

For a project-local Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Configure your analysis

1. Obtain GWAS summary statistics with appropriate permission and place them in
   `data/raw/`. All traits must have compatible genome builds and ancestries.
2. Edit only `config/analysis.yaml` to define each trait, its source columns,
   sample size and assumptions; set the reference paths and result filters; and
   write the model syntax inline. Mark exactly one model as `primary` for the
   optional SNP stage. Cohort metadata and overlap notes also live in this file.
3. Download HapMap3, ancestry-matched LD scores/weights, and the SNP allele/MAF
   reference described in the [GenomicSEM tutorials](https://github.com/GenomicSEM/GenomicSEM/wiki),
   then update their paths in `config/analysis.yaml`. Example settings are
   fictional and target European resources.

The pipeline generates its internal TSV and model-text files under `.work/`;
users do not need to edit those generated files.

## Run

Activate the Python environment created above. Install the R/GenomicSEM
dependencies once per machine or project, then validate the configuration and
run the LDSC/model/report stages:

```bash
source .venv/bin/activate
./run_pipeline.sh setup
./run_pipeline.sh check
./run_pipeline.sh validate
./run_pipeline.sh all
```

`all` prepares and munges the summary statistics, runs LDSC, fits and compares
the models, and builds the report. It does **not** run SNP associations. Review
the model fit, residual variances, and assumptions before running either GWAS
stage.

The optional pilot runs up to `results.pilot_snps` variants (500 by default) to
check SNP harmonization and the primary SNP model:

```bash
./run_pipeline.sh gwas-pilot
```

After reviewing the pilot, run the full factor GWAS:

```bash
./run_pipeline.sh gwas
```

The full stage uses all eligible variants from the configured summary-statistic
files, after harmonization and the configured INFO/MAF filters; it is not limited
to the pilot count. The example config requires real GWAS inputs and references
before it can run. Normal analysis stages do not reinstall GenomicSEM.

Final report, figures, and optional factor-GWAS outputs are stored in
`results/YYYY-MM-DD/`. Supporting QC, covariance, and model-comparison tables are
kept under `.work/results/YYYY-MM-DD/` and embedded in the HTML report.
Rebuild HTML with `./run_pipeline.sh html-report`. Set `GENOMICSEM_RUN_DATE` to
keep multiple invocations in one dated folder. Same-date runs overwrite outputs;
intermediates in `.work/` are shared. Use separate copies for independent analyses.

Snakemake reads the same `config/analysis.yaml` file. To select a different
analysis file, use `--config analysis_config=config/other_analysis.yaml`.
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
