# GenomicSEM workflow

The included example analysis prepares three sleep GWASs, reverses only the
Austin-Zimmerman short-sleep effects, runs multivariable LDSC, and compares an
unconstrained sleep-duration factor with a sensitivity model that fixes the
Dashti residual variance to zero. It can then run a constrained SNP-level
factor GWAS across all variants that pass multivariate harmonization.

For the scientific question, complete methods, results, limitations, and file
guide, read **[`docs/ANALYSIS.md`](docs/ANALYSIS.md)** first.

For adapting the workflow to different GWASs or models, read
**[`docs/REUSE.md`](docs/REUSE.md)**. Trait names, source columns, effect
directions, sample sizes, references, and SEM models are configuration-driven;
the workflow scripts contain no sleep-study names.

The project has two equivalent entry points:

- `run_pipeline.sh` for a simple sequential run;
- `Snakefile` for dependency-aware local or cluster execution.

Both call the same analysis scripts, read the same configuration, and create
the same outputs. See **[`docs/SNAKEMAKE.md`](docs/SNAKEMAKE.md)** for the
Snakemake commands.

## Clean project layout

```text
run_pipeline.sh    One command for the analysis
Snakefile          Optional Snakemake interface to the same analysis
config/            Trait, cohort, and model definitions
data/raw/          Original downloaded GWAS files
reference/         HapMap3 and European LDSC reference files
workflow/          Scripts, R library, and scheduler details
results/           Small, report-ready outputs only
.work/             Hidden intermediate files and logs
```

The original downloads in `data/raw/` are never modified. Everything large or
technical is kept under `.work/`; normal directory listings therefore show
only the files needed to understand or present the analysis.

## Run the analysis

For a fresh deployment, from the project root:

```bash
./run_pipeline.sh setup
./run_pipeline.sh check
./run_pipeline.sh validate
./run_pipeline.sh all
./run_pipeline.sh gwas-pilot
./run_pipeline.sh gwas
```

`setup` creates the project-local R library and installs GenomicSEM. `check`
records package versions and verifies the runtime. `all` intentionally stops
after LDSC, model fitting, and the model report; review those outputs before
starting either SNP-level GWAS stage.

Final report, figures (including a Manhattan plot after the full factor GWAS),
and optional factor-GWAS outputs are stored in
`results/YYYY-MM-DD/`. Supporting QC, covariance, and model-comparison tables are
kept under `.work/results/YYYY-MM-DD/` and embedded in the HTML report.
Rebuild HTML with `./run_pipeline.sh html-report`. Set `GENOMICSEM_RUN_DATE` to
keep multiple invocations in one dated folder. Same-date runs overwrite outputs;
intermediates in `.work/` are shared. Use separate copies for independent analyses.

Snakemake reads the same `config/analysis.yaml` file. To select a different
analysis file, use `--config analysis_config=config/other_analysis.yaml`.
Snakemake provides the same stages:

To rerun only one stage:

```bash
./run_pipeline.sh prepare
./run_pipeline.sh munge
./run_pipeline.sh ldsc
./run_pipeline.sh models
./run_pipeline.sh report
./run_pipeline.sh gwas-pilot
./run_pipeline.sh gwas
```

The main report-ready files are written under `results/YYYY-MM-DD/`:

```text
preparation_ldsc_qc.tsv
genetic_correlations.tsv
genetic_covariances.tsv
model_fit_comparison.tsv
model_parameter_comparison.tsv
model_comparison_summary.txt
cohort_overlap.md
genomicsem_summary.png
genomicsem_summary.pdf
sleep_duration_factor_gwas_summary.tsv
sleep_duration_factor_gwas.tsv.gz
```

The deployment copy is intentionally shipped without generated results or
intermediates. A successful fresh run repopulates `results/` and `.work/`.

## Install or verify GenomicSEM

The isolated R library is stored under `workflow/R/library/`:

```bash
workflow/bin/run-r workflow/scripts/install_genomicsem.R
workflow/bin/run-r workflow/scripts/check_environment.R
workflow/bin/run-r workflow/scripts/tutorial_smoke_test.R
```

Installation and smoke-test details are technical outputs and are stored under
`.work/setup/`.

On Wynton, load R before running the commands:

```bash
module load CBI r
```

## Configuration-driven design

`config/workflow.tsv` selects resources and QC settings. `config/traits.tsv`
maps arbitrary GWAS columns and trait-level assumptions. `config/models.tsv`
selects any number of SEM models and identifies the primary SNP model.

`run_pipeline.sh` and `Snakefile` are thin orchestration layers. Preparation,
munging, LDSC, model fitting, comparison, reporting, and factor GWAS remain
separate scripts with explicit inputs and outputs. The model table controls
standardized (`std_lv=TRUE`) versus unit-loading (`std_lv=FALSE`) identification,
and a unit-loading primary model can name its reference indicator. No analysis-
specific GWAS script is required.

Use the Bash interface for a straightforward run. Use Snakemake when dependency
tracking, selective reruns, or scheduler execution is useful on local, Wynton,
or coreHPC systems.

## Data governance

Confirm that each computing environment is approved for the data being used.
Public GWAS summary statistics may be portable, while controlled individual-
level data or restricted summary statistics may require a specifically approved
system.

## Offline HTML report

`results/YYYY-MM-DD/<analysis_id>_report.html` combines workflow steps (including which
steps call GenomicSEM), current configuration, QC tables, covariance and
correlation matrices, model syntax/fit/warnings, the embedded summary figure,
and available pilot/full GWAS results. Execution records are retained separately under `.work/runs/`.
It is self-contained and can be opened or shared without an internet connection.
The full compressed GWAS is kept separate.

`./run_pipeline.sh report` regenerates figures and HTML. To build only the HTML
from existing outputs, run `./run_pipeline.sh html-report`. Both `gwas-pilot` and
`gwas` refresh the HTML after completing. No additional packages are required.
The default Snakemake target includes HTML; `snakemake html_report --cores 1`
builds it explicitly. After manually replacing optional outputs, use
`./run_pipeline.sh html-report` to refresh the report.

Missing optional results are marked as unavailable. File presence is not a
run-completion check: existing outputs can come from earlier runs, and the
report explicitly distinguishes available artifacts from verified execution.

For execution manifests, the pinned GenomicSEM commit, and repeat-run instructions,
see [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

Generated outputs use `results/YYYY-MM-DD/`; see [dated-results instructions](docs/RESULTS_LAYOUT.md) for date selection and migration details.
