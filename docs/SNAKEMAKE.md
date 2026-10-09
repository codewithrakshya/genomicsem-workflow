# Snakemake interface

`Snakefile` and `run_pipeline.sh` are two interfaces to the same scripts,
configuration, and outputs. Snakemake does not reimplement the statistics; it
adds dependency tracking, selective reruns, parallel scheduling, and cluster
execution.

## Install Snakemake

The R/GenomicSEM environment remains project-local and is installed by the
workflow. Snakemake itself can be installed separately with the supplied
environment:

```bash
conda env create -f workflow/snakemake/environment.yml
conda activate genomicsem-snakemake
```

## Fresh local run

```bash
snakemake setup --cores 1
snakemake check --cores 1
snakemake --cores 1
```

The default target stops after the model report, matching:

```bash
./run_pipeline.sh all
```

Neither default target runs SNP associations. `gwas_pilot` uses the configured
pilot size (500 variants by default); review its output before running
`gwas`, which analyzes all eligible variants after harmonization and filters.

Review those results before running:

```bash
snakemake gwas_pilot --cores 1
snakemake gwas --cores 1
```

## Dry run

Inspect the planned DAG without executing it:

```bash
snakemake --dry-run --printshellcmds
```

Generate a DAG if Graphviz is available:

```bash
snakemake --dag | dot -Tpng > workflow_dag.png
```

## Alternate analysis configuration

Both interfaces use `config/analysis.yaml` by default. To select another config,
the Bash interface accepts `GENOMICSEM_SETTINGS`; Snakemake accepts:

```bash
snakemake --config analysis_config=config/my_analysis.yaml --cores 1
```

The selected file contains the analysis settings, traits, and model syntax.
Generated intermediate tables are written under `.work/config/`.

## HPC execution

Install the executor plugin appropriate to the cluster and use the scheduler's
normal resource policy. For a Slurm installation, for example:

```bash
snakemake --executor slurm --jobs 20
```

Package installation often requires internet access and should usually be run
once on a login or approved setup node before submitting compute jobs:

```bash
snakemake setup --cores 1
snakemake check --cores 1
```

The current rules call the same stage runner used interactively. Therefore a
given configuration uses identical preparation logic, GenomicSEM calls, model
syntax, and output paths under Bash and Snakemake.

## Targets

| Snakemake target | Bash equivalent |
|---|---|
| `setup` | `./run_pipeline.sh setup` |
| `check` | `./run_pipeline.sh check` |
| `validate` | `./run_pipeline.sh validate` |
| `prepare` | `./run_pipeline.sh prepare` |
| `munge` | `./run_pipeline.sh munge` |
| `ldsc` | `./run_pipeline.sh ldsc` |
| `models` | `./run_pipeline.sh models` |
| `report` or default | `./run_pipeline.sh report` or `all` |
| `gwas_pilot` | `./run_pipeline.sh gwas-pilot` |
| `gwas` | `./run_pipeline.sh gwas` |

## Offline HTML report

`results/YYYY-MM-DD/<analysis_id>_report.html` combines workflow steps (including which
steps call GenomicSEM), current configuration, QC tables, covariance and
correlation matrices, model syntax/fit/warnings, the embedded summary figure,
and available pilot/full GWAS results. Supporting tables and working objects are
retained under `.work/`.
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

Generated outputs use `results/YYYY-MM-DD/`; see [dated-results instructions](../docs/RESULTS_LAYOUT.md) for date selection and migration details.
