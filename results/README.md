# Generated results

This directory is intentionally clean in the deployment copy. Generated
outputs appear here only after a fresh run, so files from an older analysis
cannot be mistaken for current results.

Recommended sequence:

```bash
./run_pipeline.sh setup
./run_pipeline.sh check
./run_pipeline.sh validate
./run_pipeline.sh all
```

Review the LDSC and model outputs before continuing:

```bash
./run_pipeline.sh gwas-pilot
./run_pipeline.sh gwas
```

Expected small/report-ready outputs include genetic correlations and
covariances, model comparisons, preparation QC, the summary figure, and the
factor-GWAS summary. The complete factor-GWAS summary statistics are compressed
and ignored by Git.

`cohort_overlap.md` is input documentation rather than a generated result and
is intentionally retained between runs.

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

Generated outputs use `results/YYYY-MM-DD/`; see [dated-results instructions](../docs/RESULTS_LAYOUT.md) for date selection and migration details.
