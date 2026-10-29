# Sleep-duration GenomicSEM analysis specification

This document defines the planned analysis. It deliberately contains no
numerical findings from an earlier run. After a fresh execution, the files in
`results/` are the only active results.

## Scientific objective

Determine whether the configured sleep GWASs share a genetically coherent
latent sleep-duration factor and, if the measurement model is acceptable,
estimate SNP associations with that shared factor.

This is a multivariate summary-statistic analysis. It is not an individual-
level analysis and it is not a conventional fixed-effect meta-analysis.

## Configured traits

The authoritative trait definitions are in `config/traits.tsv`.

- Austin-Zimmerman short sleep is oriented toward longer sleep by multiplying
  its effect estimate by -1.
- Dashti sleep duration retains its reported effect direction.
- Portas accelerometer-derived sleep duration retains its reported direction.
- Trait order in `config/traits.tsv` is preserved through munging, LDSC,
  `sumstats()`, and `userGWAS()`.
- All configured files must use the same genome build and compatible ancestry
  reference resources.

The binary-trait prevalence values are explicit modeling assumptions and must
be scientifically confirmed or sensitivity-tested before publication.

## Reference resources

- `reference/w_hm3.snplist`: HapMap3 SNP filter for LDSC preparation.
- `reference/eur_w_ld_chr/`: European LD scores and regression weights for
  chromosomes 1–22.
- `reference/reference.1000G.maf.0.005.txt.gz`: GenomicSEM 1000 Genomes Phase 3
  SNP/allele/MAF table for SNP-level harmonization.

The 1000 Genomes table comes from the public UT Austin resource linked by the
GenomicSEM Common Factor GWAS tutorial:
https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v

## Models

The authoritative model selection is in `config/models.tsv`.

The baseline model is:

```text
SleepDuration =~ Austin_short_reversed + Dashti_sleep_duration + Portas_sleep_duration
```

The configured primary sensitivity model is:

```text
SleepDuration =~ Austin_short_reversed + Dashti_sleep_duration + Portas_sleep_duration
Dashti_sleep_duration ~~ 0*Dashti_sleep_duration
```

The constraint is carried into the SNP model by appending:

```text
SleepDuration ~ SNP
```

The constraint must be justified from the freshly fitted model results; it
must not be accepted solely because it removes an inadmissible estimate.

## Execution sequence

The commands below use the simple Bash interface. The same stages and outputs
are available through the optional Snakemake interface documented in
`docs/SNAKEMAKE.md`.

### 1. Bootstrap and validation

```bash
./run_pipeline.sh setup
./run_pipeline.sh check
./run_pipeline.sh validate
```

This installs the isolated R environment, records package versions, checks the
configuration, confirms model files, confirms reference paths, and verifies
that all traits use one genome build.

### 2. Preparation through model report

```bash
./run_pipeline.sh all
```

This performs:

1. configuration-driven column standardization and effect orientation;
2. HapMap3 filtering;
3. GenomicSEM `munge()`;
4. multivariable LDSC to estimate genetic covariance S and sampling covariance
   V, including covariance induced by cohort overlap;
5. every model in `config/models.tsv`;
6. generic model comparison and visualization.

Do not continue automatically to the full GWAS. First review:

- preparation row counts and allele losses;
- LDSC intercepts and heritability estimates;
- genetic correlations and their uncertainty;
- negative residual variances or loadings outside an admissible range;
- model fit and warnings;
- whether the primary factor is scientifically interpretable.

### 3. SNP-model pilot

```bash
./run_pipeline.sh gwas-pilot
```

The pilot uses a deterministic SNP subset and iterative DWLS. Confirm that the
configured primary model runs without systematic errors and that effect
direction, columns, and Q-SNP output are correct.

### 4. Full factor GWAS

```bash
./run_pipeline.sh gwas
```

The full stage prepares genome-wide inputs, harmonizes them with the configured
1000 Genomes allele/MAF reference using GenomicSEM `sumstats()`, and runs the
primary first-order factor model with analytic `userGWAS()`.

The full result contains:

- the SNP effect on the configured latent factor;
- its standard error, Z statistic, and P value;
- the Q-SNP omnibus heterogeneity statistic and P value.

## Expected active outputs

After `all`:

- `results/YYYY-MM-DD/preparation_ldsc_qc.tsv`
- `results/YYYY-MM-DD/genetic_correlations.tsv`
- `results/YYYY-MM-DD/genetic_covariances.tsv`
- `results/YYYY-MM-DD/model_fit_comparison.tsv`
- `results/YYYY-MM-DD/model_parameter_comparison.tsv`
- `results/YYYY-MM-DD/model_comparison_summary.txt`
- `results/YYYY-MM-DD/genomicsem_summary.png`
- `results/YYYY-MM-DD/genomicsem_summary.pdf`

After `gwas-pilot` and `gwas`, filenames are prefixed by the `analysis_id` in
`config/workflow.tsv`, including:

- the pilot SNP results;
- the compressed full factor-GWAS summary statistics;
- the compact factor-GWAS summary.

Technical intermediates and logs are written under `.work/`.

## Interpretation requirements

- A positive genetic correlation supports shared genetic architecture but does
  not establish that phenotypes are interchangeable.
- Approximate fit indices and chi-square must be considered together.
- A boundary constraint must be reported explicitly.
- Genome-wide significant rows are not independent loci; LD clumping or
  conditional analysis is required.
- Significant Q-SNP results indicate effects not fully explained by the common
  factor and require trait-specific inspection.
- Biological claims require locus annotation, comparison with the original
  GWASs, and ideally replication.

## Reuse

For a different set of GWASs, edit configuration and model files rather than
workflow code. See `docs/REUSE.md`.

Generated outputs use `results/YYYY-MM-DD/`; see [dated-results instructions](../docs/RESULTS_LAYOUT.md) for date selection and migration details.
