# Reusing the GenomicSEM workflow

The workflow code is study-agnostic. A new analysis should change one file,
`config/analysis.yaml`, not scripts or generated workflow inputs.

## The analysis config

`config/analysis.yaml` is the only user-edited configuration. It contains:

- `analysis_id`: short output prefix;
- `references`: HapMap3, LD-score, regression-weight, and GWAS allele/MAF paths;
- `results`: pilot size and INFO/MAF filters;
- `traits`: one entry per GWAS with its file, sample details, and source columns;
- `models`: one or more inline model definitions and the selected primary model;
- `cohort_metadata` and `cohort_overlap`: report context for study overlap.

Each trait's `columns` maps the source summary-statistic column names (including
P or -log10(P)); `effect_multiplier` can reverse effect direction. Use
`sample_prev` and `population_prev` for binary-trait prevalence assumptions.
Values in the example are fictional.

Each model has a `label`, `role`, `factor_name`, `snp_regression`, and multiline
`syntax`. Exactly one model must have `role: primary`; any number of other models
can be compared. The primary model is used for the SNP-level GWAS.

To run a different analysis file:

```bash
GENOMICSEM_SETTINGS=config/my_analysis.yaml ./run_pipeline.sh validate
GENOMICSEM_SETTINGS=config/my_analysis.yaml ./run_pipeline.sh all
```

The equivalent Snakemake command is:

```bash
snakemake --config analysis_config=config/my_analysis.yaml --cores 1
```

## Starting a different analysis

1. Copy the repository without `data/raw/`, `reference/`, `.work/`, or large
   `results/` files.
2. Place or link the new GWAS files under `data/raw/`.
3. Edit `config/analysis.yaml`: add traits and their source columns, reference
   paths, cohort details, and inline model syntax; mark one justified model
   `primary`.
4. Validate before running anything:

```bash
./run_pipeline.sh validate
```

5. Run through the model stage:

```bash
./run_pipeline.sh all
```

6. Review genetic correlations, model fit, residual variances, and cohort
   overlap. Do not run a full factor GWAS just because a model converges.
7. Test the SNP model on a small subset:

```bash
./run_pipeline.sh gwas-pilot
```

The pilot size is controlled by `results.pilot_snps` in `config/analysis.yaml`
(500 by default). It checks harmonization and execution; it is not the full GWAS.

8. Run the full factor GWAS after reviewing the pilot and model assumptions:

```bash
./run_pipeline.sh gwas
```

This stage uses all eligible variants after harmonization and configured filters.
It can take substantially longer and use more memory/storage than the pilot.
For large jobs, submit this command through an approved scheduler.

## Stable workflow stages

| Stage | Purpose | Driven by |
|---|---|---|
| `validate` | Check configuration, models, build consistency, and references | analysis config |
| `prepare` | Standardize columns, direction, P values, and HapMap3 subset | analysis config |
| `munge` | GenomicSEM/LDSC munging | analysis config |
| `ldsc` | Genetic and sampling covariance | analysis config |
| `models` | Fit every configured model and compare results | analysis config |
| `report` | Generic correlation/loading/residual figure | analysis config |
| `gwas-pilot` | Small iterative `userGWAS()` validation | analysis config |
| `gwas` | Genome-wide preparation and analytic factor GWAS | analysis config |

The Bash and Snakemake interfaces invoke these same stages. They are not two
separate implementations; see `docs/SNAKEMAKE.md` for target names and cluster
usage.

## What remains analysis-specific

Modularity cannot replace scientific decisions. Each new project must justify:

- phenotype direction and comparability;
- ancestry and genome-build compatibility;
- binary-trait population prevalence;
- sample overlap and cohort composition;
- SEM model structure and constraints;
- whether the factor is interpretable;
- appropriate follow-up for Q-SNP heterogeneity.

Those choices belong in the analysis config and analysis-specific documentation.
The software mechanics should not encode them.

The full-GWAS stage currently uses GenomicSEM's fast analytic `userGWAS()`
estimator, which is intended for first-order measurement models. A higher-order,
mediation, or other regression-heavy SEM can still use the LDSC and model-fitting
stages, but its SNP stage may require iterative `userGWAS()` or an additional
estimator option rather than the current analytic full-genome path.

## Portability

The workflow uses Bash, Python with PyYAML, R, and a project-local R library.
Paths in configuration are relative to the project root, so the same
checkout can run locally, on Wynton, or on coreHPC. Set `GENOMICSEM_CORES`,
`NSLOTS`, or `SLURM_CPUS_PER_TASK` through the scheduler rather than editing
scripts.

Raw GWASs, controlled data, references, `.work/`, and large result files should
remain outside Git. Commit generic configuration, model syntax, scripts, and documentation. Keep private configuration in config/local/ and do not commit study outputs.

Generated outputs use `results/YYYY-MM-DD/`; see [dated-results instructions](../docs/RESULTS_LAYOUT.md) for date selection and migration details.
