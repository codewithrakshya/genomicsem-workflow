# Reusing the GenomicSEM workflow

The workflow code is study-agnostic. A new analysis should change configuration
and model files, not scripts.

## Configuration layers

### `config/workflow.tsv`

Defines analysis-wide resources and settings:

- `analysis_id`: short output prefix;
- `traits_config`: trait table to use;
- `models_config`: model table to use;
- `hm3_reference`: HapMap3 list for munging/LDSC;
- `ld_reference` and `weight_reference`: LDSC directories;
- `gwas_reference`: SNP/allele/MAF table for `sumstats()`;
- `pilot_snps`: default pilot size;
- `info_filter` and `maf_filter`: QC thresholds.

To keep several analyses in one checkout, create another settings file and run:

```bash
GENOMICSEM_SETTINGS=config/my_analysis/workflow.tsv ./run_pipeline.sh validate
GENOMICSEM_SETTINGS=config/my_analysis/workflow.tsv ./run_pipeline.sh all
```

The equivalent Snakemake command is:

```bash
snakemake --config settings=config/my_analysis/workflow.tsv --cores 1
```

### `config/traits.tsv`

One row represents one GWAS. The table controls:

- internal and display names;
- original, LDSC-prepared, and genome-wide-prepared file paths;
- N and binary-trait prevalence values;
- ancestry and genome build;
- the input names of chromosome, position, rsID, alleles, allele frequency,
  beta, SE, P or -log10(P), and N columns;
- effect-direction multiplier;
- whether the effect/SE is logistic and whether the GWAS used OLS.

For example, `effect_multiplier=-1` reverses a trait's beta while leaving the
alleles unchanged. Use `1` when no reversal is required.

Copy `config/traits.example.tsv` when starting a new project. You may add or
remove trait rows; the preparation, munging, LDSC, plotting, pilot, and factor-
GWAS scripts determine the number and names of traits from this table.

### `config/models.tsv`

One row represents one lavaan/GenomicSEM model:

- `label`: safe output label;
- `model_file`: model syntax file;
- `role`: `baseline`, `candidate`, or `primary`;
- `factor_name`: latent factor tested in the SNP analysis;
- `snp_regression`: regression appended for `userGWAS()`, such as `F1 ~ SNP`.

Exactly one model must be marked `primary`. Any number of models can be fitted
and compared. The primary model is carried into the SNP-level GWAS.

## Starting a different analysis

1. Copy the repository without `data/raw/`, `reference/`, `.work/`, or large
   `results/` files.
2. Place or link the new GWAS files under `data/raw/`.
3. Copy `config/traits.example.tsv` and describe every input column accurately.
4. Create one or more model text files using the trait names from the table.
5. Update `config/models.tsv`, marking one scientifically justified model as
   `primary`.
6. Update reference paths and output prefix in `config/workflow.tsv`.
7. Validate before running anything:

```bash
./run_pipeline.sh validate
```

8. Run through the model stage:

```bash
./run_pipeline.sh all
```

9. Review genetic correlations, model fit, residual variances, and cohort
   overlap. Do not run a full factor GWAS just because a model converges.
10. Test the SNP model on a small subset:

```bash
./run_pipeline.sh gwas-pilot
```

11. Submit the full stage to a scheduler on an approved machine:

```bash
./run_pipeline.sh gwas
```

## Stable workflow stages

| Stage | Purpose | Driven by |
|---|---|---|
| `validate` | Check configuration, models, build consistency, and references | all configuration files |
| `prepare` | Standardize columns, direction, P values, and HapMap3 subset | traits + workflow |
| `munge` | GenomicSEM/LDSC munging | traits + workflow |
| `ldsc` | Genetic and sampling covariance | traits + workflow |
| `models` | Fit every configured model and compare results | models |
| `report` | Generic correlation/loading/residual figure | traits + models |
| `gwas-pilot` | Small iterative `userGWAS()` validation | traits + primary model |
| `gwas` | Genome-wide preparation and analytic factor GWAS | traits + primary model |

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

Those choices belong in the configuration, model files, and analysis-specific
documentation. The software mechanics should not encode them.

The full-GWAS stage currently uses GenomicSEM's fast analytic `userGWAS()`
estimator, which is intended for first-order measurement models. A higher-order,
mediation, or other regression-heavy SEM can still use the LDSC and model-fitting
stages, but its SNP stage may require iterative `userGWAS()` or an additional
estimator option rather than the current analytic full-genome path.

## Portability

The workflow uses Bash, Python's standard library, R, and a project-local R
library. Paths in configuration are relative to the project root, so the same
checkout can run locally, on Wynton, or on coreHPC. Set `GENOMICSEM_CORES`,
`NSLOTS`, or `SLURM_CPUS_PER_TASK` through the scheduler rather than editing
scripts.

Raw GWASs, controlled data, references, `.work/`, and large result files should
remain outside Git. Commit generic configuration, model syntax, scripts, and documentation. Keep private configuration in config/local/ and do not commit study outputs.

Generated outputs use `results/YYYY-MM-DD/`; see [dated-results instructions](../docs/RESULTS_LAYOUT.md) for date selection and migration details.
