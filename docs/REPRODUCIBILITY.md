# Reproducing an analysis

Every new `run_pipeline.sh` invocation (including those launched by Snakemake)
creates `.work/runs/<UTC-time>-<id>/`. Preserve this directory with your report.
It contains immutable copies of the configuration, model definitions and code,
a JSON manifest, and R runtime/package versions. The manifest records command,
selected environment overrides, platform, Python version, SHA-256 input hashes,
exit status and outputs changed during that invocation. Large files are hashed
in chunks; analysis-stage startup therefore adds disk-read time. Report-only
stages do not hash raw inputs. Existing results are inventoried separately and
are not claimed as products of the new invocation.

## Repeat the analysis

1. Use a clean project copy and restore code/configuration from a chosen run's
   `snapshot/`, using `snapshot_map` in its manifest for original locations.
2. Supply the same raw GWAS and reference files at the configured paths.
3. Match the R, Python, and R package versions in `r_environment.txt`. Preserve
   the original R library or provision those exact package versions for the
   target platform. The environment record is a version inventory, not an
   automatic dependency restore lockfile. Snakemake's conda YAML gives version
   ranges and is likewise not an exact environment lock.
4. `config/genomicsem_ref.txt` pins GenomicSEM to the installed commit used when
   this feature was added. Setup now checks that commit instead of silently
   accepting another installed version. An explicit `GENOMICSEM_REF` override
   is recorded; use a commit SHA. A deliberate reinstall uses
   `GENOMICSEM_REINSTALL=true ./run_pipeline.sh setup`.
5. Verify the copied code, configuration, and recorded inputs:

   ```bash
   python3 workflow/scripts/09_provenance.py verify --manifest .work/runs/RUN_ID/manifest.json
   ```

6. Run the recorded stages/settings and environment overrides. A complete
   sequence is `check`, `all`, review the model, `gwas-pilot`, review the pilot,
   then `gwas`. Do not reuse generated outputs from a different configuration.

The pilot chooses a deterministic subset; it does not randomly sample SNPs.
The existing R launcher limits BLAS/OpenMP threads to one. Exact floating-point
identity across operating systems or R/BLAS builds is not guaranteed. Successful
execution does not establish scientific validity of the model or its assumptions.
Historical outputs cannot acquire reliable original provenance after the fact.

## HTML report

`./run_pipeline.sh html-report` rebuilds the report without rerunning statistics.
The report focuses on GenomicSEM methods, findings, current output tables, and
figures. Execution manifests, runtime details, and hashes are retained separately
in `.work/runs/` and are not embedded in the scientific report.
Optional missing outputs are labeled unavailable. Keep reports and run archives
for each analysis before reusing this output directory for another dataset.

Generated outputs use `results/YYYY-MM-DD/`; see [dated-results instructions](../docs/RESULTS_LAYOUT.md) for date selection and migration details.
