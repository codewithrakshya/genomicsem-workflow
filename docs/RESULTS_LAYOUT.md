# Dated results

Final outputs are stored in `results/YYYY-MM-DD/`: the HTML report, summary
figures, the full-GWAS Manhattan plot, and optional pilot/full factor-GWAS
outputs. Supporting QC, genetic
covariance/correlation, and model-comparison tables are stored in
`.work/results/YYYY-MM-DD/` and embedded in the HTML report. Prepared and munged
sumstats and LDSC/model objects are also kept in `.work/`.
The default date is the invocation start date in America/Los_Angeles.
Set the same explicit date for a multi-day sequence:

```bash
export GENOMICSEM_RUN_DATE=2026-10-08
./run_pipeline.sh all
./run_pipeline.sh gwas-pilot
./run_pipeline.sh gwas
./run_pipeline.sh html-report
```

Snakemake: `snakemake --config run_date=2026-10-08 --cores 1`.
The selected date is exported to every stage and recorded in provenance.
Reruns on the same date replace that date's outputs; choose a different date to
keep separate daily result collections. Technical intermediates in `.work/`
remain shared: dated result folders are not isolated computational workspaces.
For different datasets/configurations, use separate project copies.

Existing flat outputs were moved to `results/2026-10-08/` during the layout
migration, preserving original modification timestamps. This folder date is
not a claim that the historical statistical analyses ran on that day.
Cohort-overlap input documentation now lives in `config/cohort_overlap.md`.
