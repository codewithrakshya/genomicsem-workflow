# Results layout

Report-ready outputs are stored in `results/YYYY-MM-DD/`. The date defaults to
America/Los_Angeles at invocation start. For a multi-day sequence, export a
single date before running stages, for example `GENOMICSEM_RUN_DATE=2026-10-08`.
Snakemake accepts `--config run_date=2026-10-08`.

Same-date runs replace existing results. `.work/` holds shared intermediates;
dated output folders do not isolate analyses. Use separate project copies for
different datasets or models. All generated results and intermediates are ignored
by Git. No results are distributed with this template.
