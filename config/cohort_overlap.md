# Cohort-overlap review

## Studies

- **Austin-Zimmerman 2023:** European short-sleep meta-analysis of UK
  Biobank and Million Veteran Program participants (47,054 short-sleep cases
  and 382,950 normal-sleep controls). Short sleep was self-reported and
  defined as 5 hours or less; normal sleep was 7-8 hours.
- **Dashti 2019:** self-reported habitual sleep duration in 446,118 European
  ancestry UK Biobank participants.
- **Portas 2026:** device-measured night-time sleep duration in 80,013 UK
  Biobank participants, selected from the wrist-accelerometer subset.

## Expected overlap

- Dashti and Portas have substantial overlap because both use UK Biobank and
  Portas is based on an accelerometer subset of that cohort.
- Austin partially overlaps both studies through its UK Biobank component;
  its MVP component is independent of Dashti and Portas.
- Exact participant counts shared across the three GWAS cannot be recovered
  from summary statistics alone.
- The non-zero LDSC cross-trait intercepts are consistent with overlap. The
  GenomicSEM sampling-covariance matrix accounts for this, so the studies
  should not be treated as independent in a conventional meta-analysis.

## Interpretation for the current factor model

The overlap is real and expected, not evidence that the summary-statistics
files were duplicated. GenomicSEM's multivariable LDSC step models correlated
sampling error from overlap. The unrestricted Heywood case is driven by the
observed genetic-correlation pattern and the three-indicator factor structure;
overlap alone should not be named as its cause.

## Sources

- Austin-Zimmerman et al. (2023):
  https://www.nature.com/articles/s41467-023-41249-y
- Dashti et al. (2019):
  https://www.nature.com/articles/s41467-019-08917-4
- Portas et al. (2026):
  https://doi.org/10.1038/s41467-026-71252-y
