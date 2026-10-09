# Reproducibility

GenomicSEM is pinned in `config/genomicsem_ref.txt`. Setup checks this commit.
An explicit commit override uses `GENOMICSEM_REF`; deliberate reinstall uses
`GENOMICSEM_REINSTALL=true ./run_pipeline.sh setup`.

The workflow does not create run manifests or snapshot source files. Preserve
the analysis config, source data, references, and installed package details
separately when you need to reproduce a run.

R dependencies are recorded, not automatically restored from an exact lockfile.
Match the saved versions and preserve the R environment when reproducing work.
The conda YAML also uses version ranges. Identical results across different
R/BLAS/operating-system builds are not guaranteed. The pilot subset is deterministic.
