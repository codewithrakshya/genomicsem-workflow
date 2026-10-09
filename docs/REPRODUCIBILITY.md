# Reproducibility

Each stage invocation records its command, date, configuration/code snapshots,
SHA-256 hashes, Python/platform details, R/package versions, and exit status in
`.work/runs/`. Analysis stages hash configured raw inputs and references; reporting
stages do not. Preserve the run directory privately with the corresponding data.
Verify a manifest with:

```bash
python3 workflow/scripts/09_provenance.py verify --manifest .work/runs/RUN_ID/manifest.json
```

GenomicSEM is pinned in `config/genomicsem_ref.txt`. Setup checks this commit.
An explicit commit override uses `GENOMICSEM_REF`; deliberate reinstall uses
`GENOMICSEM_REINSTALL=true ./run_pipeline.sh setup`.

R dependencies are recorded, not automatically restored from an exact lockfile.
Match the saved versions and preserve the R environment when reproducing work.
The conda YAML also uses version ranges. Identical results across different
R/BLAS/operating-system builds are not guaranteed. The pilot subset is deterministic.
