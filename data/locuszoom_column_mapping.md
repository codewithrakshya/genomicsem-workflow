# LocusZoom column mapping

All three datasets use GRCh37 coordinates. The upload file must be tab-delimited and sorted by chromosome and genomic position.

| Dataset | Chromosome | Position | Reference allele | Alternate allele | Variant ID | P-value | −log10(P-value) | Effect size (ALT allele) | Standard error | ALT allele frequency |
|---|---|---:|---|---|---|---|---|---|---|---|
| Austin Zimmerman 2023 — Short Sleep | `CHR` | `POS` | `REF` | `ALT` | `DBSNP_ID` | — | `LP` | `ES` | `SE` | `AF` |
| Dashti 2019 — Sleep Duration | `CHROM` | `POS` | `REF` | `ALT` | `DBSNP_ID` | `P` | — | `BETA` | `SE` | `AF` |
| Portas 2026 — Sleep Duration | `CHR` | `POS` | `REF` | `ALT` | `DBSNP_ID` | — | `LP` | `ES` | `SE` | `AF` |

## Canonical LocusZoom columns

For a prepared P-value upload:

```tsv
CHROM	POS	REF	ALT	PVAL	BETA	SE	AF	SNP
```

For a prepared −log10(P-value) upload:

```tsv
CHROM	POS	REF	ALT	LOG10P	BETA	SE	AF	SNP
```

## Dataset-specific notes

- Austin and Portas provide `LP`, so use `LP` as the LocusZoom `LOG10P`/−log10(P) field rather than applying a log transformation again.
- Dashti provides raw P-values in `P`, so map that column to LocusZoom `PVAL`, or generate a `LOG10P` column from `-log10(P)`.
- Dashti metadata lines beginning with `##` must be removed from an upload file if they are not already stripped during preparation; retain the TSV header beginning with `ID`.
- `REF` and `ALT` are required for LD annotation and must match the selected GRCh37 reference build.
- Every allele-frequency value must consistently refer to the ALT allele.
- LocusZoom's coordinate sorting requirement applies to the final upload file, not the source files as stored here.
