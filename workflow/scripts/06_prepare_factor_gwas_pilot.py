#!/usr/bin/env python3
"""Build a deterministic, trait-agnostic pilot subset for userGWAS."""

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path, mode="rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else open(path, mode, newline="")


def read_table(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--traits", default=".work/config/traits.tsv")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--snps", type=int, default=500)
    args = parser.parse_args()
    if args.snps < 1:
        raise ValueError("--snps must be positive")

    root = Path(args.project_root).resolve()
    traits = read_table(root / args.traits)
    output_dir = root / ".work" / "factor_gwas" / "pilot"
    output_dir.mkdir(parents=True, exist_ok=True)

    first = root / traits[0]["file"]
    with open_text(first) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        wanted = []
        for row in reader:
            wanted.append(row["SNP"])
            if len(wanted) == args.snps:
                break
    wanted_set = set(wanted)
    (output_dir / "snps.txt").write_text("\n".join(wanted) + "\n")

    counts = []
    for trait in traits:
        source = root / trait["file"]
        destination = output_dir / f"{trait['trait']}.tsv.gz"
        kept = 0
        with open_text(source) as input_handle, gzip.open(destination, "wt", newline="") as output_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            writer = csv.DictWriter(output_handle, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                if row["SNP"] in wanted_set:
                    writer.writerow(row)
                    kept += 1
        counts.append((trait["trait"], kept))

    reference = root / args.reference
    with open_text(reference) as input_handle, open(output_dir / "reference.tsv", "w", newline="") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        for line_number, line in enumerate(input_handle):
            fields = line.split()
            if line_number == 0:
                writer.writerow(fields)
            elif fields and fields[0] in wanted_set:
                writer.writerow(fields)

    with open(output_dir / "pilot_counts.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", "data_rows"])
        writer.writerows(counts)
    print(f"Pilot inputs prepared in {output_dir}")


if __name__ == "__main__":
    main()
