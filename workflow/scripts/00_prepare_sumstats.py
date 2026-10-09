#!/usr/bin/env python3
"""Prepare arbitrary GWAS files from a trait configuration table."""

import argparse
import csv
import gzip
import math
import os
import re
from pathlib import Path
from paths import results_dir


RSID = re.compile(r"^rs[0-9]+$")
AUTOSOME = {str(i) for i in range(1, 23)}


def open_text(path, mode="rt"):
    if str(path).endswith(".gz") or str(path).endswith(".gz.tmp"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def read_traits(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "trait", "raw_file", "file", "full_file", "N", "chr_col", "bp_col",
        "snp_col", "a1_col", "a2_col", "af_col", "beta_col", "se_col",
        "p_col", "lp_col", "n_col", "effect_multiplier",
    }
    missing = required - set(rows[0]) if rows else required
    if missing:
        raise ValueError("Missing trait configuration columns: " + ", ".join(sorted(missing)))
    return rows


def read_hm3(path):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        return {row[0] for row in reader if row and RSID.match(row[0])}


def number(value):
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def find_header(handle):
    for line in handle:
        if not line.startswith("##"):
            return next(csv.reader([line], delimiter="\t"))
    return None


def prepare_trait(root, trait, mode, hm3, stats):
    source = root / trait["raw_file"]
    destination = root / (trait["file"] if mode == "ldsc" else trait["full_file"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(destination) + f".tmp.{os.getpid()}.gz")
    multiplier = float(trait["effect_multiplier"])
    seen = set()
    counts = {"input": 0, "kept": 0, "duplicates": 0, "invalid": 0}

    if not source.exists():
        raise FileNotFoundError(f"Raw GWAS does not exist: {source}")

    with open_text(source) as input_handle, gzip.open(temporary, "wt", newline="") as output_handle:
        header = find_header(input_handle)
        if not header:
            raise ValueError(f"No header found in {source}")
        required_source = [
            trait["snp_col"], trait["chr_col"], trait["bp_col"],
            trait["a1_col"], trait["a2_col"], trait["beta_col"], trait["se_col"],
        ]
        required_source += [x for x in (trait["p_col"], trait["lp_col"]) if x]
        missing = [column for column in required_source if column not in header]
        if missing:
            raise ValueError(f"{trait['trait']} is missing columns: {', '.join(missing)}")

        reader = csv.DictReader(input_handle, fieldnames=header, delimiter="\t")
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "CHR", "BP", "A1", "A2", "MAF", "BETA", "SE", "P", "N"])

        for row in reader:
            counts["input"] += 1
            snp = row.get(trait["snp_col"], "").strip()
            if mode == "ldsc" and snp not in hm3:
                continue
            if snp in seen:
                counts["duplicates"] += 1
                continue

            chrom = row.get(trait["chr_col"], "")
            if chrom.lower().startswith("chr"):
                chrom = chrom[3:]
            bp = row.get(trait["bp_col"], "")
            a1 = row.get(trait["a1_col"], "").upper()
            a2 = row.get(trait["a2_col"], "").upper()
            beta = number(row.get(trait["beta_col"], ""))
            se = number(row.get(trait["se_col"], ""))
            if not (RSID.match(snp) and chrom in AUTOSOME and bp.isdigit() and
                    a1 in "ACGT" and len(a1) == 1 and a2 in "ACGT" and len(a2) == 1 and
                    beta is not None and se is not None and se > 0):
                counts["invalid"] += 1
                continue

            if trait["p_col"]:
                p = number(row.get(trait["p_col"], ""))
            else:
                lp = number(row.get(trait["lp_col"], ""))
                p = None if lp is None or lp < 0 else (1e-300 if lp >= 300 else 10 ** (-lp))
            if p is None or not 0 < p <= 1:
                counts["invalid"] += 1
                continue

            af = number(row.get(trait["af_col"], "")) if trait["af_col"] else None
            maf = "NA" if af is None or not 0 <= af <= 1 else format(min(af, 1 - af), ".15g")
            n_value = trait["N"].strip()
            if n_value.upper() == "NA":
                n_value = ""
            if not n_value and trait["n_col"]:
                n_value = row.get(trait["n_col"], "")
            if not n_value or n_value == ".":
                n_value = "NA"

            writer.writerow([
                snp, chrom, bp, a1, a2, maf, format(beta * multiplier, ".15g"),
                format(se, ".15g"), format(p, ".15g"), n_value,
            ])
            seen.add(snp)
            counts["kept"] += 1

    temporary.replace(destination)
    stats.append({
        "trait": trait["trait"], **counts, "effect_multiplier": multiplier,
        "output": str(destination.relative_to(root)),
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--traits", default="config/traits.tsv")
    parser.add_argument("--mode", choices=("ldsc", "full"), required=True)
    parser.add_argument("--hm3", default="reference/w_hm3.snplist")
    parser.add_argument("--trait", action="append", help="Prepare only the named trait; may be repeated")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    traits = read_traits(root / args.traits)
    if args.trait:
        requested = set(args.trait)
        traits = [trait for trait in traits if trait["trait"] in requested]
        missing = requested - {trait["trait"] for trait in traits}
        if missing:
            raise ValueError("Unknown trait(s): " + ", ".join(sorted(missing)))
    hm3 = read_hm3(root / args.hm3) if args.mode == "ldsc" else None
    stats = []
    for trait in traits:
        print(f"Preparing {trait['trait']} ({args.mode})", flush=True)
        prepare_trait(root, trait, args.mode, hm3, stats)

    result = results_dir(root) / f"preparation_{args.mode}_qc.tsv"
    result.parent.mkdir(parents=True, exist_ok=True)
    with open(result, "w", newline="") as handle:
        fields = ["trait", "input", "kept", "duplicates", "invalid", "effect_multiplier", "output"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(stats)
    print(f"Preparation summary: {result}")


if __name__ == "__main__":
    main()
