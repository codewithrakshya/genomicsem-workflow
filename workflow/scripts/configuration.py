#!/usr/bin/env python3
"""Load the single YAML analysis config and materialize workflow inputs."""

import argparse
import csv
from pathlib import Path

import yaml


TRAIT_FIELDS = [
    "trait", "display_name", "raw_file", "file", "full_file", "N",
    "sample_prev", "population_prev", "type", "ancestry", "build",
    "chr_col", "bp_col", "snp_col", "a1_col", "a2_col", "af_col",
    "beta_col", "se_col", "p_col", "lp_col", "n_col",
    "effect_multiplier", "se_logit", "OLS",
]
MODEL_FIELDS = ["label", "model_file", "role", "factor_name", "snp_regression"]
COHORT_FIELDS = [
    "trait", "study", "cohorts", "phenotype", "measurement", "ancestry",
    "sample_size", "overlap_note",
]


def require_mapping(value, name):
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a mapping")
    return value


def require_rows(value, name):
    if not isinstance(value, list) or not value:
        raise ValueError(f"{name} must be a non-empty list")
    for index, row in enumerate(value):
        require_mapping(row, f"{name}[{index}]")
    return value


def text(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return str(value)


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def materialize_config(project_root, config_path):
    root = Path(project_root).resolve()
    source = Path(config_path)
    if not source.is_absolute():
        source = root / source
    with source.open() as handle:
        config = require_mapping(yaml.safe_load(handle), "analysis config")

    analysis_id = text(config.get("analysis_id", ""))
    if not analysis_id or Path(analysis_id).name != analysis_id or analysis_id in (".", ".."):
        raise ValueError("analysis_id must be a non-empty filename component")
    references = require_mapping(config.get("references"), "references")
    results = require_mapping(config.get("results"), "results")
    traits = require_rows(config.get("traits"), "traits")
    models = require_rows(config.get("models"), "models")
    cohort_metadata = config.get("cohort_metadata", [])
    if not isinstance(cohort_metadata, list):
        raise ValueError("cohort_metadata must be a list")

    output = root / ".work/config"
    trait_rows = []
    trait_ids = set()
    for index, trait in enumerate(traits):
        trait_id = text(trait.get("id", ""))
        if not trait_id or trait_id in trait_ids or not trait_id.replace("_", "").replace("-", "").isalnum():
            raise ValueError(f"traits[{index}].id must be unique and use letters, numbers, _ or -")
        trait_ids.add(trait_id)
        columns = require_mapping(trait.get("columns"), f"traits[{index}].columns")
        required = ("raw_file", "N", "type", "ancestry", "build")
        missing = [key for key in required if key not in trait]
        if missing:
            raise ValueError(f"traits[{index}] is missing: {', '.join(missing)}")
        for key in ("chr", "bp", "snp", "a1", "a2", "beta", "se"):
            if not columns.get(key):
                raise ValueError(f"traits[{index}].columns.{key} is required")
        if not columns.get("p") and not columns.get("lp"):
            raise ValueError(f"traits[{index}] needs columns.p or columns.lp")
        trait_rows.append({
            "trait": trait_id,
            "display_name": text(trait.get("display_name", trait_id)),
            "raw_file": text(trait["raw_file"]),
            "file": f".work/prepared/{trait_id}.hm3.tsv.gz",
            "full_file": f".work/factor_gwas/full/{trait_id}.tsv.gz",
            "N": text(trait["N"]),
            "sample_prev": text(trait.get("sample_prev", "NA")),
            "population_prev": text(trait.get("population_prev", "NA")),
            "type": text(trait["type"]),
            "ancestry": text(trait["ancestry"]),
            "build": text(trait["build"]),
            "chr_col": text(columns.get("chr")),
            "bp_col": text(columns.get("bp")),
            "snp_col": text(columns.get("snp")),
            "a1_col": text(columns.get("a1")),
            "a2_col": text(columns.get("a2")),
            "af_col": text(columns.get("af")),
            "beta_col": text(columns.get("beta")),
            "se_col": text(columns.get("se")),
            "p_col": text(columns.get("p")),
            "lp_col": text(columns.get("lp")),
            "n_col": text(columns.get("n")),
            "effect_multiplier": text(trait.get("effect_multiplier", 1)),
            "se_logit": text(trait.get("se_logit", False)),
            "OLS": text(trait.get("OLS", True)),
        })

    model_rows = []
    model_labels = set()
    primary_count = 0
    models_dir = output / "models"
    for index, model in enumerate(models):
        label = text(model.get("label", ""))
        if not label or label in model_labels or not label.replace("_", "").replace("-", "").isalnum():
            raise ValueError(f"models[{index}].label must be unique and use letters, numbers, _ or -")
        model_labels.add(label)
        missing = [key for key in ("role", "factor_name", "snp_regression", "syntax") if not model.get(key)]
        if missing:
            raise ValueError(f"models[{index}] is missing: {', '.join(missing)}")
        if model["role"] == "primary":
            primary_count += 1
        model_path = models_dir / f"{label}.txt"
        model_path.parent.mkdir(parents=True, exist_ok=True)
        model_path.write_text(text(model["syntax"]).rstrip() + "\n")
        model_rows.append({
            "label": label,
            "model_file": str(model_path.relative_to(root)),
            "role": text(model["role"]),
            "factor_name": text(model["factor_name"]),
            "snp_regression": text(model["snp_regression"]),
        })
    if primary_count != 1:
        raise ValueError("Exactly one model must have role=primary")

    for index, cohort in enumerate(cohort_metadata):
        require_mapping(cohort, f"cohort_metadata[{index}]")
    write_tsv(output / "traits.tsv", TRAIT_FIELDS, trait_rows)
    write_tsv(output / "models.tsv", MODEL_FIELDS, model_rows)
    write_tsv(output / "cohort_metadata.tsv", COHORT_FIELDS, cohort_metadata)
    overlap = config.get("cohort_overlap", "")
    if not isinstance(overlap, str):
        raise ValueError("cohort_overlap must be text")
    (output / "cohort_overlap.md").write_text(overlap.rstrip() + "\n")

    settings = {
        "analysis_id": analysis_id,
        "traits_config": ".work/config/traits.tsv",
        "models_config": ".work/config/models.tsv",
        "hm3_reference": text(references.get("hm3", "")),
        "ld_reference": text(references.get("ld_scores", "")),
        "weight_reference": text(references.get("ld_weights", "")),
        "gwas_reference": text(references.get("gwas", "")),
        "pilot_snps": text(results.get("pilot_snps", 500)),
        "info_filter": text(results.get("info_filter", 0.9)),
        "maf_filter": text(results.get("maf_filter", 0.01)),
        "cohort_metadata": ".work/config/cohort_metadata.tsv",
        "cohort_overlap": ".work/config/cohort_overlap.md",
    }
    write_tsv(output / "workflow.tsv", ["key", "value"], [
        {"key": key, "value": value} for key, value in settings.items()
    ])
    return settings


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("materialize", choices=["materialize"])
    parser.add_argument("--config", default="config/analysis.yaml")
    parser.add_argument("--project-root", default=".")
    args = parser.parse_args()
    materialize_config(args.project_root, args.config)