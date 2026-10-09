"""Snakemake interface for the same stages exposed by run_pipeline.sh."""

import csv
import os
import sys
from pathlib import Path


ROOT = Path(workflow.basedir)
sys.path.insert(0, str(ROOT / "workflow/scripts"))
from configuration import materialize_config
from paths import run_date
if "run_date" in config:
    os.environ["GENOMICSEM_RUN_DATE"] = str(config["run_date"])
os.environ["GENOMICSEM_RUN_DATE"] = run_date()
RESULTS = "results/" + os.environ["GENOMICSEM_RUN_DATE"]
INTERMEDIATE_RESULTS = ".work/results/" + os.environ["GENOMICSEM_RUN_DATE"]
CONFIG_FILE = Path(config.get("settings", "config/analysis.yaml"))
materialize_config(ROOT, CONFIG_FILE)
SETTINGS_FILE = Path(".work/config/workflow.tsv")
os.environ["GENOMICSEM_SETTINGS"] = str(CONFIG_FILE)


def read_key_value(path):
    with open(ROOT / path, newline="") as handle:
        return {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def read_rows(path):
    with open(ROOT / path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


SETTINGS = read_key_value(SETTINGS_FILE)
TRAITS_FILE = Path(SETTINGS["traits_config"])
MODELS_FILE = Path(SETTINGS["models_config"])
TRAITS = read_rows(TRAITS_FILE)
MODELS = read_rows(MODELS_FILE)
PRIMARY = next(model for model in MODELS if model["role"] == "primary")
ANALYSIS_ID = SETTINGS["analysis_id"]

RAW_FILES = [trait["raw_file"] for trait in TRAITS]
PREPARED_FILES = [trait["file"] for trait in TRAITS]
FULL_FILES = [trait["full_file"] for trait in TRAITS]
MUNGED_FILES = [f".work/munged/{trait['trait']}.sumstats.gz" for trait in TRAITS]
MODEL_FILES = [model["model_file"] for model in MODELS]
MODEL_FITS = [f".work/models/{model['label']}/model_fit.rds" for model in MODELS]
PRIMARY_FIT = f".work/models/{PRIMARY['label']}/model_fit.rds"
FACTOR_GWAS = f"{RESULTS}/{ANALYSIS_ID}_factor_gwas.tsv.gz"
FACTOR_SUMMARY = f"{RESULTS}/{ANALYSIS_ID}_factor_gwas_summary.tsv"
PILOT_RESULT = f"{RESULTS}/{ANALYSIS_ID}_factor_gwas_pilot.tsv"
LD_FILES = [
    f"{SETTINGS['ld_reference']}/{chromosome}.l2.{suffix}"
    for chromosome in range(1, 23)
    for suffix in ("ldscore.gz", "M_5_50")
]
WEIGHT_FILES = [
    f"{SETTINGS['weight_reference']}/{chromosome}.l2.{suffix}"
    for chromosome in range(1, 23)
    for suffix in ("ldscore.gz", "M_5_50")
]
REFERENCE_LD_FILES = sorted(set(LD_FILES + WEIGHT_FILES))


rule all:
    input:
        f"{RESULTS}/{ANALYSIS_ID}_report.html",
        f"{RESULTS}/genomicsem_summary.png",
        f"{RESULTS}/genomicsem_summary.pdf",


rule setup:
    output:
        "workflow/R/library/GenomicSEM/DESCRIPTION",
    shell:
        "./run_pipeline.sh setup"


rule check:
    input:
        rules.setup.output,
    output:
        ".work/setup/environment/sessionInfo.txt",
    shell:
        "./run_pipeline.sh check"


rule validate:
    input:
        str(CONFIG_FILE),
        str(SETTINGS_FILE),
        str(TRAITS_FILE),
        str(MODELS_FILE),
        MODEL_FILES,
        SETTINGS["hm3_reference"],
        SETTINGS["gwas_reference"],
    output:
        touch(".work/checkpoints/config.valid"),
    shell:
        "./run_pipeline.sh validate"


rule prepare:
    input:
        rules.validate.output,
        RAW_FILES,
        SETTINGS["hm3_reference"],
    output:
        PREPARED_FILES,
        f"{INTERMEDIATE_RESULTS}/preparation_ldsc_qc.tsv",
    shell:
        "./run_pipeline.sh prepare"


rule munge:
    input:
        PREPARED_FILES,
        SETTINGS["hm3_reference"],
    output:
        MUNGED_FILES,
    shell:
        "./run_pipeline.sh munge"


rule ldsc:
    input:
        MUNGED_FILES,
        REFERENCE_LD_FILES,
    output:
        ".work/ldsc/ldsc_output.rds",
        f"{INTERMEDIATE_RESULTS}/genetic_correlations.tsv",
        f"{INTERMEDIATE_RESULTS}/genetic_covariances.tsv",
    params:
        ld=SETTINGS["ld_reference"],
        weights=SETTINGS["weight_reference"],
    shell:
        "./run_pipeline.sh ldsc"


rule models:
    input:
        ".work/ldsc/ldsc_output.rds",
        MODEL_FILES,
    output:
        MODEL_FITS,
        f"{INTERMEDIATE_RESULTS}/model_fit_comparison.tsv",
        f"{INTERMEDIATE_RESULTS}/model_parameter_comparison.tsv",
        f"{INTERMEDIATE_RESULTS}/model_comparison_summary.txt",
    shell:
        "./run_pipeline.sh models"


rule report:
    input:
        f"{INTERMEDIATE_RESULTS}/genetic_correlations.tsv",
        f"{INTERMEDIATE_RESULTS}/model_fit_comparison.tsv",
        f"{INTERMEDIATE_RESULTS}/model_parameter_comparison.tsv",
    output:
        f"{RESULTS}/genomicsem_summary.png",
        f"{RESULTS}/genomicsem_summary.pdf",
    shell:
        "./run_pipeline.sh report"


rule gwas_pilot:
    input:
        PREPARED_FILES,
        ".work/ldsc/ldsc_output.rds",
        PRIMARY_FIT,
        SETTINGS["gwas_reference"],
    output:
        PILOT_RESULT,
    shell:
        "./run_pipeline.sh gwas-pilot"


rule gwas:
    input:
        RAW_FILES,
        ".work/ldsc/ldsc_output.rds",
        PRIMARY_FIT,
        SETTINGS["gwas_reference"],
    output:
        FULL_FILES,
        FACTOR_GWAS,
        FACTOR_SUMMARY,
    shell:
        "./run_pipeline.sh gwas"


rule html_report:
    input:
        f"{RESULTS}/genomicsem_summary.png",
        f"{INTERMEDIATE_RESULTS}/genetic_covariances.tsv",
        f"{INTERMEDIATE_RESULTS}/genetic_correlations.tsv",
        f"{INTERMEDIATE_RESULTS}/model_fit_comparison.tsv",
        f"{INTERMEDIATE_RESULTS}/model_parameter_comparison.tsv",
        f"{INTERMEDIATE_RESULTS}/model_comparison_summary.txt",
        f"{INTERMEDIATE_RESULTS}/preparation_ldsc_qc.tsv",
        str(CONFIG_FILE),
        str(TRAITS_FILE),
        str(MODELS_FILE),
        MODEL_FILES,
        "workflow/scripts/08_render_report.py",
        "workflow/scripts/paths.py",
        SETTINGS["cohort_overlap"],
    output:
        f"{RESULTS}/{ANALYSIS_ID}_report.html",
    shell:
        "./run_pipeline.sh html-report"
