#!/usr/bin/env bash
set -euo pipefail

project_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "${project_root}"

stage=${1:-all}
settings_file=${GENOMICSEM_SETTINGS:-config/workflow.tsv}
export GENOMICSEM_RUN_DATE=$(python3 workflow/scripts/paths.py)
results_dir="results/${GENOMICSEM_RUN_DATE}"
mkdir -p "${results_dir}"

setting() {
  awk -F '\t' -v key="$1" 'NR > 1 && $1 == key { print $2; exit }' "${settings_file}"
}

traits_config=$(setting traits_config)
models_config=$(setting models_config)
hm3=$(setting hm3_reference)
ld_reference=$(setting ld_reference)
weight_reference=$(setting weight_reference)
gwas_reference=$(setting gwas_reference)
export GENOMICSEM_INFO_FILTER=$(setting info_filter)
export GENOMICSEM_MAF_FILTER=$(setting maf_filter)

validate() {
  workflow/bin/run-r workflow/scripts/validate_config.R "${settings_file}"
}

setup() {
  workflow/bin/run-r workflow/scripts/install_genomicsem.R
}

check() {
  workflow/bin/run-r workflow/scripts/check_environment.R
}

prepare() {
  python3 workflow/scripts/00_prepare_sumstats.py \
    --project-root . --traits "${traits_config}" --mode ldsc --hm3 "${hm3}"
}

munge() {
  workflow/bin/run-r workflow/scripts/01_munge.R "${traits_config}" "${hm3}"
}

ldsc() {
  workflow/bin/run-r workflow/scripts/02_ldsc.R "${traits_config}" "${ld_reference}" "${weight_reference}"
}

models() {
  while IFS=$'\t' read -r label model_file role factor_name snp_regression
  do
    [[ "${label}" == "label" ]] && continue
    workflow/bin/run-r workflow/scripts/03_fit_model.R \
      .work/ldsc/ldsc_output.rds "${model_file}" "${label}"
  done < "${models_config}"
  workflow/bin/run-r workflow/scripts/04_compare_models.R "${models_config}" "${results_dir}"
}

html_report() {
  python3 workflow/scripts/08_render_report.py --project-root . --settings "${settings_file}"
}

report() {
  workflow/bin/run-r workflow/scripts/05_visualize_results.R . "${results_dir}" "${settings_file}"
  html_report
}

gwas_pilot() {
  python3 workflow/scripts/06_prepare_factor_gwas_pilot.py \
    --project-root . --traits "${traits_config}" --reference "${gwas_reference}" \
    --snps "${GENOMICSEM_PILOT_SNPS:-$(setting pilot_snps)}"
  workflow/bin/run-r workflow/scripts/07_factor_gwas.R pilot "${settings_file}"
  html_report
}

gwas() {
  python3 workflow/scripts/00_prepare_sumstats.py \
    --project-root . --traits "${traits_config}" --mode full --hm3 "${hm3}"
  workflow/bin/run-r workflow/scripts/07_factor_gwas.R full "${settings_file}"
  html_report
}

case "${stage}" in
  setup|check|all|validate|prepare|munge|ldsc|models|report|html-report|gwas-pilot|gwas) ;;
  *) echo "Unknown stage: ${stage}" >&2; exit 2 ;;
esac

run_manifest=$(python3 workflow/scripts/09_provenance.py start --settings "${settings_file}" --stage "${stage}")
export GENOMICSEM_RUN_MANIFEST="${run_manifest}"
finish_run() {
  run_exit=$?
  trap - EXIT
  python3 workflow/scripts/09_provenance.py finish --manifest "${run_manifest}" --exit-code "${run_exit}" || run_exit=1
  if [[ "${run_exit}" == 0 ]]; then
    case "${stage}" in
      all|report|html-report|gwas-pilot|gwas) html_report || run_exit=$? ;;
    esac
  fi
  exit "${run_exit}"
}
trap finish_run EXIT

case "${stage}" in
  setup) setup ;;
  check) check ;;
  all)
    validate
    prepare
    munge
    ldsc
    models
    report
    ;;
  validate) validate ;;
  prepare) prepare ;;
  munge) munge ;;
  ldsc) ldsc ;;
  models) models ;;
  report) report ;;
  html-report) html_report ;;
  gwas-pilot) gwas_pilot ;;
  gwas) gwas ;;
  *)
    echo "Usage: $0 [setup|check|validate|all|prepare|munge|ldsc|models|report|html-report|gwas-pilot|gwas]" >&2
    exit 2
    ;;
esac

echo "Done. Review the report-ready files in ${project_root}/${results_dir}"
