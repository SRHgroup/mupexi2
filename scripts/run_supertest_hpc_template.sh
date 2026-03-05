#!/usr/bin/env bash
set -euo pipefail

# MuPeXi2 SUPERTEST (HPC template)
# Copy this file to scripts/run_supertest_hpc_local.sh and fill your real paths.

# Example:
# module load tools ngs anaconda3/2025.06-1 netmhcpan/4.0a perl/5.36.1 ensembl-tools/90

REPO_DIR="/path/to/mupexi2"
export PYTHONPATH="$REPO_DIR/src:${PYTHONPATH:-}"

VCF="/path/to/test_subset.vcf.gz"
EXPR="/path/to/expression.tsv"
FUSION_TSV="/path/to/fusions.tsv"
CONFIG="/path/to/config.computerome.ini"
ALLELES="HLA-A02:01"
OUT_BASE="/path/to/output/supertest"

mkdir -p "$OUT_BASE"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="$OUT_BASE/$RUN_ID"
mkdir -p "$RUN_DIR"

SUMMARY_TSV="$RUN_DIR/summary.tsv"
SUMMARY_MD="$RUN_DIR/summary.md"

echo -e "case_id\tkind\texpected_rc\tactual_rc\tstatus\tprefix\tlog" > "$SUMMARY_TSV"

# case format: case_id|kind|expected_rc|extra_args
# kind in {snv,fusion,mixed}
cases=(
  "merged_super_on|snv|0|--vcf-type merged --germlines true --rna-edit true --superpeptides true --tumor-sample TUMOR --normal-sample DNA_NORMAL --parallel-k true -l 9-11"
  "merged_super_off|snv|0|--vcf-type merged --germlines true --rna-edit true --superpeptides false --tumor-sample TUMOR --normal-sample DNA_NORMAL --parallel-k false -l 9-11"
  "merged_no_rna|snv|0|--vcf-type merged --germlines true --rna-edit false --superpeptides true --tumor-sample TUMOR --normal-sample DNA_NORMAL -l 9-11"
  "somatic_only_mutect2|snv|0|--vcf-type mutect2 --germlines false --rna-edit false --tumor-sample TUMOR --normal-sample DNA_NORMAL -l 9-11"
  "known_only|snv|0|--vcf-type merged --germlines true --rna-edit true --rnaedit-known-only --tumor-sample TUMOR --normal-sample DNA_NORMAL -l 9-11"
  "allow_novel|snv|0|--vcf-type merged --germlines true --rna-edit true --rnaedit-allow-novel --tumor-sample TUMOR --normal-sample DNA_NORMAL -l 9-11"
  "context_qc_relaxed|snv|0|--vcf-type merged --germlines true --rna-edit true --superpeptides true --tumor-sample TUMOR --normal-sample DNA_NORMAL --context-rna-min-depth 5 --context-rna-min-alt-count 2 --context-rna-min-vaf 0.02 --context-somatic-min-depth 8 --context-somatic-min-alt-count 3 --context-somatic-min-vaf 0.03 -l 9-11"
  "context_qc_strict|snv|0|--vcf-type merged --germlines true --rna-edit true --superpeptides true --tumor-sample TUMOR --normal-sample DNA_NORMAL --context-rna-min-depth 20 --context-rna-min-alt-count 6 --context-rna-min-vaf 0.10 --context-somatic-min-depth 20 --context-somatic-min-alt-count 6 --context-somatic-min-vaf 0.10 -l 9-11"
  "strict_phasing_expected_fail|snv|1|--vcf-type merged --germlines true --rna-edit true --phasing-mode strict --tumor-sample TUMOR --normal-sample DNA_NORMAL -l 9-11"
  "fusion_only|fusion|0|-l 9-11"
  "snv_plus_fusion|mixed|0|--vcf-type merged --germlines true --rna-edit true --tumor-sample TUMOR --normal-sample DNA_NORMAL -l 9-11"
)

ok_count=0
unexpected_count=0
skip_count=0
case_idx=0
total_cases="${#cases[@]}"

run_case() {
  local case_id="$1"
  local kind="$2"
  local expected_rc="$3"
  local extra_args="$4"

  local prefix="supertest_${case_id}"
  local case_out="$RUN_DIR/$case_id"
  local log_file="$RUN_DIR/${case_id}.log"
  mkdir -p "$case_out"
  case_idx=$((case_idx + 1))

  local -a cmd=(python3 -m mupexi2.cli)

  if [[ "$kind" == "snv" || "$kind" == "mixed" ]]; then
    if [[ ! -f "$VCF" ]]; then
      echo -e "${case_id}\t${kind}\t${expected_rc}\tNA\tSKIP\t${prefix}\t${log_file}" >> "$SUMMARY_TSV"
      echo "[$(date +'%F %T')] CASE=$case_id SKIP missing VCF $VCF" | tee "$log_file"
      skip_count=$((skip_count + 1))
      return
    fi
    cmd+=(-v "$VCF")
  fi

  if [[ "$kind" == "fusion" || "$kind" == "mixed" ]]; then
    if [[ ! -f "$FUSION_TSV" ]]; then
      echo -e "${case_id}\t${kind}\t${expected_rc}\tNA\tSKIP\t${prefix}\t${log_file}" >> "$SUMMARY_TSV"
      echo "[$(date +'%F %T')] CASE=$case_id SKIP missing FUSION_TSV $FUSION_TSV" | tee "$log_file"
      skip_count=$((skip_count + 1))
      return
    fi
    cmd+=(-z "$FUSION_TSV")
  fi

  cmd+=(-a "$ALLELES" -t -f -n -c "$CONFIG" -p "$prefix" -d "$case_out")
  if [[ "$kind" != "fusion" && -f "$EXPR" ]]; then
    cmd+=(-e "$EXPR")
  fi

  # shellcheck disable=SC2206
  local extra=( $extra_args )
  cmd+=("${extra[@]}")

  echo
  echo "============================================================"
  echo "CASE ${case_idx}/${total_cases}: $case_id"
  echo "kind=$kind expected_rc=$expected_rc"
  echo "log=$log_file"
  echo "CMD: ${cmd[*]}"
  echo "============================================================"
  {
    echo "[$(date +'%F %T')] CASE=$case_id START kind=$kind expected_rc=$expected_rc"
    echo "CMD: ${cmd[*]}"
  } > "$log_file"

  set +e
  "${cmd[@]}" 2>&1 | tee -a "$log_file"
  local rc=${PIPESTATUS[0]}
  set -e

  local artifact_ok=1
  if [[ "$rc" -eq 0 && "$expected_rc" -eq 0 ]]; then
    if [[ "$kind" == "snv" ]]; then
      [[ -f "$case_out/${prefix}_snv.mupexi" ]] || artifact_ok=0
    elif [[ "$kind" == "fusion" ]]; then
      [[ -f "$case_out/${prefix}_fus.mupexi" ]] || artifact_ok=0
    else
      [[ -f "$case_out/${prefix}_snv.mupexi" ]] || artifact_ok=0
      [[ -f "$case_out/${prefix}_fus.mupexi" ]] || artifact_ok=0
    fi
  fi

  local status
  if [[ "$rc" -eq "$expected_rc" && "$artifact_ok" -eq 1 ]]; then
    status="OK"
    ok_count=$((ok_count + 1))
  else
    status="UNEXPECTED"
    unexpected_count=$((unexpected_count + 1))
  fi

  echo "[$(date +'%F %T')] CASE=$case_id END rc=$rc expected=$expected_rc artifact_ok=$artifact_ok status=$status" | tee -a "$log_file"
  echo -e "${case_id}\t${kind}\t${expected_rc}\t${rc}\t${status}\t${prefix}\t${log_file}" >> "$SUMMARY_TSV"
}

for row in "${cases[@]}"; do
  IFS='|' read -r case_id kind expected_rc extra_args <<< "$row"
  run_case "$case_id" "$kind" "$expected_rc" "$extra_args"
done

{
  echo "# MuPeXi2 Supertest Summary"
  echo
  echo "- Run ID: $RUN_ID"
  echo "- Output dir: $RUN_DIR"
  echo "- Total cases: ${#cases[@]}"
  echo "- OK: $ok_count"
  echo "- Skipped: $skip_count"
  echo "- Unexpected: $unexpected_count"
  echo
  echo "| case_id | kind | expected_rc | actual_rc | status | log |"
  echo "|---|---|---:|---:|---|---|"
  tail -n +2 "$SUMMARY_TSV" | while IFS=$'\t' read -r c k e a s p l; do
    echo "| $c | $k | $e | $a | $s | $l |"
  done
} > "$SUMMARY_MD"

echo "Supertest done"
echo "Summary TSV: $SUMMARY_TSV"
echo "Summary MD : $SUMMARY_MD"

if [[ "$unexpected_count" -gt 0 ]]; then
  exit 1
fi
