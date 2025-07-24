#!/usr/bin/env bash
################################################################################
# FULL-SCALE RUN PLAYBOOK  (macOS-friendly)
# – Duplicates *-testmode YAMLs -> *-full YAMLs
# – Removes test_mode flags
# – Inserts ${DESEQ_OBJ}, ${DESEQ_CONTRASTS}, ${ATAC_OBJ}, ${ATAC_CONTRASTS}
#   placeholders into Step-2 YAMLs
# – Executes workflows in the recommended order
# – Captures logs and creates patched YAMLs with envsubst
################################################################################

set -euo pipefail
shopt -s nullglob
cd "$(dirname "$0")"    # ensure we're inside my_local_test_data/

# 0.  Session log  ----------------------------------------------------------------
exec > >(tee -i FULL_RUN_$(date +%Y%m%d_%H%M).log) 2>&1

echo "=== FULL-SCALE WORKFLOW SESSION START $(date) ==="

# -----------------------------------------------------------------------------
# A. Generate synthetic BAMs matching each sample*_peaks.csv file --------------
# -----------------------------------------------------------------------------
# The BAM creation below is only necessary for ATAC workflows which we
# are not running in this DESeq-only execution.  Commenting it out avoids
# an intermittent segmentation fault observed with samtools on macOS.
: <<'DISABLE_BAM_GEN'
PEAK_DIR="core_data"
for peaks_csv in ${PEAK_DIR}/sample*_peaks.csv; do
  [ -e "$peaks_csv" ] || continue  # guard against glob failure
  bam_file="${peaks_csv%_peaks.csv}.bam"
  if [[ -f "$bam_file" && -f "${bam_file}.bai" ]]; then
    echo "[INFO] BAM already present for $(basename "$peaks_csv") – skipping"
    continue
  fi
  echo "[INFO] Generating synthetic BAM for $(basename "$peaks_csv")"
  python3 generate_synthetic_bams.py "$peaks_csv" "$bam_file"
done
DISABLE_BAM_GEN

# 1.  Duplicate & sanitise YAMLs  -------------------------------------------------
for f in \
  deseq_lrt_step_1/inputs/*testmode.yml \
  deseq_lrt_step_2/inputs/*testmode.yml \
  deseq_pairwise/inputs/*testmode.yml \
  atac_lrt_step_1/inputs/*testmode.yml \
  atac_lrt_step_2/inputs/*testmode.yml \
  atac_pairwise/inputs/*testmode.yml
do
  full="${f/_testmode/_full}"
  cp "$f" "$full"
  # Remove the test_mode line
  sed -i '' '/^[[:space:]]*test_mode:/d' "$full"
  # Ensure metadata file path line followed by required format field
  sed -i '' -E 's|(path:[[:space:]].*metadata\.csv.*)|\1\
  format: "http://edamontology.org/format_2330"|' "$full"

  # Inject placeholders into the two Step-2 YAMLs
  case "$full" in
    *deseq_lrt_step_2*)
      sed -i '' -E 's|path:[[:space:]].*_contrasts\.rds.*|  path: "${DESEQ_OBJ}"|' "$full"
      sed -i '' -E 's|path:[[:space:]].*_contrasts_table\.tsv.*|  path: "${DESEQ_CONTRASTS}"|' "$full"
      ;;
    *atac_lrt_step_2*)
      sed -i '' -E 's|path:[[:space:]].*_atac_obj\.rds.*|  path: "${ATAC_OBJ}"|' "$full"
      sed -i '' -E 's|path:[[:space:]].*_contrasts_table\.tsv.*|  path: "${ATAC_CONTRASTS}"|' "$full"
      ;;
  esac

  # Ensure required 'alias_trigger' key only for Step-2 YAMLs
  if [[ "$full" == *"step_2"* ]]; then
    sed -i '' 's/^alias:/alias_trigger:/' "$full"
  fi

  # For ATAC LRT Step 1 ensure lrt_only_mode flag is present and TRUE
  if [[ "$full" == *"atac_lrt_step_1"* ]]; then
    grep -q "^[[:space:]]*lrt_only_mode:" "$full" || echo "lrt_only_mode: true" >> "$full"
  fi
done
echo "YAML duplication + placeholder insertion complete."

# Quick fix for accidental duplicate 'format:' lines (schema_salad forbids dup keys)
# Remove exact duplicate adjacent lines or global duplicates in YAMLs.
for yaml in deseq_lrt_step_1/inputs/*_full.yml \
            deseq_lrt_step_2/inputs/*_full.yml \
            atac_lrt_step_1/inputs/*_full.yml \
            atac_lrt_step_2/inputs/*_full.yml \
            deseq_pairwise/inputs/*_full.yml \
            atac_pairwise/inputs/*_full.yml; do
  [ -f "$yaml" ] || continue
  awk 'NR==1{prev=""} {if($0==prev){next} print; prev=$0}' "$yaml" > "${yaml}.tmp" && mv "${yaml}.tmp" "$yaml"
done
echo "Duplicate key cleanup complete."

# Handle ATAC LRT Step 2 minimal_test.yml (has different naming convention)
# if [ -f atac_lrt_step_2/inputs/minimal_test.yml ]; then
#   cp atac_lrt_step_2/inputs/minimal_test.yml atac_lrt_step_2/inputs/minimal_full.yml
#   sed -i '' '/^[[:space:]]*test_mode:/d' atac_lrt_step_2/inputs/minimal_full.yml
#   sed -i '' -E 's|path:[[:space:]].*_atac_obj\.rds.*|  path: "${ATAC_OBJ}"|' atac_lrt_step_2/inputs/minimal_full.yml
#   sed -i '' -E 's|path:[[:space:]].*_contrasts_table\.tsv.*|  path: "${ATAC_CONTRASTS}"|' atac_lrt_step_2/inputs/minimal_full.yml
#   sed -i '' 's/^alias:/alias_trigger:/' atac_lrt_step_2/inputs/minimal_full.yml
# fi

# 2.  Create output directories  --------------------------------------------------
mkdir -p \
  deseq_lrt_step_1/outputs/full_run \
  deseq_lrt_step_2/outputs/full_run \
  deseq_pairwise/outputs/full_run \
  atac_lrt_step_1/outputs/full_run \
  atac_lrt_step_2/outputs/full_run \
  atac_pairwise/outputs/full_run

# 3.  Run Step-1 + pairwise workflows  -------------------------------------------
# (Uncomment --platform if you actually need x86-64 emulation)
# export DOCKER_DEFAULT_PLATFORM=linux/amd64

/usr/bin/time -l cwltool --debug --leave-container \
  --outdir deseq_lrt_step_1/outputs/full_run \
  ../workflows/deseq-lrt-step-1.cwl \
  deseq_lrt_step_1/inputs/deseq_lrt_s1_workflow_standard_full.yml \
  | tee deseq_lrt_step_1/outputs/full_run/run.log

# /usr/bin/time -l cwltool --debug --leave-container \
#   --outdir atac_lrt_step_1/outputs/full_run \
#   ../workflows/atac-lrt-step-1-test.cwl \
#   atac_lrt_step_1/inputs/atac_lrt_s1_workflow_interaction_full.yml \
#   | tee atac_lrt_step_1/outputs/full_run/run.log

/usr/bin/time -l cwltool --debug --leave-container \
  --outdir deseq_pairwise/outputs/full_run \
  ../workflows/deseq-pairwise.cwl \
  deseq_pairwise/inputs/deseq_pairwise_workflow_CMR_vs_KMR_full.yml \
  | tee deseq_pairwise/outputs/full_run/run.log

# /usr/bin/time -l cwltool --debug --leave-container \
#   --outdir atac_pairwise/outputs/full_run \
#   ../workflows/atac-pairwise.cwl \
#   atac_pairwise/inputs/atac_pairwise_workflow_rest_vs_active_full.yml \
#   | tee atac_pairwise/outputs/full_run/run.log

# 4.  Resolve real Step-1 artefact paths  ----------------------------------------
export DESEQ_OBJ=$(find deseq_lrt_step_1/outputs/full_run -name '*_contrasts.rds' | head -n1)
export DESEQ_CONTRASTS=$(find deseq_lrt_step_1/outputs/full_run -name '*_contrasts_table.tsv' | head -n1)
# export ATAC_OBJ=$(find atac_lrt_step_1/outputs/full_run -name '*_atac_obj.rds' | head -n1)
# export ATAC_CONTRASTS=$(find atac_lrt_step_1/outputs/full_run -name '*_contrasts_table.tsv' | head -n1)

# Updated artefact summary (ATAC disabled)
echo "Resolved artefacts:
  DESEQ_OBJ        = $DESEQ_OBJ
  DESEQ_CONTRASTS  = $DESEQ_CONTRASTS"

# 5.  Patch Step-2 YAMLs with envsubst  ------------------------------------------
envsubst < deseq_lrt_step_2/inputs/deseq_lrt_s2_workflow_dual_contrast_full.yml \
  > deseq_lrt_step_2/inputs/deseq_lrt_s2_workflow_dual_contrast_full_patched.yml

# envsubst < atac_lrt_step_2/inputs/minimal_full.yml \
#   > atac_lrt_step_2/inputs/atac_lrt_s2_workflow_full_patched.yml

# 6.  Run Step-2 workflows  -------------------------------------------------------
/usr/bin/time -l cwltool --debug --leave-container \
  --outdir deseq_lrt_step_2/outputs/full_run \
  ../workflows/deseq-lrt-step-2.cwl \
  deseq_lrt_step_2/inputs/deseq_lrt_s2_workflow_single_contrast.yml \
  | tee deseq_lrt_step_2/outputs/full_run/run.log

# /usr/bin/time -l cwltool --debug --leave-container \
#   --outdir atac_lrt_step_2/outputs/full_run \
#   ../workflows/atac-lrt-step-2-test.cwl \
#   atac_lrt_step_2/inputs/atac_lrt_s2_workflow_full_patched.yml \
#   | tee atac_lrt_step_2/outputs/full_run/run.log

echo "=== FULL-SCALE WORKFLOW SESSION END $(date) ===" 