#!/usr/bin/env bash
# Orchestrator: run all ATAC workflows in full mode with asynchronous execution
# 1. ATAC LRT Step 1
# 2. ATAC Pairwise (runs in parallel with Step 1)
# 3. ATAC LRT Step 2 (runs after Step 1 completes and consumes its artefacts)
# All logs are captured to run.log in the respective output directories.

set -euo pipefail

# Ensure AMD64 containers on Apple Silicon
export DOCKER_DEFAULT_PLATFORM=linux/amd64

# Navigate to the script directory to keep paths predictable
cd "$(dirname "$0")"

# Create fresh output directories
mkdir -p \
  atac_lrt_step_1/outputs/full_run \
  atac_pairwise/outputs/full_run \
  atac_lrt_step_2/outputs/full_run

########################################
# Step 1 – ATAC LRT Step 1 (starts now)
########################################
(
  cwltool --debug --leave-container \
    --outdir atac_lrt_step_1/outputs/full_run \
    ../workflows/atac-lrt-step-1-test.cwl \
    atac_lrt_step_1/inputs/atac_lrt_s1_workflow_interaction_full.yml \
    | tee atac_lrt_step_1/outputs/full_run/run.log
) &
STEP1_PID=$!

echo "[INFO] ATAC LRT Step 1 launched (PID=$STEP1_PID)"

########################################
# Step 1b – ATAC Pairwise (parallel)
########################################
(
  cwltool --debug --leave-container \
    --outdir atac_pairwise/outputs/full_run \
    ../workflows/atac-pairwise.cwl \
    atac_pairwise/inputs/atac_pairwise_workflow_rest_vs_active_full.yml \
    | tee atac_pairwise/outputs/full_run/run.log
) &
PAIR_PID=$!

echo "[INFO] ATAC Pairwise launched (PID=$PAIR_PID)"

########################################
# Wait for Step 1 to finish
########################################
wait "$STEP1_PID" || { echo "[ERROR] ATAC LRT Step 1 failed – aborting."; exit 1; }

echo "[INFO] ATAC LRT Step 1 completed successfully. Preparing inputs for Step 2."

# Locate artefacts produced by Step 1
ATAC_OBJ=$(find atac_lrt_step_1/outputs/full_run -name '*_contrasts.rds' | head -n1 || true)
ATAC_CONTRASTS=$(find atac_lrt_step_1/outputs/full_run -name '*_contrasts_table.tsv' | head -n1 || true)

if [[ -z "$ATAC_OBJ" || -z "$ATAC_CONTRASTS" ]]; then
  echo "[ERROR] Required artefacts not found. Expected *_contrasts.rds and *_contrasts_table.tsv." >&2
  exit 1
fi

# Compose the Step 2 input YAML dynamically
cat > atac_lrt_step_2/inputs/atac_lrt_s2_full.yml <<YAML
dsq_obj_data:
  class: File
  path: "$ATAC_OBJ"

contrasts_table:
  class: File
  path: "$ATAC_CONTRASTS"

contrast_indices: "1"

# Disable test mode for full analysis
test_mode: false

output_prefix: "atac_lrt_step2_full"
YAML

echo "[INFO] Generated Step 2 job file: atac_lrt_step_2/inputs/atac_lrt_s2_full.yml"

########################################
# Step 2 – ATAC LRT Step 2
########################################
(
  cwltool --debug --leave-container \
    --outdir atac_lrt_step_2/outputs/full_run \
    ../workflows/atac-lrt-step-2-test.cwl \
    atac_lrt_step_2/inputs/atac_lrt_s2_full.yml \
    | tee atac_lrt_step_2/outputs/full_run/run.log
) &
STEP2_PID=$!

echo "[INFO] ATAC LRT Step 2 launched (PID=$STEP2_PID)"

########################################
# Wait for remaining jobs
########################################
wait "$PAIR_PID" || { echo "[ERROR] ATAC Pairwise failed."; exit 1; }
wait "$STEP2_PID" || { echo "[ERROR] ATAC LRT Step 2 failed."; exit 1; }

echo "[SUCCESS] ALL ATAC WORKFLOWS FINISHED SUCCESSFULLY" 