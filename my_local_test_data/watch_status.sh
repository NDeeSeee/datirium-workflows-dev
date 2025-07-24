#!/usr/bin/env bash
while true; do
    printf "\n$(date)\n"
    grep -E '^\[[0-9.]+ Workflow|Completed|Exception' \
        my_local_test_data/atac_pairwise/full_run_outputs/run.log | tail -n 5
    grep -E '^\[[0-9.]+ Workflow|Completed|Exception' \
        my_local_test_data/atac_lrt_step_1/outputs/full_run/run.log | tail -n 5
    sleep 30
done