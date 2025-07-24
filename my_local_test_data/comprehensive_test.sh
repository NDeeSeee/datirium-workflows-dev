#!/bin/bash
# Comprehensive test script for all 6 CWL workflows
# Must be run from my_local_test_data directory
# Usage: cd my_local_test_data && ./comprehensive_test.sh

set -e  # Exit on any error

# Use native ARM64 Docker images for optimal performance
# Performance benefits: ~3x faster execution, reduced memory usage, better battery life
echo "Using native ARM64 Docker images for optimal performance (~3x faster)"
unset DOCKER_DEFAULT_PLATFORM
# Check if we're in the right directory
if [[ ! -f "../tools/deseq-lrt-step-1.cwl" ]]; then
    echo "❌ Error: Must run from my_local_test_data directory"
    echo "Usage: cd my_local_test_data && ./comprehensive_test.sh"
    exit 1
fi

echo "=========================================="
echo "🧪 COMPREHENSIVE TEST - All 6 Workflows"
echo "=========================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test configurations (name:cwl:input:output)
workflows=(
  "DESeq LRT Step 1:../tools/deseq-lrt-step-1.cwl:deseq_lrt_step_1/inputs/deseq_lrt_s1_workflow_standard_testmode.yml:deseq_lrt_step_1/outputs/comprehensive_test"
  "DESeq LRT Step 2:../tools/deseq-lrt-step-2.cwl:deseq_lrt_step_2/inputs/deseq_lrt_s2_workflow_single_contrast.yml:deseq_lrt_step_2/outputs/comprehensive_test"
  "DESeq Pairwise:../tools/deseq-pairwise.cwl:deseq_pairwise/inputs/deseq_pairwise_workflow_CMR_vs_KMR_testmode.yml:deseq_pairwise/outputs/comprehensive_test"
  "ATAC LRT Step 1:../tools/atac-lrt-step-1.cwl:atac_lrt_step_1/inputs/atac_lrt_s1_workflow_interaction_testmode.yml:atac_lrt_step_1/outputs/comprehensive_test"
  "ATAC LRT Step 2:../tools/atac-lrt-step-2.cwl:atac_lrt_step_2/inputs/atac_lrt_s2_workflow_single_contrast.yml:atac_lrt_step_2/outputs/comprehensive_test"
  "ATAC Pairwise:../tools/atac-pairwise.cwl:atac_pairwise/inputs/atac_pairwise_workflow_rest_vs_active.yml:atac_pairwise/outputs/comprehensive_test"
)

success_count=0
total_count=6
failed_tests=()

for workflow in "${workflows[@]}"; do
  IFS=":" read -r name cwl input output <<< "$workflow"
  echo -e "\n${YELLOW}Testing: $name${NC}"
  
  # Create output directory
  mkdir -p "$output"
  
  # Run the test
  if cwltool --outdir "$output" "$cwl" "$input" > "$output/test.log" 2>&1; then
    echo -e "${GREEN}✓ PASSED: $name${NC}"
    ((success_count++))
  else
    echo -e "${RED}✗ FAILED: $name${NC}"
    echo -e "${RED}  Check log: $output/test.log${NC}"
    failed_tests+=("$name")
  fi
done

# Summary
echo -e "\n=========================================="
echo -e "${YELLOW}COMPREHENSIVE TEST SUMMARY${NC}"
echo "=========================================="
echo -e "Tests passed: ${GREEN}$success_count${NC}"
echo -e "Tests failed: ${RED}$((total_count - success_count))${NC}"

if [ $success_count -lt $total_count ]; then
    echo -e "\n${RED}Failed tests:${NC}"
    for test in "${failed_tests[@]}"; do
        echo -e "  - $test"
    done
    echo -e "\n${YELLOW}💡 Tip: Check individual log files for details${NC}"
    exit 1
else
    echo -e "\n${GREEN}🎉 All workflows working perfectly!${NC}"
    exit 0
fi