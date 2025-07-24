# CWL Workflow Testing Environment

## CRITICAL RULES FOR TESTING
⚠️ **NEVER create files outside `my_local_test_data/` directory**  
⚠️ **ALWAYS check if files exist before creating new ones**  
⚠️ **MUST rebuild Docker images after script changes**

## Testing Pipeline Overview

The testing workflow follows this sequence:
1. **Run CWL tests** → Get errors
2. **Update core scripts** in `../tools/dockerfiles/scripts/`
3. **Rebuild Docker images** (deseq: `biowardrobe2/scidap-deseq`, atacseq: `biowardrobe2/scidap-atac`)
4. **Re-run tests** → Repeat until success

## Directory Structure

```
my_local_test_data/
├── README.md                    # This file - main testing documentation
├── core_data/                   # Shared test datasets
│   ├── *.isoforms.csv          # Expression data files
│   ├── metadata.csv            # Sample metadata
│   └── batch_file.csv          # Batch correction info
├── deseq_lrt_step_1/           # DESeq2 LRT Step 1 tests
│   ├── inputs/                 # Test input files (.yml configs)
│   ├── outputs/                # Test outputs (git-ignored)
│   └── scripts/                # Test helper scripts
├── deseq_lrt_step_2/           # DESeq2 LRT Step 2 tests
├── deseq_standard/             # Standard DESeq2 tests
├── atac_lrt_step_1/            # ATAC-seq LRT Step 1 tests
├── atac_lrt_step_2/            # ATAC-seq LRT Step 2 tests
├── atac_standard/              # Standard ATAC-seq tests
├── quick_test.sh               # Fast individual test runner
└── run_all_tests.sh            # Complete test suite runner
```

## Quick Start Testing

### 1. CWL Tool Testing (Start Here)
```bash
# Validate CWL syntax first
cwltool --validate ../tools/deseq-lrt-step-1.cwl
cwltool --validate ../tools/atac-lrt-step-1.cwl

# Test individual CWL tools
cd my_local_test_data/deseq_lrt_step_1
cwltool --debug ../../tools/deseq-lrt-step-1.cwl inputs/basic_test.yml
```

### 2. Workflow Testing
```bash
# Test complete workflows
cwltool --debug ../../workflows/deseq-lrt-step-1-test.cwl inputs/basic_test.yml
```

### 3. Docker Image Management
```bash
# Check current images
docker images | grep -E "(scidap-deseq|scidap-atac)"

# Rebuild after script changes (requires repository access)
# This should be done through CI/CD pipeline
```

## Test Data Requirements

### DESeq2 Tests
- **Expression files**: `*.isoforms.csv` with columns: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand, TotalReads, Rpkm
- **Metadata**: CSV with columns: sampleID, treatment, cond
- **8 samples**: 2 treatments × 2 conditions × 2 replicates

### ATAC-seq Tests  
- **Peak files**: `*_peaks.csv` with peak coordinates
- **BAM files**: `*.bam` alignment files
- **Metadata**: CSV with sample information
- **4+ samples**: Sufficient for differential accessibility analysis

## Integrated Fixes (All Applied to Main Scripts)

### ✅ ATAC Pairwise - Missing Library Import
**Problem**: `Error: could not find function "ArgumentParser"`
**Solution**: Added `library(argparse)` to `tools/dockerfiles/scripts/functions/atac_pairwise/cli_args.R`
**Status**: ✅ **INTEGRATED**

### ✅ ATAC LRT Step 2 - S4 Object Accessor  
**Problem**: `no applicable method for '@' applied to an object of class "DESeqResults"`
**Solution**: Fixed column access order in `tools/dockerfiles/scripts/functions/atac_lrt_step_2/contrast_analysis.R`
**Status**: ✅ **INTEGRATED**

### ✅ DESeq Pairwise - Variable Shadowing
**Problem**: `Error in args$untreated_sample_names : object of type 'closure' is not subsettable`
**Solution**: Renamed parameter to avoid conflict in `tools/dockerfiles/scripts/functions/deseq/cli_args.R`
**Status**: ✅ **INTEGRATED**

## Common Issues & Solutions

### Error: "File not found"
- ✅ **Fix**: Update file paths in `.yml` input files to use absolute paths
- ✅ **Check**: Ensure test data files exist in expected locations

### Error: "Docker image not found"
- ✅ **Fix**: Rebuild Docker images or use script mounting for development
- ✅ **Check**: `docker images | grep scidap`

### Error: "R script failure"
- ✅ **Fix**: Update scripts in `../tools/dockerfiles/scripts/`
- ✅ **Must**: Rebuild Docker image after script changes
- ✅ **Test**: Mount scripts during development to avoid rebuilds

### Error: "CLI argument parsing"
- ✅ **Fix**: Update `cli_args.R` functions for proper CWL argument handling
- ✅ **Check**: Boolean flags need special handling for CWL compatibility

### Error: "YAML Syntax Errors"
**Common patterns and fixes**:
```yaml
# WRONG (missing space after colon)
hints:
- class: DockerRequirement
    dockerPull: "image:tag"

# CORRECT (proper indentation)
hints:
  - class: DockerRequirement
    dockerPull: "image:tag"
```

## Development Workflow

### Script Development Mode (Fast)
```bash
# Mount updated scripts to avoid Docker rebuilds
docker run --rm \
  -v "$(pwd)/../tools/dockerfiles/scripts:/usr/local/bin" \
  -v "$(pwd)/inputs:/data" \
  biowardrobe2/scidap-deseq:latest \
  Rscript /usr/local/bin/run_deseq_lrt_step_1.R --args...
```

### Production Mode (Complete)
1. **Update scripts** in `../tools/dockerfiles/scripts/`
2. **Commit changes** to trigger Docker rebuild
3. **Test with new image** using CWL workflows
4. **Validate outputs** are scientifically correct

## Testing Commands Reference

### Individual Tests
```bash
# Quick DESeq2 test
./quick_test.sh deseq_lrt_step_1

# Quick ATAC-seq test  
./quick_test.sh atac_lrt_step_1

# All tests
./run_all_tests.sh
```

### Debug Mode
```bash
# Full debug output
CWLTOOL_OPTS="--debug" ./run_all_tests.sh

# Single test with debug
cwltool --debug ../../tools/deseq-lrt-step-1.cwl inputs/basic_test.yml
```

## Success Criteria

### For Each Test:
- ✅ CWL validation passes
- ✅ Docker container runs without errors  
- ✅ Expected output files are generated
- ✅ Log files show no unexpected errors
- ✅ Results are scientifically plausible

### For Complete Pipeline:
- ✅ All workflow types function correctly
- ✅ `test_mode=true` reduces runtime significantly  
- ✅ Different parameter combinations work
- ✅ Batch correction scenarios function properly

## File Management Rules

### ✅ ALLOWED:
- Update existing `.yml` input files
- Create new test cases in existing structure
- Add output directories (git-ignored)
- Update helper scripts within `my_local_test_data/`

### ❌ FORBIDDEN:
- Create files outside `my_local_test_data/`
- Create new documentation files in root directory
- Modify git-tracked files without explicit approval
- Create temporary directories in root

## File Modification Priority Order (CRITICAL)
1. **R Scripts** (`../tools/dockerfiles/scripts/`) - **Highest priority**, core logic
2. **Function Libraries** (`../tools/dockerfiles/scripts/functions/`) - **High priority**, modular components
3. **CWL Tools** (`../tools/*.cwl`) - Medium priority, interface definitions
4. **CWL Workflows** (`../workflows/*.cwl`) - Lower priority, orchestration only
5. **Docker configurations** - Only when necessary for environment issues

## Scientific Validation Requirements
### Statistical Methods
- **Always verify** DESeq2 parameters are appropriate for experimental design
- **Check normalization** methods are suitable for data type
- **Validate multiple testing correction** is properly applied
- **Ensure adequate replication** for statistical power

### Result Interpretation
- **Log fold change thresholds** should be biologically meaningful
- **P-value cutoffs** should be appropriate for discovery vs validation
- **Batch effects** should be properly corrected when present
- **Interaction terms** should be interpreted correctly in LRT context

## Cost Optimization Rules
### Context Management
- **Always check .cursorignore** before starting conversations
- **Exclude test results** and temporary files from context
- **Focus on relevant files only** for current issue
- **Use parallel tool calls** when reading multiple files

### Efficient Development
- Use `test_mode=true` to reduce processing time
- Mount scripts during development to avoid Docker rebuilds
- Run validation checks before full workflow execution
- Cache Docker layers effectively

## Git Workflow Rules
### Branching Strategy
```bash
# Create feature branch for each major change
git checkout -b fix-deseq-lrt-step1-yaml-syntax

# Commit frequently with descriptive messages
git commit -m "fix: correct YAML indentation in deseq-lrt-step-1.cwl line 9"

# Merge only after complete validation
git checkout master && git merge fix-deseq-lrt-step1-yaml-syntax
```

### Commit Message Format
- `fix:` brief description of what was fixed
- `feat:` brief description of new feature  
- `refactor:` brief description of code restructuring
- `docs:` brief description of documentation changes
- `test:` brief description of testing changes

## Emergency Procedures
### If Tests Keep Failing
1. **Step back to validate basics**: CWL syntax, Docker availability, file paths
2. **Test individual components**: R scripts, functions, data formats
3. **Use minimal test data** to isolate issues
4. **Create clean feature branch** and start systematic fixes

### If Docker Issues Persist
1. **Check image availability**: `docker images | grep scidap`
2. **Try manual docker run** to test script execution
3. **Mount local scripts** to test without rebuilding
4. **Verify CI/CD pipeline** is functioning correctly

### If Context Becomes Too Large
1. **Immediately check .cursorignore** and add exclusions
2. **Focus on single file** at a time
3. **Use file_search instead of reading entire files**
4. **Break conversation and restart** with optimized context

## Current Testing Status (Last Updated: 2025-06-16 - Comprehensive Testing Complete)

### ✅ CONFIRMED WORKING (3/6):
- **DESeq2 LRT Step 1**: Complete differential expression analysis (151s)
  - Image: `biowardrobe2/scidap-deseq:v0.0.58`
  - Status: Full production functionality
  
- **DESeq2 LRT Step 2**: Multi-contrast analysis with visualizations
  - Image: `biowardrobe2/scidap-deseq:v0.0.58`
  - Status: Fixed missing input file, now functional
  
- **ATAC LRT Step 1**: ATAC-seq differential accessibility analysis (20s)
  - Image: `biowardrobe2/scidap-atac:v0.0.62-fixed`
  - Status: Test mode bypass functional

### ❌ SYSTEMATIC R SCRIPT ISSUES (3/6):
- **DESeq Pairwise**: R execution successful but missing `*summary.md` output
- **ATAC Pairwise**: R execution successful but missing `*summary.md` output  
- **ATAC LRT Step 2**: Docker stats parsing errors + output file issues

### 🔧 INFRASTRUCTURE FIXES APPLIED:
- Created missing `basic_test.yml` files in test directories
- Verified Docker image consistency with CWL tool definitions
- Updated to latest image versions (DESeq v0.0.58, ATAC v0.0.62-fixed)

### 📋 NEXT STEPS:
- Debug R output generation functions in pairwise workflows
- Investigate ATAC LRT Step 2 Docker stats parsing issues
- Focus on `tools/dockerfiles/scripts/functions/` R libraries for output file generation

---
**Remember**: Always start with CWL tool testing, then proceed to workflow testing. Update scripts → rebuild Docker → test → repeat.