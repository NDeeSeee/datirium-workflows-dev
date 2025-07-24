## Environment Setup
- You need to run conda activate cwltool_env to run CWL tests

## Bioinformatics Workflow Mastery Patterns

### Core Breakthrough: Format Detection Architecture
- Critical Discovery: The "undefined columns selected" error pattern in DiffBind workflows is a universal bioinformatics problem - not specific to one workflow
- Root cause: hardcoded format assumptions vs. actual file formats
- Golden Knowledge Extraction: ATAC Workflow MACS2 Format Resolution discovered universal auto-detection solution for peak file formats

### Advanced Debugging Framework
- Resource vs Logic Separation: Test with minimal datasets (6-9 samples) to distinguish core logic errors from resource constraints
- Component Independence Testing: Validate each workflow step separately (Step 1 → Pairwise → Step 2)
- Log-Driven Diagnosis: stderr/stdout analysis reveals exact failure points
- Progressive Isolation Strategy: Systematically test workflow components, start with minimal datasets, and use log analysis to pinpoint exact failure locations

### Docker-CWL Integration Strategies
- Permission Fallback Strategy: Docker fails → Singularity (--singularity flag)
- Version Consistency Verification: Always cross-check CWL dockerPull vs available images
- Cache Optimization: --cachedir /tmp/cwl_cache essential for large bioinformatics datasets
- Additional Docker-CWL Integration Patterns:
  - Verify Docker image versions in CWL files
  - Use Singularity for permission issues
  - Implement careful cache management for large datasets

### Systematic Fix Propagation Strategy
- Identify Root Pattern: Format assumptions, parameter mismatches, etc.
- Systematic Search: Find ALL affected workflows (not just the failing one)
- Universal Implementation: Apply fix architecture across entire workflow suite
- Incremental Validation: Fix one → validate → propagate → test all
- Multi-Workflow Fix Strategy:
  1. Identify pattern of hardcoded format assumptions
  2. Apply systematic search across related workflows
  3. Implement universal auto-detection
  4. Test incrementally (one workflow → validate → apply to others)

### Error Decode Intelligence
- "undefined columns selected" error often indicates format mismatch (CSV ≠ MACS2)
- Recommended solution: Implement auto-detection logic for peak file formats
- Comprehensive Error Message Decoding:
  - "undefined columns selected" = Format mismatch
  - "Error processing one or more read files" = Resource constraints
  - "Did not find output file with glob pattern" = Missing summary generation

### Testing and Scaling Strategy
- Subset Validation: 6-9 samples (logic verification)
- Medium Scale: 12-15 samples (performance testing)
- Production Scale: Full dataset (final validation)
- Resource Management Insights:
  - Use 6-9 samples for initial validation
  - Ensure balanced metadata representation
  - Implement progressive scaling approach

### Key Strategic Insights
- Prioritize format detection over hardcoding
- Use progressive testing approach (start small, validate logic, then scale)
- Look for universal fix patterns applicable across workflows
- Rely on log analysis for precise error identification
- Maintain Docker fallback options (e.g., Singularity)
- Implement universal peak format detection function
- Always validate file format assumptions
- Maintain comprehensive documentation of workflow states and compatibility

### Future Workflow Development Guidelines
- Implement universal peak format detection mechanism
- Create auto-detection function for peak files
- Validate with multiple file formats and sample sizes
- Document error patterns and format compatibility
- Maintain transparent logging of format detection decisions