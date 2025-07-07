cwlVersion: v1.0
class: Workflow

# -----------------------------------------------------------------------------
#                    REQUIREMENTS & GLOBAL METADATA
# -----------------------------------------------------------------------------
requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

"sd:upstream":
  atac_experiment:
    - "trim-atacseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-atacseq-pe.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-atacseq-se.cwl"

$namespaces:
  s: http://schema.org/

$schemas:
  - https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "ATAC-Seq LRT (step 1) – differential peak analysis using likelihood ratio test"
label: "ATAC-Seq LRT (step 1) – differential peak analysis using likelihood ratio test"
s:alternateName: "ATAC-Seq differential accessibility analysis based on LRT"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/atac-lrt-step-1.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

# -----------------------------------------------------------------------------
#                                 INPUTS
# -----------------------------------------------------------------------------
inputs:
  test_peak_files:
    type: File[]
    label: "Peak files"
    doc: "Merged/narrowPeak BED files for each sample"
    "sd:upstreamSource": "atac_experiment/peaks"
    "sd:localLabel": true

  peak_file_names:
    type: string[]
    label: "Peak file names"
    doc: "Aliases (sample IDs) corresponding to test_peak_files"
    "sd:upstreamSource": "atac_experiment/alias"

  bam_files:
    type: File[]
    label: "BAM files"
    doc: "Coordinate-sorted, indexed BAMs matching the peak files"
    "sd:upstreamSource": "atac_experiment/bams"

  metadata_file:
    type: File
    label: "Sample metadata table"
    doc: "CSV/TSV file describing experimental factors (time, condition, etc.)"

  design_formula:
    type: string
    label: "Design formula"
    doc: "Full design formula starting with ~ (e.g. ~ time + condition)"

  reduced_formula:
    type: string
    label: "Reduced formula"
    doc: "Reduced design formula with terms of interest removed"

  batchcorrection:
    type:
      - "null"
      - type: enum
        symbols: ["none", "combatseq", "model"]
    default: "none"
    label: "Batch correction method"

  scaling_type:
    type:
      - "null"
      - type: enum
        symbols: ["minmax", "zscore"]
    default: "zscore"
    label: "Expression data scaling"

  fdr:
    type: float?
    default: 0.1
    label: "FDR cutoff"

  lfcthreshold:
    type: float?
    default: 0.59
    label: "log2FC threshold"

  use_lfc_thresh:
    type: boolean
    default: false
    label: "Use lfcthreshold in null hypothesis"

  rpkm_cutoff:
    type: int?
    default: null
    label: "RPKM cutoff"

  cluster_method:
    type:
      - "null"
      - type: enum
        symbols: ["row", "column", "both", "none"]
    default: "none"
    label: "Hopach clustering method"

  row_distance:
    type:
      - "null"
      - type: enum
        symbols: ["cosangle", "abscosangle", "euclid", "cor", "abscor"]
    default: "cosangle"
    label: "HOPACH row distance"

  column_distance:
    type:
      - "null"
      - type: enum
        symbols: ["cosangle", "abscosangle", "euclid", "cor", "abscor"]
    default: "euclid"
    label: "HOPACH column distance"

  k_hopach:
    type: int?
    default: 3
    label: "HOPACH k"

  kmax_hopach:
    type: int?
    default: 5
    label: "HOPACH kmax"

  threads:
    type: int?
    default: 6
    label: "Threads"

  lrt_only_mode:
    type: boolean
    default: false
    label: "Run LRT only (skip contrasts)"

  test_mode:
    type: boolean
    default: false
    label: "Run only first 100 peaks for fast testing"

  output_prefix:
    type: string
    default: "atac_lrt_step_1"
    label: "Output filename prefix"

# -----------------------------------------------------------------------------
#                                 OUTPUTS
# -----------------------------------------------------------------------------
outputs:
  lrt_diff_expr:
    type: File?
    outputSource: atac_lrt_step_1/lrt_diff_expr
    "sd:visualPlugins":
      - syncfusiongrid:
          tab: "LRT results"
          Title: "Combined LRT results"

  contrasts_table:
    type: File?
    outputSource: atac_lrt_step_1/contrasts_table

  dsq_obj_data:
    type: File?
    outputSource: atac_lrt_step_1/dsq_obj_data

  mds_plots_html:
    type: File?
    outputSource: atac_lrt_step_1/mds_plots_html
    "sd:visualPlugins":
      - linkList:
          tab: "Overview"
          target: "_blank"

  mds_plots_corrected_html:
    type: File?
    outputSource: atac_lrt_step_1/mds_plots_corrected_html
    "sd:visualPlugins":
      - linkList:
          tab: "Overview"
          target: "_blank"

  read_counts_file_all:
    type: File?
    outputSource: atac_lrt_step_1/counts_all_gct

  read_counts_file_filtered:
    type: File?
    outputSource: atac_lrt_step_1/counts_filtered_gct

  heatmap_html:
    type: File?
    outputSource: morpheus_heatmap/heatmap_html
    "sd:visualPlugins":
      - linkList:
          tab: "Overview"
          target: "_blank"

  alignment_stats_barchart:
    type: File?
    outputSource: atac_lrt_step_1/alignment_stats_barchart

  stdout_log:
    type: File
    outputSource: atac_lrt_step_1/stdout_log

  stderr_log:
    type: File
    outputSource: atac_lrt_step_1/stderr_log

  morpheus_stdout_log:
    type: File
    outputSource: morpheus_heatmap/stdout_log

  morpheus_stderr_log:
    type: File
    outputSource: morpheus_heatmap/stderr_log

# -----------------------------------------------------------------------------
#                                   STEPS
# -----------------------------------------------------------------------------
steps:
  atac_lrt_step_1:
    run: ../tools/atac-lrt-step-1.cwl
    in:
      test_peak_files: test_peak_files
      peak_file_names: peak_file_names
      bam_files: bam_files
      metadata_file: metadata_file
      design_formula: design_formula
      reduced_formula: reduced_formula
      batchcorrection: batchcorrection
      scaling_type: scaling_type
      fdr: fdr
      lfcthreshold: lfcthreshold
      use_lfc_thresh: use_lfc_thresh
      rpkm_cutoff: rpkm_cutoff
      cluster_method: cluster_method
      row_distance: row_distance
      column_distance: column_distance
      k_hopach: k_hopach
      kmax_hopach: kmax_hopach
      output_prefix: output_prefix
      threads: threads
      lrt_only_mode: lrt_only_mode
      test_mode: test_mode
    out:
      - lrt_diff_expr
      - contrasts_table
      - dsq_obj_data
      - mds_plots_html
      - mds_plots_corrected_html
      - counts_all_gct
      - counts_filtered_gct
      - alignment_stats_barchart
      - stdout_log
      - stderr_log

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
      read_counts_gct: atac_lrt_step_1/counts_filtered_gct
    out:
      - heatmap_html
      - stdout_log
      - stderr_log

# -----------------------------------------------------------------------------
#                                 DOCUMENTATION
# -----------------------------------------------------------------------------

doc: |
  Runs ATAC-Seq differential accessibility analysis using DESeq2 LRT
  ==================================================================
  The workflow wraps the `atac-lrt-step-1.cwl` CommandLineTool and adds a
  Morpheus heat-map for exploratory visualisation, bringing feature-parity with
  the DESeq RNA-seq LRT workflow.

  * At least two biological replicates are required per experimental group.
  * The metadata file must include sample IDs matching **Peak file names**. 