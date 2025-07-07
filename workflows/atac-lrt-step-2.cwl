cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement
  - class: ScatterFeatureRequirement

"sd:upstream":
  atac_lrt_step_1:
    - "atac-lrt-step-1.cwl"
    - "https://github.com/datirium/workflows/workflows/atac-lrt-step-1.cwl"

inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  contrast_indices:
    type: string
    label: "Comma-separated list of integers representing contrast indices"
    default: "1,2,5"

  fdr_cutoff:
    type: float?
    default: 0.1
    label: "FDR cutoff for significance"
    'sd:layout':
      advanced: true

  lfcthreshold:
    type: float?
    default: 0.59
    label: "Log2 Fold Change threshold"
    'sd:layout':
      advanced: true

  use_lfc_thresh:
    type: boolean
    default: false
    label: "Use lfcthreshold in null hypothesis"
    'sd:layout':
      advanced: true

  cluster_method:
    type:
      - "null"
      - type: enum
        symbols: ["row", "column", "both", "none"]
    default: "none"
    label: "Hopach clustering method"
    'sd:layout':
      advanced: true

  scaling_type:
    type:
      - "null"
      - type: enum
        symbols: ["minmax", "zscore"]
    default: "zscore"
    label: "Scaling method"
    'sd:layout':
      advanced: true

  row_distance:
    type:
      - "null"
      - type: enum
        symbols: ["cosangle", "abscosangle", "euclid", "cor", "abscor"]
    default: "cosangle"
    label: "Row distance metric"
    'sd:layout':
      advanced: true

  column_distance:
    type:
      - "null"
      - type: enum
        symbols: ["cosangle", "abscosangle", "euclid", "cor", "abscor"]
    default: "euclid"
    label: "Column distance metric"
    'sd:layout':
      advanced: true

  k_hopach:
    type: int?
    default: 3
    label: "Hopach k"
    'sd:layout':
      advanced: true

  kmax_hopach:
    type: int?
    default: 5
    label: "Hopach kmax"
    'sd:layout':
      advanced: true

  regulation:
    type:
      - "null"
      - type: enum
        symbols: ["both", "up", "down"]
    default: "both"
    label: "Direction of differential accessibility"
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 6
    label: "Threads"
    'sd:layout':
      advanced: true

  test_mode:
    type: boolean
    default: false
    label: "Run test mode (top 100 peaks)"
    'sd:layout':
      advanced: true

  # Inputs sourced from step 1
  dsq_obj_data:
    type: File?
    label: "RDS contrasts object"
    "sd:upstreamSource": "atac_lrt_step_1/dsq_obj_data"

  contrasts_table:
    type: File?
    label: "Contrasts table TSV"
    "sd:upstreamSource": "atac_lrt_step_1/contrasts_table"

  batchcorrection:
    type:
      - "null"
      - type: enum
        symbols: ["none", "combatseq", "model"]
    default: "none"
    label: "Batch correction method"
    "sd:upstreamSource": "atac_lrt_step_1/batchcorrection"

outputs:

  diff_expr_files:
    type: File[]
    outputSource: atac_step_2/diff_expr_files
    "sd:visualPlugins":
      - syncfusiongrid:
          tab: 'Differential Accessibility Analysis'
          Title: 'Combined ATAC-Seq results'

  read_counts_file_all:
    type: File
    outputSource: atac_step_2/counts_all_gct

  read_counts_file_filtered:
    type: File
    outputSource: atac_step_2/counts_filtered_gct

  mds_plots_html:
    type: File
    outputSource: atac_step_2/mds_plots_html
    "sd:visualPlugins":
      - linkList:
          tab: 'Overview'
          target: '_blank'

  volcano_plots_html:
    type: File[]
    outputSource: make_volcano_plot/html_file
    "sd:visualPlugins":
      - linkList:
          tab: 'Overview'
          target: '_blank'

  heatmap_html:
    type: File
    outputSource: morpheus_heatmap/heatmap_html
    "sd:visualPlugins":
      - linkList:
          tab: 'Overview'
          target: '_blank'

  atac_stdout_log:
    type: File
    outputSource: atac_step_2/stdout_log

  atac_stderr_log:
    type: File
    outputSource: atac_step_2/stderr_log

  morpheus_stdout_log:
    type: File
    outputSource: morpheus_heatmap/stdout_log

  morpheus_stderr_log:
    type: File
    outputSource: morpheus_heatmap/stderr_log

steps:

  atac_step_2:
    run: ../tools/atac-lrt-step-2.cwl
    in:
      dsq_obj_data: dsq_obj_data
      contrasts_table: contrasts_table
      batchcorrection: batchcorrection
      contrast_indices: contrast_indices
      fdr_cutoff: fdr_cutoff
      lfcthreshold: lfcthreshold
      use_lfc_thresh: use_lfc_thresh
      row_distance: row_distance
      column_distance: column_distance
      cluster_method: cluster_method
      scaling_type: scaling_type
      k_hopach: k_hopach
      kmax_hopach: kmax_hopach
      regulation: regulation
      output_prefix: alias
      threads: threads
      test_mode: test_mode
    out:
      - diff_expr_files
      - mds_plots_html
      - counts_all_gct
      - counts_filtered_gct
      - stdout_log
      - stderr_log

  make_volcano_plot:
    run: ../tools/volcano-plot.cwl
    scatterMethod: dotproduct
    scatter:
      - diff_expr_file
    in:
      diff_expr_file: atac_step_2/diff_expr_files
      output_filename:
        valueFrom: $(inputs.diff_expr_file.basename.replace(/\.tsv$/, '.html'))
      x_axis_column:
        default: "log2FoldChange"
      y_axis_column:
        default: "padj"
      label_column:
        default: "GeneId"
    out:
      - html_file
      - html_data

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
      read_counts_gct: atac_step_2/counts_filtered_gct
    out:
      - heatmap_html
      - stdout_log
      - stderr_log

$namespaces:
  s: http://schema.org/

$schemas:
  - https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "ATAC-Seq LRT (step 2) - Differential accessibility analysis using contrasts"
label: "ATAC-Seq LRT (step 2) - Differential accessibility analysis using contrasts"
s:alternateName: "ATAC-Seq differential accessibility analysis based on LRT contrasts"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/atac-lrt-step-2.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
  - class: s:Organization
    s:legalName: "Cincinnati Children's Hospital Medical Center"
    s:location:
      - class: s:PostalAddress
        s:addressCountry: "USA"
        s:addressLocality: "Cincinnati"
        s:addressRegion: "OH"
        s:postalCode: "45229"
        s:streetAddress: "3333 Burnet Ave"
        s:telephone: "+1(513)636-4200"
    s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
    s:department:
      - class: s:Organization
        s:legalName: "Allergy and Immunology"
        s:department:
          - class: s:Organization
            s:legalName: "Barski Research Lab"
            s:member:
              - class: s:Person
                s:name: Michael Kotliar
                s:email: mailto:misha.kotliar@gmail.com
                s:sameAs:
                  - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Runs ATAC-Seq differential accessibility analysis using contrasts derived from LRT step 1.
  Produces per-contrast differential tables, volcano plots, and combined heatmap. 