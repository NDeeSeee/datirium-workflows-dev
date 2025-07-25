cwlVersion: v1.0
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement
sd:upstream:
  genelists_for_A:
  - filter-deseq-for-heatmap.cwl
  - filter-diffbind-for-heatmap.cwl
  genelists_for_B:
  - filter-deseq-for-heatmap.cwl
  - filter-diffbind-for-heatmap.cwl
inputs:
  alias:
    type: string
    label: Experiment short name/Alias
    sd:preview:
      position: 1
  set_operator:
    type:
    - 'null'
    - type: enum
      name: Set operation user selection
      symbols:
      - Intersection
      - Union
      - Relative_Complement
    label: Select set operation
    sd:preview:
      position: 3
    doc: "Set Examples (only scores from list A are reported):\n \tlist A = {1, 2, 3, 4};\n \tlist B = {3, 4, 5, 6};\n\n - Intersection: A ∩ B = {3, 4}\n - Union: A ∪ B = {1, 2, 3, 4, 5, 6}\n - Relative_Complement: A / B = {1, 2}"
    sd:localLabel: true
  filtered_list_A:
    type: File
    format: http://edamontology.org/format_3475
    label: 'Filtered genelist A (important to choose a specific list for relative complement):'
    doc: Select list A from filtered differential genelists from DESeq or diffbind pipelines
    sd:upstreamSource: genelists_for_A/filtered_file
    sd:localLabel: true
    sd:preview:
      position: 8
  filtered_list_B_group:
    type: File[]
    format: http://edamontology.org/format_3475
    label: 'Filtered genelist group B:'
    doc: Select 1 or more lists from filtered differential genelists from DESeq or diffbind pipelines. These will be used to generate intersections and unions among all lists (including list A), and will be applied against list A for the relative complement operation.
    sd:upstreamSource: genelists_for_B/filtered_file
    sd:localLabel: true
    sd:preview:
      position: 9
outputs:
  filtered_file:
    type: File
    format: http://edamontology.org/format_3003
    label: Filtered differentially expressed genes
    doc: Regions of interest formatted as headerless BED file with [chrom start end name score strand]
    outputSource: set_operation/genelist_filtered_set
  filtered_file_w_header:
    type: File
    format: http://edamontology.org/format_3003
    label: Filtered differentially expressed genes
    doc: Regions of interest formatted as headerless BED file with [chrom start end name score strand]
    outputSource: set_operation/genelist_filtered_set_with_header
    sd:visualPlugins:
    - syncfusiongrid:
        tab: Set results
        Title: Set table
  filtering_stdout_log_file:
    type: File
    format: http://edamontology.org/format_2330
    label: Filtering stdout log
    doc: Filtering stdout log
    outputSource: set_operation/log_file_stdout
  filtering_stderr_log_file:
    type: File
    format: http://edamontology.org/format_2330
    label: Filtering stderr log
    doc: Filtering stderr log
    outputSource: set_operation/log_file_stderr
steps:
  set_operation:
    run: ../tools/genelists-sets.cwl
    in:
      filtered_list_A: filtered_list_A
      filtered_list_B_group: filtered_list_B_group
      set_operation:
        source: set_operator
        valueFrom: $(self)
    out:
    - genelist_filtered_set
    - genelist_filtered_set_with_header
    - log_file_stdout
    - log_file_stderr
label: Set Operations for filtered genelists
doc: |-
  # Set Operations for filtered gene lists

  This workflow takes as input multiple filtered genelists samples and performs the user-selected set operation on them.
  There is one input for list A from which "scores" will be taken (these are fold change values from deseq or diffbind) and used in the output set list.
  The second genelist input is for 1+ genelists, that will be aggregated and used for intersection and union directly, and be applied against list A for the relative complement operation.
  The output is a single filtered gene list in the same format as the input files (headerless BED file with [chrom start end name score strand]).
  The returned score value (column 5) is always derived from file A.
sd:version: 100
