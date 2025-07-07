cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement

inputs:
  untreated_files:
    type: File[]
    label: "Untreated files"

  treated_files:
    type: File[]
    label: "Treated files"

  treated_sample_names:
    type:
      - "null"
      - string[]
    label: "Treated sample aliases"

  untreated_name:
    type: string
    default: "Control"
    label: "Untreated name"

  untreated_sample_names:
    type:
      - "null"
      - string[]
    label: "Untreated sample aliases"

  treated_name:
    type: string
    default: "Treatment"
    label: "Treated name"

  test_mode:
    type: boolean
    default: true

  output_prefix:
    type: string
    default: "atac_pairwise"

outputs:
  diff_expr_file:
    type: File?
    outputSource: atac_pairwise/diff_expr_file

  stdout_log:
    type: File
    outputSource: atac_pairwise/stdout_log

steps:
  atac_pairwise:
    run: ../tools/atac-pairwise.cwl
    in:
      untreated_files: untreated_files
      treated_files: treated_files
      untreated_name: untreated_name
      treated_name: treated_name
      untreated_sample_names: untreated_sample_names
      treated_sample_names: treated_sample_names
      test_mode: test_mode
      output_prefix: output_prefix
    out:
      - diff_expr_file
      - stdout_log
