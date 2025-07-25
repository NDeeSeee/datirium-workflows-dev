cwlVersion: v1.0
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var split_features = function(line) { function get_unique(value, index, self) { return self.indexOf(value) === index && value != ""; } var splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null; return (splitted_line && !!splitted_line.length)?splitted_line:null; };
inputs:
  alias:
    type: string
    label: Experiment short name/alias
    sd:preview:
      position: 1
  srr_id:
    type: string
    label: Comma or space separated list of SRR Identifiers
    doc: |
      Comma or space separated list of SRR Identifiers
  splitby:
    type:
    - 'null'
    - type: enum
      symbols:
      - Split into all available files
      - 3-way splitting for mate-pairs
      - Do not split
    default: 3-way splitting for mate-pairs
    label: Split reads by
    doc: |-
      Split into all available files. Write reads into separate files. Read number will be suffixed to the file name. In cases where not all spots have the same number of reads, this option will produce files that WILL CAUSE ERRORS in most programs which process split pair fastq files.
      3-way splitting for mate-pairs. For each spot, if there are two biological reads satisfying filter conditions, the first is placed in the `*_1.fastq` file, and the second is placed in the `*_2.fastq` file. If there is only one biological read satisfying the filter conditions, it is placed in the `*.fastq` file. All other reads in the spot are ignored.
      Do not split. Output all reads into as a single FASTQ file
  http_proxy:
    type: string?
    label: Optional HTTP proxy settings
    doc: |
      Optional HTTP proxy settings
    sd:layout:
      advanced: true
  https_proxy:
    type: string?
    label: Optional HTTPS proxy settings
    doc: |
      Optional HTTPS proxy settings
    sd:layout:
      advanced: true
outputs:
  fastq_files:
    type:
    - 'null'
    - type: array
      items: File
    outputSource: fastq_dump/fastq_files
    label: Gzip-compressed FASTQ files
    doc: |
      Gzip-compressed FASTQ files
  report_md:
    type: File
    outputSource: fastq_dump/report_md
    label: Collected report for downloaded FASTQ files
    doc: |
      Collected report for downloaded FASTQ files
      in Markdown format
    sd:visualPlugins:
    - markdownView:
        tab: Overview
  metadata_xml:
    type:
    - 'null'
    - type: array
      items: File
    outputSource: fastq_dump/metadata_xml
    label: SRR metadata files in XML format
    doc: |
      SRR metadata files in XML format
  collected_metadata:
    type: File
    outputSource: fastq_dump/collected_metadata
    label: Collected metadata in TSV format
    doc: |
      Collected metadata in TSV format
  run_acc:
    type:
    - 'null'
    - type: array
      items: string
    outputSource: fastq_dump/run_acc
    label: Collected Run identifiers
    doc: |
      Collected Run identifiers
  experiment_acc:
    type:
    - 'null'
    - type: array
      items: string
    outputSource: fastq_dump/experiment_acc
    label: Collected Experiment identifiers
    doc: |
      Collected Experiment identifiers
  study_acc:
    type:
    - 'null'
    - type: array
      items: string
    outputSource: fastq_dump/study_acc
    label: Collected SRA Study identifiers
    doc: |
      Collected SRA Study identifiers
  biosample:
    type:
    - 'null'
    - type: array
      items: string
    outputSource: fastq_dump/biosample
    label: Collected BioSample identifiers
    doc: |
      Collected BioSample identifiers
  bioproject:
    type:
    - 'null'
    - type: array
      items: string
    outputSource: fastq_dump/bioproject
    label: Collected BioProject identifiers
    doc: |
      Collected BioProject identifiers
  fastq_dump_stdout:
    type: File
    outputSource: fastq_dump/log_stdout
    label: stdout log generated by fastq_dump step
    doc: |
      stdout log generated by fastq_dump step
  fastq_dump_stderr_log:
    type: File
    outputSource: fastq_dump/log_stderr
    label: stderr log generated by fastq_dump step
    doc: |
      stderr log generated by fastq_dump step
steps:
  fastq_dump:
    run: ../tools/fastq-dump.cwl
    in:
      srr_id:
        source: srr_id
        valueFrom: $(split_features(self))
      split_files:
        source: splitby
        valueFrom: $(self=="Split into all available files"?true:null)
      split_3:
        source: splitby
        valueFrom: $(self=="3-way splitting for mate-pairs"?true:null)
      http_proxy:
        source: http_proxy
        valueFrom: $(self==""?null:self)
      https_proxy:
        source: https_proxy
        valueFrom: $(self==""?null:self)
    out:
    - fastq_files
    - metadata_xml
    - collected_metadata
    - run_acc
    - experiment_acc
    - study_acc
    - biosample
    - bioproject
    - report_md
    - log_stdout
    - log_stderr
label: FASTQ Download
doc: |
  FASTQ Download

  Assists in downloading problematic single-cell sequencing
  data from Sequence Read Archive (SRA)
sd:version: 100
