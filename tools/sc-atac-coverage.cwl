cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.41
inputs:
  query_data_rds:
    type: File
    inputBinding:
      prefix: --query
    doc: |
      Path to the RDS file to load Seurat object from. This file
      should include chromatin accessibility information stored
      in the ATAC assay with a proper seqinfo data.
  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    inputBinding:
      prefix: --fragments
    doc: |
      Count and barcode information for every ATAC fragment used in the
      loaded Seurat object. File should be saved in TSV format and to be
      tbi-indexed.
  splitby:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: --splitby
    doc: |
      Column from the Seurat object metadata to split cells into groups.
      May be one of the columns added with --metadata or --barcodes
      parameters. Default: split by dataset
  datasets_metadata:
    type: File?
    inputBinding:
      prefix: --metadata
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat object metadata with
      categorical values using samples identities. First column - 'library_id'
      should correspond to all unique values from the 'new.ident' column of the
      loaded Seurat object. If any of the provided in this file columns are already
      present in the Seurat object metadata, they will be overwritten. When combined
      with --barcodes parameter, first the metadata will be extended, then barcode
      filtering will be applied. Default: no extra metadata is added
  barcodes_data:
    type: File?
    inputBinding:
      prefix: --barcodes
    doc: |
      Path to the TSV/CSV file to optionally prefilter and extend Seurat object
      metadata be selected barcodes. First column should be named as 'barcode'.
      If file includes any other columns they will be added to the Seurat object
      metadata ovewriting the existing ones if those are present.
      Default: all cells used, no extra metadata is added
  flank_distance:
    type: int?
    inputBinding:
      prefix: --flank
    doc: |
      Distance in bp to flank both start and end of the each fragment in both
      direction to generate cut sites coverage. Default: 5
  verbose:
    type: boolean?
    inputBinding:
      prefix: --verbose
    doc: |
      Print debug information.
      Default: false
  export_html_report:
    type: boolean?
    default: false
    doc: |
      Export tehcnical report. HTML format.
      Note, stdout will be less informative.
      Default: false
  output_prefix:
    type: string?
    inputBinding:
      prefix: --output
    doc: |
      Output prefix.
      Default: ./sc
  parallel_memory_limit:
    type: int?
    inputBinding:
      prefix: --memory
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32
  vector_memory_limit:
    type: int?
    default: 128
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Default: 128
  threads:
    type: int?
    inputBinding:
      prefix: --cpus
    doc: |
      Number of cores/cpus to use.
      Default: 1
  seed:
    type: int?
    inputBinding:
      prefix: --seed
    doc: |
      Seed number for random values.
      Default: 42
outputs:
  peaks_bigbed_file:
    type: File
    outputBinding:
      glob: '*_peaks.bigBed'
    doc: |
      Locations of open-chromatin regions ("peaks")
      in bigBed format.
  cut_sites_bigwig_file:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_cut_cov.bigWig'
    doc: |
      Genome coverage calculated for Tn5 cut sites
      in bigWig format.
  fragments_bigwig_file:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_frg_cov.bigWig'
    doc: |
      Genome coverage calculated for ATAC fragments
      in bigWig format.
  sc_report_html_file:
    type: File?
    outputBinding:
      glob: sc_report.html
    doc: |
      Tehcnical report.
      HTML format.
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- Rscript
arguments:
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_atac_coverage.R"]:"/usr/local/bin/sc_atac_coverage.R")
stdout: sc_atac_coverage_stdout.log
stderr: sc_atac_coverage_stderr.log
label: Single-Cell ATAC-Seq Genome Coverage
doc: |
  Single-Cell ATAC-Seq Genome Coverage

  Creates genome coverage bigWig files from the provided ATAC fragments file
  and selected grouping parameters.

  --tmpdir parameter is not exposed as input.
