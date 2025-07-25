cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      var listing = [
        {
          "entry": inputs.fastq_file_r1,
          "entryname": "sample_S1_L001_R1_001.fastq",
          "writable": true
        },
        {
          "entry": inputs.fastq_file_r2,
          "entryname": "sample_S1_L001_R2_001.fastq",
          "writable": true
        }
      ];
      if (inputs.fastq_file_i1){
        listing.push(
          {
            "entry": inputs.fastq_file_i1,
            "entryname": "sample_S1_L001_I1_001.fastq",
            "writable": true
          }
        );
      };
      return listing;
    }
hints:
- class: DockerRequirement
  dockerPull: cumulusprod/cellranger:8.0.1
inputs:
  fastq_file_r1:
    type: File
    doc: |
      FASTQ read 1 file.
      It will be staged into workdir
      as sample_S1_L001_R1_001.fastq.
  fastq_file_r2:
    type: File
    doc: |
      FASTQ read 2 file.
      It will be staged into workdir
      as sample_S1_L001_R2_001.fastq.
  fastq_file_i1:
    type: File?
    doc: |
      FASTQ index file.
      If provided, it will be staged
      into workdir as
      sample_S1_L001_I1_001.fastq.
  indices_folder:
    type: Directory
    inputBinding:
      position: 5
      prefix: --transcriptome
    doc: |
      Path of folder containing 10x-compatible
      transcriptome reference. These indices
      should be generated by "cellranger mkref"
      command.
  r1_length:
    type: int?
    inputBinding:
      position: 6
      prefix: --r1-length
    doc: |
      Limit the length of the input Read 1 sequence
      of Gene Expression library to the first N bases,
      where N is a user-supplied value. Note that the
      length includes the 10x Barcode and UMI sequences
      so do not set this below 26 for Single Cell 3′ v2
      or Single Cell 5′. This and --r2-length are useful
      options for determining the optimal read length
      for sequencing.
  r2_length:
    type: int?
    inputBinding:
      position: 7
      prefix: --r2-length
    doc: |
      Limit the length of the input R2 sequence to the
      first N bases, where N is a user-supplied value.
      Trimming occurs before sequencing metrics are
      computed and therefore, limiting R2 read length
      may affect Q30 scores.
  expect_cells:
    type: int?
    inputBinding:
      position: 8
      prefix: --expect-cells
    doc: |
      Expected number of recovered cells.
      Starting in Cell Ranger 7.0, the expected number
      of cells can be either auto-estimated or specified
      with --expect-cells. To replicate an old cellranger
      count analysis, set this parameter to 3,000 cells.
  force_cells:
    type: int?
    inputBinding:
      position: 9
      prefix: --force-cells
    doc: |
      Force pipeline to use this number of cells, bypassing
      the cell detection algorithm. Use this if the number
      of cells estimated by Cell Ranger is not consistent
      with the barcode rank plot.
  no_bam:
    type: boolean?
    default: false
    inputBinding:
      prefix: --create-bam=
      valueFrom: $(self?"false":"true")
      separate: false
      position: 10
    doc: |
      Enable or disable BAM file generation. Setting
      --create-bam=false reduces the total computation
      time and the size of the output directory (BAM
      file not generated). We recommend setting
      --create-bam=true if unsure. See
      https://10xgen.com/create-bam for additional
      guidance [possible values: true, false]
  exclude_introns:
    type: boolean?
    default: false
    inputBinding:
      prefix: --include-introns=
      valueFrom: $(self?"false":"true")
      separate: false
      position: 11
    doc: |
      Starting from the Cell Ranger v7.0 the intronic reads
      are counted by default for whole transcriptome gene
      expression data, except when --target-panel is used.
      Therefore, here we provide a flag to disable this
      default behavior.
  threads:
    type: int?
    inputBinding:
      position: 12
      prefix: --localcores
    doc: |
      Set max cores the pipeline may request at one time.
      Default: all available
  memory_limit:
    type: int?
    inputBinding:
      position: 13
      prefix: --localmem
    doc: |
      Set max GB the pipeline may request at one time
      Default: all available
  virt_memory_limit:
    type: int?
    inputBinding:
      position: 14
      prefix: --localvmem
    doc: |
      Set max virtual address space in GB for the pipeline
      Default: all available
outputs:
  web_summary_report:
    type: File
    outputBinding:
      glob: sample/outs/web_summary.html
    doc: |
      Run summary metrics and charts in HTML format.
  metrics_summary_report:
    type: File
    outputBinding:
      glob: sample/outs/metrics_summary.csv
    doc: |
      Run summary metrics in CSV format.
  possorted_genome_bam_bai:
    type: File?
    outputBinding:
      glob: sample/outs/possorted_genome_bam.bam
    secondaryFiles:
    - .bai
    doc: |
      Indexed reads aligned to the genome
      and transcriptome annotated with
      barcode information.
  filtered_feature_bc_matrix_folder:
    type: Directory
    outputBinding:
      glob: sample/outs/filtered_feature_bc_matrix
    doc: |
      Folder with filtered feature-barcode matrices
      containing only cellular barcodes in MEX format.
  filtered_feature_bc_matrix_h5:
    type: File
    outputBinding:
      glob: sample/outs/filtered_feature_bc_matrix.h5
    doc: |
      Filtered feature-barcode matrices containing
      only cellular barcodes in HDF5 format.
  raw_feature_bc_matrices_folder:
    type: Directory
    outputBinding:
      glob: sample/outs/raw_feature_bc_matrix
    doc: |
      Folder with unfiltered feature-barcode matrices
      containing all barcodes in MEX format.
  raw_feature_bc_matrices_h5:
    type: File
    outputBinding:
      glob: sample/outs/raw_feature_bc_matrix.h5
    doc: |
      Unfiltered feature-barcode matrices containing
      all barcodes in HDF5 format.
  secondary_analysis_report_folder:
    type: Directory
    outputBinding:
      glob: sample/outs/analysis
    doc: |
      Folder with secondary analysis results
      including dimensionality reduction, cell
      clustering, and differential expression.
  molecule_info_h5:
    type: File
    outputBinding:
      glob: sample/outs/molecule_info.h5
    doc: |
      Molecule-level information used by cellranger
      aggr to aggregate samples into larger datasets.
  loupe_browser_track:
    type: File
    outputBinding:
      glob: sample/outs/cloupe.cloupe
    doc: |
      Loupe Browser visualization and
      analysis file.
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- cellranger
- count
- --disable-ui
- --fastqs
- .
- --sample
- sample
- --id
- sample
stdout: cellranger_count_stdout.log
stderr: cellranger_count_stderr.log
label: Cell Ranger Count (RNA)
doc: |
  Cell Ranger Count (RNA)

  Quantifies single-cell gene expression of the sequencing
  data from a single 10x Genomics library.

  Starting from the Cell Ranger v7.0 the intronic reads are
  counted by default for whole transcriptome gene expression
  data. For more details see
  https://support.10xgenomics.com/docs/intron-mode-rec

  Input parameters for Feature Barcode, Targeted Gene Expression
  analyses are not implemented, therefore the correspondent outputs
  are also excluded.

  Parameters set by default:
  --disable-ui - no need in any UI when running in Docker container
  --id         - can be hardcoded as we rename input files anyway
  --fastqs     - points to the current directory, because input
                 FASTQ files are staged there
  --sample     - hardcoded to sample as we stage input fastq files
                 with the hardcoded names

  Not implemented parameters:
  --description                   - not needed for now
  --project                       - no needed to select input files by folder
  --lanes                         - not needed for now
  --libraries                     - needed only for Gene expression + Feature Barcode analysis
  --feature-ref                   - needed only for Feature Barcode analysis
  --target-panel                  - needed only for Targeted Gene Expression analysis
  --nosecondary                   - no reason to disable it
  --chemistry                     - cell ranger will autodetect the library by default
  --no-libraries                  - used only in Feature Barcode analysis
  --check-library-compatibility   - no reason to disable it
  --min-crispr-umi                - needed only for Protospacer calling (for pooled CRISPR screens)
  --no-target-umi-filter          - needed only for Targeted Gene Expression analysis
  --dry                           - not applicable to our use case
  --jobmode                       - we use default local mode
  --mempercore                    - not used for local mode
  --maxjobs                       - not used for local mode
  --jobinterval                   - not used for local mode
  --overrides                     - not needed for now
  --uiport                        - we disabled UI
  --noexit                        - we disabled UI
  --nopreflight                   - no reason to skip preflight checks
  --output-dir                    - will be saved to the sample/outs by default
