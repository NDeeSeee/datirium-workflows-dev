cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_label = function(input_array, i) { var rootname = input_array[i].basename.split('.').slice(0,-1).join('.'); rootname = (rootname=="")?input_array[i].basename:rootname; return inputs.gem_well_labels?inputs.gem_well_labels[i].replace(/\t|\s|\[|\]|\>|\<|,|\./g, "_"):rootname; };
- class: InitialWorkDirRequirement
  listing: |
    ${
      var grouping = "library_id\tcondition\n"
      if (inputs.molecule_info_h5 != null){
        var entry = "sample_id,molecule_h5\n"
        for (var i=0; i < inputs.molecule_info_h5.length; i++){
          entry += get_label(inputs.molecule_info_h5, i) + "," + inputs.molecule_info_h5[i].path + "\n"
          grouping += get_label(inputs.molecule_info_h5, i) + "\t" + get_label(inputs.molecule_info_h5, i) + "\n"
        }
      } else if (inputs.filtered_data_folder != null){
        var entry = "sample_id,sample_outs,donor,origin\n"
        for (var i=0; i < inputs.filtered_data_folder.length; i++){
          var donor = "donor"
          var origin = "origin"
          if (inputs.clonotype_grouping == "same_donor_different_origins"){
            origin = "origin_" + i
          } else if (inputs.clonotype_grouping == "different_donors"){
            donor = "donor_" + i
            origin = "origin_" + i
          }
          entry += get_label(inputs.filtered_data_folder, i) + "," + inputs.filtered_data_folder[i].path + "," + donor  + "," + origin + "\n"
          grouping += get_label(inputs.filtered_data_folder, i) + "\t" + get_label(inputs.filtered_data_folder, i) + "\n"
        }
      } else {
        var entry = "neither molecule_info_h5 nor filtered_data_folder was provided"
        var grouping = "neither molecule_info_h5 nor filtered_data_folder was provided"
      }
      return [
        {
          "entry": entry,
          "entryname": "metadata.csv",
          "writable": true
        },
        {
          "entry": grouping,
          "entryname": "grouping.tsv",
          "writable": true
        }
      ];
    }
hints:
- class: DockerRequirement
  dockerPull: cumulusprod/cellranger:8.0.1
inputs:
  molecule_info_h5:
    type:
    - 'null'
    - File[]
    doc: |
      Array of molecule-level information files in HDF5 format.
      Outputs from "cellranger count" command. Either
      molecule_info_h5 or filtered_data_folder should be
      provided. If both inputs are provided - use molecule_info_h5.
  filtered_data_folder:
    type:
    - 'null'
    - Directory[]
    doc: |
      Array of folders containing filtered data, i.e., only
      cell-associated barcodes. Outputs from "cellranger multi"
      command. Either molecule_info_h5 or filtered_data_folder should
      be provided. If both inputs are provided - use molecule_info_h5.
  gem_well_labels:
    type:
    - 'null'
    - string[]
    doc: |
      Array of GEM well identifiers to be used for labeling purposes only.
      If not provided use rootnames of files from the molecule_info_h5 or
      directories from filtered_data_folder inputs. If labels are not
      unique, cellranger will fails.
  normalization_mode:
    type:
    - 'null'
    - type: enum
      name: normalization
      symbols:
      - none
      - mapped
    inputBinding:
      position: 5
      prefix: --normalize
    doc: |
      Library depth normalization mode: mapped, none.
      Default: mapped
  clonotype_grouping:
    type:
    - 'null'
    - type: enum
      name: clonotype_grouping
      symbols:
      - same_donor_different_origins
      - same_donor_and_origin
      - different_donors
    default: different_donors
    doc: |
      When cellranger aggr is called with cellranger multi outputs, there are three
      ways it can process the datasets depending on the combination of donor and
      origin values:
      1. If two datasets come from the same donor but have different origins, Cell Ranger
         will rerun the clonotype grouping algorithm on the combined set of cells. This
         allows cells from different datasets to belong to the same clonotype.
      2. If two datasets come from the same donor and origin, then Cell Ranger performs
         additional filtering to remove certain rare artifacts. For example, Cell Ranger
         will filter expanded exact subclonotypes that are present in one library but not
         in another from the same origin, which would be highly improbable, assuming random
         draws of cells from the tube. These are believed to arise when a plasma or
         plasmablast cell breaks up during or after pipetting from the tube, and the resulting
         fragments contaminate GEMs, yielding expanded false clonotypes that are residues of
         real single plasma cells.
      3. If two cells came from different donors, then Cell Ranger will not put them in the
         same clonotype.
      Ignored if cellranger aggr is run with molecule_info_h5 inputs.
  threads:
    type: int?
    inputBinding:
      position: 6
      prefix: --localcores
    doc: |
      Set max cores the pipeline may request at one time.
      Default: all available
  memory_limit:
    type: int?
    inputBinding:
      position: 7
      prefix: --localmem
    doc: |
      Set max GB the pipeline may request at one time
      Default: all available
  virt_memory_limit:
    type: int?
    inputBinding:
      position: 8
      prefix: --localvmem
    doc: |
      Set max virtual address space in GB for the pipeline
      Default: all available
outputs:
  web_summary_report:
    type: File
    outputBinding:
      glob: aggregated/outs/web_summary.html
    doc: |
      Aggregated run summary metrics and charts in HTML format
  metrics_summary_report_json:
    type: File
    outputBinding:
      glob: aggregated/outs/count/summary.json
    doc: |
      Aggregated RNA run summary metrics in JSON format
  secondary_analysis_report_folder:
    type: Directory
    outputBinding:
      glob: aggregated/outs/count/analysis
    doc: |
      Folder with secondary analysis of RNA data including dimensionality reduction,
      cell clustering, and differential expression
  filtered_feature_bc_matrix_folder:
    type: Directory
    outputBinding:
      glob: aggregated/outs/count/filtered_feature_bc_matrix
    doc: |
      Folder with aggregated filtered feature-barcode matrices
      containing only cellular barcodes in MEX format
  filtered_feature_bc_matrix_h5:
    type: File
    outputBinding:
      glob: aggregated/outs/count/filtered_feature_bc_matrix.h5
    doc: |
      Filtered feature-barcode matrices containing only cellular
      barcodes in HDF5 format.
  aggregation_metadata:
    type: File
    outputBinding:
      glob: aggregated/outs/aggregation.csv
    doc: |
      Copy of the input aggregation CSV file
  grouping_data:
    type: File
    outputBinding:
      glob: grouping.tsv
    doc: |
      Example of TSV file to define datasets grouping
  loupe_browser_track:
    type: File
    outputBinding:
      glob: aggregated/outs/count/cloupe.cloupe
    doc: |
      Loupe Browser visualization and analysis file
  airr_rearrangement_tsv:
    type: File?
    outputBinding:
      glob: aggregated/outs/vdj_*/airr_rearrangement.tsv
    doc: |
      Annotated contigs and consensus sequences of V(D)J
      rearrangements in the AIRR format. It includes only
      viable cells identified by both V(D)J and RNA algorithms.
  clonotypes_csv:
    type: File?
    outputBinding:
      glob: aggregated/outs/vdj_*/clonotypes.csv
    doc: |
      CSV file with high-level descriptions of each clonotype
  consensus_sequences_fasta:
    type: File?
    outputBinding:
      glob: aggregated/outs/vdj_*/consensus.fasta
    doc: |
      The consensus sequence of each assembled contig.
  consensus_annotations_csv:
    type: File?
    outputBinding:
      glob: aggregated/outs/vdj_*/consensus_annotations.csv
    doc: |
      CSV file with high-level and detailed annotations of each clonotype
      consensus sequence.
  filtered_contig_annotations_csv:
    type: File?
    outputBinding:
      glob: aggregated/outs/vdj_*/filtered_contig_annotations.csv
    doc: |
      CSV file with high-level annotations of each high-confidence contig from
      cell-associated barcodes
  loupe_vdj_browser_track:
    type: File?
    outputBinding:
      glob: aggregated/outs/vdj_*/vloupe.vloupe
    doc: |
      Loupe V(D)J Browser visualization and analysis file
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- cellranger
- aggr
- --disable-ui
- --id
- aggregated
- --csv
- metadata.csv
stdout: cellranger_aggr_stdout.log
stderr: cellranger_aggr_stderr.log
label: Cell Ranger Aggregate (RNA, RNA+VDJ)
doc: |
  Cell Ranger Aggregate (RNA, RNA+VDJ)

  Aggregates outputs from multiple runs of the "Cell Ranger Count
  (RNA)" (if molecule_info_h5 input provided) or "Cell Ranger Count
  (RNA+VDJ)" experiments (if filtered_data_folder input provided).
  If both inputs are provided - uses molecule_info_h5. If neither of
  them are provided cellranger aggr will fail.

  Parameters set by default:
  --disable-ui - no need in any UI when running in Docker container
  --id - hardcoded to `aggregated` as we want to return the content
         of the outputs folder as separate outputs

  Skipped parameters:
  --nosecondary
  --dry
  --min-crispr-umi
  --noexit
  --nopreflight
  --description
  --jobmode
  --mempercore
  --maxjobs
  --jobinterval
  --overrides
  --output-dir
  --uiport

  Not supported features when aggregating RNA experiments:
  - Batch correction caused by different versions of the Single Cell Gene
    Expression chemistry is not supported as the generated metadata file
    for merging molecule_info_h5 inputs doesn't include "batch" field.
