cwlVersion: v1.0
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var split_features = function(line) { function get_unique(value, index, self) { return self.indexOf(value) === index && value != ""; } let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null; return (splitted_line && !!splitted_line.length)?splitted_line:null; };
sd:upstream:
  sc_tools_sample:
  - sc-rna-cluster.cwl
  - sc-atac-cluster.cwl
  - sc-wnn-cluster.cwl
  - sc-ctype-assign.cwl
  - sc-rna-azimuth.cwl
inputs:
  alias:
    type: string
    label: Analysis name
    sd:preview:
      position: 1
  query_data_rds:
    type: File
    label: Single-cell Analysis with Clustered RNA-Seq Datasets
    doc: |
      Analysis that includes single-cell
      multiome RNA and ATAC-Seq or just
      RNA-Seq datasets run through either
      "Single-Cell Manual Cell Type
      Assignment" (based on the RNA or WNN
      clustering results), "Single-Cell
      RNA-Seq Cluster Analysis", or
      "Single-Cell WNN Cluster Analysis"
      at any of the processing stages.
    sd:upstreamSource: sc_tools_sample/seurat_data_rds
    sd:localLabel: true
  dimensions:
    type: int?
    label: Target dimensionality
    default: 0
    doc: |
      Number of principal components to be used
      in the trajectory analysis. Accepted values
      range from 1 to 50. Will fail if used more
      dimensions than it was available in the
      selected "Single-cell Analysis with
      Clustered RNA-Seq Datasets". If set to 0,
      use all available dimensions
      Default: 0
  query_source_column:
    type: string
    label: Cells grouping
    doc: |
      Single cell metadata column to group
      cells into clusters. Usually, in a form
      of "[rna|atac|wsnn]_res.X", where X is
      the clustering resolution. If cell types
      are available, add "custom_" prefix to
      the column name.
  trajectory_start:
    type: string?
    label: Trajectory start (optional)
    doc: |
      Value from the single cell metadata
      column used for grouping cells into
      the clusters. This value will define
      the trajectory starting point.
  genes_of_interest:
    type: string?
    default: null
    label: Genes of interest
    doc: |
      Comma or space separated list of
      genes of interest to visualize
      expression.
      Default: None
  barcodes_data:
    type: File?
    label: Selected cell barcodes (optional)
    doc: |
      A TSV/CSV file to optionally prefilter
      the single cell data by including only
      the cells with the selected barcodes.
      The provided file should include at
      least one column named "barcode", with
      one cell barcode per line. All other
      columns, except for "barcode", will be
      added to the single cell metadata loaded
      from "Single-cell Analysis with Clustered
      RNA-Seq Datasets" and can be utilized in
      the current or future steps of analysis.
  export_ucsc_cb:
    type: boolean?
    default: false
    label: Show results in UCSC Cell Browser
    doc: |
      Export results into UCSC Cell Browser
      Default: false
    sd:layout:
      advanced: true
  export_loupe_data:
    type: boolean?
    default: false
    label: Save raw counts to Loupe file. I confirm that data is generated by 10x technology and accept the EULA available at https://10xgen.com/EULA
    doc: |
      Save raw counts from the RNA assay to Loupe file. By
      enabling this feature you accept the End-User License
      Agreement available at https://10xgen.com/EULA.
      Default: false
    sd:layout:
      advanced: true
  export_html_report:
    type: boolean?
    default: true
    label: Show HTML report
    doc: |
      Export tehcnical report in HTML format.
      Default: true
    sd:layout:
      advanced: true
  color_theme:
    type:
    - 'null'
    - type: enum
      symbols:
      - gray
      - bw
      - linedraw
      - light
      - dark
      - minimal
      - classic
      - void
    default: classic
    label: Plots color theme
    doc: |
      Color theme for all plots saved
      as PNG files.
      Default: classic
    sd:layout:
      advanced: true
  threads:
    type:
    - 'null'
    - type: enum
      symbols:
      - '1'
      - '2'
      - '3'
      - '4'
      - '5'
      - '6'
    default: '4'
    label: Cores/CPUs
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 4
    sd:layout:
      advanced: true
outputs:
  trjc_gr_clst_plot_png:
    type: File?
    outputSource: rna_trajectory/trjc_gr_clst_plot_png
    label: Trajectory plot, colored by cluster
    doc: |
      Trajectory plot, colored by cluster
    sd:visualPlugins:
    - image:
        tab: Trajectory
        Caption: Trajectory plot, colored by cluster
  trjc_pstm_plot_png:
    type: File?
    outputSource: rna_trajectory/trjc_pstm_plot_png
    label: Trajectory plot, colored by pseudotime
    doc: |
      Trajectory plot, colored by pseudotime
    sd:visualPlugins:
    - image:
        tab: Trajectory
        Caption: Trajectory plot, colored by pseudotime
  dndr_gr_clst_plot_png:
    type: File?
    outputSource: rna_trajectory/dndr_gr_clst_plot_png
    label: Dendrogram plot, colored by cluster
    doc: |
      Dendrogram plot, colored by cluster
    sd:visualPlugins:
    - image:
        tab: Trajectory
        Caption: Dendrogram plot, colored by cluster
  dndr_pstm_plot_png:
    type: File?
    outputSource: rna_trajectory/dndr_pstm_plot_png
    label: Dendrogram plot, colored by pseudotime
    doc: |
      Dendrogram plot, colored by pseudotime
    sd:visualPlugins:
    - image:
        tab: Trajectory
        Caption: Dendrogram plot, colored by pseudotime
  grph_gr_clst_plot_png:
    type: File?
    outputSource: rna_trajectory/grph_gr_clst_plot_png
    label: Trajectory graph, colored by cluster
    doc: |
      Trajectory graph, colored by cluster
    sd:visualPlugins:
    - image:
        tab: Topology
        Caption: Trajectory graph, colored by cluster
  grph_pstm_plot_png:
    type: File?
    outputSource: rna_trajectory/grph_pstm_plot_png
    label: Trajectory graph, colored by pseudotime
    doc: |
      Trajectory graph, colored by pseudotime
    sd:visualPlugins:
    - image:
        tab: Topology
        Caption: Trajectory graph, colored by pseudotime
  tplg_plot_png:
    type: File?
    outputSource: rna_trajectory/tplg_plot_png
    label: Topology plot
    doc: |
      Topology plot
    sd:visualPlugins:
    - image:
        tab: Topology
        Caption: Topology plot
  xpr_htmp_plot_png:
    type: File?
    outputSource: rna_trajectory/xpr_htmp_plot_png
    label: Gene expression heatmap
    doc: |
      Gene expression heatmap
    sd:visualPlugins:
    - image:
        tab: Gene expression
        Caption: Gene expression heatmap
  xpr_pstm_plot_png:
    type: File?
    outputSource: rna_trajectory/xpr_pstm_plot_png
    label: Gene expression along pseudotime
    doc: |
      Gene expression along pseudotime
    sd:visualPlugins:
    - image:
        tab: Gene expression
        Caption: Gene expression along pseudotime
  umap_rd_rnaumap_plot_png:
    type: File?
    outputSource: rna_trajectory/umap_rd_rnaumap_plot_png
    label: UMAP, colored by pseudotime, RNA
    doc: |
      UMAP, colored by pseudotime, RNA
    sd:visualPlugins:
    - image:
        tab: Pseudotime
        Caption: UMAP, colored by pseudotime, RNA
  umap_rd_atacumap_plot_png:
    type: File?
    outputSource: rna_trajectory/umap_rd_atacumap_plot_png
    label: UMAP, colored by pseudotime, ATAC
    doc: |
      UMAP, colored by pseudotime, ATAC
    sd:visualPlugins:
    - image:
        tab: Pseudotime
        Caption: UMAP, colored by pseudotime, ATAC
  umap_rd_wnnumap_plot_png:
    type: File?
    outputSource: rna_trajectory/umap_rd_wnnumap_plot_png
    label: UMAP, colored by pseudotime, WNN
    doc: |
      UMAP, colored by pseudotime, WNN
    sd:visualPlugins:
    - image:
        tab: Pseudotime
        Caption: UMAP, colored by pseudotime, WNN
  pstm_dnst_spl_idnt_plot_png:
    type: File?
    outputSource: rna_trajectory/pstm_dnst_spl_idnt_plot_png
    label: Pseudotime density, split by dataset
    doc: |
      Pseudotime density, split by dataset
    sd:visualPlugins:
    - image:
        tab: Per dataset
        Caption: Pseudotime density, split by dataset
  pstm_hist_gr_clst_spl_idnt_plot_png:
    type: File?
    outputSource: rna_trajectory/pstm_hist_gr_clst_spl_idnt_plot_png
    label: Pseudotime histogram, colored by cluster, split by dataset
    doc: |
      Pseudotime histogram,
      colored by cluster,
      split by dataset
    sd:visualPlugins:
    - image:
        tab: Per dataset
        Caption: Pseudotime histogram, colored by cluster, split by dataset
  umap_spl_idnt_rd_rnaumap_plot_png:
    type: File?
    outputSource: rna_trajectory/umap_spl_idnt_rd_rnaumap_plot_png
    label: UMAP, colored by pseudotime, split by dataset, RNA
    doc: |
      UMAP, colored by pseudotime,
      split by dataset, RNA
    sd:visualPlugins:
    - image:
        tab: Per dataset
        Caption: UMAP, colored by pseudotime, split by dataset, RNA
  umap_spl_idnt_rd_atacumap_plot_png:
    type: File?
    outputSource: rna_trajectory/umap_spl_idnt_rd_atacumap_plot_png
    label: UMAP, colored by pseudotime, split by dataset, ATAC
    doc: |
      UMAP, colored by pseudotime,
      split by dataset, ATAC
    sd:visualPlugins:
    - image:
        tab: Per dataset
        Caption: UMAP, colored by pseudotime, split by dataset, ATAC
  umap_spl_idnt_rd_wnnumap_plot_png:
    type: File?
    outputSource: rna_trajectory/umap_spl_idnt_rd_wnnumap_plot_png
    label: UMAP, colored by pseudotime, split by dataset, WNN
    doc: |
      UMAP, colored by pseudotime,
      split by dataset, WNN
    sd:visualPlugins:
    - image:
        tab: Per dataset
        Caption: UMAP, colored by pseudotime, split by dataset, WNN
  pstm_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: rna_trajectory/pstm_dnst_spl_cnd_plot_png
    label: Pseudotime density, split by grouping condition
    doc: |
      Pseudotime density, split by
      grouping condition
    sd:visualPlugins:
    - image:
        tab: Per group
        Caption: Pseudotime density, split by grouping condition
  pstm_hist_gr_clst_spl_cnd_plot_png:
    type: File?
    outputSource: rna_trajectory/pstm_hist_gr_clst_spl_cnd_plot_png
    label: Pseudotime histogram, colored by cluster, split by grouping condition
    doc: |
      Pseudotime histogram, colored by
      cluster, split by grouping condition
    sd:visualPlugins:
    - image:
        tab: Per group
        Caption: Pseudotime histogram, colored by cluster, split by grouping condition
  umap_spl_cnd_rd_rnaumap_plot_png:
    type: File?
    outputSource: rna_trajectory/umap_spl_cnd_rd_rnaumap_plot_png
    label: UMAP, colored by pseudotime, split by grouping condition, RNA
    doc: |
      UMAP, colored by pseudotime,
      split by grouping condition, RNA
    sd:visualPlugins:
    - image:
        tab: Per group
        Caption: UMAP, colored by pseudotime, split by grouping condition, RNA
  umap_spl_cnd_rd_atacumap_plot_png:
    type: File?
    outputSource: rna_trajectory/umap_spl_cnd_rd_atacumap_plot_png
    label: UMAP, colored by pseudotime, split by grouping condition, ATAC
    doc: |
      UMAP, colored by pseudotime,
      split by grouping condition, ATAC
    sd:visualPlugins:
    - image:
        tab: Per group
        Caption: UMAP, colored by pseudotime, split by grouping condition, ATAC
  umap_spl_cnd_rd_wnnumap_plot_png:
    type: File?
    outputSource: rna_trajectory/umap_spl_cnd_rd_wnnumap_plot_png
    label: UMAP, colored by pseudotime, split by grouping condition, WNN
    doc: |
      UMAP, colored by pseudotime,
      split by grouping condition, WNN
    sd:visualPlugins:
    - image:
        tab: Per group
        Caption: UMAP, colored by pseudotime, split by grouping condition, WNN
  ucsc_cb_html_data:
    type: Directory?
    outputSource: rna_trajectory/ucsc_cb_html_data
    label: UCSC Cell Browser (data)
    doc: |
      UCSC Cell Browser html data.
  ucsc_cb_html_file:
    type: File?
    outputSource: rna_trajectory/ucsc_cb_html_file
    label: UCSC Cell Browser
    doc: |
      UCSC Cell Browser html index.
    sd:visualPlugins:
    - linkList:
        tab: Overview
        target: _blank
  seurat_data_rds:
    type: File
    outputSource: rna_trajectory/seurat_data_rds
    label: Seurat object in RDS format
    doc: |
      Seurat object.
      RDS format.
  seurat_data_cloupe:
    type: File?
    outputSource: rna_trajectory/seurat_data_cloupe
    label: Seurat object in Loupe format
    doc: |
      Seurat object.
      Loupe format.
  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: Compressed folder with all PDF plots
    doc: |
      Compressed folder with all PDF plots.
  sc_report_html_file:
    type: File?
    outputSource: rna_trajectory/sc_report_html_file
    label: Analysis log
    doc: |
      Tehcnical report.
      HTML format.
    sd:visualPlugins:
    - linkList:
        tab: Overview
        target: _blank
  rna_trajectory_stdout_log:
    type: File
    outputSource: rna_trajectory/stdout_log
    label: Output log
    doc: |
      Stdout log from the rna_trajectory step.
  rna_trajectory_stderr_log:
    type: File
    outputSource: rna_trajectory/stderr_log
    label: Error log
    doc: |
      Stderr log from the rna_trajectory step.
steps:
  rna_trajectory:
    run: ../tools/sc-rna-trajectory.cwl
    in:
      query_data_rds: query_data_rds
      barcodes_data: barcodes_data
      reduction:
        default: pca
      dimensions:
        source: dimensions
        valueFrom: $(self==0?null:self)
      query_source_column: query_source_column
      trajectory_start:
        source: trajectory_start
        valueFrom: $(self==""?null:self)
      genes_of_interest:
        source: genes_of_interest
        valueFrom: $(split_features(self))
      predictive_genes:
        default: 100
      verbose:
        default: true
      export_ucsc_cb: export_ucsc_cb
      export_loupe_data: export_loupe_data
      export_pdf_plots:
        default: true
      color_theme: color_theme
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 128
      export_html_report: export_html_report
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - trjc_gr_clst_plot_png
    - trjc_gr_clst_plot_pdf
    - trjc_pstm_plot_png
    - trjc_pstm_plot_pdf
    - grph_gr_clst_plot_png
    - grph_gr_clst_plot_pdf
    - grph_pstm_plot_png
    - grph_pstm_plot_pdf
    - dndr_gr_clst_plot_png
    - dndr_gr_clst_plot_pdf
    - dndr_pstm_plot_png
    - dndr_pstm_plot_pdf
    - tplg_plot_png
    - tplg_plot_pdf
    - xpr_htmp_plot_png
    - xpr_htmp_plot_pdf
    - xpr_pstm_plot_png
    - xpr_pstm_plot_pdf
    - pstm_dnst_spl_idnt_plot_png
    - pstm_dnst_spl_idnt_plot_pdf
    - pstm_dnst_spl_cnd_plot_png
    - pstm_dnst_spl_cnd_plot_pdf
    - pstm_hist_gr_clst_spl_idnt_plot_png
    - pstm_hist_gr_clst_spl_idnt_plot_pdf
    - pstm_hist_gr_clst_spl_cnd_plot_png
    - pstm_hist_gr_clst_spl_cnd_plot_pdf
    - umap_rd_rnaumap_plot_png
    - umap_rd_rnaumap_plot_pdf
    - umap_rd_atacumap_plot_png
    - umap_rd_atacumap_plot_pdf
    - umap_rd_wnnumap_plot_png
    - umap_rd_wnnumap_plot_pdf
    - umap_spl_idnt_rd_rnaumap_plot_png
    - umap_spl_idnt_rd_rnaumap_plot_pdf
    - umap_spl_idnt_rd_atacumap_plot_png
    - umap_spl_idnt_rd_atacumap_plot_pdf
    - umap_spl_idnt_rd_wnnumap_plot_png
    - umap_spl_idnt_rd_wnnumap_plot_pdf
    - umap_spl_cnd_rd_rnaumap_plot_png
    - umap_spl_cnd_rd_rnaumap_plot_pdf
    - umap_spl_cnd_rd_atacumap_plot_png
    - umap_spl_cnd_rd_atacumap_plot_pdf
    - umap_spl_cnd_rd_wnnumap_plot_png
    - umap_spl_cnd_rd_wnnumap_plot_pdf
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - seurat_data_cloupe
    - sc_report_html_file
    - stdout_log
    - stderr_log
  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - rna_trajectory/trjc_gr_clst_plot_pdf
        - rna_trajectory/trjc_pstm_plot_pdf
        - rna_trajectory/grph_gr_clst_plot_pdf
        - rna_trajectory/grph_pstm_plot_pdf
        - rna_trajectory/dndr_gr_clst_plot_pdf
        - rna_trajectory/dndr_pstm_plot_pdf
        - rna_trajectory/tplg_plot_pdf
        - rna_trajectory/xpr_htmp_plot_pdf
        - rna_trajectory/xpr_pstm_plot_pdf
        - rna_trajectory/pstm_dnst_spl_idnt_plot_pdf
        - rna_trajectory/pstm_dnst_spl_cnd_plot_pdf
        - rna_trajectory/pstm_hist_gr_clst_spl_idnt_plot_pdf
        - rna_trajectory/pstm_hist_gr_clst_spl_cnd_plot_pdf
        - rna_trajectory/umap_rd_rnaumap_plot_pdf
        - rna_trajectory/umap_rd_atacumap_plot_pdf
        - rna_trajectory/umap_rd_wnnumap_plot_pdf
        - rna_trajectory/umap_spl_idnt_rd_rnaumap_plot_pdf
        - rna_trajectory/umap_spl_idnt_rd_atacumap_plot_pdf
        - rna_trajectory/umap_spl_idnt_rd_wnnumap_plot_pdf
        - rna_trajectory/umap_spl_cnd_rd_rnaumap_plot_pdf
        - rna_trajectory/umap_spl_cnd_rd_atacumap_plot_pdf
        - rna_trajectory/umap_spl_cnd_rd_wnnumap_plot_pdf
        valueFrom: $(self.flat().filter(n => n))
      folder_basename:
        default: pdf_plots
    out:
    - folder
  compress_pdf_plots:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: folder_pdf_plots/folder
    out:
    - compressed_folder
label: Single-Cell RNA-Seq Trajectory Analysis
doc: |-
  Single-Cell RNA-Seq Trajectory Analysis

  Infers developmental trajectories and pseudotime from
  cells clustered by similarity of gene expression data.
sd:version: 100
