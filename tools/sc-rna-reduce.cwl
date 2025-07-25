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
      Path to the RDS file to load Seurat object from. This file should include genes
      expression information stored in the RNA assay.
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
      filtering will be applied.
      Default: no extra metadata is added
  barcodes_data:
    type: File?
    inputBinding:
      prefix: --barcodes
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata be selected barcodes.
      First column should be named as 'barcode'. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present.
      Default: all cells used, no extra metadata is added
  cell_cycle_data:
    type:
    - 'null'
    - File
    - type: enum
      symbols:
      - hg19
      - hg38
      - mm10
    inputBinding:
      prefix: --cellcycle
      valueFrom: |
        ${
          if (self.class && self.class == "File"){
            return self;
          } else if (self == "hg19") {
            return "/opt/sc_tools/human_cc_genes.csv";
          } else if (self == "hg38") {
            return "/opt/sc_tools/human_cc_genes.csv";
          } else if (self == "mm10") {
            return "/opt/sc_tools/mouse_cc_genes.csv";
          } else {
            return null;
          }
        }
    doc: |
      Path to the TSV/CSV file with the information for cell cycle score assignment.
      First column - 'phase', second column 'gene_id'. If loaded Seurat object already
      includes cell cycle scores in 'S.Score', 'G2M.Score', and 'CC.Difference' metatada
      columns they will be overwritten. If a string value provided, it should be one of
      the hg19, hg38, or mm10 as we replace it with the file location from docker image.
      Default: skip cell cycle score assignment.
  normalization_method:
    type:
    - 'null'
    - type: enum
      symbols:
      - sct
      - log
      - sctglm
    inputBinding:
      prefix: --norm
    doc: |
      Normalization method applied to genes expression counts. If loaded Seurat object
      includes multiple datasets, normalization will be run independently for each of
      them, unless integration is disabled with 'none' or set to 'harmony'
      Default: sctglm
  integration_method:
    type:
    - 'null'
    - type: enum
      symbols:
      - seurat
      - harmony
      - none
    inputBinding:
      prefix: --ntgr
    doc: |
      Integration method used for joint analysis of multiple datasets. Automatically
      set to 'none' if loaded Seurat object includes only one dataset.
      Default: seurat
  integrate_by:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: --ntgrby
    doc: |
      Column(s) from the Seurat object metadata to define the variable(s) that should
      be integrated out when running multiple datasets integration with harmony. May
      include columns from the extra metadata added with --metadata parameter. Ignored
      if --ntgr is not set to harmony.
      Default: new.ident
  highly_var_genes_count:
    type: int?
    inputBinding:
      prefix: --highvargenes
    doc: |
      Number of highly variable genes used in datasets integration, scaling and
      dimensionality reduction.
      Default: 3000
  regress_mito_perc:
    type: boolean?
    inputBinding:
      prefix: --regressmt
    doc: |
      Regress the percentage of RNA reads mapped to mitochondrial genes as a
      confounding source of variation.
      Default: false
  regress_genes:
    type: string?
    inputBinding:
      prefix: --regressgenes
    doc: |
      Regex pattern to identify genes which expression should be
      regressed as a confounding source of variation. Default: none
  regress_ccycle_full:
    type: boolean?
    inputBinding:
      prefix: --regressccfull
    doc: |
      Regress all signals associated with cell cycle phase.
      Ignored if --cellcycle is not provided. Mutually exclusive
      with --regressccdiff parameter.
      Default: false
  regress_ccycle_diff:
    type: boolean?
    inputBinding:
      prefix: --regressccdiff
    doc: |
      Regress only differences in cell cycle phase among proliferating
      cells. Signals separating non-cycling and cycling cells will be
      maintained. Ignored if --cellcycle is not provided. Mutually
      exclusive with --regressccfull
      Default: false
  dimensions:
    type: int?
    inputBinding:
      prefix: --dimensions
    doc: |
      Dimensionality to use for datasets integration (if provided RDS
      file includes multiple datasets and --ntgr is not set to 'harmony')
      and UMAP projection (from 1 to 50). Set to 0 to use auto-estimated
      dimensionality.
      Default: 10
  umap_spread:
    type: float?
    inputBinding:
      prefix: --uspread
    doc: |
      The effective scale of embedded points on UMAP. In combination with '--mindist'
      it determines how clustered/clumped the embedded points are.
      Default: 1
  umap_mindist:
    type: float?
    inputBinding:
      prefix: --umindist
    doc: |
      Controls how tightly the embedding is allowed compress points together on UMAP.
      Larger values ensure embedded points are moreevenly distributed, while smaller
      values allow the algorithm to optimise more accurately with regard to local structure.
      Sensible values are in the range 0.001 to 0.5.
      Default:  0.3
  umap_neighbors:
    type: int?
    inputBinding:
      prefix: --uneighbors
    doc: |
      Determines the number of neighboring points used in UMAP. Larger values will result
      in more global structure being preserved at the loss of detailed local structure.
      In general this parameter should often be in the range 5 to 50.
      Default: 30
  umap_metric:
    type:
    - 'null'
    - type: enum
      symbols:
      - euclidean
      - manhattan
      - chebyshev
      - minkowski
      - canberra
      - braycurtis
      - mahalanobis
      - wminkowski
      - seuclidean
      - cosine
      - correlation
      - haversine
      - hamming
      - jaccard
      - dice
      - russelrao
      - kulsinski
      - ll_dirichlet
      - hellinger
      - rogerstanimoto
      - sokalmichener
      - sokalsneath
      - yule
    inputBinding:
      prefix: --umetric
    doc: |
      The metric to use to compute distances in high dimensional space for UMAP.
      Default: cosine
  umap_method:
    type:
    - 'null'
    - type: enum
      symbols:
      - uwot
      - uwot-learn
      - umap-learn
    inputBinding:
      prefix: --umethod
    doc: |
      UMAP implementation to run. If set to 'umap-learn' use --umetric 'correlation'
      Default: uwot
  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: --pdf
    doc: |
      Export plots in PDF.
      Default: false
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
    inputBinding:
      prefix: --theme
    doc: |
      Color theme for all generated plots. One of gray, bw, linedraw, light,
      dark, minimal, classic, void.
      Default: classic
  verbose:
    type: boolean?
    inputBinding:
      prefix: --verbose
    doc: |
      Print debug information.
      Default: false
  export_h5seurat_data:
    type: boolean?
    inputBinding:
      prefix: --h5seurat
    doc: |
      Save Seurat data to h5seurat file.
      Default: false
  export_h5ad_data:
    type: boolean?
    inputBinding:
      prefix: --h5ad
    doc: |
      Save raw counts from the RNA assay to h5ad file.
      Default: false
  export_loupe_data:
    type: boolean?
    inputBinding:
      prefix: --loupe
    doc: |
      Save raw counts from the RNA assay to Loupe file. By
      enabling this feature you accept the End-User License
      Agreement available at https://10xgen.com/EULA.
      Default: false
  export_scope_data:
    type: boolean?
    inputBinding:
      prefix: --scope
    doc: |
      Save Seurat data to SCope compatible loom file.
      Default: false
  export_ucsc_cb:
    type: boolean?
    inputBinding:
      prefix: --cbbuild
    doc: |
      Export results to UCSC Cell Browser. Default: false
  low_memory:
    type: boolean?
    inputBinding:
      prefix: --lowmem
    doc: |
      Attempts to minimize RAM usage when integrating multiple datasets
      with SCTransform algorithm (slows down the computation). Ignored if
      '--ntgr' is not set to 'seurat' or if '--norm' is not set to either
      'sct' or 'sctglm'.
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
  elbow_plot_png:
    type: File?
    outputBinding:
      glob: '*_elbow.png'
    doc: |
      Elbow plot.
      PNG format
  elbow_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_elbow.pdf'
    doc: |
      Elbow plot.
      PDF format
  qc_dim_corr_plot_png:
    type: File?
    outputBinding:
      glob: '*_qc_dim_corr.png'
    doc: |
      Correlation between QC metrics and principal components.
      PNG format
  qc_dim_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_qc_dim_corr.pdf'
    doc: |
      Correlation between QC metrics and principal components.
      PDF format
  umap_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_qc_mtrcs.png'
    doc: |
      UMAP, QC metrics.
      PNG format
  umap_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_qc_mtrcs.pdf'
    doc: |
      UMAP, QC metrics.
      PDF format
  umap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap.png'
    doc: |
      UMAP, colored by dataset.
      PNG format
  umap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap.pdf'
    doc: |
      UMAP, colored by dataset.
      PDF format
  umap_spl_ph_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_spl_ph.png'
    doc: |
      UMAP, colored by dataset, split by
      cell cycle phase.
      PNG format
  umap_spl_ph_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_spl_ph.pdf'
    doc: |
      UMAP, colored by dataset, split by
      cell cycle phase.
      PDF format
  ccpca_plot_png:
    type: File?
    outputBinding:
      glob: '*_ccpca.png'
    doc: |
      PCA, colored by cell cycle phase.
      PNG format
  ccpca_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_ccpca.pdf'
    doc: |
      PCA, colored by cell cycle phase.
      PDF format
  umap_spl_mito_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_spl_mito.png'
    doc: |
      UMAP, colored by dataset, split by
      mitochondrial percentage.
      PNG format
  umap_spl_mito_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_spl_mito.pdf'
    doc: |
      UMAP, colored by dataset, split by
      mitochondrial percentage.
      PDF format
  umap_spl_umi_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_spl_umi.png'
    doc: |
      UMAP, colored by dataset, split by
      RNA reads per cell.
      PNG format
  umap_spl_umi_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_spl_umi.pdf'
    doc: |
      UMAP, colored by dataset, split by
      RNA reads per cell.
      PDF format
  umap_spl_gene_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_spl_gene.png'
    doc: |
      UMAP, colored by dataset, split by
      genes per cell.
      PNG format
  umap_spl_gene_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_spl_gene.pdf'
    doc: |
      UMAP, colored by dataset, split by
      genes per cell.
      PDF format
  umap_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_spl_idnt.png'
    doc: |
      UMAP, split by dataset.
      PNG format
  umap_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_spl_idnt.pdf'
    doc: |
      UMAP, split by dataset.
      PDF format
  ccpca_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: '*_ccpca_spl_idnt.png'
    doc: |
      PCA, colored by cell cycle phase,
      split by dataset.
      PNG format
  ccpca_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_ccpca_spl_idnt.pdf'
    doc: |
      PCA, colored by cell cycle phase,
      split by dataset.
      PDF format
  umap_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_spl_cnd.png'
    doc: |
      UMAP, colored by dataset, split by
      grouping condition.
      PNG format
  umap_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_spl_cnd.pdf'
    doc: |
      UMAP, colored by dataset, split by
      grouping condition.
      PDF format
  umap_gr_cnd_spl_ph_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_gr_cnd_spl_ph.png'
    doc: |
      UMAP, colored by grouping condition,
      split by cell cycle phase.
      PNG format
  umap_gr_cnd_spl_ph_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_gr_cnd_spl_ph.pdf'
    doc: |
      UMAP, colored by grouping condition,
      split by cell cycle phase.
      PDF format
  ccpca_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: '*_ccpca_spl_cnd.png'
    doc: |
      PCA, colored by cell cycle phase,
      split by grouping condition.
      PNG format
  ccpca_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_ccpca_spl_cnd.pdf'
    doc: |
      PCA, colored by cell cycle phase,
      split by grouping condition.
      PDF format
  umap_gr_cnd_spl_mito_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_gr_cnd_spl_mito.png'
    doc: |
      UMAP, colored by grouping condition,
      split by mitochondrial percentage.
      PNG format
  umap_gr_cnd_spl_mito_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_gr_cnd_spl_mito.pdf'
    doc: |
      UMAP, colored by grouping condition,
      split by mitochondrial percentage.
      PDF format
  umap_gr_cnd_spl_umi_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_gr_cnd_spl_umi.png'
    doc: |
      UMAP, colored by grouping condition,
      split by RNA reads per cell.
      PNG format
  umap_gr_cnd_spl_umi_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_gr_cnd_spl_umi.pdf'
    doc: |
      UMAP, colored by grouping condition,
      split by RNA reads per cell.
      PDF format
  umap_gr_cnd_spl_gene_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_gr_cnd_spl_gene.png'
    doc: |
      UMAP, colored by grouping condition,
      split by genes per cell.
      PNG format
  umap_gr_cnd_spl_gene_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_gr_cnd_spl_gene.pdf'
    doc: |
      UMAP, colored by grouping condition,
      split by genes per cell.
      PDF format
  ucsc_cb_config_data:
    type: Directory?
    outputBinding:
      glob: '*_cellbrowser'
    doc: |
      UCSC Cell Browser configuration data.
  ucsc_cb_html_data:
    type: Directory?
    outputBinding:
      glob: '*_cellbrowser/html_data'
    doc: |
      UCSC Cell Browser html data.
  ucsc_cb_html_file:
    type: File?
    outputBinding:
      glob: '*_cellbrowser/html_data/index.html'
    doc: |
      UCSC Cell Browser html index.
  seurat_data_rds:
    type: File
    outputBinding:
      glob: '*_data.rds'
    doc: |
      Seurat object.
      RDS format
  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: '*_data.h5seurat'
    doc: |
      Seurat object.
      h5Seurat format
  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: '*_counts.h5ad'
    doc: |
      Seurat object.
      H5AD format
  seurat_data_cloupe:
    type: File?
    outputBinding:
      glob: '*_counts.cloupe'
    doc: |
      Seurat object.
      Loupe format
  seurat_data_scope:
    type: File?
    outputBinding:
      glob: '*_data.loom'
    doc: |
      Seurat object.
      SCope compatible.
      Loom format
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
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_rna_reduce.R"]:"/usr/local/bin/sc_rna_reduce.R")
stdout: sc_rna_reduce_stdout.log
stderr: sc_rna_reduce_stderr.log
label: Single-Cell RNA-Seq Dimensionality Reduction Analysis
doc: |
  Single-Cell RNA-Seq Dimensionality Reduction Analysis

  Integrates multiple single-cell RNA-Seq datasets, reduces dimensionality using PCA.
