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
      Path to the RDS file to load Seurat object from. This
      file should include genes expression information
      stored in the RNA assay. Additionally, rnaumap,
      and/or atacumap, and/or wnnumap dimensionality
      reductions should be present.
  datasets_metadata:
    type: File?
    inputBinding:
      prefix: --metadata
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat
      object metadata with categorical values using samples
      identities. First column - library_id should
      correspond to all unique values from the new.ident
      column of the loaded Seurat object. If any of the
      provided in this file columns are already present in
      the Seurat object metadata, they will be overwritten.
      When combined with --barcodes parameter, first the
      metadata will be extended, then barcode filtering will
      be applied. Default: no extra metadata is added
  barcodes_data:
    type: File?
    inputBinding:
      prefix: --barcodes
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata by selected barcodes.
      First column should be named as barcode. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present. Default: all cells used, no extra
      metadata is added
  groupby:
    type: string?
    inputBinding:
      prefix: --groupby
    doc: |
      Column from the Seurat object metadata to group cells
      for optional subsetting when combined with --subset
      parameter. May be one of the extra metadata columns
      added with --metadata or --barcodes parameters.
      Ignored if --subset is not set. Default: do not
      subset, include all cells into analysis.
  subset:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: --subset
    doc: |
      Values from the column set with --groupby parameter to
      subset cells before running differential expression
      analysis. Ignored if --groupby is not provided.
      Default: do not subset cells, include all of them
  splitby:
    type: string
    inputBinding:
      prefix: --splitby
    doc: |
      Column from the Seurat object metadata to split cells
      into two groups to run --second vs --first
      differential expression analysis. May be one of the
      extra metadata columns added with --metadata or
      --barcodes parameters.
  first_cond:
    type: string
    inputBinding:
      prefix: --first
    doc: |
      Value from the Seurat object metadata column set with
      --splitby parameter to define the first group of cells
      for differential expression analysis.
  second_cond:
    type: string
    inputBinding:
      prefix: --second
    doc: |
      Value from the Seurat object metadata column set with
      --splitby parameter to define the second group of
      cells for differential expression analysis.
  analysis_method:
    type:
    - 'null'
    - type: enum
      symbols:
      - wilcoxon
      - likelihood-ratio
      - t-test
      - negative-binomial
      - poisson
      - logistic-regression
      - mast
      - deseq
      - deseq-lrt
    inputBinding:
      prefix: --test
    doc: |
      Test type to use in differential expression analysis.
      If set to deseq or deseq-lrt, gene expression will be
      aggregated to the pseudobulk level per dataset. For
      deseq, the pair-wise Wald test will be used. For
      deseq-lrt, the reduced formula will look like ~1 if
      --batchby parameter is omitted or will be set to
      ~batchby to exclude the criteria if interest (defined
      by --splitby). For all other values of the --test
      parameter the FindMarkers function will be used (genes
      will be prefiltered by minimum percentage >= 0.1 and
      by minimum log2FoldChange >= 0.25 before running
      differential expression analysis). Default: use
      FindMarkers with Wilcoxon Rank Sum test.
  batchby:
    type: string?
    inputBinding:
      prefix: --batchby
    doc: |
      Column from the Seurat object metadata to group cells
      into batches. If --test is set to deseq or deseq-lrt
      the --batchby parameter will be used in the design
      formula in the following way ~splitby+batchby. If
      --test is set to negative-binomial, poisson, logistic-
      regression, or mast it will be used as a latent
      variable in the FindMarkers function. Not supported
      for --test values equal to wilcoxon, likelihood-ratio,
      or t-test. May be one of the extra metadata columns
      added with --metadata or --barcodes parameters.
      Default: do not model batch effect.
  maximum_padj:
    type: float?
    inputBinding:
      prefix: --padj
    doc: |
      In the exploratory visualization part of the analysis
      output only differentially expressed genes with
      adjusted P-value not bigger than this value.
      Default: 0.05
  minimum_pct:
    type: float?
    inputBinding:
      prefix: --minpct
    doc: |
      Include only those genes that are detected in not lower than this
      fraction of cells in either of the two tested conditions.
      Default: 0.1
  genes_of_interest:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: --genes
    doc: |
      Genes of interest to label on the generated plots.
      Default: top 10 genes with the highest and the lowest
      log2FoldChange values.
  exclude_pattern:
    type: string?
    inputBinding:
      prefix: --exclude
    doc: |
      Regex pattern to identify and exclude specific genes
      from the differential expression analysis (not case-
      sensitive). If any of such genes are provided in the
      --genes parameter, they will be excluded from there as
      well. Default: use all genes
  cluster_method:
    type:
    - 'null'
    - type: enum
      symbols:
      - row
      - column
      - both
    inputBinding:
      prefix: --cluster
    doc: |
      Hopach clustering method to be run on the normalized
      read counts for the exploratory visualization part of
      the analysis. Clustering by column is supported only
      when --test is set to deseq or deseq-lrt. Default: do
      not run clustering
  row_distance:
    type:
    - 'null'
    - type: enum
      symbols:
      - cosangle
      - abscosangle
      - euclid
      - abseuclid
      - cor
      - abscor
    inputBinding:
      prefix: --rowdist
    doc: |
      Distance metric for HOPACH row clustering. Ignored if
      --cluster is set to column or not provided.
      Default: cosangle
  column_distance:
    type:
    - 'null'
    - type: enum
      symbols:
      - cosangle
      - abscosangle
      - euclid
      - abseuclid
      - cor
      - abscor
    inputBinding:
      prefix: --columndist
    doc: |
      Distance metric for HOPACH column clustering. Ignored
      if --cluster is set to row or not provided.
      Default: euclid
  center_row:
    type: boolean?
    inputBinding:
      prefix: --center
    doc: |
      Apply mean centering for gene expression prior to
      running clustering by row. Ignored if --cluster is
      set to column or not provided. Default: do not
      center
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
      Color theme for all generated plots. One of gray, bw,
      linedraw, light, dark, minimal, classic, void.
      Default: classic
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
      Maximum memory in GB allowed to be shared between
      the workers when using multiple --cpus.
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
  umap_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_rd_rnaumap.png'
    doc: |
      UMAP with cells selected for analysis.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction rnaumap.
      PNG format.
  umap_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_rd_rnaumap.pdf'
    doc: |
      UMAP with cells selected for analysis.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction rnaumap.
      PDF format.
  umap_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_rd_atacumap.png'
    doc: |
      UMAP with cells selected for analysis.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction atacumap.
      PNG format.
  umap_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_rd_atacumap.pdf'
    doc: |
      UMAP with cells selected for analysis.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction atacumap.
      PDF format.
  umap_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_rd_wnnumap.png'
    doc: |
      UMAP with cells selected for analysis.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction wnnumap.
      PNG format.
  umap_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_rd_wnnumap.pdf'
    doc: |
      UMAP with cells selected for analysis.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction wnnumap.
      PDF format.
  mds_plot_html:
    type: File?
    outputBinding:
      glob: '*_mds_plot.html'
    doc: |
      MDS plot of pseudobulk aggregated
      normalized reads counts.
      HTML format.
  pca_1_2_plot_png:
    type: File?
    outputBinding:
      glob: '*_pca_1_2.png'
    doc: |
      Gene expression PCA (1,2).
      PNG format.
  pca_1_2_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_pca_1_2.pdf'
    doc: |
      Gene expression PCA (1,2).
      PDF format.
  pca_2_3_plot_png:
    type: File?
    outputBinding:
      glob: '*_pca_2_3.png'
    doc: |
      Gene expression PCA (2,3).
      PNG format
  pca_2_3_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_pca_2_3.pdf'
    doc: |
      Gene expression PCA (2,3).
      PDF format.
  dxpr_vlcn_plot_png:
    type: File?
    outputBinding:
      glob: '*_dxpr_vlcn.png'
    doc: |
      Differentially expressed genes.
      Volcano plot of differentially expressed genes.
      Highlighed genes are either provided by user or
      top 10 genes with the highest log2FoldChange
      values. The direction of comparison is defined
      as --second vs --first. Cells are optionally
      subsetted to the specific group and optionally
      coerced to the pseudobulk form.
      PNG format.
  dxpr_vlcn_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_dxpr_vlcn.pdf'
    doc: |
      Differentially expressed genes.
      Volcano plot of differentially expressed genes.
      Highlighed genes are either provided by user or
      top 10 genes with the highest log2FoldChange
      values. The direction of comparison is defined
      as --second vs --first. Cells are optionally
      subsetted to the specific group and optionally
      coerced to the pseudobulk form.
      PDF format.
  xpr_dnst_plot_png:
    type: File?
    outputBinding:
      glob: '*_xpr_dnst.png'
    doc: |
      Gene expression density.
      Gene expression violin plots for either user
      provided or top 10 differentially expressed
      genes with the highest log2FoldChange values.
      The direction of comparison is defined as
      --second vs --first.
      PNG format.
  xpr_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_xpr_dnst.pdf'
    doc: |
      Gene expression density.
      Gene expression violin plots for either user
      provided or top 10 differentially expressed
      genes with the highest log2FoldChange values.
      The direction of comparison is defined as
      --second vs --first.
      PDF format.
  xpr_per_cell_rd_rnaumap_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_per_cell_rd_rnaumap_*.png'
    doc: |
      UMAP colored by gene expression.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction rnaumap.
      PNG format.
  xpr_per_cell_rd_rnaumap_plot_pdf:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_per_cell_rd_rnaumap_*.pdf'
    doc: |
      UMAP colored by gene expression.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction rnaumap.
      PDF format.
  xpr_per_cell_rd_atacumap_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_per_cell_rd_atacumap_*.png'
    doc: |
      UMAP colored by gene expression.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction atacumap.
      PNG format.
  xpr_per_cell_rd_atacumap_plot_pdf:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_per_cell_rd_atacumap_*.pdf'
    doc: |
      UMAP colored by gene expression.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction atacumap.
      PDF format.
  xpr_per_cell_rd_wnnumap_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_per_cell_rd_wnnumap_*.png'
    doc: |
      UMAP colored by gene expression.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction wnnumap.
      PNG format.
  xpr_per_cell_rd_wnnumap_plot_pdf:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_per_cell_rd_wnnumap_*.pdf'
    doc: |
      UMAP colored by gene expression.
      Split by selected criteria; optionally
      subsetted to the specific group;
      reduction wnnumap.
      PDF format.
  xpr_htmp_plot_png:
    type: File?
    outputBinding:
      glob: '*_xpr_htmp.png'
    doc: |
      Gene expression heatmap.
      Filtered by adjusted p-value; optionally
      subsetted to the specific groups.
      PNG format.
  xpr_htmp_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_xpr_htmp.pdf'
    doc: |
      Gene expression heatmap.
      Filtered by adjusted p-value; optionally
      subsetted to the specific groups.
      PDF format.
  diff_expr_genes:
    type: File?
    outputBinding:
      glob: '*_de_genes.tsv'
    doc: |
      Differentially expressed genes.
      Not filtered by adjusted p-value.
      TSV format.
  bulk_read_counts_gct:
    type: File?
    outputBinding:
      glob: '*_bulk_counts.gct'
    doc: |
      GSEA compatible not filtered normalized
      reads counts aggregated to pseudobulk
      form.
      GCT format.
  bulk_phenotypes_cls:
    type: File?
    outputBinding:
      glob: '*_bulk_phntps.cls'
    doc: |
      GSEA compatible phenotypes file defined
      based on --splitby, --first, and --second
      parameters.
      CLS format.
  cell_read_counts_gct:
    type: File?
    outputBinding:
      glob: '*_cell_counts.gct'
    doc: |
      Filtered normalized reads counts per cell.
      GCT format.
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
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_rna_de_pseudobulk.R"]:"/usr/local/bin/sc_rna_de_pseudobulk.R")
stdout: sc_rna_de_pseudobulk_stdout.log
stderr: sc_rna_de_pseudobulk_stderr.log
label: Single-Cell RNA-Seq Differential Expression Analysis
doc: |
  Single-Cell RNA-Seq Differential Expression Analysis

  Identifies differentially expressed genes between any
  two groups of cells, optionally aggregating gene
  expression data from single-cell to pseudobulk form.
