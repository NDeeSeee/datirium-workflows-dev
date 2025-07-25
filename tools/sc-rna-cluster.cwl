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
      expression information stored in the RNA assay, as well as 'pca' and 'rnaumap'
      dimensionality reductions applied to that assay.
  dimensions:
    type: int?
    inputBinding:
      prefix: --dimensions
    doc: |
      Dimensionality to use when constructing nearest-neighbor
      graph before clustering (from 1 to 50). Set to 0 to use
      auto-estimated dimensionality.
      Default: 10
  cluster_metric:
    type:
    - 'null'
    - type: enum
      symbols:
      - euclidean
      - cosine
      - manhattan
      - hamming
    inputBinding:
      prefix: --ametric
    doc: |
      Distance metric used when constructing nearest-neighbor graph before clustering.
      Default: euclidean
  cluster_algorithm:
    type:
    - 'null'
    - type: enum
      symbols:
      - louvain
      - mult-louvain
      - slm
      - leiden
    inputBinding:
      prefix: --algorithm
    doc: |
      Algorithm for modularity optimization when running clustering.
      Default: louvain
  resolution:
    type:
    - 'null'
    - float
    - float[]
    inputBinding:
      prefix: --resolution
    doc: |
      Clustering resolution applied to the constructed nearest-neighbor graph.
      Can be set as an array but only the first item from the list will be used
      for cluster labels and gene markers in the UCSC Cell Browser when running
      with --cbbuild and --diffgenes parameters.
      Default: 0.3, 0.5, 1.0
  genes_of_interest:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: --genes
    doc: |
      Genes of interest to build genes expression plots.
      Default: None
  genesets_data:
    type: File?
    inputBinding:
      prefix: --genesets
    doc: |
      Path to the GMT file for calculating average expression levels
      (module scores) per gene set. This file can be downloaded from
      the Molecular Signatures Database (MSigDB) following the link
      https://www.gsea-msigdb.org/gsea/msigdb.
      Default: do not calculate gene set expression scores.
  identify_diff_genes:
    type: boolean?
    inputBinding:
      prefix: --diffgenes
    doc: |
      Identify differentially expressed genes (putative gene markers) between each
      pair of clusters for all resolutions.
      Default: false
  minimum_logfc:
    type: float?
    inputBinding:
      prefix: --logfc
    doc: |
      For putative gene markers identification include only those genes that
      on average have log fold change difference in expression between every
      tested pair of clusters not lower than this value. Ignored if '--diffgenes'
      is not set.
      Default: 0.25
  minimum_pct:
    type: float?
    inputBinding:
      prefix: --minpct
    doc: |
      For putative gene markers identification include only those genes that
      are detected in not lower than this fraction of cells in either of the
      two tested clusters. Ignored if '--diffgenes' is not set.
      Default: 0.1
  only_positive_diff_genes:
    type: boolean?
    inputBinding:
      prefix: --onlypos
    doc: |
      For putative gene markers identification return only positive markers.
      Ignored if '--diffgenes' is not set.
      Default: false
  test_to_use:
    type:
    - 'null'
    - type: enum
      symbols:
      - wilcox
      - bimod
      - roc
      - t
      - negbinom
      - poisson
      - LR
      - MAST
      - DESeq2
    inputBinding:
      prefix: --testuse
    doc: |
      Statistical test to use for putative gene markers identification.
      Ignored if '--diffgenes' is not set.
      Default: wilcox
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
  cell_cnts_gr_clst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_cell_cnts_gr_clst_res_*.png'
    doc: |
      Number of cells per cluster.
      All cells; all resolutions.
      PNG format.
  gene_umi_spl_clst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_gene_umi_spl_clst_res_*.png'
    doc: |
      Genes vs RNA reads per cell.
      Split by cluster; all cells;
      all resolutions.
      PNG format.
  umi_mito_spl_clst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_umi_mito_spl_clst_res_*.png'
    doc: |
      RNA reads vs mitochondrial % per cell.
      Split by cluster; all cells; all
      resolutions.
      PNG format.
  rnadbl_gr_clst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_rnadbl_gr_clst_res_*.png'
    doc: |
      Percentage of RNA doublets per cluster.
      All cells; all resolutions.
      PNG format.
  qc_mtrcs_dnst_gr_clst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_qc_mtrcs_dnst_gr_clst_res_*.png'
    doc: |
      Distribution of QC metrics per cell
      colored by cluster.
      All cells; all resolutions.
      PNG format.
  umap_gr_ph_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_gr_ph_spl_idnt.png'
    doc: |
      UMAP colored by cell cycle phase.
      Split by dataset; downsampled to the
      smallest dataset.
      PNG format.
  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: '*_cmp_gr_ph_spl_idnt.png'
    doc: |
      Composition plot colored by cell cycle phase.
      Split by dataset; downsampled to the smallest
      dataset.
      PNG format.
  umap_gr_ph_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_gr_ph_spl_cnd.png'
    doc: |
      UMAP colored by cell cycle phase.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.
  cmp_gr_ph_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: '*_cmp_gr_ph_spl_cnd.png'
    doc: |
      Composition plot colored by cell cycle phase.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.
  umap_gr_clst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_umap_gr_clst_res_*.png'
    doc: |
      UMAP colored by cluster.
      All cells; all resolutions.
      PNG format.
  slh_gr_clst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_slh_gr_clst_res_*.png'
    doc: |
      Silhouette scores.
      All cells; all resolutions.
      PNG format.
  umap_gr_clst_spl_idnt_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_umap_gr_clst_spl_idnt_res_*.png'
    doc: |
      UMAP colored by cluster.
      Split by dataset; downsampled to the
      smallest dataset; all resolutions.
      PNG format.
  cmp_gr_clst_spl_idnt_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_cmp_gr_clst_spl_idnt_res_*.png'
    doc: |
      Composition plot colored by cluster.
      Split by dataset; downsampled to the
      smallest dataset; all resolutions.
      PNG format.
  cmp_gr_idnt_spl_clst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_cmp_gr_idnt_spl_clst_res_*.png'
    doc: |
      Composition plot colored by dataset.
      Split by cluster; downsampled to the
      smallest dataset; all resolutions.
      PNG format.
  umap_gr_clst_spl_ph_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_umap_gr_clst_spl_ph_res_*.png'
    doc: |
      UMAP colored by cluster.
      Split by cell cycle phase; downsampled
      to the smallest dataset (if multiple
      datasets are analyzed jointly); all
      resolutions.
      PNG format.
  cmp_gr_ph_spl_clst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_cmp_gr_ph_spl_clst_res_*.png'
    doc: |
      Composition plot colored by cell cycle phase.
      Split by cluster; downsampled to the smallest
      dataset (if multiple datasets are analyzed
      jointly); all resolutions.
      PNG format.
  umap_gr_clst_spl_cnd_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_umap_gr_clst_spl_cnd_res_*.png'
    doc: |
      UMAP colored by cluster.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group; all resolutions.
      PNG format.
  cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_cmp_gr_clst_spl_cnd_res_*.png'
    doc: |
      Composition plot colored by cluster.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group; all resolutions.
      PNG format.
  cmp_gr_cnd_spl_clst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_cmp_gr_cnd_spl_clst_res_*.png'
    doc: |
      Composition plot colored by grouping condition.
      Split by cluster; first downsampled to the
      smallest dataset, then downsampled to the
      smallest group; all resolutions.
      PNG format.
  gse_per_cell_plot_png:
    type: File?
    outputBinding:
      glob: '*_gse_per_cell.png'
    doc: |
      UMAP colored by gene set expression score.
      PNG format.
  gse_avg_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_gse_avg_res_*.png'
    doc: |
      Average gene set expression score.
      All resolutions.
      PNG format.
  gse_dnst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_gse_dnst_res_*.png'
    doc: |
      Gene set expression score density.
      All resolutions.
      PNG format.
  xpr_per_cell_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_per_cell_[!sgnl_]*.png'
    doc: |
      UMAP colored by gene expression.
      All genes of interest.
      PNG format.
  xpr_per_cell_sgnl_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_per_cell_sgnl_*.png'
    doc: |
      UMAP colored by gene expression density.
      All genes of interest.
      PNG format.
  xpr_avg_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_avg_res_*.png'
    doc: |
      Average gene expression.
      All resolutions.
      PNG format.
  xpr_dnst_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_dnst_res_*.png'
    doc: |
      Gene expression density.
      All resolutions.
      PNG format.
  xpr_htmp_res_plot_png:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_htmp_res_*.png'
    doc: |
      Gene expression heatmap.
      Top gene markers; all resolutions.
      PNG format.
  all_plots_pdf:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*.pdf'
    doc: |
      All generated plots.
      PDF format.
  xpr_htmp_res_tsv:
    type:
    - 'null'
    - type: array
      items: File
    outputBinding:
      glob: '*_xpr_htmp_res_*.tsv'
    doc: |
      Gene expression heatmap.
      Top gene markers; all resolutions.
      TSV format.
  gene_markers_tsv:
    type: File?
    outputBinding:
      glob: '*_gene_markers.tsv'
    doc: |
      Gene markers.
      All resolutions.
      TSV format.
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
      RDS format.
  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: '*_data.h5seurat'
    doc: |
      Seurat object.
      h5Seurat format.
  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: '*_counts.h5ad'
    doc: |
      Seurat object.
      H5AD format.
  seurat_data_cloupe:
    type: File?
    outputBinding:
      glob: '*_counts.cloupe'
    doc: |
      Seurat object.
      Loupe format.
  seurat_data_scope:
    type: File?
    outputBinding:
      glob: '*_data.loom'
    doc: |
      Seurat object.
      SCope compatible.
      Loom format.
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
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_rna_cluster.R"]:"/usr/local/bin/sc_rna_cluster.R")
stdout: sc_rna_cluster_stdout.log
stderr: sc_rna_cluster_stderr.log
label: Single-Cell RNA-Seq Cluster Analysis
doc: |
  Single-Cell RNA-Seq Cluster Analysis

  Clusters single-cell RNA-Seq datasets, identifies gene markers.
