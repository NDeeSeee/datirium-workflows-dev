cwlVersion: v1.0
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_root = function(basename) { return basename.split('.').slice(0,1).join('.'); };
sd:metadata:
- ../metadata/chipseq-header.cwl
sd:upstream:
  genome_indices: genome-indices.cwl
  control_file: trim-atacseq-pe.cwl
inputs:
  indices_folder:
    type: Directory
    sd:upstreamSource: genome_indices/bowtie_indices
    label: Indexed genome folder (bowtie)
    doc: Path to indexed genome folder by **bowtie**
  annotation_file:
    type: File
    sd:upstreamSource: genome_indices/annotation
    label: Annotation file
    format: http://edamontology.org/format_3475
    doc: Tab-separated annotation file
  genome_size:
    type: string
    sd:upstreamSource: genome_indices/genome_size
    label: Effective genome size
    doc: 'MACS2 effective genome size: hs, mm, ce, dm or number, for example 2.7e9'
  chrom_length:
    type: File
    sd:upstreamSource: genome_indices/chrom_length
    label: Chromosomes length file
    format: http://edamontology.org/format_2330
    doc: Chromosomes length file
  control_file:
    type: File?
    default: null
    sd:upstreamSource: control_file/bambai_pair
    sd:localLabel: true
    label: Use sample as a control (optional)
    format: http://edamontology.org/format_2572
    doc: If skipping upstream (i.e. not using an upstream sample as input here), please upload a BAM format file.
  broad_peak:
    type: boolean?
    default: false
    label: Callpeak broad
    doc: Set to call broad peak for MACS2
  fastq_file_upstream:
    type:
    - File
    - type: array
      items: File
    label: FASTQ read 1 input file(s)
    format: http://edamontology.org/format_1930
    doc: Reads data in a FASTQ format, received after paired end sequencing
  fastq_file_downstream:
    type:
    - File
    - type: array
      items: File
    label: FASTQ read 2 input file(s)
    format: http://edamontology.org/format_1930
    doc: Reads data in a FASTQ format, received after paired end sequencing
  exp_fragment_size:
    type: int?
    default: 150
    sd:layout:
      advanced: true
    label: Expected fragment size
    doc: Expected fragment size for MACS2
  force_fragment_size:
    type: boolean?
    default: false
    sd:layout:
      advanced: true
    label: Force fragment size
    doc: Force MACS2 to use exp_fragment_size
  clip_3p_end:
    type: int?
    default: 0
    sd:layout:
      advanced: true
    label: Clip from 3p end
    doc: Number of bases to clip from the 3p end
  clip_5p_end:
    type: int?
    default: 0
    sd:layout:
      advanced: true
    label: Clip from 5p end
    doc: Number of bases to clip from the 5p end
  minimum_length:
    type: int?
    default: 30
    label: Mimimum allowed read length after adapter trimming
    doc: |
      Discard reads that became shorter than the
      specified length because of either quality
      or adapter trimming. A value of 0 effectively
      disables this behaviour.
      Default: 30
    sd:layout:
      advanced: true
  remove_duplicates:
    type: boolean?
    default: false
    sd:layout:
      advanced: true
    label: Remove duplicates
    doc: Calls samtools rmdup to remove duplicates from sortesd BAM file
  peak_calling_fdr:
    type: float?
    default: 0.05
    sd:layout:
      advanced: true
    label: Minimum FDR (q-value) cutoff for peak detection
    doc: |
      Minimum FDR (q-value) cutoff for peak detection. -q, and
      -p are mutually exclusive.
  exclude_chromosome:
    type: string?
    default: chrM
    sd:layout:
      advanced: true
    label: Exclude chromosomes
    doc: Space separated list of chromosomes to be excluded
  do_not_scale:
    type: boolean?
    default: false
    sd:layout:
      advanced: true
    label: Do not scale genome coverage based on mapped reads number
    doc: When converting BAM to bigWig, the scale will be set to 1 by default
  promoter_dist:
    type: int?
    default: 1000
    sd:layout:
      advanced: true
    label: Max distance from gene TSS (in both direction) overlapping which the peak will be assigned to the promoter region
    doc: Max distance from gene TSS (in both direction) overlapping which the peak will be assigned to the promoter region
  upstream_dist:
    type: int?
    default: 20000
    sd:layout:
      advanced: true
    label: Max distance from the promoter (only in upstream direction) overlapping which the peak will be assigned to the upstream region
    doc: Max distance from the promoter (only in upstream direction) overlapping which the peak will be assigned to the upstream region
  threads:
    type: int?
    default: 2
    sd:layout:
      advanced: true
    doc: Number of threads for those steps that support multithreading
    label: Number of threads
outputs:
  unaligned_fastq:
    type:
    - 'null'
    - File[]
    format: http://edamontology.org/format_1930
    label: Unaligned FASTQ file(s)
    doc: Unaligned FASTQ file(s)
    outputSource: bowtie_aligner/unaligned_fastq
  multimapped_fastq:
    type:
    - 'null'
    - File[]
    format: http://edamontology.org/format_1930
    label: Multimapped FASTQ file(s)
    doc: Multimapped FASTQ file(s)
    outputSource: bowtie_aligner/multimapped_fastq
  bigwig:
    type: File
    format: http://edamontology.org/format_3006
    label: BigWig file
    doc: Generated BigWig file
    outputSource: bam_to_bigwig/bigwig_file
    sd:visualPlugins:
    - igvbrowser:
        tab: IGV Genome Browser
        id: igvbrowser
        type: wig
        name: BigWig Track
        height: 120
  fastx_statistics_upstream:
    type: File
    label: FASTQ 1 statistics
    format: http://edamontology.org/format_2330
    doc: fastx_quality_stats generated FASTQ 1 quality statistics file
    outputSource: fastx_quality_stats_upstream/statistics_file
    sd:visualPlugins:
    - line:
        tab: QC Plots
        Title: FASTQ 1 Base frequency plot
        xAxisTitle: Nucleotide position
        yAxisTitle: Frequency
        colors:
        - '#b3de69'
        - '#888888'
        - '#fb8072'
        - '#fdc381'
        - '#99c0db'
        data:
        - $13
        - $14
        - $15
        - $16
        - $17
    - boxplot:
        tab: QC Plots
        Title: FASTQ 1 Quality Control
        xAxisTitle: Nucleotide position
        yAxisTitle: Quality score
        colors:
        - '#b3de69'
        - '#888888'
        - '#fb8072'
        - '#fdc381'
        - '#99c0db'
        data:
        - $11
        - $7
        - $8
        - $9
        - $12
  fastx_statistics_downstream:
    type: File
    label: FASTQ 2 statistics
    format: http://edamontology.org/format_2330
    doc: fastx_quality_stats generated FASTQ 2 quality statistics file
    outputSource: fastx_quality_stats_downstream/statistics_file
    sd:visualPlugins:
    - line:
        tab: QC Plots
        Title: FASTQ 2 Base frequency plot
        xAxisTitle: Nucleotide position
        yAxisTitle: Frequency
        colors:
        - '#b3de69'
        - '#888888'
        - '#fb8072'
        - '#fdc381'
        - '#99c0db'
        data:
        - $13
        - $14
        - $15
        - $16
        - $17
    - boxplot:
        tab: QC Plots
        Title: FASTQ 2 Quality Control
        xAxisTitle: Nucleotide position
        yAxisTitle: Quality score
        colors:
        - '#b3de69'
        - '#888888'
        - '#fb8072'
        - '#fdc381'
        - '#99c0db'
        data:
        - $11
        - $7
        - $8
        - $9
        - $12
  bowtie_log:
    type: File
    label: BOWTIE alignment log
    format: http://edamontology.org/format_2330
    doc: BOWTIE generated alignment log
    outputSource: bowtie_aligner/log_file
  iaintersect_log:
    type: File
    label: Island intersect log
    format: http://edamontology.org/format_3475
    doc: Iaintersect generated log
    outputSource: island_intersect/log_file
  iaintersect_result:
    type: File
    label: Island intersect results
    format: http://edamontology.org/format_3475
    doc: Iaintersect generated results
    outputSource: island_intersect/result_file
    sd:visualPlugins:
    - syncfusiongrid:
        tab: Peak Calling
        Title: Islands list
  atdp_log_file:
    type: File
    label: ATDP log
    format: http://edamontology.org/format_3475
    doc: Average Tag Density generated log
    outputSource: average_tag_density/log_file
  atdp_result:
    type: File
    label: ATDP results
    format: http://edamontology.org/format_3475
    doc: Average Tag Density generated results
    outputSource: average_tag_density/result_file
    sd:visualPlugins:
    - scatter:
        tab: QC Plots
        Title: Average Tag Density
        xAxisTitle: Distance From TSS (bases)
        yAxisTitle: Average Tag Density (per bp)
        colors:
        - '#b3de69'
        height: 500
        data:
        - $1
        - $2
        comparable: atdp
  bambai_pair:
    type: File
    format: http://edamontology.org/format_2572
    label: Coordinate sorted BAM alignment file (+index BAI)
    doc: Coordinate sorted BAM file and BAI index file
    outputSource: samtools_remove_duplicates/deduplicated_bam_bai_pair
    sd:visualPlugins:
    - igvbrowser:
        tab: IGV Genome Browser
        id: igvbrowser
        optional: true
        type: alignment
        format: bam
        name: BAM Track
        displayMode: SQUISHED
  macs2_called_peaks:
    type: File?
    label: Called peaks
    format: http://edamontology.org/format_3468
    doc: XLS file to include information about called peaks
    outputSource: macs2_callpeak/peak_xls_file
  macs2_narrow_peaks:
    type: File?
    label: Narrow peaks
    format: http://edamontology.org/format_3613
    doc: Contains the peak locations together with peak summit, pvalue and qvalue
    outputSource: macs2_callpeak/narrow_peak_file
    sd:visualPlugins:
    - igvbrowser:
        tab: IGV Genome Browser
        id: igvbrowser
        type: annotation
        name: Narrow peaks
        displayMode: COLLAPSE
        height: 40
  macs2_broad_peaks:
    type: File?
    label: Broad peaks
    format: http://edamontology.org/format_3614
    doc: Contains the peak locations together with peak summit, pvalue and qvalue
    outputSource: macs2_callpeak/broad_peak_file
    sd:visualPlugins:
    - igvbrowser:
        tab: IGV Genome Browser
        id: igvbrowser
        type: annotation
        name: Broad peaks
        displayMode: COLLAPSE
        height: 40
  macs2_peak_summits:
    type: File?
    label: Peak summits
    format: http://edamontology.org/format_3003
    doc: Contains the peak summits locations for every peaks
    outputSource: macs2_callpeak/peak_summits_file
  macs2_moder_r:
    type: File?
    label: MACS2 generated R script
    format: http://edamontology.org/format_2330
    doc: R script to produce a PDF image about the model based on your data
    outputSource: macs2_callpeak/moder_r_file
  macs2_gapped_peak:
    type: File?
    label: Gapped peaks
    format: http://edamontology.org/format_3586
    doc: Contains both the broad region and narrow peaks
    outputSource: macs2_callpeak/gapped_peak_file
    sd:visualPlugins:
    - igvbrowser:
        tab: IGV Genome Browser
        id: igvbrowser
        type: annotation
        name: Gapped peaks
        displayMode: COLLAPSE
        height: 40
  macs2_log:
    type: File?
    label: MACS2 log
    format: http://edamontology.org/format_2330
    doc: MACS2 output log
    outputSource: macs2_callpeak/macs_log
  get_stat_log:
    type: File?
    label: YAML formatted combined log
    format: http://edamontology.org/format_3750
    doc: YAML formatted combined log
    outputSource: get_stat/collected_statistics_yaml
  get_stat_markdown:
    type: File?
    label: Markdown formatted combined log
    format: http://edamontology.org/format_3835
    doc: Markdown formatted combined log
    outputSource: get_stat/collected_statistics_md
    sd:visualPlugins:
    - markdownView:
        tab: Overview
  get_stat_formatted_log:
    type: File?
    label: Bowtie & Samtools Rmdup combined formatted log
    format: http://edamontology.org/format_3475
    doc: Processed and combined Bowtie aligner and Samtools rmdup formatted log
    outputSource: get_stat/collected_statistics_tsv
    sd:visualPlugins:
    - tableView:
        vertical: true
        tab: Overview
    sd:preview:
      sd:visualPlugins:
      - pie:
          colors:
          - '#b3de69'
          - '#99c0db'
          - '#fb8072'
          - '#fdc381'
          data:
          - $2
          - $3
          - $4
          - $5
  bam_statistics_report:
    type: File
    label: BAM statistics report (original)
    format: http://edamontology.org/format_2330
    doc: BAM statistics report (right after alignment and sorting)
    outputSource: get_bam_statistics/log_file
  bam_statistics_report_after_filtering:
    type: File
    label: BAM statistics report (after filtering)
    format: http://edamontology.org/format_2330
    doc: BAM statistics report (after all filters applied)
    outputSource: get_bam_statistics_after_filtering/log_file
  insert_size_report_after_filtering:
    type: File
    label: Insert size distribution report (after filtering)
    format: http://edamontology.org/format_3475
    doc: Insert size distribution report (after all filters applied)
    outputSource: get_bam_statistics_after_filtering/ext_is_section
    sd:visualPlugins:
    - scatter:
        tab: QC Plots
        Title: Insert Size Distribution (after filtering)
        xAxisTitle: Insert size
        yAxisTitle: Pairs total
        colors:
        - '#4b78a3'
        height: 500
        data:
        - $1
        - $2
        comparable: isdp
  macs2_fragment_stat:
    type: File?
    label: FRAGMENT, FRAGMENTE, ISLANDS
    format: http://edamontology.org/format_2330
    doc: fragment, calculated fragment, islands count from MACS2 results
    outputSource: macs2_callpeak/macs2_stat_file
  trim_report_upstream:
    type: File
    label: TrimGalore report FASTQ 1
    doc: TrimGalore generated log for FASTQ 1
    outputSource: trim_fastq/report_file
  trim_report_downstream:
    type: File
    label: TrimGalore report FASTQ 2
    doc: TrimGalore generated log for FASTQ 2
    outputSource: trim_fastq/report_file_pair
  preseq_estimates_plot_data:
    type: File?
    label: Preseq estimates
    format: http://edamontology.org/format_3475
    doc: Preseq estimated results
    outputSource: preseq_plot_data/estimates_file_plot_data
    sd:visualPlugins:
    - line:
        tab: QC Plots
        Title: Distinct Read Counts Estimates
        xAxisTitle: Mapped Reads/Fragments/Tags (millions)
        yAxisTitle: Distinct Reads Count
        colors:
        - '#4b78a3'
        - '#a3514b'
        height: 500
        data:
        - $2
        - $5
  estimated_fragment_size:
    type: int
    label: Estimated fragment size
    doc: Estimated fragment size for downstream analyses
    outputSource: macs2_callpeak/macs2_fragments_calculated
  mapped_reads_number:
    type: int
    label: Mapped reads number
    doc: Mapped reads number for downstream analyses
    outputSource: get_stat/mapped_reads
steps:
  extract_fastq_upstream:
    label: Loading unmapped sequence data for read 1
    doc: "Most DNA cores and commercial NGS companies return unmapped sequence data in FASTQ format.\nThe data can be uploaded from users computer, downloaded directly from an ftp server of\nthe core facility by providing a URL or from GEO by providing SRA accession number.  \n"
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_upstream
      output_prefix:
        default: read_1
    out:
    - fastq_file
  extract_fastq_downstream:
    label: Loading unmapped sequence data for read 2
    doc: "Most DNA cores and commercial NGS companies return unmapped sequence data in FASTQ format.\nThe data can be uploaded from users computer, downloaded directly from an ftp server of\nthe core facility by providing a URL or from GEO by providing SRA accession number.  \n"
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_downstream
      output_prefix:
        default: read_2
    out:
    - fastq_file
  trim_fastq:
    label: Adapter trimming
    doc: "For libraries sequenced on the Illumina platform it’s recommended to remove adapter sequences\nfrom the reads. If adapters are not trimmed there is a high risk of reads being unmapped to a\nreference genome. This becomes particularly important when the reads are long and the fragments\nare short - resulting in sequencing adapters at the end of read. If adapter trimming will cause\nall the reads become too short (<30bp), this step will be skipped.  \n"
    run: ../tools/trimgalore.cwl
    in:
      input_file: extract_fastq_upstream/fastq_file
      input_file_pair: extract_fastq_downstream/fastq_file
      dont_gzip:
        default: true
      length: minimum_length
      trim1:
        default: true
      paired:
        default: true
    out:
    - trimmed_file
    - trimmed_file_pair
    - report_file
    - report_file_pair
  bypass_trim:
    run: ../tools/bypass-trimgalore-pe.cwl
    in:
      original_fastq_file_1: extract_fastq_upstream/fastq_file
      trimmed_fastq_file_1: trim_fastq/trimmed_file
      trimming_report_file_1: trim_fastq/report_file
      original_fastq_file_2: extract_fastq_downstream/fastq_file
      trimmed_fastq_file_2: trim_fastq/trimmed_file_pair
      trimming_report_file_2: trim_fastq/report_file_pair
      min_reads_count:
        default: 100
    out:
    - selected_fastq_file_1
    - selected_report_file_1
    - selected_fastq_file_2
    - selected_report_file_2
  rename_upstream:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_1
      target_filename:
        source: extract_fastq_upstream/fastq_file
        valueFrom: $(self.basename)
    out:
    - target_file
  rename_downstream:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_2
      target_filename:
        source: extract_fastq_downstream/fastq_file
        valueFrom: $(self.basename)
    out:
    - target_file
  fastx_quality_stats_upstream:
    label: Quality control of unmapped sequence data for read 1
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_upstream/target_file
    out:
    - statistics_file
  fastx_quality_stats_downstream:
    label: Quality control of unmapped sequence data for read 2
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_downstream/target_file
    out:
    - statistics_file
  bowtie_aligner:
    label: Alignment to reference genome
    doc: |
      Aligns reads to the reference genome.
      Reads are assumed to be mapped if they
      have less than 3 mismatches.
      sam_file output includes both mapped
      and unmapped reads.
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: rename_upstream/target_file
      downstream_filelist: rename_downstream/target_file
      indices_folder: indices_folder
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
      v:
        default: 3
      m:
        default: 1
      best:
        default: true
      strata:
        default: true
      sam:
        default: true
      unaligned_prefix:
        default: unaligned_reads
      multimapped_prefix:
        default: multimapped_reads
      threads: threads
      q:
        default: true
      X:
        default: 500
    out:
    - sam_file
    - log_file
    - unaligned_fastq
    - multimapped_fastq
  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: bowtie_aligner/sam_file
      threads: threads
    out:
    - bam_bai_pair
  samtools_mark_duplicates:
    run: ../tools/samtools-markdup.cwl
    in:
      bam_bai_pair: samtools_sort_index/bam_bai_pair
      keep_duplicates:
        default: true
      threads: threads
    out:
    - deduplicated_bam_bai_pair
  filter_bam:
    run: ../tools/samtools-filter.cwl
    in:
      bam_bai_pair: samtools_mark_duplicates/deduplicated_bam_bai_pair
      exclude_chromosome: exclude_chromosome
    out:
    - filtered_bam_bai_pair
  clean_sam_headers_for_preseq:
    run: ../tools/samtools-clean-headers.cwl
    in:
      bam_file: filter_bam/filtered_bam_bai_pair
    out:
    - preseq_bam
  preseq:
    label: Sequencing depth estimation
    doc: |
      Estimates the complexity of the sequencing library, evaluates how many reads can
      be expected from the additional sequencing of the same experiment.
    run: ../tools/preseq-lc-extrap.cwl
    in:
      bam_file: clean_sam_headers_for_preseq/preseq_bam
      pe_mode:
        default: true
      extrapolation:
        default: 1000000000
    out:
    - estimates_file
    - log_file_stdout
    - log_file_stderr
  samtools_remove_duplicates:
    run: ../tools/samtools-markdup.cwl
    in:
      bam_bai_pair: filter_bam/filtered_bam_bai_pair
      keep_duplicates:
        source: remove_duplicates
        valueFrom: $(!self)
      threads: threads
    out:
    - deduplicated_bam_bai_pair
  macs2_callpeak:
    label: Peak detection
    doc: |
      Identifies enriched with aligned reads genome areas. Those areas correspond to the
      transcription factor binding sites.
    run: ../tools/macs2-callpeak-biowardrobe-only.cwl
    in:
      treatment_file: samtools_remove_duplicates/deduplicated_bam_bai_pair
      control_file: control_file
      nolambda:
        source: control_file
        valueFrom: $(!self)
      genome_size: genome_size
      mfold:
        default: 4 40
      verbose:
        default: 3
      nomodel: force_fragment_size
      extsize: exp_fragment_size
      bw: exp_fragment_size
      broad: broad_peak
      call_summits:
        source: broad_peak
        valueFrom: $(!self)
      keep_dup:
        default: auto
      q_value: peak_calling_fdr
      format_mode:
        default: BAMPE
      buffer_size:
        default: 10000
    out:
    - peak_xls_file
    - narrow_peak_file
    - peak_summits_file
    - broad_peak_file
    - moder_r_file
    - gapped_peak_file
    - treat_pileup_bdg_file
    - control_lambda_bdg_file
    - macs_log
    - macs2_stat_file
    - macs2_fragments_calculated
  bam_to_bigwig:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_remove_duplicates/deduplicated_bam_bai_pair
      chrom_length_file: chrom_length
      mapped_reads_number:
        source:
        - do_not_scale
        - get_stat/mapped_reads
        valueFrom: $(self[0]?null:self[1])
      fragment_size:
        default: 9
    out:
    - bigwig_file
  get_bam_statistics:
    label: Quality control of aligned sequence data
    doc: |
      Calculates alignment statistics, such as reads mapped/unmapped, average
      read length and quality score, etc.
    run: ../tools/samtools-stats.cwl
    in:
      bambai_pair: samtools_mark_duplicates/deduplicated_bam_bai_pair
      output_filename:
        source: samtools_mark_duplicates/deduplicated_bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_bam_statistics_report.txt")
    out:
    - log_file
  get_bam_statistics_after_filtering:
    run: ../tools/samtools-stats.cwl
    in:
      bambai_pair: samtools_remove_duplicates/deduplicated_bam_bai_pair
      output_filename:
        source: samtools_remove_duplicates/deduplicated_bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_bam_statistics_report_after_filtering.txt")
    out:
    - log_file
    - ext_is_section
    - reads_mapped
  get_stat:
    run: ../tools/collect-statistics-chip-seq.cwl
    in:
      trimgalore_report_fastq_1: bypass_trim/selected_report_file_1
      trimgalore_report_fastq_2: bypass_trim/selected_report_file_2
      bowtie_alignment_report: bowtie_aligner/log_file
      bam_statistics_report: get_bam_statistics/log_file
      bam_statistics_after_filtering_report: get_bam_statistics_after_filtering/log_file
      macs2_called_peaks: macs2_callpeak/peak_xls_file
      atdp_results: average_tag_density/result_file
      preseq_results: preseq/estimates_file
      paired_end:
        default: true
    out:
    - collected_statistics_yaml
    - collected_statistics_tsv
    - mapped_reads
    - collected_statistics_md
  preseq_plot_data:
    label: Formats sequencing depth estimation data for plotting
    doc: |
      Formats estimates file from preseq standard output for QC plotting. This adds a new
      column that includes the actual read count point on the plot.
    run: ../tools/preseq-plot-data.cwl
    in:
      preseq_stderr_log_file: preseq/log_file_stderr
      estimates_file: preseq/estimates_file
      mapped_reads: get_stat/mapped_reads
    out:
    - estimates_file_plot_data
  island_intersect:
    label: Peak annotation
    doc: |
      Assigns nearest genes to peaks to explore the biological implication of the open
      chromatin binding sites.
    run: ../tools/iaintersect.cwl
    in:
      input_filename: macs2_callpeak/peak_xls_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
    out:
    - result_file
    - log_file
  samtools_sort_index_for_atdp:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: samtools_remove_duplicates/deduplicated_bam_bai_pair
      threads: threads
    out:
    - bam_bai_pair
  average_tag_density:
    label: Read enrichment around genes TSS
    doc: |
      Generates average tag density plot around genes TSS as a lot of cis-regulatory
      elements are close to the TSS of their targets.
    run: ../tools/atdp.cwl
    in:
      input_file: samtools_sort_index_for_atdp/bam_bai_pair
      annotation_filename: annotation_file
      fragmentsize_bp: macs2_callpeak/macs2_fragments_calculated
      avd_window_bp:
        default: 5000
      avd_smooth_bp:
        default: 50
      ignore_chr:
        default: chrM
      double_chr:
        default: chrX chrY
      avd_heat_window_bp:
        default: 200
      mapped_reads:
        source: get_bam_statistics_after_filtering/reads_mapped
        valueFrom: $(parseInt(self/2))
    out:
    - result_file
    - log_file
label: Trim Galore ATAC-Seq pipeline paired-end
doc: "This ATAC pipeline is based on original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)\n**ChIP-Seq** basic analysis workflow for a **paired-end** experiment with Trim Galore. The pipeline was adapted for\nATAC-Seq single-read data analysis by updating genome coverage step.\n### Data Analysis Steps\nSciDAP starts from the .fastq files which most DNA cores and commercial NGS companies return. Starting from raw data allows us to ensure that all experiments have been processed in the same way and simplifies the deposition of data to GEO upon publication. The data can be uploaded from users computer, downloaded directly from an ftp server of the core facility by providing a URL or from GEO by providing SRA accession number.\nOur current pipelines include the following steps:\n1. Trimming the adapters with TrimGalore. This step is particularly important when the reads are long and the fragments are short as in ATAC -resulting in sequencing adapters at the end of read. If adapter is not removed the read will not map. TrimGalore can recognize standard adapters, such as Nexterra/Tn5 adapters.\n2. QC\n3. (Optional) trimming adapters on 5' or 3' end by the specified number of bases.\n4. Mapping reads with BowTie. Only uniquely mapped reads with less than 3 mismatches are used in the downstream analysis. Results are saved as a .bam file.\n5. Reads mapping to chromosome M are removed. Since there are many copies of chromosome M in the cell and it is not protected by histones, some ATAC libraries have up to 50% of reads mapping to chrM. We recommend using OMNI-ATAC protocol that reduces chrM reads and provides better specificity. \n6.  (Optional) Removal of duplicates (reads/pairs of reads mapping to exactly same location). This step is used to remove reads overamplified in PCR. Unfortunately, it may also remove \"good\" reads. We usually do not remove duplicates unless the library is heavily duplicated. Please note that MACS2 will remove 'excessive' duplicates during peak calling ina smart way (those not supported by other nearby reads).\n7.  Peakcalling by MACS2. (Optionally), it is possible to specify read extension length for MACS2 to use if the length determined automatically is wrong. \n8.  Generation of BigWig coverage files for display on the browser. Since the cuts by the Tn5 transposome are 9bp apart, we show coverage by 9bp reads rather than fragments as in ChIP-Seq. The coverage shows the number of fragments at each base in the genome normalized to the number of millions of mapped reads. This way the peak of coverage will be located at the most accessible site. \n\n### Details\n\n\n_Trim Galore_ is a wrapper around [Cutadapt](https://github.com/marcelm/cutadapt)\nand [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to consistently\napply adapter and quality trimming to FastQ files, with extra functionality for RRBS data.\n\nA [FASTQ](http://maq.sourceforge.net/fastq.shtml) input file has to be provided.\n\n\nIn outputs it returns coordinate sorted BAM file alongside with index BAI file,\nquality statistics for both the input FASTQ files, reads coverage in a form of BigWig file,\npeaks calling data in a form of narrowPeak or broadPeak files, islands with the assigned nearest\ngenes and region type, data for average tag density plot (on the base of BAM file).\n\nWorkflow starts with running fastx_quality_stats (steps fastx_quality_stats_upstream and\nfastx_quality_stats_downstream) from FASTX-Toolkit to calculate quality statistics for both upstream\nand downstream input FASTQ files. At the same time Bowtie is used to align reads from input FASTQ\nfiles to reference genome (Step bowtie_aligner). The output of this step is unsorted SAM file which\nis being sorted and indexed by samtools sort and samtools index (Step samtools_sort_index).\nDepending on workflow’s input parameters indexed and sorted BAM file\ncould be processed by `samtools markdup` *samtools\\_remove\\_duplicates* to remove all possible read duplicates.\nRight after that macs2 callpeak performs\npeak calling (Step macs2_callpeak). On the base of returned outputs the next step\n(Step macs2_island_count) calculates the number of islands and estimated fragment size. If the last\none is less that 80 (hardcoded in a workflow) macs2 callpeak is rerun again with forced fixed\nfragment size value (Step macs2_callpeak_forced). If at the very beginning it was set in workflow\ninput parameters to force run peak calling with fixed fragment size, this step is skipped and the\noriginal peak calling results are saved. In the next step workflow again calculates the number\nof islands and estimated fragment size (Step macs2_island_count_forced) for the data obtained from\nmacs2_callpeak_forced step. If the last one was skipped the results from macs2_island_count_forced step\nare equal to the ones obtained from macs2_island_count step.\nNext step (Step macs2_stat) is used to define which of the islands and estimated fragment size should be used\nin workflow output: either from macs2_island_count step or from macs2_island_count_forced step. If input\ntrigger of this step is set to True it means that macs2_callpeak_forced step was run and it returned different\nfrom macs2_callpeak step results, so macs2_stat step should return [fragments_new, fragments_old, islands_new],\nif trigger is False the step returns [fragments_old, fragments_old, islands_old], where sufix \"old\" defines\nresults obtained from macs2_island_count step and sufix \"new\" - from macs2_island_count_forced step.\nThe following two steps (Step bamtools_stats and bam_to_bigwig) are used to calculate coverage on the base\nof input BAM file and save it in BigWig format. For that purpose bamtools stats returns the number of\nmapped reads number which is then used as scaling factor by bedtools genomecov when it performs coverage\ncalculation and saves it in BED format. The last one is then being sorted and converted to BigWig format by\nbedGraphToBigWig tool from UCSC utilities. To adapt the pipeline for ATAC-Seq data analysis we calculate genome\ncoverage using only the first 9 bp from every read. Step get_stat is used to return a text file with statistics\nin a form of [TOTAL, ALIGNED, SUPRESSED, USED] reads count. Step island_intersect assigns genes and regions\nto the islands obtained from macs2_callpeak_forced. Step average_tag_density is used to calculate data for\naverage tag density plot on the base of BAM file.\n"
sd:version: 100
