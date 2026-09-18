# Run all Downstream analysis pipeline

Uses the edgeR differentiall expression results and the formatted
Deseq2Dataset Performs volcano_plot.R which filters for the chosen
log2fc and pval_adj

## Usage

``` r
run_de_pipeline(
  dds_object,
  edgeR_results,
  deg_log2fc_thresh = 1,
  deg_pval_adj_thresh = 0.05,
  vol_plot_lower_xlim = -10,
  vol_plot_upper_xlim = 10,
  cartesian_step = 5,
  ora_convert_ids = FALSE,
  ora_pval_adj_threshold = 0.05,
  orgdb = org.Hs.eg.db,
  ensembl_dataset = "hsapiens_gene_ensembl"
)
```

## Arguments

- dds_object:

  A formatted DESeq2 object where each sample and its counts are aligned
  with the respective metadata in `@colData`.

- edgeR_results:

  dataframe of log2fc results from edge

- deg_log2fc_thresh:

  absolute value cutoff of log2fc

- deg_pval_adj_thresh:

  pvalue cutoff for initial log2fc values

- vol_plot_lower_xlim:

  Numeric value specifying the lower limit of the x-axis. Should be less
  than 0. For example, `-5`.

- vol_plot_upper_xlim:

  Numeric value specifying the upper limit of the x-axis. Should be
  greater than 0. For example, `5`.

- cartesian_step:

  Numeric value specifying the interval between x-axis ticks.

- ora_convert_ids:

  Boolean option to convert the gene IDs from Ensemble to External
  format

- ora_pval_adj_threshold:

  pvalue cutoff of the over-representation test

- orgdb:

  An organism-specific annotation database used to facilitate gene
  identifier conversion, such as `org.Hs.eg.db` for human or
  `org.Mm.eg.db` for mouse.

- ensembl_dataset:

  Character string specifying the Ensembl dataset used for annotation.
  For example, `"hsapiens_gene_ensembl"` for human genes.

## Value

A list of the following objects

- filtered_results: The initial filtering of DEGs log2fc and
  pval_adj_threshold

- volcano_plots: visualization of the filtered results on a volcano plot

- ora_up_raw: Result of the over-representation analysis on upregulated
  DEGs via clusterProfiler

- ora_down_raw: Result of the over-represenation analysis on
  downregulated DEGs via clusterProfiler

- unpacked_results_up: unpacked results table consisting of a nested
  list of ORA results on the upregulated DEGs

- unpacked_results_down: unpacked results table consisting of a nested
  list of ORA results on the down-regulated DEGs
