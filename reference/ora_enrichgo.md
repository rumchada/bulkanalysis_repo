# Over-Representation Analysis Using enrichGO

Performs Gene Ontology (GO) over-representation analysis (ORA) using
[`enrichGO`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html) on
differential expression results. Genes can be analyzed separately based
on their direction of differential expression (upregulated or
downregulated).

## Usage

``` r
ora_enrichgo(
  filtered_results,
  direction = c("up", "down"),
  log2fc_thresh = 1,
  pval_thresh = 0.05,
  top_terms = 15,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP"
)
```

## Arguments

- filtered_results:

  list of tables or a table of differential expression results

- direction:

  Which direction +log2fc -log2fc should be looked at ? e.g. `"up"`,
  `"down"`

- log2fc_thresh:

  what is log2FC cutoff you would want for the analysis

- pval_thresh:

  after running fisher exact test for ORA, what should the distribution
  of pvalues be?

- top_terms:

  How many resulting terms should be display in the visuals.

- OrgDb:

  rganism database. Default is human: `org.Hs.eg.db`. Other organism
  annotation databases are available from
  https://bioconductor.org/packages/3.23/data/annotation/.

- keyType:

  which type of geneid are being entered into the cluster

- ont:

  which ontological DB will you be using ? options `c("BP", "CC", "MF")`

## Value

A list containing:

- `enrichgo_results`: Results from the GO over-representation analysis
  returned by
  [`enrichGO`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

- `visuals`: A collection of visualizations summarizing the enriched GO
  terms, including dotplots and barplots.

## Details

Genes are first separated according to their direction of differential
expression and filtered using the supplied log2 fold-change and p-value
thresholds. Over-representation analysis is then performed using
[`enrichGO`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

## Examples

``` r
if (FALSE) { # \dontrun{
    enrich_ora_results <- ora_enrichgo(filtered_results,
                                 direction = c("up", "down"),
                                     log2fc_thresh = 1,
                                      pval_thresh = 0.05,
                                      top_terms = 15,
                                      OrgDb = org.Hs.eg.db,
                                     keyType = "ENSEMBL",
                                       ont = "BP")
} # }
```
