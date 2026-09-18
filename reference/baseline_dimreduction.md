# Baseline Dimension Reduction of RNA-seq Data

Takes a formatted `DESeqDataSet` object and performs several exploratory
analyses to evaluate sample similarity and major sources of variation
within the dataset. Performs baseline exploratory dimension reduction
and sample-level visualization of a `DESeqDataSet`. Variance-stabilized
expression values are used to calculate sample distances and principal
components. The most variable genes are used for PCA and downstream
visualization.

The analysis:

1.  applies a variance-stabilizing transformation to the expression
    matrix,

2.  identifies the top variable genes,

3.  calculates Euclidean distances between samples,

4.  performs principal component analysis (PCA) using the top variable
    genes,

5.  generates sample embedding plots,

6.  generates an elbow plot showing the variance explained by the
    leading principal components, and

7.  performs UMAP for nonlinear dimensionality reduction.

## Usage

``` r
baseline_dimreduction(bulk_dataset, top_var = 2000, grouping = NA)
```

## Arguments

- bulk_dataset:

  A formatted `DESeqDataSet` (`DDS`) object containing the RNA-seq
  expression data and associated sample metadata.

- top_var:

  Integer specifying the number of genes with the highest variance to
  retain for PCA and downstream dimension reduction. For example,
  `top_var = 2000` uses the 2,000 most variable genes.

- grouping:

  Character string specifying the sample metadata variable used to
  annotate or group samples in the distance and dimension reduction
  visualizations. The variable should correspond to a column in the
  `colData` of `bulk_dataset`.

## Value

- `distance_plot`: euclidean distance plot,

- `elbow_plot`: shows which PC has the majority of the variance within
  the data

- `PCA plot`: scatterplot of the left unitary vectors

- `UMAP_plot`: A UMAP visualization of the sample embeddings.

  results\$distance_plot results\$PCA_plot results\$elbow_plot
  results\$umap_plot

## Examples

``` r
if (FALSE)  baseline_dimreduction(dds_object, top_var = 2000, grouping = NA) # \dontrun{}

```
