# Principal Component Variation Association

Quantifies the association between sample-level metadata covariates and
factors and principal components (PCs) derived from a bulk RNA-seq
expression dataset. The function evaluates both categorical variables
(covariates) and continuous variables across the top principal
components and returns a matrix of statistical significance values.

This function is a custom implementation inspired by Principal Variance
Component Analysis (PVCA) described by Li, J., Bushel, P., et.
al.Chapter 12, 2010. It can be used both before and after batch
correction to evaluate whether experimental, biological, or technical
variables explain major sources of variation in the dataset.

## Usage

``` r
pc_var_association(bulk_ds)
```

## Arguments

- bulk_ds:

  object containing normalized or appropriately transformed expression
  data and corresponding sample metadata.

## Value

A matrix containing p-values representing the association between each
metadata variable and each principal component. The matrix has principal
components as rows and metadata variables as columns.

A heatmap visualizing the resulting p-value matrix is also generated.

## Details

INPUT:

A `DESeqDataSet` object containing normalized or appropriately
transformed expression data and corresponding sample metadata.

Metadata variables are evaluated according to their data type:

- Categorical variables are evaluated using ANOVA and an F-test.

- Continuous variables are evaluated using an appropriate linear model
  or association test.

## References

Li, J., Bushel, P., Chu, T., Wolfinger, R., Batch Effects and Noise in
Microarray Experiment, Chapter 12, 2010

## See also

[`DESeqDataSet`](https://rdrr.io/pkg/DESeq2/man/DESeqDataSet.html),
[`vst`](https://rdrr.io/pkg/DESeq2/man/vst.html),
[`rlog`](https://rdrr.io/pkg/DESeq2/man/rlog.html)
