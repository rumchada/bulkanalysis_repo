Welcome to the bulkanalysis package.

The purpose of the bulkanalysis package is to streamline the downstream
analysis of bulk RNA-seq data by providing a collection of commonly used
analysis and visualizations in a single package. The package was
developed as a wrapper around statistical and visualization methods that
I have frequently encountered during my experience analyzing bulk
RNA-seq data. Rather than introducing new statistical methods,
bulkanalysis is designed to reduce repetitive coding and provide a
consistent workflow for obtaining commonly used downstream results. This
package performs techniques such as differential expression via edgeR
with traditional outputs such as heatmaps, volcano plots. As well as
over-representation analysis via clusterProfiler using Fishcher’s Exact
Test with summarized output visuals (As of 9/17/26).

**Installation Instructions**:

``` bash
devtools::install_github("https://github.com/rumchada/bulkanalysis_repo")
```

**Currently Includes**:

[`baseline_dimreduction()`](https://rumchada.github.io/bulkanalysis_repo/reference/baseline_dimreduction.md):
ariance-stabilized expression values are used to calculate sample
distances and principal components. The most variable genes are used for
PCA, and euclidean distance visualization.

[`edgeR_diffexp()`](https://rumchada.github.io/bulkanalysis_repo/reference/edgeR_diffexp.md):
Performs EdgeR based pairwise (group A vs group B) differential
expression analysis for each experimental state in the entered colData
of the DeseqDataset input.

[`volcano_plot()`](https://rumchada.github.io/bulkanalysis_repo/reference/volcano_plot.md):
Contrsucts a volcano plot of the initial DEGs

[`heatmap_function()`](https://rumchada.github.io/bulkanalysis_repo/reference/heatmap_function.md):
Constructs a column clustered, z-scaled heatmap of DEGs

[`gene_id_converter_ver2()`](https://rumchada.github.io/bulkanalysis_repo/reference/gene_id_converter_ver2.md):
Converts genes via getBM() Ensembl dataset API pulls. If fails to
connect to API, then will use the downloaded AnnotationDbi of chosen
organism to convert available geneids.

[`enrichgo_unpack_ver3()`](https://rumchada.github.io/bulkanalysis_repo/reference/enrichgo_unpack_ver3.md):
Unpacks the key results of clusterProfiler::enrichgo() function. Outputs
are set into the set list marked by over-represented terms to view
individual DEGs for eahc term.

[`ora_enrichgo()`](https://rumchada.github.io/bulkanalysis_repo/reference/ora_enrichgo.md):
Performs over-representation test on upregulated or downregulated DEGs
and plots the top 15 over-represented terms into a dotplot.

[`pc_var_association()`](https://rumchada.github.io/bulkanalysis_repo/reference/pc_var_association.md):
Custom implementation inspired by Principal Variance Component Analysis
(PVCA) described by Li, J., Bushel, P., et. al.Chapter 12, 2010.
Quantifies the association between sample-level metadata covariates and
factors and principal components (PCs) derived from a bulk RNA-seq
expression dataset using ANOVA or a linear regression.

`run_de_pipleine()`: Runs
[`volcano_plot()`](https://rumchada.github.io/bulkanalysis_repo/reference/volcano_plot.md),
[`heatmap_function()`](https://rumchada.github.io/bulkanalysis_repo/reference/heatmap_function.md),
[`ora_enrichgo()`](https://rumchada.github.io/bulkanalysis_repo/reference/ora_enrichgo.md),
and
[`enrichgo_unpack_ver3()`](https://rumchada.github.io/bulkanalysis_repo/reference/enrichgo_unpack_ver3.md)
in succession and reports the up/down regulated DEGs for each comparison
in your initial list of results.

**Future Goals:**

Implement GSEA

Protein-Protein Interaction Network

New ORA visuals:

``` R
-Radial Upset Plot to determine which enriched term contains most DEGs

-IGraph based ORA network to determine central enriched terms 
```

**Citations:**

1.  Mark D. Robinson, Davis J. McCarthy, Gordon K. Smyth, edgeR: a
    Bioconductor package for differential expression analysis of digital
    gene expression data, Bioinformatics, Volume 26, Issue 1, January
    2010, Pages 139–140, <https://doi.org/10.1093/bioinformatics/btp616>

2.  Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan
    L, Fu x, Liu S, Bo X, Yu G (2021). “clusterProfiler 4.0: A universal
    enrichment tool for interpreting omics data.” The Innovation, 2(3),
    100141.doi:10.1016/j.xinn.2021.100141.

3.  Love MI, Huber W, Anders S (2014).“Moderated estimation of fold
    change and dispersion for RNA-seq data with DESeq2.”Genome Biology,
    15, 550.doi:10.1186/s13059-014-0550-8.

4.  Li, J., Bushel, P., Chu, T., Wolfinger, R., Batch Effects and Noise
    in Microarray Experiment, Chapter 12, 2010
