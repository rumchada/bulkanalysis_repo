# Gene ID Converter

Enter the type of conversion from_type and to_type (e.g. `"ensembl"`,
`"entrez"`, `"symbol"`).

## Usage

``` r
gene_id_converter_ver2(vector, from_type, to_type, org_db, ensembl_dataset)
```

## Arguments

- vector:

  Character vector of genes.

- from_type:

  Character string specifying the initial annotation type of the genes.
  Supported values include `"ensembl"`, `"entrez"`, and `"symbol"`.

- to_type:

  Character string specifying the desired annotation type of the genes.
  Supported values include `"ensembl"`, `"entrez"`, and `"symbol"`.

- org_db:

  An organism-specific annotation database used to facilitate gene
  identifier conversion, such as `org.Hs.eg.db` or `org.Mm.eg.db`.

- ensembl_dataset:

  A character string specifying the Ensembl dataset used for annotation.
  For example, `"hsapiens_gene_ensembl"` for human genes.

## Details

Uses `Bioconductor::getBM()` and
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) to communicate
with the Ensembl database via `getBM()` to convert a character vector of
genes into any of the three supported identifier forms.

## Examples

``` r
if (FALSE) { # \dontrun{
gene_id_conversion <- gene_id_converter_ver2(
  vector = c("ENSG00000141510", "ENSG00000171862"),
  from_type = "ensembl",
  to_type = "symbol",
  org_db = org.Hs.eg.db,
  ensembl_dataset = "hsapiens_gene_ensembl"
)
} # }
```
