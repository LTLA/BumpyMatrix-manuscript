---
title: "BumpyMatrix: a Bioconductor package to represent a 2-dimensional array of non-scalar objects"
tags:
  - R
  - Bioconductor
  - bioinformatics
authors:
  - name: Aaron T. L. Lun
    orcid: 0000-0002-3564-4813
    affiliation: 1
affiliations:
  - name: Genentech Inc., South San Francisco, USA
    index: 1
date: 17 January 2022
bibliography: ref.bib
---

# Summary

# Statement of need

# Basic usage

Let's mock up some data for a multiplexed FISH experiment.
This consists of a `DataFrame` (Bioconductor's wrapper around the base `data.frame` class) where each row contains the coordinates of a single detected transcript molecule in a particular cell.

```r
library(S4Vectors)
df <- DataFrame(
    x=rnorm(10000), y=rnorm(10000), 
    gene=paste0("GENE_", sample(100, 10000, replace=TRUE)),
    cell=paste0("CELL_", sample(20, 10000, replace=TRUE))
)
print(df)
## DataFrame with 10000 rows and 4 columns
##                x           y        gene        cell
##        <numeric>   <numeric> <character> <character>
## 1     -0.5120619   0.8008861     GENE_36      CELL_8
## 2      0.6668935  -0.0625209     GENE_22      CELL_5
## 3     -0.0075178  -0.4204114     GENE_49     CELL_10
## 4      2.3256096  -0.2246588    GENE_100      CELL_1
## 5     -0.2200939  -0.1674417     GENE_43      CELL_4
## ...          ...         ...         ...         ...
## 9996    2.235438  0.29477272     GENE_22     CELL_18
## 9997    0.801421  0.90281187     GENE_21     CELL_10
## 9998   -1.967398  0.00963577     GENE_20     CELL_16
## 9999   -0.471152  1.79703060     GENE_32     CELL_11
## 10000  -0.517331 -0.73837155     GENE_43     CELL_16
```

We use the `splitAsBumpyMatrix()` function to create a `BumpyDataFrameMatrix` where the rows are genes and the columns are cells.
Each entry of this `BumpyDataFrameMatrix` is itself a `DataFrame` that corresponds to a particular gene/cell combination,
and containing the x/y-coordinates for all transcript molecules of that gene in that cell.
Each `DataFrame` may have zero, one or multiple rows, depending on the cell's expression of the gene.

```r
library(BumpyMatrix)
mat <- splitAsBumpyMatrix(df[,c("x", "y")], row=df$gene, column=df$cell)
mat
## 100 x 20 BumpyDataFrameMatrix
## rownames: GENE_1 GENE_10 ... GENE_98 GENE_99
## colnames: CELL_1 CELL_10 ... CELL_8 CELL_9
## preview [1,1]:
##   DataFrame with 2 rows and 2 columns
##             x         y
##     <numeric> <numeric>
##   1 -0.852709  0.397673
##   2  0.562927  0.589408
```

We can now treat `mat` in the same manner as any other matrix-like R object.
For example, we can slice by row and column, equivalent to filtering our original `df` by gene or cell respectively.

```r
mat[1:5,]
## 5 x 20 BumpyDataFrameMatrix
## rownames: GENE_1 GENE_10 GENE_100 GENE_11 GENE_12
## colnames: CELL_1 CELL_10 ... CELL_8 CELL_9
## preview [1,1]:
##   DataFrame with 2 rows and 2 columns
##             x         y
##     <numeric> <numeric>
##   1 -0.852709  0.397673
##   2  0.562927  0.589408

mat[,c("CELL_5", "CELL_10", "CELL_20")]
## 100 x 3 BumpyDataFrameMatrix
## rownames: GENE_1 GENE_10 ... GENE_98 GENE_99
## colnames: CELL_5 CELL_10 CELL_20
## preview [1,1]:
##   DataFrame with 4 rows and 2 columns
##             x          y
##     <numeric>  <numeric>
##   1  1.920933 -0.7523476
##   2  1.948563  1.5275105
##   3  0.489665  0.9025059
##   4 -0.197640  0.0319262
```

The `lengths()` function will return the number of entries in each `DataFrame`.
In this context, the return value is a count matrix that quantifies the expression of each gene in each cell.

```r
counts <- lengths(mat)
counts[1:10,1:5]
##          CELL_1 CELL_10 CELL_11 CELL_12 CELL_13
## GENE_1        2       6       6       3       2
## GENE_10       7       5       6       3       3
## GENE_100      3       6       7       6      10
## GENE_11       9       8       5       6       3
## GENE_12       6       8      10       5       4
## GENE_13       4       2       3       3       9
## GENE_14       5       4       6       9       3
## GENE_15       4       4       4       5       3
## GENE_16       4       1       4       7       1
## GENE_17       8       8       4       3      11
```

## Compatibility with the `SummarizedExperiment`

Of particular interest is the ability to store `mat` inside Bioconductor's `SummarizedExperiment` data structure.
The `SummarizedExperiment` provides a holistic representation of a dataset that contains all experimental assays, feature- and sample-level annotations as well as additional metadata. 
This representation is convenient as it ensures that any operations are synchronised across all aspects of the dataset, avoiding any account-keeping errors across the lifetime of an analysis workflow.
By storing the `BumpyDataFrameMatrix` inside the `SummarizedExperiment`, we can synchronise the manipulations of the transcript coordinates to the rest of the dataset.

```r
library(SummarizedExperiment)
se <- SummarizedExperiment(list(positions = mat, counts = lengths(mat)))
se
## class: SummarizedExperiment
## dim: 100 20
## metadata(0):
## assays(2): positions counts
## rownames(100): GENE_1 GENE_10 ... GENE_98 GENE_99
## rowData names(0):
## colnames(20): CELL_1 CELL_10 ... CELL_8 CELL_9
## colData names(0):

# Limiting our analysis to the first 10 cells:
sub.se <- se[,1:10]
dim(assay(sub.se, "counts"))
## [1] 100  10

# Our BumpyDataFrameMatrix is also subject to the same operation:
assay(sub.se, "positions")
## 100 x 10 BumpyDataFrameMatrix
## rownames: GENE_1 GENE_10 ... GENE_98 GENE_99
## colnames: CELL_1 CELL_10 ... CELL_17 CELL_18
## preview [1,1]:
##   DataFrame with 2 rows and 2 columns
##             x         y
##     <numeric> <numeric>
##   1 -0.852709  0.397673
##   2  0.562927  0.589408
```

This compatibility also simplifies the consumption of `BumpyDataFrameMatrix` assays by frameworks that operate on a `SummarizedExperiment`.
For example, the **alabaster** framework provides a mechanism for serializing R/Bioconductor objects into language-agnostic file formats, for use with other programming ecosystems like Python or Javascript.
As **alabaster** already provides a mechanism to serialize `SummarizedExperiment` objects to file, it can be easily extended to handle objects containing `BumpyDataFrameMatrix` objects in its assays.

```r
library(alabaster)

# Saving the SummarizedExperiment into a staging directory.
staging <- tempfile()
dir.create(staging)
info <- stageObject(se, staging, "experiment_1")

# Loading the SummarizedExperiment recovers our BumpyDataFrameMatrix.
reloaded <- loadObject(info, staging)
assay(reloaded, "positions")
## 100 x 20 BumpyDataFrameMatrix
## rownames: GENE_1 GENE_10 ... GENE_98 GENE_99
## colnames: CELL_1 CELL_10 ... CELL_8 CELL_9
## preview [1,1]:
##   DataFrame with 2 rows and 2 columns
##             x         y
##     <numeric> <numeric>
##   1 -0.852709  0.397673
##   2  0.562927  0.589408
```

# References

