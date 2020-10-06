# 5k_pbmc_v3

## About the dataset

The **5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor (v3 chemistry)** dataset is a single cell gene expression dataset retreived from the [10X Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_v3). 

## About the R objects

The downloaded scRNA-seq data file (`5k_pbmc_v3_filtered_feature_bc_matrix.tar.gz`) is uncompressed and processed in R. The R objects are write to files:

- `pbmc_2k_v3_Seurat.rds` - A `Seurat` object of the pbmc dataset
- `pbmc_2k_v3_Seurat_Idents.rds` - A factor object containing the cluster assignments (saved from Seurat workflow)
- `pbmc_2k_v3_df.rds` - A `data.frame` object containing the normalised expression (saved from Seurat workflow)
- `pbmc_2k_v3_sce.rds` - A `SingleCellExperiment` object of the pbmc dataset (saved from scater/scran workflow)

## Creating the R objects

The R objects are created as follow:

### Seurat

```R
library(Seurat)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc5k-downsampled", min.cells = 3, min.features = 200)
pbmc
```

> An object of class Seurat<br />
> 18791 features across 4962 samples within 1 assay <br />
> Active assay: RNA (18791 features, 0 variable features)

```R
# Downsampling to 2000 cells to reduce file size
set.seed(12345)
pbmc = subset(pbmc, cells = sample(Cells(pbmc), 2000))
pbmc
```

> An object of class Seurat<br />
> 18791 features across 2000 samples within 1 assay<br />
> Active assay: RNA (18791 features, 0 variable features)

```R
# QC
pbmc[["percent.mito"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Subset data
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 25)
pbmc
```

> An object of class Seurat<br />
> 18791 features across 1865 samples within 1 assay<br />
> Active assay: RNA (18791 features, 0 variable features)

```R
# Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identifies features that are outliers
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scales and centers features
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Regress out percent.mito
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mito")

# Run a PCA dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)

# Constructs a Shared Nearest Neighbor (SNN) Graph
set.seed(12345)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)

# Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm
pbmc <- FindClusters(pbmc, resolution = 0.5, random.seed = 12345)

# Save pmbc obj
saveRDS(pbmc, "pbmc_2k_v3_Seurat.rds")

# Get object's identity classes
identity <- Idents(pbmc)

# Save identity obj
saveRDS(identity, "data/pbmc_2k_v3_Seurat_Idents.rds")
```

### data.frame (continue from Seurat)

```R
# Save expression as data.frame
pbmcNorm = FetchData(pbmc, all.genes, slot = "data")

# Transpose so that rows are genes and columns are cells
pbmcNorm = t(pbmcNorm)

# Save pbmcNorm obj
saveRDS(pbmcNorm, "pbmc_2k_v3_df.rds")
```

### scater/scran

```R
library(scater)
library(scran)

# Load the PBMC dataset
sce <- DropletUtils::read10xCounts("filtered_feature_bc_matrix", sample.names = "pbmc5k-downsampled")

# Add barcode ID to column names and gene name to row names
colnames(sce) <- colData(sce)$Barcode
rownames(sce) <- rowData(sce)$Symbol
```

> \> dim(sce)<br />
> [1] 33538  5025

```R
# Compute and add QC metrics
is.mito <- grep("^MT-", rowData(sce)$Symbol)
sce <- addPerCellQC(sce, list(MT = is.mito))
sce <- addPerFeatureQC(sce)
rowData(sce)$n_cells <- rowData(sce)$detected/100 * ncol(sce)

# Obtain identical genes & cells as the starting Seurat object by using same fileter parameters: 
# min.cells = 3, min.features = 200
sce <- sce[rowData(sce)$n_cells >= 3, colData(sce)$detected >= 200]
```

> \> dim(sce)<br />
> [1] 18791  4962

```R
# Downsampling to 2000 cells to reduce file size
set.seed(12345)
sce <- sce[,sample(colnames(sce), 2000)]
```

> \> dim(sce)<br />
> [1] 18791  2000

```R
# Subset data using same parameters as the Seurat object
sce <- subset(sce, , detected > 200 & detected < 5000 & subsets_MT_percent < 25)
```

> \> dim(sce)<br />
> [1] 18791  1865

```R
# Normalization by deconvolution
set.seed(12345)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, min.mean = 0.1, cluster = clusters)
sce <- logNormCounts(sce, size_factors = sizeFactors(sce))

# Variance modelling
dec <- modelGeneVar(sce)

# Select the top 10% of genes with the highest biological components
hvg <- getTopHVGs(stats = dec, prop = 0.1)

# Performing PCA only on the chosen HVGs
set.seed(12345)
sce <- runPCA(sce, subset_row = hvg)

# Choosing the number of PCs based on population structure
choices <- getClusteredPCs(reducedDim(sce))
```

> \> metadata(choices)$chosen<br />
> [1] 16

```R
# Subset PC matrix
reducedDim(sce, "PCA") <- reducedDim(sce, "PCA")[,1:metadata(choices)$chosen]

# Graph-based clustering
g <- buildSNNGraph(sce, k = 10, use.dimred = "PCA")
clust <- igraph::cluster_walktrap(g)$membership

# Add cluster assignments to sce
sce$label = factor(clust)

# Save sce obj
saveRDS(sce, "data/pbmc_2k_v3_sce.rds")
```
