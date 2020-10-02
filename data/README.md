# 5k_pbmc_v3

## About the dataset

The **5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor (v3 chemistry)** dataset is a single cell gene expression dataset retreived from the [10X Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_v3). 

## About the R objects

The downloaded scRNA-seq data file (`5k_pbmc_v3_filtered_feature_bc_matrix.tar.gz`) is uncompressed and processed in R. The R objects are write to files:

- `pbmc_2k_v3_Seurat.rds` - A `Seurat` object of the pbmc dataset
- `pbmc_2k_v3_df.rds` - A `data.frame` objects containing the normalised expression

## Creating the R objects

The R objects are created as follow:

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
> Active assay: RNA (18791 features, 0 variable features)<br />

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

```R
# Save expression as data.frame
pbmcNorm = FetchData(pbmc, all.genes, slot = "data")

# Transpose so that rows are genes and columns are cells
pbmcNorm = t(pbmcNorm)

# Save pbmcNorm obj
saveRDS(pbmcNorm, "pbmc_2k_v3_df.rds")
```
