# StackedVlnPlot

## Demo data

The scRNA-seq demo data (`*rds`) files available in the `data` folder of this repository.

## With `Seurat`

Stacked violin plot functionality is added to `Seurat` in version 3.2.1. 

```R
library(Seurat)
library(ggplot2)

# Load Seurat obj
pbmc <- readRDS("data/pbmc_2k_v3_Seurat.rds")

features <- c("CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3", "S100A8", "GNLY",
              "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3")

a <- VlnPlot(pbmc, features, stack = TRUE, sort = TRUE) +
        theme(legend.position = "none") + ggtitle("identities on y-axis")
        
b <- VlnPlot(pbmc, features, stack = TRUE, sort = TRUE, flip = TRUE) +
        theme(legend.position = "none") + ggtitle("identities on x-axis")

# Use patchwork to join plots
a + b
```

<img src="https://github.com/ycl6/StackedVlnPlot/raw/master/images/StackedVlnPlot_Seurat.png" width="700px" alt="Using Seurat's VlnPlot function">
