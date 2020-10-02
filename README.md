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

## With `ggplot2`

Given a `data.frame` and a vector of identity classes (cluster ID), a stacked violin plot can be created with the `ggplot2` package.

```R
library(ggplot2)
library(cowplot)
library(patchwork)

# Load data.frame obj
pbmc <- readRDS("data/pbmc_2k_v3_df.rds")
identity <- readRDS("data/pbmc_2k_v3_Seurat_Idents.rds")

features <- c("CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3", "S100A8", "GNLY",
              "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3")

# Subset data.frame
pbmc <- pbmc[features,]

# Use melt to change data.frame format
pbmc <- reshape2::melt(pbmc, varnames = c("Feat", "Cell"), value.name = "Expr")

# Merge with identity classes
pbmc <- merge(pbmc, data.frame(Cell = names(identity), Idents = identity), by = "Cell")
```

The converted *long format* `pbmc`:

```
                Cell   Feat     Expr Idents
1 AAACGCTGTGTCCGGT-1  CD79A 0.000000      1
2 AAACGCTGTGTCCGGT-1  MS4A1 0.000000      1
3 AAACGCTGTGTCCGGT-1   CD8A 0.000000      1
4 AAACGCTGTGTCCGGT-1   CD8B 0.000000      1
5 AAACGCTGTGTCCGGT-1    LYZ 5.332104      1
6 AAACGCTGTGTCCGGT-1 LGALS3 2.306096      1
```

There are different ways to show the stacked violin plot:

| Plot | X-axis | Y-axis | Facet |
| --- | --- | --- | --- |
| a | Identities | Expression | Feature |
| b | Expression | Identities | Feature |
| c | Feature | Expression | Identities |
| d | Expression | Feature | Identities |

```R
# Identities on x-axis
a <- ggplot(pbmc, aes(factor(Idents), Expr, fill = Feat)) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
                           c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
        facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
        theme_cowplot(font_size = 12) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0)) +
        ggtitle("identities on x-axis") + xlab("Identity") + ylab("Expression Level")

# Identities on y-axis
b <- ggplot(pbmc, aes(Expr, factor(Idents), fill = Feat)) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        scale_x_continuous(expand = c(0, 0), labels = function(x)
                           c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
        facet_grid(cols = vars(Feat), scales = "free")  +
        theme_cowplot(font_size = 12) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
        ggtitle("identities on y-axis") + xlab("Expression Level") + ylab("Identity")

# Use patchwork to join plots
a + b
```

<img src="https://github.com/ycl6/StackedVlnPlot/raw/master/images/StackedVlnPlot_dataframe1.png" width="900px" alt="Using ggplot2 (Plot a and b)">

```R
# Feature on x-axis
c <- ggplot(pbmc, aes(factor(Feat), Expr, fill = Feat)) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
                           c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
        facet_grid(rows = vars(Idents), scales = "free", switch = "y") +
        theme_cowplot(font_size = 12) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggtitle("Feature on x-axis") + xlab("Feature") + ylab("Expression Level")

# Feature on y-axis
d <- ggplot(pbmc, aes(Expr, factor(Feat), fill = Feat)) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        scale_x_continuous(expand = c(0, 0), labels = function(x)
                           c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
        facet_grid(cols = vars(Idents), scales = "free") +
        theme_cowplot(font_size = 12) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold")) +
        ggtitle("Feature on y-axis") + xlab("Expression Level") + ylab("Feature")

# Use patchwork to join plots
c + d
```

<img src="https://github.com/ycl6/StackedVlnPlot/raw/master/images/StackedVlnPlot_dataframe2.png" width="900px" alt="Using ggplot2  (Plot c and d)">
