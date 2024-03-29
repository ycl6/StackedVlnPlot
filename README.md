# StackedVlnPlot

## Demo data

The scRNA-seq demo data (`*rds`) files are available in the `data` folder of this repository.

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

<img src="https://github.com/ycl6/StackedVlnPlot/raw/master/images/StackedVlnPlot_Seurat.png" width="800px" alt="Using Seurat's VlnPlot function">

## With `ggplot2`

Given a `data.frame` and a vector of identity classes (cluster ID), a stacked violin plot can be created with the `ggplot2` package.

### Prepare `data.frame`

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
pbmc <- pbmc[,features]

# Add cell ID and identity classes
pbmc$Cell <- rownames(pbmc)
pbmc$Idents <- identity

# Use melt to change data.frame format
pbmc <- reshape2::melt(pbmc, id.vars = c("Cell","Idents"), measure.vars = features,
                       variable.name = "Feat", value.name = "Expr")
```

The converted *long format* `pbmc`:

```
                 Cell Idents  Feat     Expr
1  AACAACCTCACCTCTG-1      0 CD79A 0.000000
2  AGGAGGTTCGCGGACT-1      0 CD79A 1.743733
3  AGGCATTCAAGACGGT-1      1 CD79A 0.000000
4  GCAACCGCAGTTTCGA-1      1 CD79A 0.000000
5  TTTCACATCGTCCTCA-1      1 CD79A 0.000000
6  CTGCCTAAGCGTTCAT-1      1 CD79A 0.000000
7  CCTCCTCAGCGTCAGA-1      5 CD79A 3.104723
8  AACCATGAGAGCCTGA-1      0 CD79A 0.000000
9  ATGAGTCTCACATTGG-1      2 CD79A 0.000000
10 AGTCATGCACTAACCA-1      4 CD79A 2.756005
```

### Create plots

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

### Sort identity classes and features

> **Note:** Some of the codes below are taken and modified from the `Seurat` package.

Below demonstrates how to recreate the reordering of the identity classes and features seen in Seurat's stacked violin plots. 

```R
# Calculate average expression per Idents, output as wide format
avg <- sapply(X = split(x = pbmc, f = pbmc$Idents),
              FUN = function(df) { return(tapply(X = df$Expr, INDEX = df$Feat, FUN = mean)) })

# L2Norm (Euclidean norm) function
L2Norm <- function(mat, MARGIN){
        normalized <- sweep(x = mat, MARGIN = MARGIN,
                            STATS = apply(X = mat, MARGIN = MARGIN,
                                          FUN = function(x){ sqrt(x = sum(x ^ 2)) }), FUN = "/")
        normalized[!is.finite(x = normalized)] <- 0
        return(normalized)
}

# Performs hierarchical clustering
idents.order <- hclust(d = dist(t(L2Norm(mat = avg, MARGIN = 2))))$order
avg <- avg[,idents.order]
avg <- L2Norm(mat = avg, MARGIN = 1)
mat <- hclust(d = dist(avg))$merge

# Order feature clusters by position of their "rank-1 idents"
position <- apply(X = avg, MARGIN = 1, FUN = which.max)
orderings <- list()
for (i in 1:nrow(mat)) {
        x <- if (mat[i,1] < 0) -mat[i,1] else orderings[[mat[i,1]]]
        y <- if (mat[i,2] < 0) -mat[i,2] else orderings[[mat[i,2]]]
        x.pos <- min(x = position[x])
        y.pos <- min(x = position[y])
        orderings[[i]] <- if (x.pos < y.pos) { c(x, y) } else { c(y, x) }
}
features.order <- orderings[[length(orderings)]]

# Update Feature and Identity factor orders
pbmc$Idents <- factor(pbmc$Idents, levels = levels(pbmc$Idents)[idents.order])
pbmc$Feat <- factor(pbmc$Feat, levels = levels(pbmc$Feat)[features.order])

# Plot stacked violin plot with reordered identity classes and features
e <- ggplot(pbmc, aes(factor(Feat), Expr, fill = Feat)) +
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

e
```

<img src="https://github.com/ycl6/StackedVlnPlot/raw/master/images/StackedVlnPlot_dataframe3.png" width="450px" alt="Hierarchical clustering of identity classes and features">

### Add gene grouping annotation

Below demonstrates how to add gene grouping annotation to sorted stacked violin plots. 

```R
# Create grouping info
df <- data.frame(x = levels(pbmc$Feat), group = c("A","A","B","B","B","B","B","C","C","C","D","D","D"), 
                 stringsAsFactors = FALSE)
df$x <- factor(df$x, levels = levels(pbmc$Feat))
df$group <- factor(df$group)
df
```

| x | group |
| :---: | :---: |
| MS4A1 | A |
| CD79A | A |
| LGALS3 | B |
| LYZ | B |
| CST3 | B |
| S100A8 | B |
| FCER1A | B |
| KLRB1 | C |
| CD8B | C |
| CD8A | C |
| NKG7 | D |
| GNLY | D |
| FCGR3A | D |

```R
color <- c("cyan", "pink", "green", "darkorange")

# Same as plot e, but hide x-axis labels, change plot.margin to reduce spacing between plots
f <- ggplot(pbmc, aes(Feat, Expr, fill = Feat)) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
                           c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
        facet_grid(rows = vars(Idents), scales = "free", switch = "y") +
        theme_cowplot(font_size = 12) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              panel.background = element_rect(fill = NA, color = "black"),
              plot.margin = unit(c(0,0,0,0), "cm"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank()
              axis.text.x = element_blank()) +
        ggtitle("Feature on x-axis") + ylab("Expression Level")

# Use geom_tile() to add grouping colorings and geom_text() to add grouping labels
g <- ggplot(df, aes(x = x, y = 1, fill = group, label = group)) + geom_tile() +
        geom_text(fontface = "bold", size = 3) + theme_bw(base_size = 12) +
        scale_fill_manual(values = color) + scale_y_continuous(expand = c(0, 0)) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              panel.background = element_blank(), 
              panel.border = element_blank(),
              plot.background = element_blank(), 
              plot.margin = unit(c(0,0,0,0), "cm"), 
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) + xlab("Feature")

# Use patchwork to join plots
f + g + plot_layout(ncol = 1, heights = c(0.95, 0.05))
```

<img src="https://github.com/ycl6/StackedVlnPlot/raw/master/images/StackedVlnPlot_dataframe_grouping1.png" width="450px" alt="Hierarchical clustering of identity classes and features">

Legend is used to defind the grouping labels when the labels are too long to fit within the annotation bar.

```R
# Change to long names
levels(df$group) = c("long long name A", "long long name B", "long long name C", "long long name D")

# guides() is used to specify some aesthetic parameters of legend key
h <- ggplot(df, aes(x = x, y = 1, fill = group)) + geom_tile() + theme_bw(base_size = 12) +
        scale_fill_manual(values = color) + scale_y_continuous(expand = c(0, 0)) +
        guides(fill = guide_legend(direction = "vertical", label.position = "right",
                             title.theme = element_blank(), keyheight = 0.5, nrow = 2)) +
        theme(legend.position = "bottom", 
              legend.justification = "left",
              legend.margin = margin(0,0,0,0), 
              legend.box.margin = margin(-10,05,0,0),
              panel.spacing = unit(0, "lines"),
              panel.background = element_blank(),
              panel.border = element_blank(),
              plot.background = element_blank(),
              plot.margin = unit(c(0,0,0,0), "cm"),
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) + xlab("Feature")
              
# Use patchwork to join plots
f + h + plot_layout(ncol = 1, heights = c(0.95, 0.05))
```

<img src="https://github.com/ycl6/StackedVlnPlot/raw/master/images/StackedVlnPlot_dataframe_grouping2.png" width="450px" alt="Hierarchical clustering of identity classes and features">

## Given a `SingleCellExperiment` obj

The expression and cluster information can be extracted from a *processed* `SingleCellExperiment` object to create a stacked violin plot with the `ggplot2` package.

```R
library(scater)
library(cowplot)

# Load sce obj
sce <- readRDS("data/pbmc_2k_v3_sce.rds")
```

The `SingleCellExperiment` object provided in this repository contains both raw and normalised counts. The cluster assignments are stored in the `colData`.

```
> sce
class: SingleCellExperiment 
dim: 18791 1865 
metadata(1): Samples
assays(2): counts logcounts
rownames(18791): AL627309.1 AL627309.3 ... AL354822.1 AC240274.1
rowData names(6): ID Symbol ... detected n_cells
colnames(1865): AACAACCTCACCTCTG-1 AGGAGGTTCGCGGACT-1 ...
  AATGGAACAGTAGGAC-1 CCCAACTTCTCGAGTA-1
colData names(13): Sample Barcode ... total label
reducedDimNames(1): PCA
spikeNames(0):
altExpNames(0):
```

Store the required information from the `sce` object in a `data.frame`, and create stacked violin plot. 

```R
features <- c("CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3", "S100A8", "GNLY",
              "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3")

# Subset dgCMatrix
pbmc <- assay(sce, "logcounts")[features,]

# Transpose and convert to data.frame
pbmc <- as.data.frame(t(as.matrix(pbmc)))

# Add cell ID and identity classes
pbmc$Cell <- rownames(pbmc)
pbmc$Cluster <- sce$label

# Use melt to change data.frame format
pbmc <- reshape2::melt(pbmc, id.vars = c("Cell","Cluster"), measure.vars = features,
                       variable.name = "Feat", value.name = "Expr")

f <- ggplot(pbmc, aes(factor(Feat), Expr, fill = Feat)) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
                           c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
        facet_grid(rows = vars(Cluster), scales = "free", switch = "y") +
        theme_cowplot(font_size = 12) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggtitle("Feature on x-axis") + xlab("Feature") + ylab("Expression Level")
f
```

<img src="https://github.com/ycl6/StackedVlnPlot/raw/master/images/StackedVlnPlot_scater.png" width="450px" alt="Using ggplot2 on SingleCellExperiment object">
