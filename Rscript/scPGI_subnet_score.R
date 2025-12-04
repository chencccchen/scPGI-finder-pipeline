# Compute hub subnet score for epithelial cells (example: MLF2 as hub)
library(Seurat)
library(future.apply)   # If you don't need parallelization, replace future_lapply with lapply

setwd("path/data")

# Seurat object (epithelial cells)
scRNA   <- readRDS("scRNA_Epithelial_1000.rds")

# scPGI edge table; must contain columns 'genei' and 'genej'
giresef <- readRDS("Epithelial_scPGI.rds")

# 1. Define hub gene
hub_gene <- "MLF2"

# 2. Subset edges that involve the hub gene
giresef_hub <- giresef[giresef$genei == hub_gene | giresef$genej == hub_gene, ]
featuregene <- unique(c(giresef_hub$genei, giresef_hub$genej))

# 3. Extract raw counts for the hub and its partners
celldata <- GetAssayData(scRNA[featuregene, ], slot = "counts", assay = "RNA")
celldata <- as.matrix(celldata)

# 4. Discretize expression per gene:
#    0 = no expression
#    1 = low/medium expression
#    2 = high expression (> median of non-zero values)
expls <- function(mm) {
  x <- celldata[mm, ]
  temp <- rep(1L, length(x))  # default: low/medium
  
  nz <- which(x != 0)
  if (length(nz) > 0) {
    thr <- median(x[nz])
    idx_high <- which(x > thr)
    temp[idx_high] <- 2L
  }
  
  idx_zero <- which(x == 0)
  temp[idx_zero] <- 0L
  
  temp
}

cellls <- future_lapply(
  1:nrow(celldata),
  FUN = expls,
  future.seed = TRUE
)
cellls <- do.call(rbind, cellls)
mode(cellls) <- "integer"
rownames(cellls) <- rownames(celldata)
colnames(cellls) <- colnames(celldata)

# 5. For each hub–partner edge, sum the discretized states in each cell
#    2 + 2 = 4 means both genes are highly expressed in that cell
i1 <- match(giresef_hub$genei, rownames(cellls))
i2 <- match(giresef_hub$genej, rownames(cellls))

pair_sum <- cellls[i1, ] + cellls[i2, ]

# 6. For each cell, compute the hub subnet score as the fraction of edges with pair_sum == 4
#    i.e., the proportion of hub–partner pairs that are simultaneously highly expressed
subnetscore <- colMeans(pair_sum == 4L)

# 7. Store the result in Seurat meta.data
scRNA$subnetscore_MLF2 <- subnetscore

# (optional) Save the updated Seurat object
setwd("path/scPGI_finder")
saveRDS(scRNA, "scRNA_Epithelial_1000_with_MLF2_subnetscore.rds")
