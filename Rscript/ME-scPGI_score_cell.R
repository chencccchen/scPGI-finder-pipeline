## ==========================================
## Calculate ME-scPGI score for each Epithelial cell
## ==========================================

library(Seurat)
library(future.apply)

## 1. Input files -------------------------------------------------------
## scRNA: Seurat object
## ME_scPGI_pairs: data.frame with columns Gene1, Gene2, pair (and optionally pairtype)
# Example paths (modify to your own)
setwd("path/data/")
scRNA <- readRDS("scRNA_Epithelial_1000.rds")
ME_scPGI <- read.table("ME-scPGI.txt",header = T,sep = "\t")
dim(ME_scPGI)
ME_scPGI$pair<-paste(ME_scPGI$Gene1,ME_scPGI$Gene2,sep="_")

## 2. Keep only genes present in scRNA ---------------------------------

feature_genes <- intersect(
  unique(c(ME_scPGI$Gene1, ME_scPGI$Gene2)),
  rownames(scRNA)
)

ME_scPGI <- ME_scPGI[
  ME_scPGI$Gene1 %in% feature_genes &
    ME_scPGI$Gene2 %in% feature_genes,
]

## 3. Extract counts matrix for these genes ----------------------------

counts_mat <- GetAssayData(
  object = scRNA[feature_genes, ],
  slot   = "counts",
  assay  = "RNA"
)
counts_mat <- as.matrix(counts_mat)
rownames(counts_mat) <- feature_genes  # just to be explicit

## 4. Discretize gene expression (0/1/2) -------------------------------
##    0 = not expressed
##    1 = low/medium (non-zero but <= median of non-zero values)
##    2 = high (non-zero and > median of non-zero values)

expls <- function(idx) {
  x <- counts_mat[idx, ]
  tmp <- rep(1L, length(x))
  
  nz <- which(x != 0)
  if (length(nz) > 0) {
    thr <- median(x[nz])
    tmp[x > thr] <- 2L
  }
  
  tmp[x == 0] <- 0L
  tmp
}

plan(multisession, workers = 4)  # adjust workers as needed

gene_state_list <- future_lapply(
  X   = seq_len(nrow(counts_mat)),
  FUN = expls,
  future.seed = TRUE
)

gene_states <- do.call(rbind, gene_state_list)
storage.mode(gene_states) <- "integer"
rownames(gene_states) <- rownames(counts_mat)
colnames(gene_states) <- colnames(counts_mat)

## 5. Build pair-level state matrix for ME-scPGIs ----------------------

i1 <- match(ME_scPGI$Gene1, rownames(gene_states))
i2 <- match(ME_scPGI$Gene2, rownames(gene_states))

pair_states <- gene_states[i1, ] + gene_states[i2, ]
mode(pair_states) <- "integer"
rownames(pair_states) <- ME_scPGI$pair
# rows: ME-scPGI pairs
# cols: cells
# values: 0,1,2,3,4

## 6. Define ME-scPGI score per cell -----------------------------------
score_fun <- function(x) {
  sum(x == 4L) / length(x)
}
ME_scPGI_score <- apply(pair_states, 2, score_fun)

## 7. Attach score to Seurat object & save -----------------------------

scRNA$ME_scPGI_score <- ME_scPGI_score
saveRDS(
  scRNA,
  file = "scRNA_MescPGI.rds"
)



