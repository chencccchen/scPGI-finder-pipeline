#!/usr/bin/env Rscript
##
## Epithelial scPGI pipeline:
## 1) Cell-ratio test
## 2) Proliferation dependence
## 3) CS-CORE co-expression
## 4) Functional similarity (GOSemSim)
##

## =========================
## 0. Libraries & options
## =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(CSCORE)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(sctransform)
  library(future)
  library(future.apply)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(clusterProfiler)  # bitr
  library(reshape2)         # melt()
  library(AUCell)
})

plan("multisession", workers = 24)
options(future.globals.maxSize = +Inf)
set.seed(12345)

## =========================
## Path settings (modify as needed)
## =========================
## You can also replace them with absolute paths, e.g. /storage/...
DATA_DIR <- "path/data"
OUT_DIR  <- "path/scPGI_finder"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

scRNA_file <- file.path(
  DATA_DIR,
  "scRNA_Epithelial_1000.rds"
)
## List of cell-type–specific up-regulated genes
edif_file  <- file.path(DATA_DIR, "Epithelial-specific gene.rds")
## Example of how feature_genes can be obtained from a full object:
## feature_genes <- rownames(
##   FindMarkers(
##     scRNA_all,
##     ident.1 = "Epithelial",
##     group.by = "manual_marker.labels",
##     logfc.threshold = 0.25,
##     test.use = "wilcox",
##     min.pct = 0.1,
##     only.pos = TRUE
##   )
## )

## =========================
## Helper function: FS calculation (GOSemSim)
## =========================
compute_FS_epithelial <- function(edif_geneid, out_dir = NULL) {
  ## edif_geneid: data.frame that contains SYMBOL and ENTREZID
  edif_geneid$ENTREZID <- as.character(edif_geneid$ENTREZID)
  
  ## 1. Build GO semantic data for BP / CC / MF
  bp <- godata('org.Hs.eg.db', ont = "BP", computeIC = FALSE)
  cc <- godata('org.Hs.eg.db', ont = "CC", computeIC = FALSE)
  mf <- godata('org.Hs.eg.db', ont = "MF", computeIC = FALSE)
  
  ## 2. Compute gene–gene semantic similarity matrices (ENTREZID level)
  simbp <- mgeneSim(
    edif_geneid$ENTREZID,
    semData = bp,
    measure = "Wang",
    drop = NULL,
    combine = "BMA"
  )
  simcc <- mgeneSim(
    edif_geneid$ENTREZID,
    semData = cc,
    measure = "Wang",
    drop = NULL,
    combine = "BMA"
  )
  simmf <- mgeneSim(
    edif_geneid$ENTREZID,
    semData = mf,
    measure = "Wang",
    drop = NULL,
    combine = "BMA"
  )
  
  ## 3. Convert lower triangle of each matrix into a “gene pair + similarity” table
  ## ---- BP ----
  simbp1 <- matrix(0, ncol(simbp), ncol(simbp))
  simbp1[lower.tri(simbp1)] <- simbp[lower.tri(simbp)]
  rownames(simbp1) <- rownames(simbp)  # ENTREZID
  colnames(simbp1) <- colnames(simbp)
  
  simbp2 <- melt(simbp1)
  simbp2 <- simbp2[simbp2$value != 0, ]
  symbol1 <- edif_geneid[match(simbp2$Var1, edif_geneid$ENTREZID), "SYMBOL"]
  symbol2 <- edif_geneid[match(simbp2$Var2, edif_geneid$ENTREZID), "SYMBOL"]
  simbp2$symbol1 <- symbol1
  simbp2$symbol2 <- symbol2
  simbp2$pair    <- paste(simbp2$symbol2, simbp2$symbol1, sep = "_")
  
  ## ---- CC ----
  simcc1 <- matrix(0, ncol(simcc), ncol(simcc))
  simcc1[lower.tri(simcc1)] <- simcc[lower.tri(simcc)]
  rownames(simcc1) <- rownames(simcc)
  colnames(simcc1) <- colnames(simcc)
  
  simcc2 <- melt(simcc1)
  simcc2 <- simcc2[simcc2$value != 0, ]
  symbol1 <- edif_geneid[match(simcc2$Var1, edif_geneid$ENTREZID), "SYMBOL"]
  symbol2 <- edif_geneid[match(simcc2$Var2, edif_geneid$ENTREZID), "SYMBOL"]
  simcc2$symbol1 <- symbol1
  simcc2$symbol2 <- symbol2
  simcc2$pair    <- paste(simcc2$symbol2, simcc2$symbol1, sep = "_")
  
  ## ---- MF ----
  simmf1 <- matrix(0, ncol(simmf), ncol(simmf))
  simmf1[lower.tri(simmf1)] <- simmf[lower.tri(simmf)]
  rownames(simmf1) <- rownames(simmf)
  colnames(simmf1) <- colnames(simmf)
  
  simmf2 <- melt(simmf1)
  simmf2 <- simmf2[simmf2$value != 0, ]
  symbol1 <- edif_geneid[match(simmf2$Var1, edif_geneid$ENTREZID), "SYMBOL"]
  symbol2 <- edif_geneid[match(simmf2$Var2, edif_geneid$ENTREZID), "SYMBOL"]
  simmf2$symbol1 <- symbol1
  simmf2$symbol2 <- symbol2
  simmf2$pair    <- paste(simmf2$symbol2, simmf2$symbol1, sep = "_")
  
  ## 4. Merge BP/CC/MF similarity and reorder columns
  simres <- full_join(full_join(simbp2, simcc2, by = "pair"), simmf2, by = "pair")
  simres <- simres[, c(6, 2, 1, 5, 4, 3, 9, 14)]  ## keep same order as original script
  colnames(simres) <- c(
    "pair",
    "ENTREZID1", "ENTREZID2",
    "symbol1", "symbol2",
    "simbp", "simcc", "simmf"
  )
  
  ## 5. Geometric mean of BP/CC/MF similarities
  simres$fsim <- (simres$simbp * simres$simcc * simres$simmf)^(1 / 3)
  
  ## 6. Optionally save intermediate results
  if (!is.null(out_dir)) {
    saveRDS(simres, file.path(out_dir, "Epithelial_simres_merge.rds"))
  }
  
  simres
}

## =========================
## 1. Load data & compute cell proliferation scores
## =========================
scRNA    <- readRDS(scRNA_file)          ## Seurat object for epithelial cells
feature_genes <- readRDS(edif_file)      ## List of up-regulated genes (feature genes)

expr_counts <- GetAssayData(
  scRNA[feature_genes, ],
  slot  = "counts",
  assay = "RNA"
)
expr_counts <- as.matrix(expr_counts)
rownames(expr_counts) <- feature_genes

## Compute proliferation scores using AUCell
genelist <- readRDS(
  file.path(DATA_DIR, "cell_proliferation_genelist.rds")
)
cells_rankings <- AUCell_buildRankings(scRNA@assays$RNA$counts)
cells_AUC <- AUCell_calcAUC(genelist, cells_rankings)
genescore <- data.frame(t(cells_AUC@assays@data$AUC))
rownames(genescore) == colnames(scRNA)
scRNA@meta.data <- cbind(scRNA@meta.data, genescore)
drcc <- scRNA$drcc   ## proliferation score per cell, derived from AUCell
## The cell proliferation score was evaluated using the enrichment
## of cell-cycle and DNA-replication signatures with AUCell (v1.22.0).

## =========================
## 3. Build all gene pairs & discretize expression
## =========================
## 3.1 Build all unique gene pairs (upper triangle)
build_gene_pairs <- function(genes) {
  n <- length(genes)
  genei <- rep(genes[1:(n - 1)], (n - 1):1)
  genej_list <- lapply(2:n, function(j) genes[j:n])
  genej <- do.call(c, genej_list)
  data.frame(
    genei = genei,
    genej = genej,
    stringsAsFactors = FALSE
  )
}

pair_df <- build_gene_pairs(feature_genes)
genei   <- pair_df$genei
genej   <- pair_df$genej
i1      <- match(genei, rownames(expr_counts))
i2      <- match(genej, rownames(expr_counts))

## 3.2 Discretize expression:
##     0 = no expression
##     1 = low / intermediate expression
##     2 = high expression (> median of non-zero values)
expls <- function(mm) {
  x <- expr_counts[mm, ]
  tmp <- rep(1L, length(x))
  nz  <- which(x != 0)
  if (length(nz) > 0) {
    thr <- median(x[nz])
    idx_high <- which(x > thr)
    tmp[idx_high] <- 2L
  }
  idx_zero <- which(x == 0)
  tmp[idx_zero] <- 0L
  tmp
}

cellls_list <- future_lapply(
  1:nrow(expr_counts),
  FUN = expls,
  future.seed = TRUE
)
cellls <- do.call(rbind, cellls_list)
storage.mode(cellls) <- "integer"
rownames(cellls) <- rownames(expr_counts)
colnames(cellls) <- colnames(expr_counts)

n_cells <- ncol(cellls)

## =========================
## 4. Step 1 + 2:
##    Cell-ratio test + Proliferation test
##    (computed in one loop)
## =========================

## 4.1 Per-gene counts of 0/1/2 states, used for random co-activation probability
pct0 <- rowSums(cellls == 0)
pct1 <- rowSums(cellls == 1)
pct2 <- rowSums(cellls == 2)

g1 <- pct0[i1]
g2 <- pct0[i2]
l1 <- pct1[i1]
l2 <- pct1[i2]
h1 <- pct2[i1]
h2 <- pct2[i2]

prandom_co <- h1 * h2 / (n_cells * n_cells)
prandom_co[prandom_co == 0] <- 1e-7
prandom_co[prandom_co == 1] <- 1 - 1e-7

ptj <- data.frame(
  g1 = g1, g2 = g2,
  l1 = l1, l2 = l2,
  h1 = h1, h2 = h2,
  preal_co = 0,
  prandom_co = prandom_co
)

## 4.2 For each gene pair:
##     - Step 1: cell-ratio test, does the observed high–high frequency
##       exceed the random expectation? (prop.test)
##     - Step 2: proliferation dependence, high–high vs non-high / zero
##       in terms of drcc (Wilcoxon test)

pair_stat_fun <- function(mm,
                          min_cells_hh   = 3,  ## minimum cell count for high–high
                          min_cells_rest = 3,  ## minimum cell count for non-high/non-zero
                          min_cells_zero = 3)  ## minimum cell count for zero group
{
  genedata <- cellls[i1[mm], ] + cellls[i2[mm], ]
  
  idx_hh   <- which(genedata == 4)  ## 2 + 2
  idx_00   <- which(genedata == 0)  ## 0 + 0
  idx_rest <- setdiff(seq_along(genedata), c(idx_00, idx_hh))
  
  zero_num  <- length(idx_00)
  nzero_num <- length(idx_rest)
  hh_num    <- length(idx_hh)
  
  ## ---- Step 1: cell-ratio test (high–high vs random expectation) ----
  co_grep <- prop.test(
    hh_num,
    n_cells,
    ptj$prandom_co[mm],
    alternative = "greater",
    correct = TRUE
  )$p.value
  
  ## ---- Step 2a: proliferation, high–high vs non-high/non-zero ----
  if (length(idx_hh) >= min_cells_hh && length(idx_rest) >= min_cells_rest) {
    mean_hh   <- mean(drcc[idx_hh],   na.rm = TRUE)
    mean_rest <- mean(drcc[idx_rest], na.rm = TRUE)
    
    fd1 <- mean_hh / mean_rest
    
    p1 <- as.numeric(
      wilcox.test(
        drcc[idx_hh],
        drcc[idx_rest],
        alternative = "greater",
        exact = FALSE
      )$p.value
    )
    
    drccres1 <- c(mean_hh, mean_rest, fd1, p1)
  } else {
    drccres1 <- c(NA, NA, NA, NA)
  }
  
  ## ---- Step 2b: proliferation, high–high vs zero-expression ----
  if (length(idx_hh) >= min_cells_hh && length(idx_00) >= min_cells_zero) {
    mean_hh   <- mean(drcc[idx_hh], na.rm = TRUE)
    mean_zero <- mean(drcc[idx_00], na.rm = TRUE)
    
    fd2 <- mean_hh / mean_zero
    
    p2 <- as.numeric(
      wilcox.test(
        drcc[idx_hh],
        drcc[idx_00],
        alternative = "greater",
        exact = FALSE
      )$p.value
    )
    
    drccres2 <- c(mean_hh, mean_zero, fd2, p2)
  } else {
    drccres2 <- c(NA, NA, NA, NA)
  }
  
  c(
    zero_num, nzero_num, hh_num,
    co_grep,
    drccres1,
    drccres2
  )
}


t1 <- Sys.time()
pair_stat_list <- future_lapply(
  seq_along(i1),
  FUN = pair_stat_fun,
  future.seed = TRUE
)
pair_stat_mat <- do.call(rbind, pair_stat_list)
t2 <- Sys.time()
print(t2 - t1)

colnames(pair_stat_mat) <- c(
  "zero_num", "nzero_num", "hh_num",
  "co_grep",
  "hhexp_mean1", "nzexp_mean", "co_fd1", "co_p1",
  "hhexp_mean2", "zeroexp_mean", "co_fd2", "co_p2"
)

pair_stat_df <- as.data.frame(pair_stat_mat)
preal_co_p   <- pair_stat_df$hh_num / n_cells

gires_step12 <- cbind(
  genei = genei,
  genej = genej,
  ptj,
  preal_co_p = preal_co_p,
  pair_stat_df
)

## FDR correction for cell-ratio and proliferation tests
gires_step12$co_grefdr    <- p.adjust(gires_step12$co_grep, method = "BH")
gires_step12$co_drcc1fdr  <- p.adjust(gires_step12$co_p1,   method = "BH")
gires_step12$co_drcc2fdr  <- p.adjust(gires_step12$co_p2,   method = "BH")

saveRDS(
  gires_step12,
  file.path(OUT_DIR, "gires_Epithelial_step12_cellratio_prolif.rds")
)

## =========================
## 5. Step 3: CS-CORE co-expression network
## =========================
expr_counts <- GetAssayData(
  scRNA[feature_genes, ],
  assay = "RNA",
  layer = "counts"   ## Seurat v5: use layer for counts
)
expr_counts <- as.matrix(expr_counts)
X <- t(expr_counts)           ## n cells × p genes
seq_depth <- rowSums(X)
CSCORE_Epithelial <- CSCORE_IRLS_cpp(
  X         = X,
  seq_depth = seq_depth
)
CSCORE_coexp <- CSCORE_Epithelial$est
CSCORE_p     <- CSCORE_Epithelial$p_value
CSCORE_r_vec <- CSCORE_coexp[lower.tri(CSCORE_coexp)]
CSCORE_p_vec <- CSCORE_p[lower.tri(CSCORE_p)]

## Attach CS-CORE statistics to the pair table
gires_step3 <- gires_step12
gires_step3$CSCORE_r <- CSCORE_r_vec
gires_step3$CSCORE_p <- CSCORE_p_vec
gires_step3$co_CSCOREfdr <- p.adjust(gires_step3$CSCORE_p, method = "BH")
gires_step3$rate_cop <- gires_step3$preal_co_p / gires_step3$prandom_co

rownames(gires_step3) <- paste(gires_step3$genei, gires_step3$genej, sep = "_")
gires_step3$pair <- rownames(gires_step3)

saveRDS(
  gires_step3,
  file.path(OUT_DIR, "gires_Epithelial_step3_CSCORE.rds")
)

## =========================
## 6. Step 4: Functional similarity (FS, GOSemSim)
## =========================

## Map feature_genes to ENTREZID for GO semantic similarity
edif_geneid <- bitr(
  feature_genes,
  fromType = "SYMBOL",
  toType   = c("ENTREZID"),
  OrgDb    = "org.Hs.eg.db"
)

## Compute FS and save
simres <- compute_FS_epithelial(
  edif_geneid,
  out_dir = OUT_DIR
)

## Keep pair + BP/CC/MF + fsim
simres1 <- simres[, c("pair", "simbp", "simcc", "simmf", "fsim")]

## Merge FS into gires_step3
gires_step4 <- dplyr::left_join(
  gires_step3,
  simres1,
  by = "pair"
)

saveRDS(
  gires_step4,
  file.path(OUT_DIR, "gires_Epithelial_step4_with_FS.rds")
)

## =========================
## (Optional) Final scPGI filtering example
## =========================
## You may change the thresholds according to the manuscript:
##   - co_grefdr    < 0.01  (cell-ratio test)
##   - co_drcc1fdr  < 0.01  (proliferation test 1)
##   - co_drcc2fdr  < 0.01  (proliferation test 2)
##   - co_CSCOREfdr < 0.01  (CS-CORE co-expression)
##   - CSCORE_r     > 0.1   (CS-CORE effect size)
##   - fsim         > 0.5   (functional similarity)

final_scPGI <- subset(
  gires_step4,
  co_grefdr    < 0.01 &
    co_drcc1fdr  < 0.01 &
    co_drcc1fdr  < 0.01 &
    co_CSCOREfdr < 0.01 &
    CSCORE_r     > 0.1 &
    fsim         > 0.5
)

saveRDS(
  final_scPGI,
  file.path(OUT_DIR, "Epithelial_final_scPGI.rds")
)
