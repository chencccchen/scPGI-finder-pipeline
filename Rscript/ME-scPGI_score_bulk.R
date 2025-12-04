## =======================================================================
## Calculate ME-scPGI score for each cell lines or bulk samples
## =======================================================================
setwd("path/data/")

## -----------------------------------------------
## 1. Read ME-scPGI gene pairs and expression matrix
## -----------------------------------------------
ME_scPGI <- read.table("ME-scPGI.txt", header = TRUE, sep = "\t")
ME_scPGI$pair <- paste(ME_scPGI$Gene1, ME_scPGI$Gene2, sep = "_")

## -----------------------------------------------
## 2. Function to compute me-scPGI score
##    expr: gene x sample expression matrix
##    pairs_df: data.frame with columns Gene1, Gene2, pair
## -----------------------------------------------
compute_me_scPGI_score <- function(expr, pairs_df) {
  expr <- as.matrix(expr)
  
  ## Keep only pairs whose both genes are in the expression matrix
  genes_in_expr <- rownames(expr)
  keep_idx <- which(pairs_df$Gene1 %in% genes_in_expr &
                      pairs_df$Gene2 %in% genes_in_expr)
  pairs_use <- pairs_df[keep_idx, ]
  
  if (nrow(pairs_use) == 0) {
    stop("No ME-scPGI pairs found in the expression matrix.")
  }
  
  ## 0/1 discretization for each gene:
  ## res[j] = 1 if expr_ij > 2/3 quantile (per gene), otherwise 0
  expls <- function(x) {
    x <- as.numeric(x)
    res <- rep(0, length(x))
    thr <- stats::quantile(x, 2/3, na.rm = TRUE)
    res[which(x > thr)] <- 1
    res
  }
  
  genels <- t(apply(expr, 1, expls))  ## gene x sample (0/1)
  colnames(genels) <- colnames(expr)
  
  ## Match gene1/gene2 to rows
  i1 <- match(pairs_use$Gene1, rownames(genels))
  i2 <- match(pairs_use$Gene2, rownames(genels))
  
  ## Sum of two genes: 0,1,2
  gp <- genels[i1, ] + genels[i2, ]
  rownames(gp) <- pairs_use$pair
  
  ## Convert to 0/1 for “both high” (2 -> 1, others -> 0)
  gp[gp == 1] <- 0
  gp[gp == 2] <- 1
  
  ## Remove pairs with NA (genes missing in expression)
  if (length(which(is.na(gp[, 1]))) != 0) {
    gptemp <- gp[-which(is.na(gp[, 1])), , drop = FALSE]
  } else {
    gptemp <- gp
  }
  
  ## me-scPGI score per sample = mean of pair scores
  me_scPGI_score <- apply(gptemp, 2, sum) / nrow(gptemp)
  return(me_scPGI_score)
}

## -----------------------------------------------
## 3. Example: bulk and cell-line me-scPGI scores
## -----------------------------------------------
## bulkf      : gene x bulk-sample matrix
## cellf      : gene x cell-line-sample matrix


## bulkf  <- readRDS("bulk_expression.rds")      # or read.table(...)
## cellf  <- readRDS("cellline_expression.rds")
bulkf<-readRDS("TCGA_LUAD_Expression_matrix.rds")

# bulk_me_scPGI  <- compute_me_scPGI_score(bulkf,  ME_scPGI)
# cell_me_scPGI  <- compute_me_scPGI_score(cellf,  ME_scPGI)
bulk_me_scPGI  <- compute_me_scPGI_score(bulkf,  ME_scPGI)

## 结果：
## bulk_me_scPGI  and cell_me_scPGI are named numeric vector，
## ME-scPGI score per bulk smaple and celllines
