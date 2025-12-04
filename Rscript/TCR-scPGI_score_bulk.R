#!/usr/bin/env Rscript
##
## TCR-scPGI score for bulk samples or cell lines
## 1) Gene-wise discretization (2/3 quantile -> 0/1)
## 2) Build pair-level 0/1 matrix
## 3) Apply lasso coefficients (TCRTGI) to get TCR-scPGI score
##

## =======================================================================
## 0. Set working directory and load data
## =======================================================================
setwd("path/data/")   # <-- change to your folder


## Expression matrix: gene x sample
## (e.g. TCGA LUAD bulk expression)
bulk_expr <- readRDS("TCGA_LUAD_Expression_matrix.rds")

## Lasso model for TCR-scPGI (TCRTGI)
## rownames: pair IDs + last row as intercept / meta
## column "b": coefficients (including intercept as last value)
TCRTGI <- readRDS("TCR-scPGI.rds")
TCRTGI_cut <- 1.022312    # decision threshold used in your previous code

## =======================================================================
## 1. Helper: build 0/1 pair matrix from expression and gene pairs
##    - expr:      gene x sample expression matrix
##    - pairs_df:  data.frame with columns Gene1, Gene2, pair
##    Output:
##      gp: pair x sample 0/1 matrix for "both genes high"
## =======================================================================
build_pair_matrix <- function(expr, pairs_df) {
  expr <- as.matrix(expr)
  
  ## keep only pairs whose both genes are in the expression matrix
  genes_in_expr <- rownames(expr)
  keep_idx <- which(pairs_df$gene1 %in% genes_in_expr &
                      pairs_df$gene2 %in% genes_in_expr)
  pairs_use <- pairs_df[keep_idx, ]
  
  if (nrow(pairs_use) == 0) {
    stop("No gene pairs found in the expression matrix.")
  }
  
  ## 0/1 discretization per gene:
  ## res[j] = 1 if expr_ij > 2/3 quantile (per gene), otherwise 0
  expls <- function(x) {
    x   <- as.numeric(x)
    res <- rep(0, length(x))
    thr <- stats::quantile(x, 2/3, na.rm = TRUE)
    res[which(x > thr)] <- 1
    res
  }
  
  genels <- t(apply(expr, 1, expls))   # gene x sample (0/1)
  colnames(genels) <- colnames(expr)
  
  ## match Gene1 / Gene2 to rows
  i1 <- match(pairs_use$gene1, rownames(genels))
  i2 <- match(pairs_use$gene2, rownames(genels))
  
  ## sum of the two genes: 0,1,2
  gp <- genels[i1, ] + genels[i2, ]
  rownames(gp) <- pairs_use$pair
  
  ## convert to 0/1 for "both high"
  gp[gp == 1] <- 0
  gp[gp == 2] <- 1
  
  ## remove pairs with NA (genes missing)
  if (any(is.na(gp[, 1]))) {
    gp <- gp[!is.na(gp[, 1]), , drop = FALSE]
  }
  rownames(gp)<-paste0(rownames(genels[i1, ]),"_",rownames(genels[i2, ]))
  gp
}

## =======================================================================
## 2. TCR-scPGI score
##    - expr:      gene x sample expression matrix
##    - pairs_df:  data.frame with Gene1, Gene2, pair
##    - coef_df:   lasso model (e.g. TCRTGI), with:
##                  * rownames: pair IDs + last row for intercept
##                  * column "b": coefficients
##    - cutoff:    decision threshold (optional); if provided,
##                 return both score and 0/1 prediction
## =======================================================================
compute_TCR_scPGI_score <- function(expr, pairs_df, coef_df) {
  
  ## 1) build pair 0/1 matrix: pair x sample
  gp_all <- build_pair_matrix(expr, pairs_df)
  
  ## 2) extract coefficients
  b <- as.numeric(coef_df$b)
  intercept <- b[length(b)]              # last value = intercept
  w <- b[1:(length(b) - 1)]              # first (n-1) = pair coefficients
  
  model_pairs <- rownames(coef_df)[1:(nrow(coef_df) - 1)]
  
  ## 3) match model pairs to pair matrix rows
  common_pairs <- intersect(model_pairs, rownames(gp_all))
  
  if (length(common_pairs) == 0) {
    stop("No overlap between model gene pairs and pair matrix.")
  }
  
  ## keep the same order as in coef_df
  idx_model  <- match(common_pairs, model_pairs)
  idx_matrix <- match(common_pairs, rownames(gp_all))
  w_use      <- w[idx_model]
  gp_use     <- gp_all[idx_matrix, , drop = FALSE]
  
  ## 4) linear score: w' * gp + intercept
  ##    score is length = n_sample
  score <- as.numeric(w_use %*% gp_use + intercept)
  names(score) <- colnames(gp_use)
  
  return(score)
}


## =======================================================================
## 3. Example: TCR-scPGI scores for TCGA LUAD bulk
## =======================================================================

TCR_scPGI_score <- compute_TCR_scPGI_score(
  expr     = bulk_expr,
  pairs_df = TCRTGI,   # or a TCR-specific pair list if you have one
  coef_df  = TCRTGI
)

## Results:
##   TCR_scPGI_score : numeric TCR-scPGI score for each sample

