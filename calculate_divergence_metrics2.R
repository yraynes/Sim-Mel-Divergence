library(dplyr)

calculate_divergence_metrics2 <- function(df, tissue_name = "Tissue", n_boot = 1000) {
  # Ensure required columns exist
  required_cols <- c("log2FoldChange.SS", "log2FoldChange.OO", "padj.SS", "padj.OO")
  if (!all(required_cols %in% names(df))) {
    stop("Missing required columns in dataframe")
  }
  
  sex <- strsplit(tissue_name, split = " ")[[1]][1]
  tissue <- strsplit(tissue_name, split = " ")[[1]][2]
  
  # compute metrics
  compute_metrics <- function(x, y) {
    list(
      Spearman   = cor(x, y, method = "spearman", use = "complete.obs"),
      Pearson    = cor(x, y, method = "pearson",  use = "complete.obs"),
      Euclidean  = sqrt(sum((x - y)^2, na.rm = TRUE))
    )
  }
  
  #  bootstrap CIs
  bootstrap_ci <- function(x, y, n = n_boot) {
    n_genes <- length(x)
    boot_res <- replicate(n, {
      idx <- sample(seq_len(n_genes), replace = TRUE)
      metrics <- compute_metrics(x[idx], y[idx])
      unlist(metrics)
    })
    
    boot_res <- t(boot_res)
    ci <- apply(boot_res, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    ci <- as.data.frame(t(ci))
    names(ci) <- c("low", "high")
    ci
  }
  
  # Subsets
  all_genes <- df
  lfc_thresh <- df %>%
    dplyr::filter(abs(log2FoldChange.SS) > 0.58 | abs(log2FoldChange.OO) > 0.58)
  padj_sig <- df %>%
    dplyr::filter(padj.SS < 0.05 | padj.OO < 0.05)
  
  # Compute metrics + CIs for each subset
  process_subset <- function(subset_df) {
    metrics <- compute_metrics(subset_df$log2FoldChange.SS, subset_df$log2FoldChange.OO)
    ci <- bootstrap_ci(subset_df$log2FoldChange.SS, subset_df$log2FoldChange.OO)
    
    tibble::tibble(
      Spearman      = metrics$Spearman,
      Spearman_low  = ci["Spearman", "low"],
      Spearman_high = ci["Spearman", "high"],
      OneMinusRho      = 1 - metrics$Spearman,
      OneMinusRho_low  = 1 - ci["Spearman", "high"], # flip order
      OneMinusRho_high = 1 - ci["Spearman", "low"],
      Pearson       = metrics$Pearson,
      Pearson_low   = ci["Pearson", "low"],
      Pearson_high  = ci["Pearson", "high"],
      Euclidean     = metrics$Euclidean,
      Euclidean_low = ci["Euclidean", "low"],
      Euclidean_high= ci["Euclidean", "high"]
    )
  }
  
  # Apply
  metrics_all <- process_subset(all_genes)
  metrics_lfc <- process_subset(lfc_thresh)
  metrics_sig <- process_subset(padj_sig)
  
  # Final tibble
  dplyr::bind_rows(
    dplyr::mutate(metrics_all, Tissue = tissue, Sex = sex, Subset = "All Genes"),
    dplyr::mutate(metrics_lfc, Tissue = tissue, Sex = sex, Subset = "|LFC| > 0.58"),
    dplyr::mutate(metrics_sig, Tissue = tissue, Sex = sex, Subset = "padj < 0.05")
  )
}
