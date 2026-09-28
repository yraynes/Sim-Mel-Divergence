library(dplyr)

calculate_divergence_metrics_NES <- function(df, tissue_name = "Tissue", n_boot = 1000) {
  # Ensure required columns exist
  required_cols <- c("NES_SS", "NES_OO")
  if (!all(required_cols %in% names(df))) {
    stop("Missing required columns in dataframe")
  }
  
  sex <- strsplit(tissue_name, split = " ")[[1]][1]
  tissue <- strsplit(tissue_name, split = " ")[[1]][2]
  
  # Helper: compute metrics
  compute_metrics <- function(x, y) {
    list(
      Spearman   = cor(x, y, method = "spearman", use = "complete.obs"),
      Pearson    = cor(x, y, method = "pearson",  use = "complete.obs"),
      Euclidean  = sqrt(sum((x - y)^2, na.rm = TRUE))
    )
  }
  
  # bootstrap CIs
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

  process_subset <- function(subset_df) {
    metrics <- compute_metrics(subset_df$NES_SS, subset_df$NES_OO)
    ci <- bootstrap_ci(subset_df$NES_SS, subset_df$NES_OO)
    
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

  process_subset(df) %>%
    dplyr::mutate(Tissue = tissue, Sex = sex)

}
