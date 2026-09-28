compare_tissues_metrics <- function(df_A, df_B, nperm = 1000) {
  # df_A and df_B should have the same required columns: log2FoldChange.SS, log2FoldChange.OO
  
  # Observed metrics
  obs_spearman_A <- cor(df_A$log2FoldChange.SS, df_A$log2FoldChange.OO, method = "spearman")
  obs_spearman_B <- cor(df_B$log2FoldChange.SS, df_B$log2FoldChange.OO, method = "spearman")
  obs_spearman_diff <- obs_spearman_A - obs_spearman_B
  
  obs_euc_A <- sqrt(sum((df_A$log2FoldChange.SS - df_A$log2FoldChange.OO)^2, na.rm = TRUE))
  obs_euc_B <- sqrt(sum((df_B$log2FoldChange.SS - df_B$log2FoldChange.OO)^2, na.rm = TRUE))
  obs_euc_diff <- obs_euc_A - obs_euc_B
  
  # Prepare permutation
  nA <- nrow(df_A)
  nB <- nrow(df_B)
  all_genes <- rbind(df_A, df_B)
  
  perm_spearman_diff <- numeric(nperm)
  perm_euc_diff <- numeric(nperm)
  
  for (i in 1:nperm) {
    perm_idx <- sample(rep(c("A", "B"), times = c(nA, nB)))
    df_Ap <- all_genes[perm_idx == "A", ]
    df_Bp <- all_genes[perm_idx == "B", ]
    
    # Spearman
    rho_A <- cor(df_Ap$log2FoldChange.SS, df_Ap$log2FoldChange.OO, method = "spearman")
    rho_B <- cor(df_Bp$log2FoldChange.SS, df_Bp$log2FoldChange.OO, method = "spearman")
    perm_spearman_diff[i] <- rho_A - rho_B
    
    # Euclidean
    euc_A <- sqrt(sum((df_Ap$log2FoldChange.SS - df_Ap$log2FoldChange.OO)^2, na.rm = TRUE))
    euc_B <- sqrt(sum((df_Bp$log2FoldChange.SS - df_Bp$log2FoldChange.OO)^2, na.rm = TRUE))
    perm_euc_diff[i] <- euc_A - euc_B
  }
  
  # p-values (two-sided)
  pval_spearman <- mean(abs(perm_spearman_diff) >= abs(obs_spearman_diff))
  pval_euc <- mean(abs(perm_euc_diff) >= abs(obs_euc_diff))
  
  list(
    spearman = list(
      obs_A = obs_spearman_A,
      obs_B = obs_spearman_B,
      obs_diff = obs_spearman_diff,
      pval = pval_spearman,
      perm_dist = perm_spearman_diff
    ),
    euclidean = list(
      obs_A = obs_euc_A,
      obs_B = obs_euc_B,
      obs_diff = obs_euc_diff,
      pval = pval_euc,
      perm_dist = perm_euc_diff
    )
  )
}
