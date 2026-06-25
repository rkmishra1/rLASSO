# Script to parse simulation results and generate Markdown tables for README.md

load("results/simulation_results.RData")

cat("\n--- PARSING SIMULATION RESULTS ---\n")

methods <- c("Lasso", "ALasso", "BLasso", "Horseshoe", "BayesA", "BayesB", "BayesC", "rLASSO (S5)")

scenarios_names <- names(results_mse)

for (scen_name in scenarios_names) {
  cat("\n###", scen_name, "\n\n")
  
  mse_mat <- results_mse[[scen_name]]
  bar_mat <- results_bar[[scen_name]]
  
  mean_mse <- colMeans(mse_mat)
  sd_mse <- apply(mse_mat, 2, sd)
  sem_mse <- sd_mse / sqrt(nrow(mse_mat))
  
  mean_bar <- colMeans(bar_mat)
  sd_bar <- apply(bar_mat, 2, sd)
  sem_bar <- sd_bar / sqrt(nrow(bar_mat))
  
  # Format as Mean (SD) or Mean (SEM)
  # Let's use Mean (SD)
  tbl <- data.frame(
    Method = methods,
    MSE = sprintf("%.3f (%.3f)", mean_mse, sd_mse),
    BAR = sprintf("%.3f (%.3f)", mean_bar, sd_bar)
  )
  
  # Print as Markdown Table
  cat("| Method | MSE (SD) | BAR (SD) |\n")
  cat("| :--- | :--- | :--- |\n")
  for (i in 1:nrow(tbl)) {
    cat(sprintf("| %s | %s | %s |\n", tbl$Method[i], tbl$MSE[i], tbl$BAR[i]))
  }
}

# Save a CSV file for general use
res_list <- list()
for (scen_name in scenarios_names) {
  mse_mat <- results_mse[[scen_name]]
  bar_mat <- results_bar[[scen_name]]
  
  mean_mse <- colMeans(mse_mat)
  sd_mse <- apply(mse_mat, 2, sd)
  mean_bar <- colMeans(bar_mat)
  sd_bar <- apply(bar_mat, 2, sd)
  
  df_scen <- data.frame(
    Scenario = scen_name,
    Method = methods,
    Mean_MSE = mean_mse,
    SD_MSE = sd_mse,
    Mean_BAR = mean_bar,
    SD_BAR = sd_bar
  )
  res_list[[scen_name]] <- df_scen
}
final_df <- do.call(rbind, res_list)
write.csv(final_df, "results/simulation_summary_table.csv", row.names = FALSE)
cat("\nSummary table saved to results/simulation_summary_table.csv\n")
