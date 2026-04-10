## ---- load libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(grid)
library(scales)
library(cowplot)  # for get_legend()
library(corrplot)
library(patchwork)
library(vegan)
library(ggnewscale)
library(stringr)
library(softImpute)
library(ecodist)
library(rstatix)
library(ggrepel)
library(pROC)
library(PRROC)

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 18, color = "black"),
              axis.title = element_text(size = 18, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- parameters ----
Study_id_chose <- "1_Bartomeus"
n_sim <- 100
prop_ones_to_remove <- 0.20


set.seed(123)


## ---- load data ----
# read data for the Study_id_chose
Interaction_data <- read.csv("Interaction_edges.csv", encoding = "latin1", stringsAsFactors = FALSE)
one_study_interactions <- Interaction_data %>% 
  filter(Study_id == Study_id_chose)

# map the layers 
layer_names <- unique(one_study_interactions$layer_from)
layer_map <- data.frame(
  original = layer_names,
  new_id = paste0("layer_", seq_along(layer_names))
)

# apply mapped layer names to keep layer labels consistent across the pipeline
one_study_interactions <- one_study_interactions %>%
  left_join(layer_map, by = c("layer_from" = "original")) %>%
  mutate(layer_from = new_id) %>%
  select(-new_id) %>%
  left_join(layer_map, by = c("layer_to" = "original")) %>%
  mutate(layer_to = new_id) %>%
  select(-new_id)


## ---- functions ----

# building matrices from edge list
build_interaction_matrix <- function(data, layers_to_filter) {
  # Step 1: Filter rows based on specified layers
  layers <- paste0("layer_", layers_to_filter)
  filtered_data <- subset(data, layer_from %in% layers)
  
  # Step 2: Aggregate weights for identical species pairs
  aggregated_data <- filtered_data %>%
    group_by(node_from, node_to) %>%
    summarise(weight = sum(weight), .groups = 'drop')
  
  # Step 3: Create the matrix with specific row and column species
  species_from <- unique(aggregated_data$node_from)  # Columns
  species_to <- unique(aggregated_data$node_to)      # Rows
  
  # Initialize an empty matrix
  interaction_matrix <- matrix(0, nrow = length(species_to), ncol = length(species_from),
                               dimnames = list(species_to, species_from))
  
  # Fill the matrix with aggregated weights
  for (i in 1:nrow(aggregated_data)) {
    row <- aggregated_data$node_to[i]    # Rows represent pollinators species
    col <- aggregated_data$node_from[i]  # Columns represent plant species
    interaction_matrix[row, col] <- aggregated_data$weight[i]
  }
  
  return(interaction_matrix)
}

# transforming raw predictions for binary evaluation
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# cutoff01 <- function(x) {
#   case_when(
#     x < 0 ~ 0,
#     x > 1 ~ 1,
#     TRUE ~ x
#   )
# }


# predict links using SVD - softImpute
implement_impute <- function(C, k, lambda) {
  # Apply softImpute
  
  fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 600)

  
  # Reconstruct the matrix
  C_reconstructed <- softImpute::complete(C, fit)
  
  # Extract the reconstructed P matrix from C_reconstructed
  P_reconstructed <- C_reconstructed[rownames(P), colnames(P)]
  
  # Combine indices of removed ones and zeros
  if (is.null(dim(remove_indices))) { # handles when remove_indices has only one row
    test_indices <- rbind(
      data.frame(row = remove_indices["row"], col = remove_indices["col"], label = 1),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = rep(0, nrow(zeros_to_remove_indices)))
    )
  } else {
    test_indices <- rbind(
      data.frame(row = remove_indices[, "row"], col = remove_indices[, "col"], label = rep(1, nrow(remove_indices))),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = rep(0, nrow(zeros_to_remove_indices)))
    )
  }
  
  # Get the row and column names of the test links
  test_rows <- rownames(P)[test_indices$row]  # These are the "node_to"
  test_cols <- colnames(P)[test_indices$col]  # These are the "node_from"
  
  # Actual labels and predictions
  original_links <- P_original[cbind(test_rows, test_cols)]
  predicted_values <- P_reconstructed[cbind(test_rows, test_cols)]
  
  # Store the results with node information
  results <- data.frame(k = k,
                        lambda = lambda,
                        original_links = original_links,
                        predicted_values = predicted_values,
                        node_to = test_rows,
                        node_from = test_cols,
                        removed = 1 # mark these links as removed
  )
  
  # ---- (A) List all edges in P
  all_edges <- expand.grid(
    node_to = rownames(P),
    node_from = colnames(P),
    k = k,
    lambda = lambda,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Fill in the original link values from P_original
  all_edges$original_links <- mapply(
    function(r, c) P_original[r, c],
    all_edges$node_to,
    all_edges$node_from
  )
  
  # Helper data frame of removed edges
  removed_edges_idx <- data.frame(
    node_to = test_rows,
    node_from = test_cols,
    stringsAsFactors = FALSE
  )
  
  # ---- (B) Subset edges NOT removed
  not_removed <- all_edges[
    !paste(all_edges$node_to, all_edges$node_from) %in%
      paste(removed_edges_idx$node_to, removed_edges_idx$node_from), # all node pairs that are not in th removed list
  ]
  not_removed$removed <- 0
  not_removed$predicted_values <- NA
  
  # ---- (C) Combine removed + not removed
  not_removed$k <- k
  not_removed$lambda <- lambda
  
  list_results <- list(results=results, not_removed=not_removed)
  
  return(list_results)
}



## ---- 1. prediction ----
### ---- Prepare data ----

aggregated_df_weighted <- one_study_interactions %>%
  filter(layer_from == layer_to) %>%
  select(layer_from, node_from, layer_to, node_to, weight, type)

# Convert to binary (presence/absence)
aggregated_df_binary <- aggregated_df_weighted %>%
  mutate(weight = if_else(weight > 0, 1, 0))

# Print summary
cat("Study:", Study_id_chose, "\n")
cat("Number of interactions:", nrow(aggregated_df_weighted), "\n")
cat("Number of unique layers:", length(unique(aggregated_df_weighted$layer_from)), "\n")

# Get layer information
num_layers <- length(unique(aggregated_df_weighted$layer_from))
layer_names <- unique(aggregated_df_weighted$layer_from)

### ---- Set up file paths and check for existing results ----
results_file <- paste0("predictions_", Study_id_chose, "_weighted_1.rds")

# Initialize combined results
combined_results <- data.frame()

if (file.exists(results_file)) {
  cat("Reading existing results file...\n")
  combined_results <- readRDS(results_file)
  cat("Predictions loaded successfully\n")
  
} else {
  cat("Generating new predictions...\n")
  
  # Loop through each layer (diagonal: train_layer == test_layer)
  for (layer in 1:num_layers) {
    cat("\n** Processing layer:", layer_names[layer], "**\n")
    
    # Build interaction matrix for this layer
    P <- build_interaction_matrix(data = aggregated_df_weighted, layers_to_filter = layer)
    P_original <- P  # Save original P
    
    # Identify links and non-links (zeros) to withhold
    ones_in_P <- which(P > 0, arr.ind = TRUE)
    zeros_in_P <- which(P == 0, arr.ind = TRUE)
    
    # Calculate number of links and zeros to withhold (20% of links, equal number of zeros)
    num_1_to_remove <- floor(sum(P > 0, na.rm = TRUE) * prop_ones_to_remove)
    num_0_to_remove <- num_1_to_remove
    
    cat("Links to withhold:", num_1_to_remove, "out of", nrow(ones_in_P), "\n")
    cat("Non-links to withhold (0):", num_0_to_remove, "out of", nrow(zeros_in_P), "\n")
    
    # Initialize bootstrap results and stratified sampling pools
    bootstrapping_results <- NULL
    
    ### ---- Stratified Sampling: Create pools ensuring all edges withheld at least once ----
    # Pool for links
    ones_pool_size <- num_1_to_remove * n_sim
    ones_pool_indices <- rep(1:nrow(ones_in_P), length.out = ones_pool_size)
    ones_pool_shuffled <- sample(ones_pool_indices)  # Shuffle for randomness
    
    # Pool for zeros (non-links)
    zeros_pool_size <- num_0_to_remove * n_sim
    zeros_pool_indices <- rep(1:nrow(zeros_in_P), length.out = zeros_pool_size)
    zeros_pool_shuffled <- sample(zeros_pool_indices)  # Shuffle for randomness
    
    ### ---- Bootstrap loop: iterate n_sim times ----
    for (i in 1:n_sim) {
      # Calculate indices for this iteration from the stratified pools
      chunk_start <- (i - 1) * num_1_to_remove + 1
      chunk_end <- i * num_1_to_remove
      
      # Get indices for links and zeros from shuffled pools
      ones_indices_this_iter <- ones_pool_shuffled[chunk_start:chunk_end]
      zeros_indices_this_iter <- zeros_pool_shuffled[chunk_start:chunk_end]
      
      # Select the actual edges to withhold
      remove_indices <- ones_in_P[ones_indices_this_iter, ]
      P[remove_indices] <- NA
      
      # Randomly select non-links to withhold
      zeros_to_remove_indices <- zeros_in_P[zeros_indices_this_iter, ]
      P[zeros_to_remove_indices] <- NA
      
      # Create combined matrix C
      all_row_ids <- rownames(P)
      all_col_ids <- colnames(P)
      C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                  dimnames = list(all_row_ids, all_col_ids))
      C[rownames(P), colnames(P)] <- P
      
      # Apply biScale normalization
      C <- biScale(C, row.center = TRUE, col.center = TRUE, row.scale = FALSE, col.scale = FALSE)
      
      ### ---- Test different k and lambda combinations ----
      k_values <- c(2, 5, 10)
      lam0 <- lambda0(C)
      lambda_values <- c(1, 5, 50, 100, lam0)
      
      results <- data.frame(k = integer(),
                            lambda = numeric(),
                            original_links = numeric(),
                            predicted_values = numeric(),
                            input_lambda = numeric())
      not_removed_all <- NULL
      
      # Loop over k and lambda combinations
      for (k in k_values) {
        for (lambda in lambda_values) {
          r <- implement_impute(C, k, lambda)
          r$results$input_lambda <- lambda
          r$not_removed$input_lambda <- lambda
          results <- rbind(results, r$results)
          not_removed_all <- rbind(not_removed_all, r$not_removed)
        }
      }
      
      # Combine results and add iteration number
      complete_edges_all <- rbind(results, not_removed_all)
      complete_edges_all$itr <- i
      bootstrapping_results <- rbind(bootstrapping_results, complete_edges_all)
      
      # Reset P for next iteration
      P <- P_original
    }
    
    # Add layer metadata and combine with overall results
    combined_results <- rbind(
      combined_results,
      cbind(
        data.frame(
          study_id = Study_id_chose,
          train_layer = layer,
          test_layer = layer,
          prop_ones_removed = prop_ones_to_remove,
          amount_of_removed_1 = num_1_to_remove,
          amount_of_removed_0 = num_0_to_remove
        ),
        bootstrapping_results
      )
    )
  }
  
  ### ---- Coverage check ----
  cat("\n=== Coverage Check ===\n")
  for (layer in 1:num_layers) {
    P_check <- build_interaction_matrix(data = aggregated_df_weighted, layers_to_filter = layer)
    
    # Create all possible edges
    all_edges <- expand.grid(
      node_to = rownames(P_check),
      node_from = colnames(P_check),
      stringsAsFactors = FALSE
    )
    all_edges$original <- mapply(
      function(r, c) P_check[r, c],
      all_edges$node_to,
      all_edges$node_from,
      USE.NAMES = FALSE
    )
    
    # Get withheld edges from results
    layer_results <- combined_results %>%
      filter(train_layer == layer & test_layer == layer & removed == 1) %>%
      distinct(node_to, node_from)
    
    # Find missing edges
    withheld_keys <- paste(layer_results$node_to, layer_results$node_from)
    all_keys <- paste(all_edges$node_to, all_edges$node_from)
    missing_edges <- all_edges[!all_keys %in% withheld_keys, ]
    
    missing_ones <- sum(missing_edges$original == 1)
    missing_zeros <- sum(missing_edges$original == 0)
    
    cat(sprintf("Layer %s: total edges = %d; never withheld = %d (ones=%d, zeros=%d)\n",
                layer_names[layer], nrow(all_edges), nrow(missing_edges), missing_ones, missing_zeros))
  }
  
  ### ---- Iteration check ----
  cat("\n=== Iteration Count Check ===\n")
  for (layer in 1:num_layers) {
    layer_results <- combined_results %>%
      filter(train_layer == layer & test_layer == layer)
    
    observed_itrs <- sort(unique(layer_results$itr))
    missing_itrs <- setdiff(seq_len(n_sim), observed_itrs)
    
    cat(sprintf("Layer %s: expected = %d; observed = %d\n",
                layer_names[layer], n_sim, length(observed_itrs)))
    
    if (length(missing_itrs) > 0) {
      warning(sprintf("Layer %s missing iterations: %s",
                      layer_names[layer], paste(missing_itrs, collapse = ", ")))
    }
  }
  
  # Save results
  saveRDS(combined_results, file = results_file)
  cat("\nPrediction results saved to", results_file, "\n")
} 




## ---- 2. Analysis and Evaluation ----

### ---- Prepare data ----
# Filter to best model and ensure valid predictions
combined_results <- combined_results %>% 
  filter(train_layer == test_layer) %>%
  filter(k == 2) %>% 
  filter(!(input_lambda %in% c(1, 5, 50, 100)))

# Clean negative predictions
df <- combined_results %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

cat("Summary of filtered predictions:\n")
print(summary(df))

### ---- Network composition ----
plant_species <- unique(df$node_from)
pollinator_species <- unique(df$node_to)
cat("Plant species:", length(plant_species), "\n")
cat("Pollinator species:", length(pollinator_species), "\n")

### ---- Threshold ----
# select the threshold for classifying links as 1s or 0s based on max f0.5

# a) set an array of thresholds
thresholds <- seq(0, 1, by = 0.1)

# b) filter & prep
df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  )

# c) expand to one row per threshold
df_thresh <- df_prepped %>%
  tidyr::expand_grid(threshold = thresholds) %>%  
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by(study_id, train_layer, test_layer, itr, threshold) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    specificity = TN / (TN + FP),
    precision   = TP / (TP + FP),
    recall      = TP / (TP + FN),
    f05_score   = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    balanced_accuracy= (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) /
      sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)),
    mse  = mean((predicted_values - original_links)^2, na.rm = TRUE),
    rmse = sqrt(mse)#,
    #.groups = "drop"
  ) %>%
  ungroup() %>%
  group_by(study_id, train_layer, test_layer, threshold) %>%
  summarise(
    TP = mean(TP, na.rm = TRUE),
    FN = mean(FN, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f05_score = mean(f05_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE),
    mse = mean(mse, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE)
  ) %>%
  ungroup() 

# d) average across emln_id/layer combos and pivot long
df_avg <- df_thresh %>%
  group_by(threshold) %>%
  summarise(across(
    c(specificity, precision, recall, 
      f05_score, balanced_accuracy, mcc),
    \(x) mean(x, na.rm = TRUE)
  )) %>%
  pivot_longer(-threshold,
               names_to  = "metric",
               values_to = "value")

df_avg_plot <- df_avg %>% mutate(metric = recode(metric,
                                                 specificity       = "Specificity",
                                                 precision         = "Precision",
                                                 recall            = "Recall",
                                                 f05_score         = "F0.5 score",
                                                 balanced_accuracy = "Balanced accuracy",
                                                 mcc               = "MCC"
))

# e) plot
optimal_threshold <- ggplot(df_avg_plot, aes(threshold, value, color = metric)) +
  geom_line(size = 1) +
  labs(
    x     = "Probability threshold",
    y     = "Average metric",
    color = "Metric"
  ) +
  scale_color_brewer(palette = "Pastel2") +
  tme

optimal_threshold

pdf(
  file   = "optimal_threshold_Bartomeus_1.pdf",
  width  = 7,    # inches
  height = 5,
  family = "Helvetica"   # or another installed font
)
print(optimal_threshold)
dev.off()     # close the file


df_eval <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob  = sigmoid(predicted_values),
    original_binary = if_else(original_links > 0, 1L, 0L)
  ) %>%
  group_by(study_id, train_layer, test_layer, itr) %>%
  summarise(
    # ROC-AUC (coerce to numeric!)
    auc_roc = tryCatch({
      roc_obj <- roc(response = original_binary,
                     predictor = predicted_prob,
                     quiet = TRUE, na.rm = TRUE,
                     levels = c(0,1), direction = "<")
      as.numeric(auc(roc_obj)) 
    }, error = function(e) NA_real_),
    
    # PR-AUC (guard against all-one-class cases)
    auc_pr = tryCatch({
      pos <- predicted_prob[original_binary == 1]
      neg <- predicted_prob[original_binary == 0]
      if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
      pr_obj <- pr.curve(scores.class0 = pos, scores.class1 = neg, curve = FALSE)
      pr_obj$auc.integral
    }, error = function(e) NA_real_)
  ) %>%
  ungroup()

df_eval_summary <- df_eval %>%
  group_by(study_id, train_layer, test_layer) %>%
  summarise(
    auc_roc_mean = mean(auc_roc, na.rm = TRUE),
    auc_roc_sd   = sd(auc_roc,   na.rm = TRUE),
    auc_pr_mean  = mean(auc_pr,  na.rm = TRUE),
    auc_pr_sd    = sd(auc_pr,    na.rm = TRUE),
    .groups = "drop"
  )


# 1) pivot to wide so F0.5 and balanced_accuracy are columns
df_wide <- df_avg %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(threshold)

# 2) find the threshold with the optimal f05 score
best_discrete <- df_wide %>%
  slice_max(f05_score, n = 1)

# results:
best_discrete_threshold <- best_discrete$threshold
best_discrete_threshold




### ---- Evaluation  ----

df_prepped <- df_prepped %>%
  mutate(predicted_bin_sigm = if_else(predicted_prob > best_discrete_threshold, 1, 0)) 


result_summary <- df_prepped %>%
  group_by(study_id, train_layer, test_layer, itr) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin_sigm == 1),
    FN = sum(original_binary == 1 & predicted_bin_sigm == 0),
    TN = sum(original_binary == 0 & predicted_bin_sigm == 0),
    FP = sum(original_binary == 0 & predicted_bin_sigm == 1),
    specificity = TN / (TN + FP),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    f05_score = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    balanced_accuracy = (recall + specificity) / 2,
    nse  = 1 - sum((predicted_values - original_links)^2, na.rm = TRUE) /
      sum((original_links   - mean(original_links, na.rm = TRUE))^2, na.rm = TRUE),
    nnse = 1 / (2 - nse)
  ) %>%
  ungroup() %>%
  group_by(study_id, train_layer, test_layer) %>%
  summarise(
    TP = mean(TP, na.rm = TRUE),
    FN = mean(FN, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f05_score = mean(f05_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    nse  = mean(nse,  na.rm = TRUE),
    nnse = mean(nnse, na.rm = TRUE)
  ) %>%
  ungroup()


head(result_summary)
# result_summary includes evaluation results across all iterations for each layer prediction
summary(result_summary) 


### ---- Box plot: TN, TP, FN, FP by layer  ----
# data per-iteration confusion matrix values
confusion_per_iter <- df_prepped %>%
  group_by(study_id, train_layer, test_layer, itr) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin_sigm == 1),
    FN = sum(original_binary == 1 & predicted_bin_sigm == 0),
    TN = sum(original_binary == 0 & predicted_bin_sigm == 0),
    FP = sum(original_binary == 0 & predicted_bin_sigm == 1),
    .groups = "drop"
  )

confusion_long <- confusion_per_iter %>%
  pivot_longer(
    cols = c(TP, FN, TN, FP),
    names_to = "confusion_category",
    values_to = "count"
  ) %>%
  mutate(
    confusion_category = factor(confusion_category, 
                               levels = c("TP", "FP", "TN", "FN")),
    layer = factor(train_layer)
  )

# Create box plot with all layers 
boxplot_confusion <- ggplot(confusion_long, 
                            aes(x = layer, 
                                y = count, 
                                fill = confusion_category)) +
  geom_boxplot(color = "black", alpha = 0.6, outlier.shape = 21, outlier.size = 1) +
  scale_fill_manual(
    values = c(
      "TP" = "#66C2A5",  
      "TN" = "#B2DF8A",  
      "FP" = "#FC8D62",  
      "FN" = "#E78AC3"
    ),
    name = "Confusion Matrix"
  ) +
  labs(
    # title = "Confusion Matrix Values Distribution Across All Layers",
    # subtitle = "Box plots showing TP, FP, TN, FN for each layer across all iterations",
    x = "Layer",
    y = "Count"
  ) +
  theme_bw(base_size = 10) +
  tme

print(boxplot_confusion)

#Save the plot
ggsave(
  filename = "boxplot_confusion_matrix_Bartomeus_1.png",
  plot = boxplot_confusion,
  width = 14,
  height = 7,
  dpi = 300
)



#Bar chart: metrics per layer 
# Prepare data for bar chart - showing key metrics per layer for best threshold 
df_metrics_per_layer <- result_summary %>%
  select(train_layer, train_layer, f05_score, balanced_accuracy, precision, recall) %>%
  pivot_longer(
    cols = c(f05_score, balanced_accuracy, precision, recall),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = recode(metric,
                         f05_score = "F0.5 score",
                         balanced_accuracy = "Balanced accuracy",
                         precision = "Precision",
                         recall = "Recall"))

# Create bar chart
metrics_per_layer_plot <- ggplot(df_metrics_per_layer, 
                                 aes(x = reorder(train_layer, train_layer), 
                                     y = value, 
                                     fill = metric)) +
  geom_col(position = "dodge", color = "black", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Layer",
    y = "Metric value",
    fill = "Metric",
    title = "Performance metrics across layers (diagonal analysis only)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  ylim(0, 1) +
  tme

print(metrics_per_layer_plot)




### ---- Distribution of predicted_values_sigmoid ----


# Histogram of predicted values
dist_plot_hist <- ggplot(df_prepped, aes(x = predicted_prob)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  labs(
    title = "Distribution of Predicted Values sigmoid",
    x = "Predicted Value (0-1)",
    y = "Count"
  ) +
  theme_bw() +
  tme



# Display both plots
print(dist_plot_hist)

## ---- 3. Validation Part ----
### ------ Looping through all layers for validation ------

# Initialize a list to store results for each layer
all_layers_validation <- list()
all_layers_categorized <- list()
all_layers_final_tables <- list()

# Get all layers for validation loop
all_layers_to_validate <- layer_map %>%
  pull(new_id)

# Loop through each layer as target_layer
for (target_layer_idx in seq_along(all_layers_to_validate)) {
  target_layer <- sub("layer_", "", all_layers_to_validate[target_layer_idx])
  target_layer_name <- all_layers_to_validate[target_layer_idx]
  
  cat(sprintf("\n\n===== PROCESSING LAYER: %s =====\n", target_layer_name))
  
  # Step 1: Get all other layers (excluding target_layer)
  other_layers <- layer_map %>%
    filter(new_id != target_layer_name) %>%
    pull(new_id)
  
  
  cat(sprintf("Target layer: %s\n", target_layer_name))
  cat(sprintf("Other layers: %s\n", paste(other_layers, collapse = ", ")))
  
  # Step 2: Aggregate edges from all other layers into one edge list
  other_layers_edges <- aggregated_df_weighted %>%
    filter(layer_from %in% other_layers) %>%
    # Aggregate across all other layers (same interaction may appear in multiple layers)
    group_by(node_from, node_to) %>%
    summarise(
      n_layers_observed = n(),                    # Number of layers where this interaction was observed
      total_weight = sum(weight),                 # Total weight across layers
      layers_list = paste(unique(layer_from), collapse = ", "),  # Which layers it appeared in
      .groups = "drop"
    ) %>%
    mutate(
      observed_in_other = 1,                      # Flag: observed in at least one other layer
      interaction_id = paste0(node_from, " -> ", node_to)
    )
  
  cat(sprintf("Number of unique interactions in other layers: %d\n", nrow(other_layers_edges)))
  cat(sprintf("Interactions appearing in multiple layers: %d\n", sum(other_layers_edges$n_layers_observed > 1)))
  
  # Step 3: Prepare target layer predictions for this layer
  df_self <- df_prepped %>%
    filter(train_layer == as.numeric(target_layer) & test_layer == as.numeric(target_layer)) %>%
    mutate(
      predicted_bin_sigm = if_else(predicted_prob > best_discrete_threshold, 1, 0),
      interaction_id = paste0(node_from, " -> ", node_to),
      original_binary = if_else(original_links > 0, 1, 0)  # Create original_binary here
    )
  
  if (nrow(df_self) == 0) {
    cat(sprintf("WARNING: No predictions found for layer %s\n", target_layer_name))
    next
  }
  
  # Step 4: Count observations per interaction across all layers (target + other)
  # First, get observation counts from the target layer
  target_layer_obs <- df_self %>%
    distinct(interaction_id, original_binary) %>%
    group_by(interaction_id) %>%
    summarise(obs_in_target = max(original_binary), .groups = "drop")  # 1 if observed in target, 0 otherwise
  
  # Combine with other layers observation info
  obs_counts <- target_layer_obs %>%
    left_join(
      other_layers_edges %>% select(interaction_id, n_layers_observed, observed_in_other),
      by = "interaction_id"
    ) %>%
    mutate(
      n_layers_observed = replace_na(n_layers_observed, 0),
      observed_in_other = replace_na(observed_in_other, 0),
      
      # Total observation count: target layer + other layers
      n_obs_total = obs_in_target + n_layers_observed
    )
  
  # Step 5: Flag interactions based on observation patterns
  df_self_flagged <- df_self %>%
    left_join(obs_counts, by = "interaction_id") %>%
    mutate(
      is_all_zero   = (n_obs_total == 0),                              # Never observed anywhere
      is_unique     = (n_obs_total == 1),                              # Observed in exactly 1 layer
      is_shared     = (n_obs_total >= 2),                              # Observed in 2+ layers
      obs_elsewhere = (observed_in_other == 1)                         # Observed in at least one other layer
    ) %>%
    select(target_layer = train_layer, interaction_id, itr, removed, original_binary, predicted_bin_sigm,
           is_all_zero, is_unique, is_shared, obs_elsewhere, n_obs_total)
  
  # Summary of flagged interactions
  cat("\n===== Interaction Observation Summary for", target_layer_name, "=====\n")
  cat(sprintf("Total unique interactions: %d\n", n_distinct(df_self_flagged$interaction_id)))
  cat(sprintf("All zero (never observed): %d\n", sum(obs_counts$is_all_zero <- obs_counts$n_obs_total == 0)))
  cat(sprintf("Unique (observed in 1 layer only): %d\n", sum(obs_counts$n_obs_total == 1)))
  cat(sprintf("Shared (observed in 2+ layers): %d\n", sum(obs_counts$n_obs_total >= 2)))
  cat(sprintf("Observed in other layers: %d\n", sum(obs_counts$observed_in_other == 1)))
  
  ### ---- Produce Results Table by Iteration ----
  final_table_by_iter <- df_self_flagged %>%
    group_by(itr) %>%
    summarise(
      TP = sum(original_binary == 1 & predicted_bin_sigm == 1, na.rm = TRUE),
      FP = sum(original_binary == 0 & predicted_bin_sigm == 1, na.rm = TRUE),
      TN = sum(original_binary == 0 & predicted_bin_sigm == 0, na.rm = TRUE),
      FN = sum(original_binary == 1 & predicted_bin_sigm == 0, na.rm = TRUE),
      
      # Unique interactions (observed in exactly 1 layer)
      locally_unique_links = sum(is_unique & original_binary == 1 & predicted_bin_sigm == 1, na.rm = TRUE),
      unsupported_links    = sum(is_unique & original_binary == 1 & predicted_bin_sigm == 0, na.rm = TRUE),
      
      # Shared interactions (observed in 2+ layers)
      recurrent_links      = sum(is_shared & original_binary == 1 & predicted_bin_sigm == 1, na.rm = TRUE),
      model_elusive_links        = sum(is_shared & original_binary == 1 & predicted_bin_sigm == 0, na.rm = TRUE),
      
      # Never observed anywhere
      possibly_forbidden_links     = sum(is_all_zero & original_binary == 0 & predicted_bin_sigm == 0, na.rm = TRUE),
      uncomfirmed_links       = sum(is_all_zero & original_binary == 0 & predicted_bin_sigm == 1, na.rm = TRUE),
      
      # Observed elsewhere but not in target layer
      locally_absent_links       = sum(original_binary == 0 & obs_elsewhere & predicted_bin_sigm == 0, na.rm = TRUE),
      possibly_missing_links = sum(original_binary == 0 & obs_elsewhere & predicted_bin_sigm == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = -itr,
      names_to = "Category", 
      values_to = "Count"
    ) %>%
    # Add target_layer AFTER pivot (so it doesn't get included in Count column)
    mutate(target_layer = target_layer_name) %>%
    arrange(itr, Category)
  
  ### ---- Summary table: sum all then turn to proportions ----
  final_table_summary <- final_table_by_iter %>%
    group_by(Category) %>%
    summarise(
      total = sum(Count, na.rm = TRUE),
      min   = min(Count, na.rm = TRUE),
      max   = max(Count, na.rm = TRUE),
      sd    = sd(Count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Category)
  
  cat("\nSummary for", target_layer_name, ":\n")
  print(final_table_summary)
  
  ### ---- Categorize each interaction ----
  df_categorized <- df_self_flagged %>%
    mutate(
      link_category = case_when(
        is_unique   & original_binary == 1 & predicted_bin_sigm == 1 ~ "locally_unique_links",
        is_unique   & original_binary == 1 & predicted_bin_sigm == 0 ~ "unsupported_links",
        is_shared   & original_binary == 1 & predicted_bin_sigm == 1 ~ "recurrent_links",
        is_shared   & original_binary == 1 & predicted_bin_sigm == 0 ~ "model_elusive_links",
        is_all_zero & original_binary == 0 & predicted_bin_sigm == 0 ~ "possibly_forbidden_links",
        is_all_zero & original_binary == 0 & predicted_bin_sigm == 1 ~ "uncomfirmed_links",
        obs_elsewhere & original_binary == 0 & predicted_bin_sigm == 0 ~ "locally_absent_links",
        obs_elsewhere & original_binary == 0 & predicted_bin_sigm == 1 ~ "possibly_missing_links",
        TRUE ~ "unclassified"
      )
    )
  
  # Check category distribution
  cat("\n===== Link Category Distribution for", target_layer_name, "=====\n")
  print(table(df_categorized$link_category))
  
  # Store results in list for later use
  all_layers_validation[[target_layer_name]] <- list(
    final_table_by_iter = final_table_by_iter,
    final_table_summary = final_table_summary,
    df_self_flagged = df_self_flagged
  )
  
  all_layers_categorized[[target_layer_name]] <- df_categorized
  all_layers_final_tables[[target_layer_name]] <- final_table_summary
  
} 



## --- 4. Aluvial plot for all layers ----

# Step 1: Combine all final tables from all layers
all_tables_combined <- bind_rows(all_layers_final_tables)


# Step 2: Calculate totals per category across all layers and all iterations
cat_totals_all_layers <- all_tables_combined %>%
  group_by(Category) %>%
  summarise(
    total = sum(total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tibble::deframe()

get_total <- function(nm) if (nm %in% names(cat_totals_all_layers)) unname(cat_totals_all_layers[[nm]]) else 0

# Step 3: Build flows for alluvial plot with all categories
flows <- tibble::tribble(
  ~L1, ~L2,                   ~L3,                     ~value,
  "TP","Validated elsewhere", "recurrent_links",        get_total("recurrent_links"),
  "TP","Not validated",       "locally_unique_links",   get_total("locally_unique_links"),
  "FN","Validated elsewhere", "model_elusive_links",    get_total("model_elusive_links"),
  "FN","Not validated",       "unsupported_links",      get_total("unsupported_links"),
  "FP","Validated elsewhere", "possibly_missing_links", get_total("possibly_missing_links"),
  "FP","Not validated",       "uncomfirmed_links",      get_total("uncomfirmed_links"),
  "TN","Validated elsewhere", "locally_absent_links",   get_total("locally_absent_links"),
  "TN","Not validated",       "possibly_forbidden_links", get_total("possibly_forbidden_links")
) %>%
  mutate(
    total_confusion = get_total("TP") + get_total("FP") + get_total("TN") + get_total("FN"),
    prop = ifelse(total_confusion > 0, value / total_confusion, 0),
    alluvium_id = paste(L1, L3, sep = "⟂")
  )

# Aesthetics
col_confusion <- c(
  "TP"="lightsteelblue","FP"="lightsteelblue2","TN"="rosybrown","FN"="rosybrown2"
)
col_validation <- c(
  "Validated elsewhere"="sandybrown","Not validated"="thistle3"
)
col_subcats <- c(
  "recurrent_links"="coral3","locally_unique_links"="thistle3",
  "model_elusive_links"="coral","unsupported_links"="thistle1",
  "possibly_missing_links"="coral2","uncomfirmed_links"="thistle",
  "locally_absent_links"="coral1","possibly_forbidden_links"="thistle2"
)
stratum_fill <- c(col_confusion, col_validation, col_subcats)

flow_alpha <- 0.7
flow_colour <- NA
stratum_label_size <- 4

# Choose whether to plot mean counts or proportions
metric <- "prop"  

# Custom orders for axes
order_confusion  <- c("TP", "FP", "TN", "FN")                 # LEFT column order
order_validation <- c("Not validated", "Validated elsewhere") # MIDDLE column order
order_subtypes   <- c(                                       # RIGHT column order
  "locally_unique_links", "uncomfirmed_links",
  "possibly_forbidden_links", "unsupported_links",
  "recurrent_links", "possibly_missing_links",
  "locally_absent_links", "model_elusive_links"
)

# Rebuild flows_long with custom orders applied
flows_long <- flows %>%
  select(L1, L2, L3, value, prop, alluvium_id) %>%
  tidyr::pivot_longer(c(L1, L2, L3), names_to = "axis", values_to = "stratum") %>%
  dplyr::mutate(
    axis = dplyr::recode(axis, L1 = "Confusion", L2 = "Validation", L3 = "Subtype"),
    # make stratum a factor with the order you chose, depending on axis
    stratum = dplyr::case_when(
      axis == "Confusion"  ~ factor(stratum, levels = order_confusion),
      axis == "Validation" ~ factor(stratum, levels = order_validation),
      axis == "Subtype"    ~ factor(stratum, levels = order_subtypes),
      TRUE ~ factor(stratum)
    ),
    # also lock the axis order (left -> middle -> right)
    axis = factor(axis, levels = c("Confusion", "Validation", "Subtype"))
  )

# Pretty labels for strata
stratum_labeller <- function(x) {
  x %>% str_replace_all("_", " ") %>% str_to_sentence()
}

# Background color for masking
bg_col <- "white"

# Build the alluvial plot
gg <- ggplot(
  flows_long,
  aes(x = axis,
      stratum = stratum,
      alluvium = alluvium_id,
      y = prop,
      fill = stratum)
) +
  # Wide "mask" strata: same width as flows, fill = background,
  # so they trim the flows exactly at the axis
  geom_stratum(
    width  = 0.02,
    color  = NA,
    fill   = bg_col,
    alpha  = 1
  ) +
  
  # Narrow visible strata on top
  geom_stratum(
    width = 0.03,
    color = "white"
  ) +
  # Flows first
  geom_alluvium(
    color = flow_colour,
    alpha = flow_alpha,
    width = 0.25,
    knot.pos = 0.2
  ) +
  
  scale_fill_manual(values = stratum_fill, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1)
  ) +
  labs(
    # title    = "Alluvial plot across all layers with all iterations",
    # subtitle = "Confusion classes → Validation group → Subcategories",
    x = NULL, y = "Proportion of interactions"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 12, face = "bold"),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(color = "grey30")
  )

# Totals per stratum (for % text)
stratum_totals <- flows_long %>%
  dplyr::distinct(axis, stratum, prop) %>%
  dplyr::group_by(axis, stratum) %>%
  dplyr::summarise(val = sum(prop), .groups = "drop") %>%
  dplyr::mutate(stratum_chr = as.character(stratum))

# Build one hidden stratum layer to grab box geometry (same width as above: 0.25)
tmp_build <- ggplot_build(
  ggplot(flows_long,
         aes(x = axis, stratum = stratum, alluvium = alluvium_id, y = prop)) +
    geom_stratum(width = 0.03, color = NA)
)

geo <- as.data.frame(tmp_build$data[[1]])
ax_levels <- levels(flows_long$axis)

# Midpoints and axis name for each stratum box
label_geom <- geo %>%
  dplyr::transmute(
    x_mid       = (xmin + xmax)/2,
    y_mid       = (ymin + ymax)/2,
    stratum_chr = as.character(stratum),
    axis        = ax_levels[pmax(1, pmin(length(ax_levels), round((xmin + xmax)/2)))]
  )

# Join values + geometry, build the two-line label strings
label_df <- dplyr::left_join(
  stratum_totals,
  label_geom,
  by = c("axis","stratum_chr")
) %>%
  dplyr::mutate(
    name_txt = stratum_chr %>% stringr::str_replace_all("_"," ") %>% stringr::str_to_sentence(),
    pct_txt  = scales::percent(val, accuracy = 1)
  )

# Tweakable label settings
label_nudge_x   <- 0.03                                 # push labels to the right of each box
label_nudge_y   <- 0.05  
y_off           <- 0.5 * diff(range(flows_long$prop, na.rm = TRUE))  # vertical gap for % line
name_size       <- 4.0
pct_size        <- 3.0
font_family     <- ""                                     # "" = default device font
label_color_map <- c("Confusion"="mistyrose4","Validation"="mistyrose4","Subtype"="mistyrose4")

label_df <- label_df %>%
  mutate(
    label_final = paste0(name_txt, " (", pct_txt, ")")
  )

gg <- gg +
  geom_text(
    data = label_df,
    inherit.aes = FALSE,
    aes(x = x_mid + label_nudge_x, y = y_mid, label = label_final, color = axis),
    fontface = "bold",
    size = name_size,
    family = font_family,
    hjust = 0
  ) +
  scale_color_manual(values = label_color_map, guide = "none") +
  coord_cartesian(clip = "off")

print(gg)

gg <- gg +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),     # remove tick labels
    axis.ticks = element_blank(),    # remove tick marks
    axis.title = element_blank(),    # remove axis titles
    axis.line = element_blank(),     # remove axis lines
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    plot.margin = margin(20, 20, 20, 20)  # keep right margin for labels
  )

gg
# Save the plot
ggsave(filename = paste0("alluvial_", metric, "_all_layers_Bartomeus_1_notitle.png"),
       plot = gg, width = 12, height = 8, dpi = 300)



## ---- 5. Species interactions heatmaps for each layer ----

### ---- Prepare data ----
all_layers_for_heatmap <- bind_rows(
  imap(all_layers_categorized, ~ {
    .x %>%
      filter(removed == 1) %>%  # ONLY include predictions with removed == 1
      mutate(target_layer = .y)
  })
) %>%
  separate(interaction_id, into = c("node_from", "node_to"), sep = " -> ", remove = FALSE) %>%
  # Collapse across iterations: keep first occurrence per plant-pollinator pair
  group_by(target_layer, node_from, node_to) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    link_display = recode(
      link_category,
      "locally_unique_links"   = "Locally unique",
      "unsupported_links"      = "Unsupported",
      "recurrent_links"        = "Recurrent",
      "model_elusive_links"    = "Model-elusive",
      "possibly_forbidden_links" = "Possibly forbidden",
      "uncomfirmed_links"      = "Uncomfirmed",
      "locally_absent_links"   = "Locally absent",
      "possibly_missing_links" = "Probably missing",
      "unclassified"           = "Unclassified"
    )
  )

### ---- Calculate degree per species for each layer ----
# Calculate plant degree (node_from) - number of unique pollinators each plant interacts with
plant_degree_by_layer <- all_layers_for_heatmap %>%
  filter(original_binary == 1) %>%  # Only count actual interactions
  group_by(target_layer, node_from) %>%
  summarise(plant_degree = n_distinct(node_to), .groups = "drop")

# Calculate pollinator degree (node_to) - number of unique plants each pollinator interacts with
poll_degree_by_layer <- all_layers_for_heatmap %>%
  filter(original_binary == 1) %>%  # Only count actual interactions
  group_by(target_layer, node_to) %>%
  summarise(poll_degree = n_distinct(node_from), .groups = "drop")

# Join degree information back to heatmap data
all_layers_for_heatmap <- all_layers_for_heatmap %>%
  left_join(plant_degree_by_layer, by = c("target_layer", "node_from")) %>%
  left_join(poll_degree_by_layer, by = c("target_layer", "node_to")) %>%
  mutate(
    plant_degree = replace_na(plant_degree, 0),
    poll_degree = replace_na(poll_degree, 0)
  )

# Define colors for link categories
link_colors_final <- c(
  "Locally unique"    = "#66C2A5",
  "Unsupported"       = "#FC8D62",
  "Recurrent"         = "#8DA0CB",
  "Model-elusive"     = "#E78AC3",
  "Possibly forbidden" = "#A6D854",
  "Uncomfirmed"       = "#FFD92F",
  "Locally absent"    = "#E5C494",
  "Probably missing"  = "#B3B3B3",
  "Unclassified"      = "#999999"
)

### ---- Create heatmap function with degree-based sorting ----
create_degree_sorted_heatmap <- function(data, layers_to_plot) {
  heatmap_data <- data %>%
    filter(target_layer %in% layers_to_plot)
  
  # For each layer, create order by degree (descending = higher degree first)
  for (layer in layers_to_plot) {
    layer_data <- heatmap_data %>% filter(target_layer == layer)
    
    # Get unique species and their degrees
    plant_order <- layer_data %>%
      distinct(node_from, plant_degree) %>%
      arrange(desc(plant_degree)) %>%
      pull(node_from)
    
    poll_order <- layer_data %>%
      distinct(node_to, poll_degree) %>%
      arrange(desc(poll_degree)) %>%
      pull(node_to)
    
    # Set factors for this layer (we'll need to handle this per-facet)
    heatmap_data <- heatmap_data %>%
      mutate(
        node_from = factor(node_from, levels = plant_order),
        node_to = factor(node_to, levels = poll_order)
      )
  }
  
  return(heatmap_data)
}

# Create faceted heatmap for layers 1-6
heatmap_layers_1_6 <- all_layers_for_heatmap %>%
  filter(target_layer %in% c("layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6"))

# Apply degree-based sorting per layer
heatmap_layers_1_6_sorted <- heatmap_layers_1_6 %>%
  nest(data = -target_layer) %>%
  mutate(
    data = map(data, ~ {
      layer_data <- .x
      plant_order <- layer_data %>%
        distinct(node_from, plant_degree) %>%
        arrange(desc(plant_degree)) %>%
        pull(node_from)
      poll_order <- layer_data %>%
        distinct(node_to, poll_degree) %>%
        arrange(desc(poll_degree)) %>%
        pull(node_to)
      
      layer_data %>%
        mutate(
          node_from = factor(node_from, levels = plant_order),
          node_to = factor(node_to, levels = poll_order)
        )
    })
  ) %>%
  unnest(data)

faceted_heatmap_1_6 <- ggplot(
  heatmap_layers_1_6_sorted,
  aes(x = node_to, y = node_from, fill = link_display)
) +
  geom_tile(color = "white", size = 0.3) +
  facet_wrap(~target_layer, scales = "free", ncol = 3) +
  scale_fill_manual(values = link_colors_final, name = "Link category", na.value = "white") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 8, face = "italic"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "right",
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Network Link Categories - Layers 1 to 6",
    subtitle = "Plants (rows, left→right: high→low degree) × Pollinators (columns, top→bottom: high→low degree)",
    x = "Pollinator (sorted by degree)",
    y = "Plant (sorted by degree)"
  ) +
  scale_y_discrete(labels = function(x) gsub("_", " ", substr(x, 1, 30))) +
  scale_x_discrete(labels = function(x) gsub("_", " ", substr(x, 1, 30)))

print(faceted_heatmap_1_6)

# # Save the heatmap for layers 1-6
# ggsave(
#   filename = "heatmap_layers_1_6_validation_Bartomeus_1.png",
#   plot = faceted_heatmap_1_6,
#   width = 16,
#   height = 10,
#   dpi = 300
# )


# Create faceted heatmap for layers 7-12
heatmap_layers_7_12 <- all_layers_for_heatmap %>%
  filter(target_layer %in% c("layer_7", "layer_8", "layer_9", "layer_10", "layer_11", "layer_12"))

# Apply degree-based sorting per layer
heatmap_layers_7_12_sorted <- heatmap_layers_7_12 %>%
  nest(data = -target_layer) %>%
  mutate(
    data = map(data, ~ {
      layer_data <- .x
      plant_order <- layer_data %>%
        distinct(node_from, plant_degree) %>%
        arrange(desc(plant_degree)) %>%
        pull(node_from)
      poll_order <- layer_data %>%
        distinct(node_to, poll_degree) %>%
        arrange(desc(poll_degree)) %>%
        pull(node_to)
      
      layer_data %>%
        mutate(
          node_from = factor(node_from, levels = plant_order),
          node_to = factor(node_to, levels = poll_order)
        )
    })
  ) %>%
  unnest(data)

faceted_heatmap_7_12 <- ggplot(
  heatmap_layers_7_12_sorted,
  aes(x = node_to, y = node_from, fill = link_display)
) +
  geom_tile(color = "white", size = 0.3) +
  facet_wrap(~target_layer, scales = "free", ncol = 3) +
  scale_fill_manual(values = link_colors_final, name = "Link category", na.value = "white") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 8, face = "italic"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "right",
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Network Link Categories - Layers 7 to 12",
    subtitle = "Plants (rows, left→right: high→low degree) × Pollinators (columns, top→bottom: high→low degree)",
    x = "Pollinator (sorted by degree)",
    y = "Plant (sorted by degree)"
  ) +
  scale_y_discrete(labels = function(x) gsub("_", " ", substr(x, 1, 30))) +
  scale_x_discrete(labels = function(x) gsub("_", " ", substr(x, 1, 30)))

print(faceted_heatmap_7_12)

# #Save the heatmap for layers 7-12
# ggsave(
#   filename = "heatmap_layers_7_12_validation_Bartomeus_1.png",
#   plot = faceted_heatmap_7_12,
#   width = 16,
#   height = 10,
#   dpi = 300
# )




