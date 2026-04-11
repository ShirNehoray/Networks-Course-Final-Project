# ---- check k and lambda influence on f0.5 ----
# here we perform a sensitivity analysis of predictive performance according to the maximal number of dimensions allowed in the SVD-based prediction and differend values of lambda.

library(tidyverse)
library(scales)

tme <-  theme(axis.text = element_text(size = 18, color = "black"),
              axis.title = element_text(size = 18, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

#sigmoid function
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# 0) set an array of thresholds
thresholds <- seq(0, 10, by = 1)
thresholds <- thresholds/10 

# set lambda to filter
las <- c(1, 5, 50, 100)

# read data
alll <- readRDS("predictions_1_Bartomeus_weighted_1.rds")
all_ks <- alll %>% filter(input_lambda %in% las) %>%
  mutate(input_lambda = factor(input_lambda, levels = las))

# convert negatives to zeros
df <- all_ks %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

# select the threshold for classifying links as 1s or 0s based on max F0.5
# 1) filter & prep
df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  )

# 2) expand to one row per threshold
df_thresh <- df_prepped %>%
  tidyr::expand_grid(threshold = thresholds) %>%  # <-- switch here
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by(train_layer, test_layer, k, itr, input_lambda, threshold) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    specificity      = TN / (TN + FP),
    precision        = TP / (TP + FP),
    recall           = TP / (TP + FN),
    f05_score   = (1.25) * (precision * recall) / ((0.25 * precision) + recall)
  ) %>%
  ungroup() 

# plot histograms of f05_score as a function of K and threshold
df_f05 <- df_thresh %>%
  select(k, itr,  input_lambda, threshold, f05_score)


default_threshold <- 0.6 #based on the threshold analysis

plot_df <- df_f05 %>%
  filter(threshold == default_threshold) %>%
  mutate(
    k = factor(k),
    input_lambda = factor(input_lambda)  # keeps facet order stable
  )

# --- 2) per-facet omnibus test across k (Kruskal–Wallis)
facet_p <- plot_df %>%
  group_by(input_lambda) %>%
  summarise(
    p = summary(aov(f05_score ~ k))[[1]][["Pr(>F)"]][1],
    .groups = "drop"
  ) %>%
  mutate(
    p_lab = paste0("p = ", scales::pvalue(p, accuracy = 0.001)),
    facet_lab = paste0("\u03BB = ", input_lambda, "\n", p_lab)
  )

# --- 3) build a named labeller for facet titles
lab_map <- setNames(facet_p$facet_lab, facet_p$input_lambda)

k_sensitivity_swap <- ggplot(plot_df, aes(x = k, y = f05_score, fill = k)) +
  geom_boxplot(color = "grey25", linewidth = 0.6, outlier_alpha = 0.35) +
  facet_wrap(
    ~ input_lambda,
    labeller = labeller(input_lambda = lab_map)
  ) +
  scale_fill_brewer(palette = "Pastel2", name = "k") +  # pastel fills + legend title
  labs(
    x = "Number of dimensions (k)",
    y = expression(F[0.5]~"score")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_text(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "grey35", linewidth = 0.6),
    strip.text = element_text(color = "grey10", face = "bold")  # not gray
  ) + tme


pdf(
  file   = "k_sensitivity_swap_Bartomeus_1.pdf",
  width  = 10,    # inches
  height = 8,
  family = "Helvetica"   # or another installed font
)
print(k_sensitivity_swap)
dev.off()     # close the file

# ---- Fig. S8: k sensitivity analysis ----
# overall k difference
pval <- kruskal.test(f05_score ~ k, data = plot_df)$p.value

p_text <- paste0("Kruskal–Wallis test: p = ",
                 formatC(pval, format = "e", digits = 2))

k_overall_sensitivity <- ggplot(plot_df, aes(x = k, y = f05_score, fill = k)) +
  geom_boxplot(color = "grey25", linewidth = 0.7, outlier_alpha = 0.35, notch = TRUE) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    x = "Number of dimensions (k)",
    y = expression(F[0.5]~"score"),
    title = p_text
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",          # ← removes legend
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "plain", size = 12),
    axis.title.x = element_text(face = "plain", size = 14),
    axis.title.y = element_text(face = "plain", size = 14)
  ) +
  tme

pdf(
  file   = "k_overall_sensitivity_Bartomeus_1.pdf",
  width  = 5,    # inches
  height = 4,
  family = "Helvetica"   # or another installed font
)
print(k_overall_sensitivity)
dev.off()     
