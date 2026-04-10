# Networks-Course-Final-Project
This repository contains the code and analysis for the final project in networks in biology and ecology course. The project examining the use of cross-network validation for evaluating interaction link prediction accuracy in ecological networks.

## Overview

Ecological interaction networks are inherently incomplete. This project applies a matrix completion approach (soft-impute SVD) to predict missing plant–pollinator interactions across 12 network layers from the same study system (Bartomeus, sourced from the [EuPPollNet database](http://dx.doi.org/10.1111/geb.70000)). Rather than evaluating predictions within each network, a cross-network validation framework is used to assess predictions against observations from other layers in the same system. All analysis was conducted in **R**. 

## Repository Structure

```

├── Bartomeus_1_full_code.R                     # Main analysis script
├── k_sensitivity_analysis.R                    # Sensitivity analysis for model parameters
├── Interaction_edges.csv                       # Network interaction data
├── predictions_1_Bartomeus_weighted_1.rds      # Model predictions output
└── plots_output_Bartomeus_1/                   # Results plots
    ├── k_overall_sensitivity_Bartomeus_1.pdf
    ├── k_sensitivity_swap_Bartomeus_1.pdf
    ├── optimal_threshold_Bartomeus_1.pdf
    ├── boxplot_confusion_matrix_Bartomeus_1.png
    ├── alluvial_prop_all_layers_Bartomeus_1.png
    ├── heatmap_layers_1_6_validation_Bartomeus_1.png
    └── heatmap_layers_7_12_validation_Bartomeus_1.png
```

## Pipeline

The script `Bartomeus_1_full_code.R` runs the full analysis in five sequential steps:

1. **Prediction** — For each of the 12 network layers, a bipartite interaction matrix is constructed, and 20% of interactions and an equal number of non-interactions are withheld per bootstrap iteration (n = 100). Missing links are predicted using soft-thresholded SVD via the softImpute package.

2. **Threshold optimisation** — Classification metrics (precision, recall, F0.5, balanced accuracy, MCC, specificity) are evaluated across probability thresholds to identify the optimal threshold.

3. **Within-network evaluation** — Model performance is assessed using confusion matrix components (TP, FP, TN, FN), averaged across bootstrap iterations for each layer.

4. **Cross-network validation** — Predicted interactions in each target layer are cross-referenced against observations from all other layers, and classified into eight ecologically interpretable categories (e.g., recurrent, locally unique, model-elusive, probably missing, possibly forbidden).

5. **Visualisation** — Results are visualised as threshold curves, boxplots and an alluvial plot.

## Key Parameters

| Parameter | Value | Description |
|---|---|---|
| `Study_id_chose` | `1_Bartomeus` | Study selected from EuPPollNet |
| `n_sim` | 100 | Number of bootstrap iterations |
| `prop_ones_to_remove` | 0.20 | Proportion of links withheld per iteration |
| `k` | 2 | Rank for SVD matrix completion |
| `lambda` | λ₀ | Regularisation parameter (upper bound) |
| `seed` | 123 | Random seed for reproducibility |

## Data

Input data (`Interaction_edges.csv`) were derived from 'Interaction_data.csv' sourced from the the **EuPPollNet** database (Lanuza et al., 2025), a European plant–pollinator interaction database containing over 1.1 million interactions across 1,864 networks. The data can be accessed at:

> Lanuza et al. (2025). EuPPollNet: A European database of plant–pollinator networks. *Global Ecology and Biogeography*.  https://doi.org/10.1111/geb.70000


