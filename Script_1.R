
rm(list = ls()) # Clear all
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## Libraries
# Load required packages

library(openxlsx)
library(ggpubr)
library(gridExtra)
library(analogue)
library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(corrplot)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(compositions)
library(caret)
library(ggExtra)
library(randomForest)  
library(reshape2)
library(truncnorm) 
library(rcarbon)
library(tidyr)
library(mgcv)
library(sf)
library(blockCV)
library(spdep)
library(grid)

## 1) SUMMARY STATISTICS & EXPLORATORY PLOTS ####


# 1. --- Read the data ---
data <- read.xlsx("Dataset_1.xlsx", sheet = "Dataset_R")


# 2. --- Define compound groups ---
tetramethylated <- c("fIa", "fIb", "fIc")
pentamethylated <- c("fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_")
hexamethylated  <- c("fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_")

# 3. --- Sum groups ---
data$Tetramethylated <- rowSums(data[, tetramethylated], na.rm = TRUE)
data$Pentamethylated <- rowSums(data[, pentamethylated], na.rm = TRUE)
data$Hexamethylated  <- rowSums(data[, hexamethylated], na.rm = TRUE)

# 4. --- Compute CBT ---
# Using only fIIa-c and fIIIa-c 
data$CBT <- log10((data$fIc + data$fIIa_+data$fIIb_+data$fIIc_+data$fIIIa_+data$fIIIb_+data$fIIIc_) / 
                    (data$fIa + data$fIIa + data$fIIIa))

# 5. --- Compute MBT (sum of all Ia-IIIa compounds) ---
all_compounds <- c("fIa", "fIb", "fIc",
                   "fIIa", "fIIb", "fIIc",
                   "fIIIa")

data$MBT <- (data$fIa + data$fIb + data$fIc) / rowSums(data[, all_compounds], na.rm = TRUE)

# 6. --- Select variables for Spearman correlations, including CBT and MBT ---
selected_vars <- c("NPP","Tetramethylated", "Pentamethylated", "Hexamethylated", "MBT" , "CBT",
                   "PET","BIO01", "BIO12", "pH", 
                    "Bulk", "Clay", "CEC", "Nitrogen", "Phosphorus", "SOC_0-5", "SOC_5-15",
                   "Elevation" )


# 7. --- Subset numeric ---
cor_data <- data[, selected_vars]
cor_data <- as.data.frame(lapply(cor_data, as.numeric))

# 8. --- Compute correlation matrix ---
cor_matrix <- cor(cor_data, use = "complete.obs", method = "spearman")

# Change variable names for better visualization
colnames(cor_matrix)[colnames(cor_matrix) == "NPP"] <- "NPP"
colnames(cor_matrix)[colnames(cor_matrix) == "BIO01"] <- "MAT"
colnames(cor_matrix)[colnames(cor_matrix) == "BIO12"] <- "MAP"

rownames(cor_matrix) <- colnames(cor_matrix)  

# --- Plot ---

corrplot(cor_matrix,
         method = "circle",
         type = "upper",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         number.cex = 0.8,
         cl.pos = "r")


# Plot specific individual compounds and NPP

# --- Define brGDGT compound groups (in case it was not done before) ---
tetramethylated <- c("fIa", "fIb", "fIc")
pentamethylated <- c("fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_")
hexamethylated  <- c("fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_")

# --- Select columns by compound names
tet_cols  <- intersect(tetramethylated, names(data))
penta_cols <- intersect(pentamethylated, names(data))
hexa_cols <- intersect(hexamethylated, names(data))

# --- Ensure numeric type ---
for(col in c(tet_cols, penta_cols, hexa_cols)){
  data[[col]] <- as.numeric(as.character(data[[col]]))
}

# --- Create summed groups ---
data2 <- data %>%
  mutate(
    tetramethylated = rowSums(across(all_of(tet_cols)), na.rm = TRUE),
    pentamethylated = rowSums(across(all_of(penta_cols)), na.rm = TRUE),
    hexamethylated  = rowSums(across(all_of(hexa_cols)), na.rm = TRUE)
  )

# --- Pivot to long format for plotting ---
all_cols <- c("tetramethylated", tetramethylated,
              "pentamethylated", pentamethylated,
              "hexamethylated", hexamethylated)

data_long <- data2 %>%
  pivot_longer(cols = all_of(all_cols), names_to = "Fraction", values_to = "Value") %>%
  mutate(
    Group = case_when(
      Fraction %in% c("tetramethylated", tetramethylated) ~ "Tetra",
      Fraction %in% c("pentamethylated", pentamethylated) ~ "Penta",
      Fraction %in% c("hexamethylated", hexamethylated)  ~ "Hexa"
    ),
    is_sum = Fraction %in% c("tetramethylated","pentamethylated","hexamethylated"),
    Fraction_plot = gsub("_", "'", Fraction)
  )

# --- Make sure Sampletype is a factor ---
data_long$Sampletype <- factor(data_long$Sampletype)

# --- Color palette for sample types ---
sample_colors <- brewer.pal(n = length(levels(data_long$Sampletype)), "Set1")

# --- Function to compute  statistics obtained from Spearman test, linear, quadrating and GAM models ---

compute_all_stats <- function(datf) {
  if(nrow(datf) < 3) 
    return(list(
      n = nrow(datf),
      spearman_rho = NA, spearman_p = NA,
      linear_r2 = NA, linear_p = NA,
      quad_r2 = NA, quad_p = NA,
      gam_dev_expl = NA, gam_p = NA, gam_edf = NA
    ))
  
  
  # Initialise with NA values
  results <- list(
    n = nrow(datf),
    spearman_rho = NA, spearman_p = NA,
    linear_r2 = NA, linear_p = NA,
    quad_r2 = NA, quad_p = NA,
    gam_dev_expl = NA, gam_p = NA, gam_edf = NA
  )
  
  # 1. Spearman correlation
  tryCatch({
    sp <- suppressWarnings(cor.test(datf$NPP, datf$Value, method = "spearman"))
    results$spearman_rho <- sp$estimate
    results$spearman_p <- sp$p.value
  }, error = function(e) {})
  
  # 2. Linear model
  tryCatch({
    lm_fit <- lm(Value ~ NPP, data = datf)
    lm_summary <- summary(lm_fit)
    results$linear_r2 <- lm_summary$r.squared
    results$linear_p <- ifelse(nrow(coef(lm_summary)) >= 2, 
                               coef(lm_summary)[2, 4], NA)
  }, error = function(e) {})
  
  # 3. Quadratic model (U shape)
  tryCatch({
    quad_fit <- lm(Value ~ NPP + I(NPP^2), data = datf)
    quad_summary <- summary(quad_fit)
    results$quad_r2 <- quad_summary$r.squared
    results$quad_p <- ifelse(nrow(coef(quad_summary)) >= 3, 
                             coef(quad_summary)[3, 4], NA)
  }, error = function(e) {})
  
  # 4. GAM model
  tryCatch({
    k_val <- min(8, nrow(datf) - 1)
    if(k_val >= 3) {
      gam_fit <- gam(Value ~ s(NPP, k = k_val), data = datf, method = "REML") # Checked with other options in 'method', main results do not change.
      gam_summary <- summary(gam_fit)
      results$gam_dev_expl <- gam_summary$dev.expl
      results$gam_p <- ifelse(!is.null(gam_summary$s.table) && nrow(gam_summary$s.table) > 0,
                              gam_summary$s.table[1, "p-value"], NA)
      results$gam_edf <- ifelse(!is.null(gam_summary$s.table) && nrow(gam_summary$s.table) > 0,
                                gam_summary$s.table[1, "edf"], NA)
    }
  }, error = function(e) {})
  
  return(results)
}

# --- Compute statistics for all fractions ---
fractions <- unique(data_long$Fraction)
stats_all <- list()

for(frac in fractions) {
  datf <- data_long %>% 
    filter(Fraction == frac) %>%
    drop_na(NPP, Value)
  
  stats_all[[frac]] <- compute_all_stats(datf)
}

# --- Function to format p-values ---
format_p <- function(p) {
  if(is.na(p)) return("NA")
  if(p < 0.001) return("<0.001")
  return(sprintf("%.3f", p))
}



# --- Figure 1. Relationship between NPP and specific brGDGT compounds and correlation ---


# Function:
make_comparison_plot <- function(frac) {
  
  # Filter data for the selected fraction and remove missing values
  datf <- data_long %>% 
    filter(Fraction == frac) %>%
    drop_na(NPP, Value)
  
  stats <- stats_all[[frac]]
  
  # Check that statistics exist
  if (is.null(stats) || length(stats) == 0) {
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "No data") +
      theme_void()
    return(p)
  }
  
  # Base scatter plot
  p <- ggplot(datf, aes(x = NPP, y = Value, color = Sampletype)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(
      title = gsub("_", "'", frac),
      x = expression("NPP (g C m"^-2*" yr"^-1*")"),
      y = "Percentage (%)"
    ) +
    scale_color_manual(values = sample_colors, drop = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
  
  pred_data <- datf %>% 
    arrange(NPP)
  
  # 1. Linear model
  if (!is.na(stats$linear_r2)) {
    tryCatch({
      lm_fit <- lm(Value ~ NPP, data = datf)
      pred_data$linear <- predict(lm_fit, newdata = pred_data)
      
      p <- p + geom_line(
        data = pred_data,
        aes(x = NPP, y = linear),
        color = "black",
        size = 1.5,
        linetype = "dashed",
        alpha = 0.9
      )
    }, error = function(e) {})
  }
  
  # 2. Quadratic model
  if (!is.na(stats$quad_r2)) {
    tryCatch({
      quad_fit <- lm(Value ~ NPP + I(NPP^2), data = datf)
      pred_data$quadratic <- predict(quad_fit, newdata = pred_data)
      
      p <- p + geom_line(
        data = pred_data,
        aes(x = NPP, y = quadratic),
        color = "black",
        size = 1.5,
        linetype = "dotted"
      )
    }, error = function(e) {})
  }
  
  # 3. GAM model 
  if (!is.na(stats$gam_dev_expl) && stats$n >= 8) {
    tryCatch({
      k_val <- min(8, nrow(datf) - 1)
      
      if (k_val >= 3) {
        gam_fit <- gam(Value ~ s(NPP, k = k_val),
                       data = datf,
                       method = "REML") # Main results do not change when other options are used in "method".
        
        pred_data$gam <- predict(gam_fit, newdata = pred_data)
        
        p <- p + geom_line(
          data = pred_data,
          aes(x = NPP, y = gam),
          color = "black",
          size = 1.5,
          linetype = "solid",
          alpha = 0.6
        )
      }
    }, error = function(e) {})
  }
  
  # Build statistics annotation text // Just to show them on plots
  stats_text <- paste0("n = ", stats$n)
  
  if (!is.na(stats$spearman_rho)) {
    stats_text <- paste0(
      stats_text, "\n",
      "Spearman: ρ = ", sprintf("%.2f", stats$spearman_rho),
      ", p = ", format_p(stats$spearman_p)
    )
  }
  
  if (!is.na(stats$linear_r2)) {
    stats_text <- paste0(
      stats_text, "\n",
      "Linear: R² = ", sprintf("%.2f", stats$linear_r2),
      ", p = ", format_p(stats$linear_p)
    )
  }
  
  if (!is.na(stats$quad_r2)) {
    stats_text <- paste0(
      stats_text, "\n",
      "Quadratic: R² = ", sprintf("%.2f", stats$quad_r2),
      ", p = ", format_p(stats$quad_p)
    )
  }
  
  if (!is.na(stats$gam_dev_expl)) {
    stats_text <- paste0(
      stats_text, "\n",
      "GAM: R² = ", sprintf("%.2f", stats$gam_dev_expl),
      ", p = ", format_p(stats$gam_p),
      ", EDF = ",
      ifelse(!is.na(stats$gam_edf),
             sprintf("%.1f", stats$gam_edf),
             "NA")
    )
  }
  
  return(p)
}



# --- Build plots ---
plots_tetra <- lapply(c("tetramethylated", tetramethylated), make_comparison_plot)
plots_penta <- lapply(c("pentamethylated", pentamethylated), make_comparison_plot)
plots_hexa  <- lapply(c("hexamethylated", hexamethylated), make_comparison_plot)


# --- Optional: Save to file w --- 
# ggsave(7.png", 
 #      plot = plots_tetra[[7]], 
  #     width = 3, height = 3, dpi = 300)


# --- Create a summary table ---

summary_table <- data.frame(
  Fraction = fractions,
  n = sapply(stats_all, function(x) ifelse(is.null(x$n), NA, x$n)),
  Spearman_rho = sapply(stats_all, function(x) ifelse(is.null(x$spearman_rho), NA, x$spearman_rho)),
  Spearman_p = sapply(stats_all, function(x) ifelse(is.null(x$spearman_p), NA, x$spearman_p)),
  Linear_R2 = sapply(stats_all, function(x) ifelse(is.null(x$linear_r2), NA, x$linear_r2)),
  Linear_p = sapply(stats_all, function(x) ifelse(is.null(x$linear_p), NA, x$linear_p)),
  Quadratic_R2 = sapply(stats_all, function(x) ifelse(is.null(x$quad_r2), NA, x$quad_r2)),
  GAM_R2 = sapply(stats_all, function(x) ifelse(is.null(x$gam_dev_expl), NA, x$gam_dev_expl))
)


print(summary_table)



#### 2) RANDOM FOREST MODEL WITH NORMALISED PREDICTORS ####

# --- Load data ---
data <- read.xlsx("Dataset_1.xlsx", sheet = "Dataset_R")
#Check
head(data)

# --- Select compounds ---
vars <- colnames(data)[11:25]

# ---  and include sample type, NPP & coordinates --- 
vars_extra <- c("Sampletype", "NPP", "Latitude", "Longitude")

# --- Select columns with specific compounds ---
comp_data <- data[, vars]  

# --- Adjust percentage in a 0-1 range ---
comp_data <- comp_data / 100  

# --- Avoid 0 values for the CLR transformation ---
comp_data <- comp_data + 1e-6  

# --- CLR transformation to avoid multicollinearity ---
clr_data <- clr(comp_data)  

# --- Convert into a dataframe and rename cols ---
clr_df <- as.data.frame(clr_data)
colnames(clr_df) <- paste0("clr_", colnames(clr_df))

# --- Merge ---
model_data <- data.frame(clr_data, data[, vars_extra])

# --- Remove missing values ---
model_data <- model_data[complete.cases(model_data), ]

# --- Combine aquatic_SPM categories SPM :

# It is important to bear in mind that there are only 16 samples where 
# Sample type == "Lacustrine Meso/Microcosm", so
# the training dataset does not include this variable in the testing set
# because of the random partitioning of the data. This can lead to issues when
# making predictions with the model, as it encounters a sample type during testing
# that it has not seen during training. This results in an error because the model
# cannot assign a prediction to an unseen category. To address this, we need to ensure
# that all sample types present in the testing set are also represented in the training set.
# Therefore, we collapse into the same cateogiry
# all aquatic SPM measurements. To differentiate it from the original
# sample types used in the dataset, this is named "Sampletype_fixed":

model_data <- model_data %>%
  mutate(
    Sampletype_fixed = case_when(
      # All aquatic SPM small (<200 samples) categories are combined under the 
      # new ""Aquatic_SPM_combined" category.
      Sampletype %in% c("Lacustrine SPM", 
                        "Low DO Lacustrine SPM",
                        "Riverine Sediment and SPM",
                        "Lacustrine Meso/Microcosm") ~ "Aquatic_SPM_combined",
      TRUE ~ as.character(Sampletype)
    ),
    Sampletype_fixed = as.factor(Sampletype_fixed)
  )


# --- Formula ---
formula <- as.formula(paste("NPP ~", paste(c(colnames(clr_data), "Sampletype_fixed"), collapse = " + ")))
# Check
formula


# --- RF model with all data --- ####
# All data is used for training and for testing

model_all_data <- train(
  formula,
  data = model_data,
  method = "rf",
  trControl = trainControl(method = "none"),  # No CV
  importance = TRUE,
  ntree = 1000,
  tuneGrid = data.frame(mtry = 11)  
)

print(model_all_data)
# Compute r^2
predictions_all <- predict(model_all_data, model_data)
r2_all <- cor(model_data$NPP, predictions_all)^2
r2_all
#Save model
saveRDS(model_all_data, file = "model_brGDGT-NPP_all1.rds")
### --- Load model ---
model_all_data <- readRDS("model_brGDGT-NPP_all1.rds")


# --- RF model with random cross validation --- ####
# 5-fold random cross validation, each time, ca. 20% of the
# dataset is excluded from the training as used for testing
# model predictions

# --- Reproducibility ---
set.seed(2026) 

# --- Cross-validation ---
train_control <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final",
  returnResamp = "final"
)

# --- RF model with random cross-validation ---
model_cv <- train(
  formula,
  data = model_data,
  method = "rf",
  trControl = train_control,
  importance = TRUE,
  ntree = 1000, #500 # Model was re-run with ntrees=500 and ntrees=1000; outcomes are identical
  tuneLength = 3
)

print(model_cv)


### --- Save model ---
  # saveRDS(model_cv, file = "model_brGDGT-NPP_cv1.rds")

### --- Load model ---
model_cv <- readRDS("model_brGDGT-NPP_cv1.rds")


# --- RF model with nested spatial cross validation --- ####
# 5-fold nested spatial cross-validation:
# The study area is partitioned into 5 spatial block types (Block 1, 2, 3, 4, 5), 
# which are randomly distributed across the map. 
# For each fold, all samples located within a given block type (e.g., Block 1) 
# are excluded from model training and used exclusively for testing. 
# This procedure is repeated for each of the five block types.
# The entire analysis is conducted four times, using spatial blocks of 
# 50 km, 100 km, 500 km, and 1000 km to assess the effect of block size on model performance.


nested_spatial_cv_complete <- function(model_data_sf, block_size_km, k_outer = 5, k_inner = 5, 
                                       save_models = FALSE, model_name = NULL) {
  # Empty data frame to store results from each model
  results <- data.frame()
  model_df <- st_drop_geometry(model_data_sf) # transforms model_data_sf into dataframe
  
  # List to save results
  all_fold_predictions <- list()
  all_final_models <- list()  # Save model
  
  # Blocks
  outer_blocks <- cv_spatial(model_data_sf, k = k_outer, size = block_size_km)
  
  for(outer_fold in 1:k_outer) {
    cat(sprintf("\n--- Outer Fold %d/%d ---\n", outer_fold, k_outer))
    
    # Train/test indices
    train_idx <- outer_blocks$folds_list[[outer_fold]][[1]]
    test_idx <- outer_blocks$folds_list[[outer_fold]][[2]]
    
    train_df <- model_df[train_idx, ]
    test_df <- model_df[test_idx, ]
    
    # inner cv for nested cv
    train_sf <- model_data_sf[train_idx, ]
    inner_blocks <- cv_spatial(train_sf, k = k_inner, size = block_size_km)
    
    inner_train_ids <- lapply(inner_blocks$folds_list, function(x) x[[1]])
    inner_test_ids <- lapply(inner_blocks$folds_list, function(x) x[[2]])
    
    tr_control_inner <- trainControl(
      method = "cv",
      index = inner_train_ids,
      indexOut = inner_test_ids
    )
    
    # Train // Bear in mind that "Sampletype_fixed" can be removed from the predictive model
    inner_model <- train(
      NPP ~ Sampletype_fixed + fIa + fIb + fIc + fIIa + fIIa_ + fIIb + fIIb_ + 
        fIIc + fIIc_ + fIIIa + fIIIa_ + fIIIb + fIIIb_ + fIIIc + fIIIc_,
      data = train_df,
      method = "rf",
      trControl = tr_control_inner,
      tuneGrid = expand.grid(mtry = c(3, 7, 11, 15)),
      ntree = 1000 
    )
    
    # Final// Check whether you want to include "Sampletype_fixed" or not.
    final_model <- train(
      NPP ~ Sampletype_fixed + fIa + fIb + fIc + fIIa + fIIa_ + fIIb + fIIb_ + 
        fIIc + fIIc_ + fIIIa + fIIIa_ + fIIIb + fIIIb_ + fIIIc + fIIIc_,
      data = train_df,
      method = "rf",
      tuneGrid = inner_model$bestTune,
      trControl = trainControl(method = "none"),
      ntree = 1000 # Checked with 500 trees, results are identical
    )
    
    # Save model
    if(save_models) {
      all_final_models[[outer_fold]] <- final_model
    }
    
    # Test model performance
    predictions <- predict(final_model, test_df)
    residuals <- test_df$NPP - predictions  
    r = cor(test_df$NPP, predictions)
    r2 <- cor(test_df$NPP, predictions)^2
    rmse <- sqrt(mean((test_df$NPP - predictions)^2))
    mae <- mean(abs(test_df$NPP - predictions))
    
    # Save predictions and residuals for this fold 
    all_fold_predictions[[outer_fold]] <- data.frame(
      Observed = test_df$NPP,
      Predicted = predictions,
      Fold = outer_fold,
      BlockSize_km = block_size_km/1000,
      Residuals = residuals,  
      Sampletype_fixed = test_df$Sampletype_fixed,
      X = st_coordinates(model_data_sf[test_idx, ])[, 1],
      Y = st_coordinates(model_data_sf[test_idx, ])[, 2]
    )
    
    # Store results for this fold
    results <- rbind(results, data.frame(
      outer_fold = outer_fold,
      r = r,
      r2 = r2,
      rmse = rmse,
      mae = mae,
      best_mtry = inner_model$bestTune$mtry,
      n_train = nrow(train_df),
      n_test = nrow(test_df)
    ))
    
    cat(sprintf("  r² = %.3f, RMSE = %.1f, MAE = %.1f, mtry = %d\n", 
                r2, rmse, mae, inner_model$bestTune$mtry))
  }
  
  # Summary statistics across folds
  cat(sprintf("\n=== Summary: r² = %.3f ± %.3f ===\n",
              mean(results$r2), sd(results$r2)))
  
  # Merge all predictions into a single data frame
  all_predictions_df <- do.call(rbind, all_fold_predictions)
  
  # Create list with result and predictions
  return_list <- list(
    results = results,
    predictions = all_predictions_df
  )
  
  # Add models if saved and save to file
  if(save_models && length(all_final_models) > 0) {
    return_list$models <- all_final_models
    
    # Save
    if(!is.null(model_name)) {
      saveRDS(return_list, file = paste0(model_name, ".rds"))
      cat(sprintf("Models saved: %s.rds\n", model_name))
    }
  }
  
  return(return_list)
}



# RUN FUNCTION FOR EACH SPATIAL BLOCK SIZE
model_data_sf_fixed <- st_as_sf(
  model_data,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)
result_caret_50km_complete <- nested_spatial_cv_complete(
  model_data_sf = model_data_sf_fixed,
  block_size_km = 50000, # meters
  save_models = TRUE,
  model_name = "spatial_cv_50km_models1"
)


result_caret_100km_complete <- nested_spatial_cv_complete(
  model_data_sf = model_data_sf_fixed,
  block_size_km = 100000, # meters
  save_models = TRUE,
  model_name = "spatial_cv_100km_models1"
)

result_caret_500km_complete <- nested_spatial_cv_complete(
  model_data_sf = model_data_sf_fixed,
  block_size_km = 500000, # meters
  save_models = TRUE,
  model_name = "spatial_cv_500km_models1"
)

result_caret_1000km_complete <- nested_spatial_cv_complete(
  model_data_sf = model_data_sf_fixed,
  block_size_km = 1000000, # meters
  save_models = TRUE,
  model_name = "spatial_cv_1000km_models1"
)



# LOAD MODELS FROM FILES
load_spatial_cv_models <- function(block_size_km) {
  file_name <- paste0("spatial_cv_", block_size_km, "km_models1.rds")
  if (file.exists(file_name)) {
    cat(sprintf("Loading %d km models...\n", block_size_km))
    return(readRDS(file_name))
  } else {
    cat(sprintf("File %s not found\n", file_name))
    return(NULL)
  }
}


model_spatial_cv_50<-load_spatial_cv_models(50) # Here the Block size is in kilometers
model_spatial_cv_100<-load_spatial_cv_models(100) # Here the Block size is in kilometers
model_spatial_cv_500<-load_spatial_cv_models(500) # Here the Block size is in kilometers
model_spatial_cv_1000<-load_spatial_cv_models(1000) # Here the Block size is in kilometers



# Get predictions for best mtry value from random CV model
best_preds <- model_cv$pred[model_cv$pred$mtry == model_cv$bestTune$mtry, ]

# Check performance metrics for each block size
# First, specify the block size (model_spatial_cv_50, model_spatial_cv_100, 
# model_spatial_cv_500 or model_spatial_cv_1000):

bs<- model_spatial_cv_1000

results_by_type <- bs$predictions %>%
  bind_rows(
    model_spatial_cv_50$predictions %>% 
      mutate(Sampletype_fixed = "All Samples")
  ) %>%
  group_by(Sampletype_fixed) %>%
  summarise(
    n = n(),
    r = cor(Observed, Predicted),
    r2 = r^2,
    rmse = sqrt(mean((Observed - Predicted)^2))
  ) %>%
  arrange(desc(n))

print(results_by_type)



# --- Plot RF models' validation--- ####

# For the model trained with all data:

obs_in <- model_data$NPP
pred_in <- predict(model_all_data, model_data) # model trained with all data
r2_in <- round((cor(obs_in, pred_in)^2),2)


#For the models trained through cross-validation:

# For random CV:
cv_predictions <- model_cv$pred
obs_cv <- cv_predictions$obs
pred_cv <- cv_predictions$pred
r2_cv <- round( (cor(obs_cv, pred_cv)^2), 2)



# For nested saptial cross-validation it is included withing the plot function (Figure 3):

plot_spatial_cv <- function(result_obj, block_name, color) {
  if(!is.null(result_obj$predictions)) {
    p <- ggplot(result_obj$predictions, aes(x = Observed, y = Predicted)) +
      geom_point(color = color, alpha = 0.4, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, 
                  color = color, fill = adjustcolor(color, alpha.f = 0.1)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkgray") +
      labs(title = paste("Spatial CV:", block_name),
           subtitle = sprintf("R² = %.3f ± %.3f", 
                              mean( round((result_obj$results$r2), 2)  ),
                              sd(result_obj$results$r2)),
           x = "Observed NPP", y = "Predicted NPP") +
      theme_minimal()
    return(p)
  }
  return(NULL)
}

p50 <- plot_spatial_cv(model_spatial_cv_50, "50 km", "orange")
p100 <- plot_spatial_cv(model_spatial_cv_100, "100 km", "darkgreen")
p500 <- plot_spatial_cv(model_spatial_cv_500, "500 km", "purple")
p1000 <- plot_spatial_cv(model_spatial_cv_1000, "1000 km", "brown")
p50
p100
p500
p1000

# Same plot for the model trained with all data and with random CV

# First, data prepared as previously done with the training dataset:
all_fracs <- colnames(model_data)[grepl("^clr_", colnames(model_data))] 
newdata_in <- model_data  # aldready contains CRL+ Sampletype
colnames(model_data)

obs_in <- newdata_in$NPP
pred_in <- predict(model_cv, newdata = newdata_in)
r2_in <- round((1 - sum((obs_in - pred_in)^2) / sum((obs_in - mean(obs_in))^2)), 2) 

# --- Cross-validated predictions ---
cv_pred <- model_cv$pred
cv_pred <- subset(cv_pred, mtry == model_cv$bestTune$mtry) 
obs_cv <- cv_pred$obs

pred_cv <- cv_pred$pred
r2_cv <- 1 - sum((obs_cv - pred_cv)^2) / sum((obs_cv - mean(obs_cv))^2)

# --- Dataframes for plot ---
df_cv <- data.frame(Observed = obs_cv, Predicted = pred_cv, Type = "Cross-validated")
df_in <- data.frame(Observed = obs_in, Predicted = pred_in, Type = "In-sample")




p1 <- ggplot(df_in, aes(x = Observed, y = Predicted)) +
  geom_point(color = "darkblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkblue", fill = "lightblue", alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "Random Forest: In-sample Predictions",
       subtitle = "Model trained and tested on all data",
       x = "Observed NPP",
       y = "Predicted NPP") +
  annotate("text", x = min(obs_in), y = max(pred_in), 
           label = paste("R² =", round(r2_in, 2)), 
           hjust = 0, vjust = 1.5, color = "black", size = 5, fontface = "bold") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"))

print(p1)

# PLOT 2: PREDICTIONS OBTAINED FROM MODEL WITH RANDOM CV 
p2 <- ggplot(df_cv, aes(x = Observed, y = Predicted)) +
  geom_point(color = "darkred", alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkred", fill = "pink", alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "Random Forest: Random Cross-Validation",
       subtitle = "5-fold CV assuming independent samples",
       x = "Observed NPP",
       y = "Predicted NPP") +
  annotate("text", x = min(obs_cv), y = max(pred_cv), 
           label = paste("R² =", round(r2_cv, 2)), 
           hjust = 0, vjust = 1.5, color = "black", size = 5, fontface = "bold") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"))

print(p2)



# --- Plot RF prediction residuals--- ####

# 1. Function to compute the 95% and 75% of the predictions range 
calc_auto_tolerance <- function(residuals, percentiles = c(95, 75)) {
  # Residual values
  abs_res <- abs(residuals)
  # Range for each percentil
  thresholds <- quantile(abs_res, probs = percentiles/100)
  
  return(thresholds)
}

# 2. Function to plot ranges
auto_residual_plot <- function(predictions_df, title, color, 
                               percentiles = c(95, 75)) {
  # Compute residuals
  if(!"Residuals" %in% names(predictions_df)) {
    predictions_df$Residuals <- predictions_df$Observed - predictions_df$Predicted
  }
  
  # Ranges
  tolerances <- calc_auto_tolerance(predictions_df$Residuals, percentiles)
  tol_95 <- tolerances[1]
  tol_75 <- tolerances[2]
  
  # % Within each range
  pct_95 <- mean(abs(predictions_df$Residuals) <= tol_95) * 100
  pct_75 <- mean(abs(predictions_df$Residuals) <= tol_75) * 100
  
  # Plot
  p <- ggplot(predictions_df, aes(x = Predicted, y = Residuals)) +
    # Band 95%
    annotate("rect", xmin = -Inf, xmax = Inf, 
             ymin = -tol_95, ymax = tol_95,
             fill = "lightgreen", alpha = 0.05) +
    # Band 75%
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -tol_75, ymax = tol_75,
             fill = "green", alpha = 0.1) +
    # Points
    geom_point(color = color, alpha = 0.4, size = 1.2) +
    # Lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = c(-tol_95, tol_95), 
               linetype = "dotted", color = "darkgreen", alpha = 0.5) +
    geom_hline(yintercept = c(-tol_75, tol_75),
               linetype = "dashed", color = "darkgreen", alpha = 0.7) +
    # Labels 
    labs(title = paste("Residuals:", title),
         subtitle = sprintf("%.0f%% within ±%.0f | %.0f%% within ±%.0f", 
                            pct_75, tol_75, pct_95, tol_95),
         x = "Predicted NPP", y = "Residuals (Obs - Pred)") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8))+
    ylim(-500, 500)+
    scale_x_continuous(limits = c(300, 1300), breaks = seq(300, 1200, 300)) 
  
  p_with_marginal <- ggMarginal(p, type = "histogram", margins = "y", 
                                groupColour = FALSE, groupFill = FALSE,
                                fill = color, alpha = 0.6)
  
  return(list(plot = p_with_marginal, tol_95 = tol_95, tol_75 = tol_75))
}


# 3. Plots

res_in <- obs_in - pred_in
auto_r1 <- auto_residual_plot(
  data.frame(Predicted = pred_in, Residuals = res_in, Observed = obs_in),
  "In-sample", "darkblue"
)


# Random CV
res_cv <- obs_cv - pred_cv
auto_r2 <- auto_residual_plot(
  data.frame(Predicted = pred_cv, Residuals = res_cv, Observed = obs_cv),
  "Random CV", "darkred"
)

# Spatial CVs
auto_r50 <- auto_residual_plot(model_spatial_cv_50$predictions, "Spatial CV 50 km", "orange")
auto_r100 <- auto_residual_plot(model_spatial_cv_100$predictions, "Spatial CV 100 km", "darkgreen")
auto_r500 <- auto_residual_plot(model_spatial_cv_500$predictions, "Spatial CV 500 km", "purple")
auto_r1000 <- auto_residual_plot(model_spatial_cv_1000$predictions, "Spatial CV 1000 km", "brown")


# Final plot (Figure 3)
grid.arrange(p1,auto_r1$plot, p2,auto_r2$plot, p50,auto_r50$plot, p100,auto_r100$plot,
             p500,auto_r500$plot, p1000, auto_r1000$plot, ncol=4)






##########   3) TEST RF-MODELS AGAINST INDEPENDENT IN-FIELD MEASURED NPP DATA ####

# --- Load data ---
m_NPP_df <- read.xlsx("Dataset_1.xlsx", sheet = "Measured_NPP_brGDGTs")

# --- Define fractions ---
tetramethylated <- c("fIa", "fIb", "fIc")
pentamethylated <- c("fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_")
hexamethylated  <- c("fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_")
all_fracs <- c(tetramethylated, pentamethylated, hexamethylated)

# --- Transform propotions and avoid 0 ---
frac_data <- m_NPP_df %>%
  select(all_of(all_fracs)) %>%
  mutate(across(everything(), ~ .x / 100 + 1e-6))

# --- Apply CLR ---
clr_transformed <- clr(frac_data)
clr_df <- as.data.frame(clr_transformed)
colnames(clr_df) <- paste0( colnames(clr_df))

# --- Prepare new dataset for prediction ---
newdata <- bind_cols(
  Sampletype = factor(m_NPP_df$Sampletype, levels = levels(model_data$Sampletype)),
  clr_df
)
newdata$Sampletype_fixed<- m_NPP_df$Sampletype

# --- Predict NPP with all the RF models created so far---

# Check all models are loaded:

model_all_data
model_cv
model_spatial_cv_50
model_spatial_cv_100
model_spatial_cv_500

####################   All data / Predicted_NPP_all
predicted_npp <- predict(model_all_data, newdata = newdata)
# Check:
predicted_npp
# --- Insert predictions into the original dataframe ---
m_NPP_df$Predicted_NPP_all<- predicted_npp

####################   Random cross validation / Predicted_NPP_cv
predicted_npp <- predict(model_cv, newdata = newdata)
# Check:
predicted_npp
# --- Insert predictions into the original dataframe ---
m_NPP_df$Predicted_NPP_cv<- predicted_npp

####################   Nested spatial cross validation / Predicted_NPP_50, Predicted_NPP_100, Predicted_NPP_500, Predicted_NPP_1000

# 50 km
predicted_npp <- predict(model_spatial_cv_50$models, newdata = newdata)
# Check:
predicted_npp
# --- Insert predictions into the original dataframe ---
# Transform list into matrx
pred_matrix <- do.call(cbind, predicted_npp)
# Compute mean
pred_mean <- rowMeans(pred_matrix)
m_NPP_df$Predicted_NPP_50<- pred_mean

#100 km
predicted_npp <- predict(model_spatial_cv_100$models, newdata = newdata)
# Check:
predicted_npp
# --- Insert predictions into the original dataframe ---
# Transform list into matrx
pred_matrix <- do.call(cbind, predicted_npp)
# Compute mean
pred_mean <- rowMeans(pred_matrix)
m_NPP_df$Predicted_NPP_100<- pred_mean

#500 km
predicted_npp <- predict(model_spatial_cv_500$models, newdata = newdata)
# Check:
predicted_npp
# --- Insert predictions into the original dataframe ---
# Transform list into matrx
pred_matrix <- do.call(cbind, predicted_npp)
# Compute mean
pred_mean <- rowMeans(pred_matrix)
m_NPP_df$Predicted_NPP_500<- pred_mean


#1000 km
predicted_npp <- predict(model_spatial_cv_1000$models, newdata = newdata)
# Check:
predicted_npp
# --- Insert predictions into the original dataframe ---
# Transform list into matrx
pred_matrix <- do.call(cbind, predicted_npp)
# Compute mean
pred_mean <- rowMeans(pred_matrix)
m_NPP_df$Predicted_NPP_1000<- pred_mean





# --- Compare observed vs predicted NPP between all models and MODIS17A NPP ---

data <- m_NPP_df

# Function to compute statistic metrics
metrics <- function(obs, pred) {
  cor_val <- cor(obs, pred, use = "complete.obs")
  lm_model <- lm(obs ~ pred)
  r2_val <- summary(lm_model)$r.squared
  rmse_val <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  list(correlation = cor_val, R2 = r2_val, RMSE = rmse_val)
}

# Compute metrics for each model
results <- tibble(
  model = c( "Predicted_NPP_all", "Predicted_NPP_cv", "Predicted_NPP_50", "Predicted_NPP_100",
            "Predicted_NPP_500", "Predicted_NPP_1000",  "MOD17A3_NPP"),
  correlation = c(
    metrics(data$NPP, data$Predicted_NPP_all)$correlation,
    metrics(data$NPP, data$Predicted_NPP_cv)$correlation,
    metrics(data$NPP, data$Predicted_NPP_50)$correlation,
    metrics(data$NPP, data$Predicted_NPP_100)$correlation,
    metrics(data$NPP, data$Predicted_NPP_500)$correlation,
    metrics(data$NPP, data$Predicted_NPP_1000)$correlation,
    metrics(data$NPP, data$MOD17A3_NPP)$correlation
  ),
  R2 = c(
    metrics(data$NPP, data$Predicted_NPP_all)$R2,
    metrics(data$NPP, data$Predicted_NPP_cv)$R2,
    metrics(data$NPP, data$Predicted_NPP_50)$R2,
    metrics(data$NPP, data$Predicted_NPP_100)$R2,
    metrics(data$NPP, data$Predicted_NPP_500)$R2,
    metrics(data$NPP, data$Predicted_NPP_1000)$R2,
    metrics(data$NPP, data$MOD17A3_NPP)$R2
  ),
  RMSE = c(
    metrics(data$NPP, data$Predicted_NPP_all)$RMSE,
    metrics(data$NPP, data$Predicted_NPP_cv)$RMSE,
    metrics(data$NPP, data$Predicted_NPP_50)$RMSE,
    metrics(data$NPP, data$Predicted_NPP_100)$RMSE,
    metrics(data$NPP, data$Predicted_NPP_500)$RMSE,
    metrics(data$NPP, data$Predicted_NPP_1000)$RMSE,
    metrics(data$NPP, data$MOD17A3_NPP)$RMSE
  )
)

print(results)

# --- Plot results --- 

plot_data <- data %>%
  select(
    NPP,
    Predicted_NPP_all,
    Predicted_NPP_cv,
    Predicted_NPP_50,
    Predicted_NPP_100,
    Predicted_NPP_500,
    Predicted_NPP_1000,
    MOD17A3_NPP
  ) %>%
  pivot_longer(
    cols = -NPP,
    names_to = "model",
    values_to = "predicted"
  )


labels <- results %>%
  mutate(
    label = paste0(
      "r = ", round(correlation, 2), "\n",
      "R² = ", round(R2, 2), "\n",
      "RMSE = ", round(RMSE, 0)
    )
  )


plot_obs_pred <- function(obs, pred, model_name, stats) {
  
  df <- data.frame(
    Observed = obs,
    Predicted = pred
  )
  m <- stats %>% filter(model == model_name)
  subtitle_txt <- sprintf(
    "r = %.3f | r² = %.3f | RMSE = %.0f",
    m$correlation,
    m$R2,
    m$RMSE
  )
  
  ggplot(df, aes(x = Observed, y = Predicted)) +
    geom_point(color = "black", alpha = 0.4, size = 1.5) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "darkred",
      fill = adjustcolor("gray50")
    ) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "darkgray"
    ) +
    labs(
      title = model_name,
      subtitle = subtitle_txt,
      x = "Observed NPP",
      y = "Predicted NPP"
    ) +
    theme_minimal()
}


p1 <- plot_obs_pred(data$NPP, data$Predicted_NPP_all,    "Predicted_NPP_all",    results)
p2 <- plot_obs_pred(data$NPP, data$Predicted_NPP_cv,     "Predicted_NPP_cv",     results)
p3 <- plot_obs_pred(data$NPP, data$Predicted_NPP_50,     "Predicted_NPP_50",     results)
p4 <- plot_obs_pred(data$NPP, data$Predicted_NPP_100,    "Predicted_NPP_100",    results)
p5 <- plot_obs_pred(data$NPP, data$Predicted_NPP_500,    "Predicted_NPP_500",    results)
p6 <- plot_obs_pred(data$NPP, data$Predicted_NPP_1000,   "Predicted_NPP_1000",   results)
p7 <- plot_obs_pred(data$NPP, data$MOD17A3_NPP,          "MOD17A3_NPP",          results)


# The empty space of this plot (Figure 4) is for the NPP map
grid.arrange(
  p7, nullGrob(), nullGrob(),
  p1, p2, p3,
  p4, p5, p6,
  ncol = 3
)


#### 4) PREDICT NPP IN PADUL ####

# --- Load data ---

Padul_df <- read.xlsx("Dataset_1.xlsx", sheet = "Padul_GDGTs")
head(Padul_df)

# --- Select compositional columns ---
# Columns used in model training
vars <- colnames(Padul_df)[4:18]  
# Check you should have the same number of columns as the model training dataset (model_data): 15 variables 
vars

# Extra variable required by the model if trained with Sample type
vars_extra <- c("Sampletype")  

comp_new <- Padul_df[, vars]


# --- Convert percentages to proportions ---
# This only necessary in case the specific percentages are expressed in a 0-100 range
# comp_new <- comp_new / 100 

# ---  Avoid zeros for CLR transformation --- 
# comp_new <- comp_new + 1e-6


# --- Apply CLR transformation ---

clr_new <- clr(comp_new)
clr_new_df <- as.data.frame(clr_new)

# --- PREPARE DATA FOR PREDICTION ---

# 1. Combine CLR-transformed data with required variables
newdata_model <- cbind(clr_new_df, 
                       Sampletype = Padul_df$Sampletype,
                       Sampletype_fixed = "Lacustrine Sediment")  # Assuming all are lacustrine / Changed in a sensitivity test (Supplementary Fig. 3)

# 2. Ensure column names match exactly with model training
# Get the model variables from the trained model
rf_model <- model_all_data$finalModel  
model_vars <- rownames(rf_model$importance)

# 3. Create dummy variables for Sampletype_fixed (as done before)
newdata_model$`Sampletype_fixedLacustrine Sediment` <- 0
newdata_model$`Sampletype_fixedPeat` <- 0
newdata_model$`Sampletype_fixedSoil` <- 0
newdata_model$`Aquatic_SPM_combined` <- 0

# Assign 1 to Lacustrine Sediment category (assuming all Padul samples are lacustrine; this was modified in a senstivity test)
for (i in 1:nrow(newdata_model)) {
  newdata_model$`Sampletype_fixedLacustrine Sediment`[i] <- 1
}

# 4. Select and order columns exactly as the model expects
newdata_ready <- newdata_model[, model_vars]


# --- MAKE PREDICTIONS WITH ALL TREES ---

# 5. Get predictions from all 1000 trees
rf_preds_all <- predict(rf_model, 
                        newdata = newdata_ready,
                        predict.all = TRUE)


# --- CALCULATE 9% CONFIDENCE INTERVALS ---

# 6. Calculate percentiles for 95% CI 
ci_95_lower <- apply(rf_preds_all$individual, 1, quantile, probs = 0.025, na.rm = TRUE)

ci_95_upper <- apply(rf_preds_all$individual, 1, quantile, probs = 0.975, na.rm = TRUE)

# 7. Create results dataframe
results_complete <- data.frame(
  SampleID = 1:nrow(Padul_df),
  Age_calBP = Padul_df$`Age.(yr.cal.BP)`,
  NPP_pred_mean = rf_preds_all$aggregate,  # Mean prediction
  NPP_CI95_lower = ci_95_lower,
  NPP_CI95_upper = ci_95_upper,
  NPP_range = ci_95_upper - ci_95_lower  # Range of of 95% CI
)

head(results_complete)

# Optional, save predictions to file
# write.csv(results_complete, "Padul_NPP_predictions_95CI.csv", row.names = FALSE)

# --- PREPARE DATA FOR PLOTTING ---

tree_indices <- sample(1:ncol(rf_preds_all$individual), 1000)

# Create dataframe with individual tree predictions
tree_data <- data.frame(Age_calBP = Padul_df$`Age.(yr.cal.BP)`)

for (i in 1:1000) {
  tree_data[[paste0("Tree_", i)]] <- rf_preds_all$individual[, tree_indices[i]]
}

# Convert to long format for ggplot
tree_long <- melt(tree_data, 
                  id.vars = "Age_calBP",
                  variable.name = "Tree",
                  value.name = "NPP")

# --- CREATE INITIAL EXPLORATORY THE PLOT ---

Padul_NPP_plot_simple <- ggplot() +
  geom_line(data = tree_long,
            aes(x = Age_calBP, y = NPP, group = Tree),
            alpha = 0.2, color = "gray", linewidth = 0.3) +
  geom_line(data = results_complete,
            aes(x = Age_calBP, y = NPP_pred_mean),
            color = "black", linewidth = 1) +
  geom_line(data = results_complete,
            aes(x = Age_calBP, y = NPP_CI95_lower),
            color = "orange", linewidth = 0.8) +
  geom_line(data = results_complete,
            aes(x = Age_calBP, y = NPP_CI95_upper),
            color = "orange", linewidth = 0.8) +
  scale_x_reverse(limits = c(27000, 5000),
                  breaks = seq(27000, 5000, by = -2000)) +
  labs(x = "Age (yr cal BP)", y = "NPP",
       title = "brGDGT-NPP") +
  theme_minimal()



print(Padul_NPP_plot_simple)




# Check whether predictions change substantially when Sample type changes (Supplementary Fig. 3):
# --- Load data ---
Padul_df <- read.xlsx("Dataset_1.xlsx", sheet = "Padul_GDGTs")

# --- Select compositional columns used in model ---
vars <- colnames(Padul_df)[4:18]  
comp_new <- Padul_df[, vars]

# --- CLR transformation ---
clr_new <- clr(comp_new)
clr_new_df <- as.data.frame(clr_new)

# --- Define all Sampletype_fixed categories exactly as in model ---
sampletypes_fixed <- c("Lacustrine Sediment", "Peat", "Soil", "Aquatic_SPM")

# --- Prepare list to store results for each sample type ---
results_list <- list()
tree_long_list <- list()

rf_model <- model_all_data$finalModel  
model_vars <- rownames(rf_model$importance)

# --- Loop over each Sampletype_fixed ---
for (stype in sampletypes_fixed) {
  
  # 1. Prepare data for model
  newdata_model <- cbind(
    clr_new_df,
    Sampletype = Padul_df$Sampletype,
    Sampletype_fixed = stype
  )
  
  # 2. Create dummy variables for Sampletype_fixed
  for (st in sampletypes_fixed) {
    colname <- paste0("Sampletype_fixed", st)
    newdata_model[[colname]] <- ifelse(st == stype, 1, 0)
  }
  
  # 3. Select columns exactly as in the model
  newdata_ready <- newdata_model[, model_vars]
  
  # --- Make predictions with all trees ---
  rf_preds_all <- predict(rf_model, newdata = newdata_ready, predict.all = TRUE)
  
  # --- Calculate 95% CI ---
  ci_95_lower <- apply(rf_preds_all$individual, 1, quantile, probs = 0.025, na.rm = TRUE)
  ci_95_upper <- apply(rf_preds_all$individual, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  results_complete <- data.frame(
    SampleID = 1:nrow(Padul_df),
    Age_calBP = Padul_df$`Age.(yr.cal.BP)`,
    NPP_pred_mean = rf_preds_all$aggregate,
    NPP_CI95_lower = ci_95_lower,
    NPP_CI95_upper = ci_95_upper,
    NPP_range = ci_95_upper - ci_95_lower,
    Sampletype_fixed = stype
  )
  
  results_list[[stype]] <- results_complete
  
  # --- Prepare individual tree data for plotting ---
  tree_indices <- sample(1:ncol(rf_preds_all$individual), 1000)
  tree_data <- data.frame(Age_calBP = Padul_df$`Age.(yr.cal.BP)`)
  
  for (i in 1:1000) {
    tree_data[[paste0("Tree_", i)]] <- rf_preds_all$individual[, tree_indices[i]]
  }
  
  tree_data$Sampletype_fixed <- stype
  tree_long <- melt(tree_data, id.vars = c("Age_calBP", "Sampletype_fixed"), 
                    variable.name = "Tree", value.name = "NPP")
  
  tree_long_list[[stype]] <- tree_long
}

# --- Combine results for plotting ---
all_results <- bind_rows(results_list)
all_tree_long <- bind_rows(tree_long_list)


NPP_plot_combined <- ggplot(data = all_results, 
                            aes(x = Age_calBP, y = NPP_pred_mean, color = Sampletype_fixed)) +
  geom_line(linewidth = 0.8) +   # Only mean lines
  labs(x = "Age (yr cal BP)",
       y = "NPP (g/m^2/yr)",
       title = "brGDGT-NPP by Sample Type",
       color = "Sample Type") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5))+
  scale_x_reverse(limits=c(27000, 5000))

print(NPP_plot_combined)


# Plot of Total Organic Content (%) in the Padul-15-05 core
TOC_df<- read.xlsx("Dataset_1.xlsx", sheet ="Padul_GDGTs")

TOC_df<- subset(TOC_df, TOC_df$`Age.(yr.cal.BP)`>4000) # remove the first three rows because they correspond to present-day 

TOC_plot<- ggplot(data=TOC_df,
                  aes(x=TOC_df$`Age.(yr.cal.BP)`, y=as.numeric(TOC))) +
  geom_point()+
  geom_line()+
  scale_x_reverse()+
  theme_minimal()+
  xlab("Age (kyr BP)")+
  ylab("Total Organic Content (%)")

TOC_plot

plot_grid(
  NPP_plot_combined, TOC_plot,
  ncol = 1,
  align = "v"
)


### Climate-derived (i.e., potential) VS brGDGT (i.e., actual) NPP ####


# --- Load data ---
# To use the pollen-derived data:
Padul_Miami_climate <- read.xlsx("Dataset_1.xlsx", sheet = "Padul_climate")
head(Padul_Miami_climate)

# To use the other paleoclimate data used in the sensitivity test:
# Padul_Miami_climate <- read.xlsx("Dataset_1.xlsx", sheet = "Regional_Paleoclimate")
# head(Padul_Miami_climate)

# --- Propagate uncertainty from MAT and MAP 95% CI to NPP 95% CI --- 
# Parameters
n_iter <- 1000
n_samples <- nrow(Padul_df)

# Store NPP outputs
npp_miami_mc <- matrix(NA, nrow = n_samples, ncol = n_iter) 


# Loop Monte Carlo
for(i in 1:n_iter){
  
  # 1) Sample MAT within 95% CI 
  mat_sim <- rtruncnorm(
    n = nrow(Padul_Miami_climate),
    a = Padul_Miami_climate$Lower95_MAT,
    b = Padul_Miami_climate$Upper95_MAT,
    mean = Padul_Miami_climate$Predicted_MAT,
    sd = (Padul_Miami_climate$Upper95_MAT - Padul_Miami_climate$Predicted_MAT)/1.96
  )
  
  # 2) Sample MAP within 95% CI
  map_sim <- rtruncnorm(
    n = nrow(Padul_Miami_climate),
    a = Padul_Miami_climate$Lower95_MAP,
    b = Padul_Miami_climate$Upper95_MAP,
    mean = Padul_Miami_climate$Predicted_MAP,
    sd = (Padul_Miami_climate$Upper95_MAP - Padul_Miami_climate$Predicted_MAP)/1.96
  )
  
  # 3) Compute NPP (Miami model) 
  npp_miami_mat <- 3000 / (1 + exp(1.315 - 0.119 * mat_sim))
  npp_miami_map <- 3000 * (1 - exp(-0.000664 * map_sim))
  npp_miami_full <- pmin(npp_miami_mat, npp_miami_map)
  
  # 4) Same chronology
  npp_miami_interp <- approx(
    x = Padul_Miami_climate$Age,
    y = npp_miami_full,
    xout = Padul_df$`Age.(yr.cal.BP)`,
    rule = 2  # extiende valores fuera del rango
  )$y
  
  # 5) Store
  npp_miami_mc[, i] <- npp_miami_interp
}

tree_long
# --- Summary statistics ---

npp_miami_mean <- apply(npp_miami_mc, 1, mean, na.rm = TRUE)
npp_miami_low  <- apply(npp_miami_mc, 1, quantile, probs = 0.025, na.rm = TRUE)
npp_miami_high <- apply(npp_miami_mc, 1, quantile, probs = 0.975, na.rm = TRUE)

# --- DF for plot

df_plot <- data.frame(
  Age = Padul_df$`Age.(yr.cal.BP)`,
  NPP_mean = npp_miami_mean,
  NPP_low  = npp_miami_low,
  NPP_high = npp_miami_high
)

# --- Plot ---

ggplot(df_plot, aes(x = Age)) +
  geom_line(aes(y = NPP_mean), color = "black", size = 1) +
  geom_ribbon(aes(ymin = NPP_low, ymax = NPP_high), fill = "black", alpha = 0.1) +
  geom_line(data = results_complete,
            aes(x = Age_calBP, y = NPP_pred_mean),
            color = "darkgreen", linewidth = 1) +
  geom_ribbon(data = results_complete,
              aes(x = Age_calBP, ymin = NPP_CI95_lower, ymax = NPP_CI95_upper), fill = "darkgreen", alpha = 0.1) +
  
  scale_x_reverse(limits = c(27000, 5000), breaks = seq(35000, 2000, by = -5000)) +
  labs(
    x = "Age (yr cal BP)",
    y = "Potential NPP"
  ) +

  theme_minimal()


# --- Deficit of NPP (DNPP), commonly known as Human Appropriation of NPP (HANPP)  --- 

Padul_Miami_climate <- read.xlsx("Dataset_1.xlsx", sheet = "Regional_Paleoclimate")
head(Padul_Miami_climate)

# 1. Calculate mean brGDGT-NPP for each age (from all 1000 trees)
rf_mean <- tree_long %>%
  group_by(Age_calBP) %>%
  summarise(NPP_rf = mean(NPP, na.rm = TRUE)) %>%
  arrange(desc(Age_calBP))


# 2. Interpolate Miami NPP to match brGDGT ages
miami_interp <- approx(x = Padul_Miami_climate$Age,
                       y = Padul_Miami_climate$NPP_Miami,
                       xout = rf_mean$Age_calBP)$y

# 3. Calculate DNPP (HANPP) (Miami - brGDGT) for each tree
hanpp_data <- tree_long %>%
  mutate(NPP_Miami = approx(Padul_Miami_climate$Age,
                            Padul_Miami_climate$NPP_Miami,
                            xout = Age_calBP)$y) %>%
  mutate(HANPP =   NPP_Miami -NPP  ) %>%  
  arrange(desc(Age_calBP))

hanpp_data

# 4. Calculate mean HANPP across trees
hanpp_mean <- hanpp_data %>%
  group_by(Age_calBP) %>%
  summarise(HANPP_mean = mean(HANPP, na.rm = TRUE),
            HANPP_sd = sd(HANPP, na.rm = TRUE))

head(hanpp_mean)


# 5. Plot DNPP/HANPP with individual trees and mean
DNPP_plot2 <- ggplot() +
  geom_line(data = hanpp_data,
            aes(x = Age_calBP, y = HANPP, group = Tree),
            alpha = 0.1, color = "darkgrey", linewidth = 0.3) +
   geom_line(data = hanpp_mean,
            aes(x = Age_calBP, y = HANPP_mean),
           color = "blue", linewidth = 0.5) +
  scale_x_reverse(limits = c(27000, 5000),
                  breaks = seq(35000, 5000, by = -5000)) +
  #geom_hline(yintercept = 743.8, linetype = "dashed", color = "blue", linewidth = 0.8) + # Present-day DNPP
  labs(x = "Age (cal yr BP)",
       y = expression(paste("NPP Deficit (g C m"^{-2}, " yr"^{-1}, ")")),
       title = " ") +
  theme_minimal()


DNPP_plot2

grid.arrange(DNPP_plot,
             DNPP_plot2,
             DNPP_plot3)




Padul_NPP_plot <- ggplot() +
  geom_line(data = tree_long,
            aes(x = Age_calBP, y = NPP, group = Tree),
            alpha = 0.2, color = "gray", linewidth = 0.3) +
  geom_line(data = results_complete,
            aes(x = Age_calBP, y = NPP_pred_mean),
            color = "black", linewidth = 1) +
  geom_line(data = results_complete,
            aes(x = Age_calBP, y = NPP_CI95_lower),
            color = "orange", linewidth = 0.8) +
  geom_line(data = results_complete,
            aes(x = Age_calBP, y = NPP_CI95_upper),
            color = "orange", linewidth = 0.8) +
 geom_line(data = Padul_Miami_climate, aes(x = Age, y = NPP_Miami),
           color = "darkblue", linewidth = 0.8, linetype = "dashed" ) +
  scale_x_reverse(limits = c(27000, 5000),
                  breaks = seq(35000, 5000, by = -5000)) +
  labs(x = "Age (yr cal BP)", y = "NPP",
       title = "brGDGT-NPP") +
  theme_minimal()



#### 5)  SUMMED PROBABILITY DISTRIBUTIONS ####

# ---- SPD by culture ----

# Delta values for marine samples
deltaR_val  <- 94
deltaR_sd   <- 61

# List of cultures
cultures <- c("Aurignacian", "Gravettian", "Solutrean", "Magdalenian", "Epipalaeolithic", "Neolithic")


# Load data
SPD_df <- read.xlsx("Dataset_1.xlsx", sheet = "SPD", detectDates = FALSE)
nrow(SPD_df)
# Select only filtered dates:
SPD_df<- subset(SPD_df, SPD_df$Included=="y")
nrow(SPD_df)
# Dataframe to store outputs
Regions_SPD <- data.frame()

# Check the number of sites, layers and dates per culture
for(cult in cultures){
  
  # Filter by culture and non-missing data
  cult_df <- SPD_df %>%
    filter(Included == "y", Culture == cult) %>%
    filter(!is.na(Age), !is.na(s.dev.), !is.na(Curve), !is.na(Labref), !is.na(layer))
  
  if(nrow(cult_df) == 0){
    message("No dates for culture: ", cult)
    next
  }
  
  # Count unique layers for this culture
  n_layers <- length(unique(cult_df$layer))
  message("Culture: ", cult, " - Number of layers: ", n_layers, " (", nrow(cult_df), " total dates)")
  
}



# Initialize empty dataframe before the next loop
Regions_SPD <- data.frame()

# 1) Callibrate, 2) Bin, 3) SPD:
for(cult in cultures){
  
  # Filter by culture and non-missing data
  cult_df <- SPD_df %>%
    filter(Included == "y", Culture == cult) %>%
    filter(!is.na(Age), !is.na(s.dev.), !is.na(Curve), !is.na(Labref), !is.na(layer))
  
  if(nrow(cult_df) == 0){
    message("No dates for culture: ", cult)
    next
  }
  
  # Calibration curves and deltaR
  curves <- ifelse(cult_df$Curve == "Terrestrial", "intcal20", "marine20")
  deltaR_vals <- ifelse(cult_df$Curve == "Marine", deltaR_val, 0)
  deltaR_sds  <- ifelse(cult_df$Curve == "Marine", deltaR_sd, 0)
  
  # Calibrate individual dates with unique IDs (= code)
  caldates_ind <- calibrate(
    x = cult_df$Age,
    errors = cult_df$s.dev.,
    calCurves = curves,
    ids = paste0(cult_df$Labref, "_", seq_len(nrow(cult_df))),
    delta.R = deltaR_vals,
    delta.STD = deltaR_sds,
    calMatrix = TRUE
  )
  
  # Combine dates by archaeological layer
  combined <- cult_df %>%
    group_by(layer) %>%
    summarise(
      combined_rads = list(
        combine(caldates_ind[which(cult_df$layer == unique(layer))], fixIDs = TRUE)$date
      ),
      .groups = "drop"
    )
  
  # Prepare vectors for SPD
  ages <- unlist(lapply(combined$combined_rads, function(x) x$age))
  sdev <- unlist(lapply(combined$combined_rads, function(x) x$error))
  
  # 300-year bins
  bins <- binPrep(
    sites = rep(combined$layer, times = sapply(combined$combined_rads, length)),
    ages  = ages,
    h     = 300
  )
  
  # Compute SPD
  spd_cult <- spd(caldates_ind, timeRange = c(35000, 2000), bins = bins)
  

  # Append results for this culture to Regional_SPD
  Regions_SPD <- rbind(
    Regions_SPD,
    data.frame(
      PrDens  = spd_cult$grid$PrDens,
      Age     = spd_cult$grid$calBP,
      Culture = cult
    )
  )
  
  # Optional: print progress
  message("Added ", nrow(spd_cult$grid), " rows for ", cult)
}


# Now check the results - this will show ALL cultures
head(Regions_SPD)
table(Regions_SPD$Culture)  # This will show counts for each culture


# Plot SPD
SPD_l3 <- ggplot(Regions_SPD, aes(x = Age, y = PrDens, fill = Culture)) +
  geom_area(alpha = 0.5, position = "identity", color = "black") +
  scale_x_reverse(limits = c(27000, 5000),
                  breaks = seq(35000, 5000, by = -5000)) +
  theme_classic() +
  labs(x = "cal BP", y = "Probability Density", title = "Summed Probability Distribution by Culture") +
  scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "none")

SPD_l3



#NGRIP data for plot

NGRIP_df <- read.xlsx("Dataset_1.xlsx", sheet = "NGRIP", detectDates = FALSE)
head(NGRIP_df)

# Create df with events (for visual purposes)
events <- data.frame(
  Event = c( "H2", "LGM", "H1",  "Younger Dryas", "8.2"),
  start = c( 24500, 21500, 17000,  12900,  8400),
  end   = c( 23000, 19000, 15000, 11700,  8100)
)


# Plot
NGRIP_plot <- ggplot() +
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4) +
  geom_line(data = NGRIP_df, aes(x = Age * 1000, y = d18O),
            color = "blue", size = 1) +
  scale_x_reverse(limits = c(27000, 5000),
                  breaks = seq(35000, 5000, by = -5000)) +
  labs(
    x = " ", # Age (yr BP)
    y = expression(delta^{18}*O~("\u2030")),
    title = " "
  ) +
  theme_minimal(base_size = 14)

NGRIP_plot


# ---- Pollen plot ----
# Load data
pollen_df <- read.xlsx("Dataset_1.xlsx", sheet = "Pollen", detectDates = FALSE)
head(pollen_df)

# % of Mediterranean pollen taxa and Precipitation Index combined in the same plot:
lp_combined <- ggplot() +
  geom_line(data = pollen_df, aes(x = Age, y = Mediterranean), color = "darkgreen") +
  geom_line(data = pollen_df, aes(x = Age, y = P.Index), color = "black") +
  scale_y_continuous(
    name = "Mediterranean",
    sec.axis = sec_axis(~., name = "P.Index")
  ) +
  scale_x_reverse(limits = c(27000, 5000),
                  breaks = seq(35000, 5000, by = -5000)) +
  labs(x = "Age (years BP)") +
  theme_minimal()

lp_combined


# Charcoal concentration (n charcoals/cm3) throughout the sequence:
Charcoal <- ggplot() +
  geom_point(data = pollen_df, aes(x = Age, y = Charcoals_cm3, color = "Charcoals_cm3")) +
  scale_x_reverse(limits = c(27000, 5000),
                  breaks = seq(35000, 5000, by = -5000)) +
  labs(x = "Age (years BP)", color = "Charcoals_cm3") +
  scale_y_continuous(limits = c(0, 300000)) +
  theme_minimal() +
  theme(legend.position = "none")  # <--- Corrección aquí

print(Charcoal)



plot_grid(
  NGRIP_plot,
  lp_combined,
  Padul_NPP_plot,
  #Padul_NPP_plot_simple,
  DNPP_plot, 
  Charcoal,
  SPD_l3,
  ncol = 1,
  align = "v",
  rel_heights = c(0.8, 0.7, 1.3, 1.5, 0.9, 1)  # Adjust these numbers as needed
)



# Check data for summary statstics
head(Regions_SPD)
head(tree_long)
head(hanpp_data)

# Create a dataframe that combines PrDens, NPP_brGDGT and DNPP 

if(exists("tree_long") && exists("hanpp_data")) {
  
  npp_by_age <- tree_long %>%
    group_by(Age_calBP) %>%
    summarise(NPP_brGDGT = mean(NPP, na.rm = TRUE))
  
  hanpp_by_age <- hanpp_data %>%
    group_by(Age_calBP) %>%
    summarise(HANPP = mean(HANPP, na.rm = TRUE))
  
  all_data <- Regions_SPD %>%
    left_join(npp_by_age, by = c("Age" = "Age_calBP")) %>%
    left_join(hanpp_by_age, by = c("Age" = "Age_calBP"))
  
} else {
  print("Check names in dataframes")
}

# Check
head(all_data)

# Filter PrDens >= 0.001
all_data_fixed_thresh <- all_data %>%
  filter(PrDens >= 0.001)

# Correlation between PrDens and Defcit of NPP (DNPP/HANPP)
DNPP_Dens_plot<- ggplot(all_data_fixed_thresh, aes(x = log10(PrDens), y = log10(HANPP))) +
  geom_point(aes(color = Culture), size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "darkred") +
  stat_cor(method = "pearson", 
           label.x.npc = "left", 
           label.y.npc = "top",
           size = 3) +
  labs(
    x = "PrDens",
    y = "HANPP",
    title = " "
  ) +
  theme_minimal()



DNPP_Dens_plot


Pollen_df<- read.xlsx("Dataset_1.xlsx", sheet="Pollen")

colnames(Pollen_df)

#  Same ages for DNPP and NPP_rf
hanpp_interp <- approx(
  x = hanpp_data$Age_calBP,
  y = hanpp_data$HANPP,
  xout = rf_mean$Age_calBP,
  rule = 2 
)$y

# The same for Pollen_df (Pollen_cm3, Carbones_cm3, PAR) 
pollen_interp <- data.frame(
  Pollen_cm3    = approx(Pollen_df$Age, Pollen_df$Pollen_cm3, xout = rf_mean$Age_calBP, rule = 2)$y,
  Carbones_cm3  = approx(Pollen_df$Age, Pollen_df$Charcoals_cm3, xout = rf_mean$Age_calBP, rule = 2)$y,
  PAR           = approx(Pollen_df$Age, Pollen_df$PAR, xout = rf_mean$Age_calBP, rule = 2)$y,
  Total_f           = approx(Pollen_df$Age, Pollen_df$Total_f, xout = rf_mean$Age_calBP, rule = 2)$y,
  NPP_climate    = approx(Pollen_df$Age, Pollen_df$NPP_Miami, xout = rf_mean$Age_calBP, rule = 2)$y,
  CAR           = approx(Pollen_df$Age, Pollen_df$CAR, xout = rf_mean$Age_calBP, rule = 2)$y,
  Charcoals           = approx(Pollen_df$Age, Pollen_df$Charcoals, xout = rf_mean$Age_calBP, rule = 2)$y
  )

#  Combined dataframe
combined_df <- rf_mean %>%
  mutate(
    HANPP = hanpp_interp
  ) %>%
  bind_cols(pollen_interp)

# Check
head(combined_df)



combined_df<- subset(combined_df, combined_df$Age_calBP> 4999 & combined_df$Age_calBP<27000)

# Correlation potential vs actual NPP
NPP_actvs_pot<-ggplot(combined_df, aes(x = NPP_rf, y = NPP_climate)) +
  # Points colored by Age
  geom_point(aes(color = Age_calBP ), size = 3, alpha = 0.8) +
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "darkgrey", alpha = 0.2) +
  # Correlation text
  stat_cor(method = "spearman", 
           label.x.npc = "left", 
           label.y.npc = "top",
           size = 5) +
  # Color gradient for Age
  scale_color_gradient(low = "blue", high = "red", name = "Age (yr BP)") +
  
  labs(
    x = "NPP brGDGT",
    y = "NPP Climate"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

NPP_actvs_pot


# Actual NPP vs Pollen Accumulation Rates
NPP_PAR<-ggplot(combined_df, aes(x = log10(NPP_rf), y = log10(PAR))) +
  # Points colored by Age
  geom_point(aes(color = Age_calBP ), size = 3, alpha = 0.8) +
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "darkgrey", alpha = 0.2) +
  # Correlation text
  stat_cor(method = "spearman", 
           label.x.npc = "left", 
           label.y.npc = "top",
           size = 5) +
  # Color gradient for Age
  scale_color_gradient(low = "blue", high = "red", name = "Age (yr BP)") +
  
  labs(
    x = "NPP brGDGT",
    y = "PAR"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

NPP_PAR


colnames(combined_df)

# Actual NPP vs Total pollen count (grains/cm3)
NPP_PollenTotal<-ggplot(combined_df, aes(x = log10(NPP_rf), y = log10(Total_f))) +
  # Points colored by Age
  geom_point(aes(color = Age_calBP ), size = 3, alpha = 0.8) +
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "darkgrey", alpha = 0.2) +
  # Correlation text
  stat_cor(method = "spearman", 
           label.x.npc = "left", 
           label.y.npc = "top",
           size = 5) +
  # Color gradient for Age
  scale_color_gradient(low = "blue", high = "red", name = "Age (yr BP)") +
  
  labs(
    x = "NPP brGDGT",
    y = "Pollen count"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

NPP_PollenTotal


grid.arrange(NPP_actvs_pot, NPP_PollenTotal, NPP_PAR, ncol =1 )



head(combined_df)


# Now check association between charcoal concentration (charcoals/cm3) and
# deficit of NPP for the whole sequence, and since the Magdalenian:

combined_df_P<- subset(combined_df, combined_df$Age_calBP<16000)

Pot_Carchoal_16<-ggplot(combined_df_P, aes(x = combined_df_P$Carbones_cm3, y = HANPP)) +
  # Points colored by Age
  geom_point(aes(color = Age_calBP ), size = 3, alpha = 0.8) +
  # Regression line
   geom_smooth(method = "lm", se = TRUE, color = "black", fill = "darkgrey", alpha = 0.2) +
  # Correlation text
  stat_cor(method = "pearson", 
        label.x.npc = "left", 
      label.y.npc = "top",
      size = 5) +
  # Color gradient for Age
  scale_color_gradient(low = "blue", high = "red", name = "Age (yr BP)") +
  
  labs(
    x = "Charcoal",
    y = "HANPP"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )
Pot_Carchoal_All
Act_Carchoal_All
Act_Carchoal_16
Pot_Carchoal_16


grid.arrange(Pot_Carchoal_All, Pot_Carchoal_16, Act_Carchoal_All, Act_Carchoal_16, ncol=2)



colnames(combined_df)

# The same with charcoal accumulation rates
DNPP_CAR<-ggplot(combined_df, aes(x = log10(HANPP), y = log10(CAR))) +
  # Points colored by Age
  geom_point(aes(color = Age_calBP ), size = 3, alpha = 0.8) +
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "darkgrey", alpha = 0.2) +
  # Correlation text
  stat_cor(method = "pearson", 
           label.x.npc = "left", 
           label.y.npc = "top",
           size = 5) +
  # Color gradient for Age
  scale_color_gradient(low = "blue", high = "red", name = "Age (yr BP)") +
  
  labs(
    x = "NPP brGDGT",
    y = "CAR"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

NPP_CAR


DNPP_Charcoal_cm3<-ggplot(combined_df, aes(x = log10(Carbones_cm3) , y =log10(HANPP)    )) +
  # Points colored by Age
  geom_point(aes(color = Age_calBP ), size = 3, alpha = 0.8) +
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "darkgrey", alpha = 0.2) +
  # Correlation text
  stat_cor(method = "pearson", 
           label.x.npc = "left", 
           label.y.npc = "top",
           size = 5) +
  # Color gradient for Age
  scale_color_gradient(low = "blue", high = "red", name = "Age (yr BP)") +
  
  labs(
    x = "DNPP",
    y = "Total Charcoal"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

DNPP_Charcoal_cm3
NPP_CAR




# CORRELATIONS WTH MEAN VALUES FOR EACH CULTURE ####

spd_ages <- Regions_SPD %>%
  filter(PrDens >= 0.001) %>%
  select(Age, Culture, PrDens) %>%
  distinct()  

npp_interp <- approx(
  x = tree_long$Age_calBP,
  y = tree_long$NPP,
  xout = spd_ages$Age,
  rule = 2  
)

hanpp_by_age <- hanpp_data %>%
  group_by(Age_calBP) %>%
  summarise(HANPP_mean = mean(HANPP, na.rm = TRUE))

hanpp_interp <- approx(
  x = hanpp_by_age$Age_calBP,
  y = hanpp_by_age$HANPP_mean,
  xout = spd_ages$Age,
  rule = 2
)

if("NPP_Miami" %in% names(Pollen_df)) {
  miami_interp <- approx(
    x = Pollen_df$Age,
    y = Pollen_df$NPP_Miami,
    xout = spd_ages$Age,
    rule = 2
  )
} else {
  stop("NPP_Miami not found: ", 
       paste(names(Pollen_df), collapse = ", "))
}

# Create dataframe
all_data_complete <- spd_ages %>%
  mutate(
    NPP_brGDGT = npp_interp$y,
    HANPP = hanpp_interp$y,
    NPP_Miami = miami_interp$y
  )


summary(all_data_complete)

#### 6) SUMMARY STATISTICS AND CORRELATIONS --------------------------------

stats_by_culture <- all_data_complete %>%
  group_by(Culture) %>%
  summarise(
    # NPP brGDGT
    NPP_brGDGT_mean = mean(NPP_brGDGT, na.rm = TRUE),
    NPP_brGDGT_sd = sd(NPP_brGDGT, na.rm = TRUE),
    NPP_brGDGT_min = min(NPP_brGDGT, na.rm = TRUE),
    NPP_brGDGT_max = max(NPP_brGDGT, na.rm = TRUE),
    
    # NPP climate
    NPP_Miami_mean = mean(NPP_Miami, na.rm = TRUE),
    NPP_Miami_sd = sd(NPP_Miami, na.rm = TRUE),
    NPP_Miami_min = min(NPP_Miami, na.rm = TRUE),
    NPP_Miami_max = max(NPP_Miami, na.rm = TRUE),
    
    # DNPP
    HANPP_mean = mean(HANPP, na.rm = TRUE),
    HANPP_sd = sd(HANPP, na.rm = TRUE),
    HANPP_min = min(HANPP, na.rm = TRUE),
    HANPP_max = max(HANPP, na.rm = TRUE),
    
    # SPD
    PrDens_mean = mean(PrDens, na.rm = TRUE),
    PrDens_sd = sd(PrDens, na.rm = TRUE),
    
    # Range
    Age_min = min(Age, na.rm = TRUE),
    Age_max = max(Age, na.rm = TRUE),
    Age_span = Age_max - Age_min,
    n_edades = n()
  ) %>%
  arrange(desc(Age_min))  # Order from older to recent

# Results
print(stats_by_culture)


#  Exploratory plots ------------------------------------------------
all_data_complete<- subset(all_data_complete, all_data_complete$Culture != "Aurignacian") # Remove the Aurignaciain because it is not within the scope of this study

# 4.1 Boxplot NPP_brGDGT / Culture
p1 <- ggplot(all_data_complete, aes(x = reorder(Culture, -Age), y = NPP_brGDGT, fill = Culture)) +
  geom_boxplot(alpha = 0.7,   outlier.color = "gray",  
               outlier.alpha = 0.5, size = 0.7) +
  stat_summary(fun = mean, geom = "point",  color = "black") +
  labs(
    title = "",
    x = "",
    y = "NPP brGDGT (g C m⁻² yr⁻¹)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


p1


# 4.2 Boxplot NPP_POT / Culture
p2 <- ggplot(all_data_complete, aes(x = reorder(Culture, -Age), y = NPP_Miami, fill = Culture)) +
  geom_boxplot(alpha = 0.7,   outlier.color = "gray",  
               outlier.alpha = 0.5, size = 0.7) +
  stat_summary(fun = mean, geom = "point",  color = "black") +
  labs(
    title = "",
    x = " ",
    y = "NPP climate (g C m⁻² yr⁻¹)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
p2

# 4.3 Boxplot DNPP / Culture
p3 <- ggplot(all_data_complete, aes(x = reorder(Culture, -Age), y = HANPP, fill = Culture)) +
  geom_boxplot(alpha = 0.7,   outlier.color = "gray",  # Color gris para los outliers
               outlier.alpha = 0.5, size = 0.7) +
  stat_summary(fun = mean, geom = "point",  color = "black") +
  labs(
    title = " ",
    x = " ",
    y = "Deficit NPP (g C m⁻² yr⁻¹)"
  ) +
  geom_hline(yintercept = 743.8, linetype = "dashed", color = "darkred", linewidth = 0.8) + # Present-day DNPP
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


grid.arrange(p1, p2,p3 ,ncol = 3, nrow = 1)


#  Summary Table ---------------------------------


article_table <- stats_by_culture %>%
  select(Culture, Age_min, Age_max, 
         NPP_brGDGT_mean, NPP_brGDGT_sd,
         NPP_Miami_mean, NPP_Miami_sd,
         HANPP_mean, HANPP_sd,
         PrDens_mean) %>%
  mutate(
    Age_range = paste(round(Age_min/1000, 1), "-", round(Age_max/1000, 1), "ka BP"),
    NPP_brGDGT = paste0(round(NPP_brGDGT_mean, 1), " ± ", round(NPP_brGDGT_sd, 1)),
    NPP_Miami = paste0(round(NPP_Miami_mean, 1), " ± ", round(NPP_Miami_sd, 1)),
    HANPP = paste0(round(HANPP_mean, 1), " ± ", round(HANPP_sd, 1)),
    Population = round(PrDens_mean * 10000, 2) 
  ) %>%
  select(Culture, Age_range, NPP_brGDGT, NPP_Miami, HANPP, Population)

print(article_table)

# Save if needed
# write.csv(stats_by_culture, "Results.csv", row.names = FALSE)


# CORRELATIONS BETWEEN PRDENS AND DNPP

# SPD
spd_ages <- Regions_SPD %>%
  filter(PrDens >= 0.001) %>%
  select(Age, Culture, PrDens) %>%
  distinct()

# Same ages
hanpp_by_age <- hanpp_data %>%
  group_by(Age_calBP) %>%
  summarise(HANPP_mean = mean(HANPP, na.rm = TRUE))

hanpp_interp <- approx(
  x = hanpp_by_age$Age_calBP,
  y = hanpp_by_age$HANPP_mean,
  xout = spd_ages$Age,
  rule = 2
)

# Create DF
all_data_complete <- spd_ages %>%
  mutate(HANPP = hanpp_interp$y) %>%
  filter(!is.na(HANPP) & !is.na(PrDens))

# Additional summary --------------------------------

culture_stats <- all_data_complete %>%
  group_by(Culture) %>%
  summarise(
    mean_PrDens = mean(PrDens, na.rm = TRUE),
    sd_PrDens = sd(PrDens, na.rm = TRUE),
    mean_HANPP = mean(HANPP, na.rm = TRUE),
    sd_HANPP = sd(HANPP, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()


print(culture_stats)

# Correlation tests and plots -------------------------------

# Correlación usando las medias por cultura
cor_test <- cor.test(culture_stats$mean_PrDens, culture_stats$mean_HANPP, method = "pearson")
r_value <- round(cor_test$estimate, 3)
r2_value <- round(r_value^2, 3)
p_value <- round(cor_test$p.value, 4)

print(paste(" r =", r_value, "r² =", r2_value, "p =", p_value))



culture_colors <- c(
  "Aurignacian" = "#440154",
  "Gravetian" = "#3b528b",
  "Solutrean" = "#21918c",
  "Magdalenian" = "#5ec962",
  "Epipaleolithic" = "#fde725",
  "Neolithic" = "#f98e09"
)


head(culture_stats)



culture_stats$Layers<- c(3, 37, 11, 31, 42, 21) # number of layers
culture_stats$Dates<- c(6, 41, 28, 50, 96, 43) # number of dates

p <- ggplot(culture_stats, aes(x = mean_PrDens, y = mean_HANPP, color = Culture)) +
  geom_errorbar(aes(ymin = mean_HANPP - sd_HANPP, ymax = mean_HANPP + sd_HANPP), 
                width = 0.05 * max(culture_stats$mean_PrDens, na.rm = TRUE), 
                size = 0.8, alpha = 0.7) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", alpha = 0.2) +
  annotate("text", 
           x = max(culture_stats$mean_PrDens) * 0.8, 
           y = max(culture_stats$mean_HANPP) * 0.95,
           label = paste("r =", r_value, "\nr² =", r2_value, "\np =", p_value),
           size = 2.5, hjust = 0, vjust = 1) +
  geom_text(aes(label = Culture), 
            vjust = -1.5, hjust = 0.5, size = 2.5) +
  scale_color_manual(values = culture_colors) +
  labs(
    title = "",
    x = "Probability Density",
    y = "DNPP (g C m⁻² yr⁻¹)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey70"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )


print(p)

culture_stats$Layers<- c(3, 37, 11, 31, 42, 21) # number of layers
culture_stats$Dates<- c(6, 41, 28, 50, 96, 43) # number of dates

cor_test <- cor.test(culture_stats$Layers, culture_stats$mean_HANPP, method = "pearson")
r_value <- round(cor_test$estimate, 3)
r2_value <- round(r_value^2, 3)
p_value <- round(cor_test$p.value, 4)

colnames(culture_stats)
p_l <- ggplot(culture_stats, aes(x = Layers, y = mean_HANPP, color = Culture)) +
  geom_errorbar(aes(ymin = mean_HANPP - sd_HANPP, ymax = mean_HANPP + sd_HANPP), 
                width = 0.05 * max(culture_stats$mean_PrDens, na.rm = TRUE), 
                size = 0.8, alpha = 0.7) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", alpha = 0.2) +
  annotate("text", 
           x = max(culture_stats$mean_PrDens) * 0.8, 
           y = max(culture_stats$mean_HANPP) * 0.95,
           label = paste("r =", r_value, "\nr² =", r2_value, "\np =", p_value),
           size = 2.5, hjust = 0, vjust = 1) +

  geom_text(aes(label = Culture), 
            vjust = -1.5, hjust = 0.5, size = 2.5) +
  scale_color_manual(values = culture_colors) +
  labs(
    title = "",
    x = "Number of layers",
    y = "DNPP (g C m⁻² yr⁻¹)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey70"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )


print(p_l)




cor_test <- cor.test(culture_stats$Dates, culture_stats$mean_HANPP, method = "pearson")
r_value <- round(cor_test$estimate, 3)
r2_value <- round(r_value^2, 3)
p_value <- round(cor_test$p.value, 4)

colnames(culture_stats)
p_d <- ggplot(culture_stats, aes(x = Dates, y = mean_HANPP, color = Culture)) +
  geom_errorbar(aes(ymin = mean_HANPP - sd_HANPP, ymax = mean_HANPP + sd_HANPP), 
                width = 0.05 * max(culture_stats$mean_PrDens, na.rm = TRUE), 
                size = 0.8, alpha = 0.7) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", alpha = 0.2) +
  annotate("text", 
           x = max(culture_stats$mean_PrDens) * 0.8, 
           y = max(culture_stats$mean_HANPP) * 0.95,
           label = paste("r =", r_value, "\nr² =", r2_value, "\np =", p_value),
           size = 2.5, hjust = 0, vjust = 1) +
  geom_text(aes(label = Culture), 
            vjust = -1.5, hjust = 0.5, size = 2.5) +
  scale_color_manual(values = culture_colors) +
  labs(
    title = "",
    x = "Number of dates",
    y = "DNPP (g C m⁻² yr⁻¹)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey70"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )


print(p_d)

grid.arrange(p, p_l)
grid.arrange(p, p_l, DNPP_Dens_plot, ncol = 3, nrow = 1)
grid.arrange(p1, p2,p4, ncol = 3, nrow = 1)
grid.arrange(p, p_l, p_d,
             ncol = 3, nrow = 1)
grid.arrange(DNPP_Dens_plot,DNPP_Charcoal_cm3,
             ncol = 2, nrow = 1)


############# ADDITIONAL INFORMATION ------------------------------

# The rest of the code is available in the Script_2.R file. It includes all the modelling of
# herbivore biomass and carrying capapcity of secondary consumers.
#


