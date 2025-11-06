
rm(list = ls()) # Clear all
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## Libraries

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

## a) SUMMARY STATISTICS & EXPLORATORY PLOTS ####


# 1. --- Read the data ---
data <- read.xlsx("Dataset.xlsx", sheet = "Dataset_R")

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

# Cambiar nombres de columnas y filas en la matriz
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


# --- Define brGDGT compound groups ---
tetramethylated <- c("fIa", "fIb", "fIc")
pentamethylated <- c("fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_")
hexamethylated  <- c("fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_")

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

# --- Compute Spearman correlation per fraction ---
fractions <- unique(data_long$Fraction)
stats_df <- lapply(fractions, function(frac){
  datf <- data_long %>% filter(Fraction == frac)
  if(nrow(datf) < 3) {
    return(tibble(Fraction = frac, spearman_rho = NA, spearman_p = NA))
  }
  sp <- suppressWarnings(cor.test(datf$NPP, datf$Value, method = "spearman"))
  tibble(Fraction = frac, spearman_rho = sp$estimate, spearman_p = sp$p.value)
}) %>% bind_rows()

# --- Function to build a single fraction plot ---
make_frac_plot <- function(frac){
  datf <- data_long %>% filter(Fraction == frac)
  st <- stats_df %>% filter(Fraction == frac)
  
  ggplot(datf, aes(x = NPP, y = Value, color = Sampletype)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "loess", color = "black", se = FALSE, size = 1.1) +
    labs(
      title = gsub("_", "'", frac),
      x = expression("NPP (g C m"^-2*" yr"^-1*")"),
      y = "Percentage (%)",
      color = "Sample type"
    ) +
    scale_color_manual(values = sample_colors, drop = FALSE) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    annotate("text",
             x = min(datf$NPP, na.rm = TRUE),
             y = max(datf$Value, na.rm = TRUE),
             label = paste0("ρ=", round(st$spearman_rho, 2),
                            ", p=", signif(st$spearman_p, 2)),
             hjust = -0.05, vjust = 1, size = 3)
}

# --- Build plots for each group ---
plots_tetra <- lapply(c("tetramethylated", tetramethylated), make_frac_plot)
plots_penta <- lapply(c("pentamethylated", pentamethylated), make_frac_plot)
plots_hexa  <- lapply(c("hexamethylated", hexamethylated), make_frac_plot)

# --- Pad plot lists to equal length (empty placeholders first) ---
max_cols <- max(length(plots_tetra), length(plots_penta), length(plots_hexa))
pad_plots <- function(plist, n = max_cols){
  if(length(plist) < n){
    empties <- replicate(n - length(plist), ggplot() + theme_void(), simplify = FALSE)
    plist <- c(empties, plist)
  }
  return(plist)
}
plots_tetra <- pad_plots(plots_tetra)
plots_penta <- pad_plots(plots_penta)
plots_hexa  <- pad_plots(plots_hexa)

# --- Flatten into one list for the final grid ---
all_plots <- list(plots_tetra, plots_penta, plots_hexa) %>% do.call(c, .)

# --- Remove legends from individual subplots ---
plots_noleg <- lapply(all_plots, function(p) {
  if(inherits(p, "ggplot")) p + theme(legend.position = "none") else p
})

# --- Build a unified legend ---
legend_df <- data_long %>% distinct(Sampletype) %>% mutate(dummy = 1)
legend_plot <- ggplot(legend_df, aes(x = dummy, y = Sampletype, color = Sampletype)) +
  geom_point() +
  scale_color_manual(values = sample_colors) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank())

legend_grob <- get_legend(legend_plot)

# --- Assemble final figure with cowplot ---
plot_matrix <- plot_grid(plotlist = plots_noleg, ncol = max_cols, align = "hv")
final_plot <- plot_grid(plot_matrix, legend_grob, ncol = 1, rel_heights = c(1, 0.08))

# --- Display final plot ---
print(final_plot)








#### b) RANDOM FOREST MODEL WITH NORMALISED PREDICTORS ####

# --- Load data if necessary ---
data <- read.xlsx("Dataset.xlsx", sheet = "Dataset_R")
#Check
head(data)

# --- Select compounds ---
vars <- colnames(data)[11:25]
# and include sample type and NPP
vars_extra <- c("Sampletype", "NPP")

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

# --- Formula ---
formula <- as.formula(paste("NPP ~", paste(c(colnames(clr_data), "Sampletype"), collapse = " + ")))

# --- Reproducibility ---
set.seed(123) 

# --- Cross-validation ---
train_control <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final",
  returnResamp = "final"
)

# --- RF model with cross-validation ---
model_cv <- train(
  formula,
  data = model_data,
  method = "rf",
  trControl = train_control,
  importance = TRUE,
  ntree = 1000, #500 # Model was re-run with ntrees=500 and ntrees=1000; outcomes identical
  tuneLength = 3
)


### --- Save model ---
#  saveRDS(model_cv, file = "model_brGDGT-NPP_cv.rds")

### --- Load model ---
model_cv <- readRDS("model_brGDGT-NPP_cv.rds")

print(model_cv)

# Variable importance
importance <- varImp(model_cv, scale = FALSE)

plot(importance)


# --- In-sample predictions ---
# Now, the model predictions are evaluated with the whole dataset
# First, data prepared as previously done with the training dataset:
all_fracs <- colnames(model_data)[grepl("^clr_", colnames(model_data))] 
newdata_in <- model_data  # aldready contains CRL+ Sampletype

obs_in <- newdata_in$NPP
pred_in <- predict(model_cv, newdata = newdata_in)
r2_in <- 1 - sum((obs_in - pred_in)^2) / sum((obs_in - mean(obs_in))^2)

# --- Cross-validated predictions ---
cv_pred <- model_cv$pred
cv_pred <- subset(cv_pred, mtry == model_cv$bestTune$mtry) 
obs_cv <- cv_pred$obs
pred_cv <- cv_pred$pred
r2_cv <- 1 - sum((obs_cv - pred_cv)^2) / sum((obs_cv - mean(obs_cv))^2)

# --- Dataframes for plot ---
df_cv <- data.frame(Observed = obs_cv, Predicted = pred_cv, Type = "Cross-validated")
df_in <- data.frame(Observed = obs_in, Predicted = pred_in, Type = "In-sample")



# --- Plot ---
ggplot() +
  geom_point(data = df_cv, aes(x = Observed, y = Predicted), color = "red", alpha = 0.5, size = 2) +
  geom_smooth(data = df_cv, aes(x = Observed, y = Predicted), method = "lm", se = TRUE,
              color = "red", fill = "pink", alpha = 0.2) +
  geom_point(data = df_in, aes(x = Observed, y = Predicted), color = "blue", alpha = 0.6, size = 2) +
  geom_smooth(data = df_in, aes(x = Observed, y = Predicted), method = "lm", se = TRUE,
              color = "blue", fill = "lightblue", alpha = 0.3) +
  labs(title = "Random Forest: Observed vs Predicted NPP",
       x = "Observed NPP",
       y = "Predicted NPP") +
  annotate("text", x = min(obs_in), y = max(pred_in), 
           label = paste("r² of the entire dataset=", round(r2_in, 2)), 
           hjust = 0, vjust = 1.5, color = "blue", size = 4) +
  annotate("text", x = min(obs_in), y = max(pred_in)*0.95, 
           label = paste("r² of the cross-validation=", round(r2_cv, 2)), 
           hjust = 0, vjust = 1.5, color = "red", size = 4) +
  theme_minimal()




# --- Calculate residuals ---
summary(model_data$NPP)
# Set a reasonable tolerance
tol <- 100  # ±100 units, roughly 15%-16% of median NPP

# Calculate percentages within tolerance
res_in <- obs_in - pred_in
res_cv <- obs_cv - pred_cv

pct_in <- mean(abs(res_in) <= tol) * 100
pct_cv <- mean(abs(res_cv) <= tol) * 100

cat("All:", round(pct_in, 1), "% within ±", tol, "\n")
cat("Cross-validated:", round(pct_cv, 1), "% within ±", tol, "\n")




# Percentage of data within range
pct_in <- mean(abs(res_in) <= tol) * 100
pct_cv <- mean(abs(res_cv) <= tol) * 100

# Compute residuals
df_in$Residuals <- df_in$Observed - df_in$Predicted
df_cv$Residuals <- df_cv$Observed - df_cv$Predicted

# Combine in-sample and cross-validated residuals
df_res_all <- rbind(df_in, df_cv)


# Base scatter plot 
p <- ggplot(df_res_all, aes(x = Predicted, y = Residuals, color = Type)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -tol, ymax = tol),
            fill = "lightgreen", alpha = 0.01, inherit.aes = FALSE) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("In-sample" = "blue", "Cross-validated" = "red")) +
  labs(title = "Residuals vs Predicted NPP",
       x = "Predicted NPP",
       y = "Residuals (Observed - Predicted)") +
  theme_minimal()

# With histogram
ggMarginal(p, type = "histogram", margins = "y", groupColour = TRUE, groupFill = TRUE)




#### c) TEST RF-MODEL AGAINST MEASURED NPP DATA ####

# --- Load data ---
m_NPP_df <- read.xlsx("Dataset.xlsx", sheet = "Measured_NPP_brGDGTs")

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
newdata$Sampletype<- m_NPP_df$Sampletype
predictors(model_cv)

# --- Predict NPP ---
predicted_npp <- predict(model_cv, newdata = newdata)
# Check:
predicted_npp
# --- Insert predictions into the original dataframe ---
m_NPP_df$Predicted_NPP <- predicted_npp
head(m_NPP_df)



# --- Miami NPP model---

T <- m_NPP_df$wc2.1_5m_bio_1.tif    # °C
P <- m_NPP_df$wc2.1_5m_bio_12.tif       # mm

# --- Compute climate-potentia) NPP with Miami model --- 
NPP_T <- 3000 / (1 + exp(1.315 - 0.119 * T))
NPP_P <- 3000 * (1 - exp(-0.000664 * P))

m_NPP_df$Miami_NPP <- pmin(NPP_T, NPP_P)

# Check rows
head(m_NPP_df[, c("NPP", "Miami_NPP")])


# --- Compare observed vs predicted NPP between three models: Miami, MODIS and brGDGT-NPP ---

data <- m_NPP_df

# Function to compute different metrics
metrics <- function(obs, pred) {
  cor_val <- cor(obs, pred, use = "complete.obs")
  lm_model <- lm(obs ~ pred)
  r2_val <- summary(lm_model)$r.squared
  rmse_val <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  list(correlation = cor_val, R2 = r2_val, RMSE = rmse_val)
}

# Compute metrics for each model
results <- tibble(
  model = c( "Predicted_NPP", "MOD17A3_NPP", "Miami_NPP"),
  correlation = c(
    metrics(data$NPP, data$Predicted_NPP)$correlation,
    metrics(data$NPP, data$MOD17A3_NPP)$correlation,
    metrics(data$NPP, data$Miami_NPP)$correlation
  ),
  R2 = c(
    metrics(data$NPP, data$Predicted_NPP)$R2,
    metrics(data$NPP, data$MOD17A3_NPP)$R2,
    metrics(data$NPP, data$Miami_NPP)$R2
  ),
  RMSE = c(
    metrics(data$NPP, data$Predicted_NPP)$RMSE,
    metrics(data$NPP, data$MOD17A3_NPP)$RMSE,
    metrics(data$NPP, data$Miami_NPP)$RMSE
  )
)

print(results)

# --- Plot results --- 

plot_pred_obs_serio <- function(obs, pred, model_name) {
  df <- data.frame(Observed = obs, Predicted = pred)
  lm_fit <- lm(Predicted ~ Observed, data = df)
  r_val <- cor(obs, pred, use = "complete.obs")
  r2_val <- summary(lm_fit)$r.squared
  rmse_val <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  
  ggplot(df, aes(x = Observed, y = Predicted)) +
    geom_point(color = "steelblue", alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", color = "black", fill = "grey", alpha = 0.4, se = TRUE) +
    labs(
      title = model_name,
      subtitle = sprintf("r = %.3f | r² = %.3f | RMSE = %.2f", r_val, r2_val, rmse_val),
      x = "Observed NPP",
      y = "Predicted NPP"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(face = "bold")
    )
}

p1 <- plot_pred_obs_serio(data$NPP, data$Predicted_NPP, "brGDGT")
p2 <- plot_pred_obs_serio(data$NPP, data$MOD17A3_NPP, "MOD17A3")
p3 <- plot_pred_obs_serio(data$NPP, data$Miami_NPP, "Miami")

grid.arrange(p1, p2, p3, ncol = 1)





### d) PREDICT NPP IN PADUL ####

# --- Load data ---

Padul_df <- read.xlsx("Dataset.xlsx", sheet = "Padul_GDGTs")
head(Padul_df)


# --- Select compositional columns ---
# Columns used in model training
vars <- colnames(Padul_df)[3:17]  
# Check
vars
# Extra variables required by the model
vars_extra <- c("Sampletype")  

comp_new <- Padul_df[, vars]


# --- Convert percentages to proportions ---
# This only necessary in case the specific percentages are extressed in a 0-100 range
# comp_new <- comp_new / 100 

# ---  Avoid zeros for CLR transformation --- 
# comp_new <- comp_new + 1e-6


# --- Apply CLR transformation ---

clr_new <- clr(comp_new)
clr_new_df <- as.data.frame(clr_new)

# --- Add extra variables ---

clr_new_df$Sampletype <- "Lacustrine Sediment"
# Convert Sampletype to factor with same levels as training
clr_new_df$Sampletype <- factor(clr_new_df$Sampletype,
                                levels = c(levels(model_cv$trainingData$Sampletype), 
                                           "Lacustrine Sediment"))

# --- Build dataset for prediction ---

newdata_model <- clr_new_df

# Keep only CLR columns for the model
newdata_model <- clr_new_df[, 1:15]  

# Create dummy variables for Sampletype
sampletypes <- c("Lacustrine Sediment", "Lacustrine SPM", 
                 "Low DO Lacustrine SPM", "Peat", 
                 "Riverine Sediment and SPM", "Soil")

for (stype in sampletypes) {
  colname <- paste0("Sampletype", stype)
  newdata_model[[colname]] <- ifelse(clr_new_df$Sampletype == stype, 1, 0)
}


# --- Predict NPP using all trees (100) ---

rf_preds <- predict(model_cv$finalModel, newdata = newdata_model, predict.all = TRUE)


# --- Step 9: Calculate mean, SD, and 95% CI ---

Padul_df$NPP_predicted <- apply(rf_preds$individual, 1, mean)
Padul_df$NPP_sd        <- apply(rf_preds$individual, 1, sd)
Padul_df$NPP_low       <- Padul_df$NPP_predicted - 1.96 * Padul_df$NPP_sd
Padul_df$NPP_high      <- Padul_df$NPP_predicted + 1.96 * Padul_df$NPP_sd


# --- Prepare data for plotting

rf_df <- as.data.frame(rf_preds$individual)
rf_df$Age <- Padul_df$`Age.(yr.cal.BP)`  # Add age column
rf_long <- melt(rf_df, id.vars = "Age")  # Convert to long format for ggplot


# --- Plot individual trees, mean, and 95% CI ---

Padul_NPP_brGDGTs<-ggplot() +
  # Individual tree lines
  geom_line(data = rf_long, aes(x = Age, y = value, group = variable),
            alpha = 0.2, color = "gray") +
  # Mean prediction line
  geom_line(data = Padul_df, aes(x = `Age.(yr.cal.BP)`, y = NPP_predicted),
            color = "black", linewidth = 1) +
  # 95% confidence interval
  geom_line(data = Padul_df, aes(x = `Age.(yr.cal.BP)`, y = NPP_low),
            color = "orange", linewidth = 0.8) +
  geom_line(data = Padul_df, aes(x = `Age.(yr.cal.BP)`, y = NPP_high),
            color = "orange", linewidth = 0.8) +
  scale_x_reverse(limits = c(35000, 8000), 
                  breaks = seq(35000, 8000, by = -2000)) +
  labs(x = "Age (yr cal BP)", y = "NPP",
       title = "brGDGT-NPP") +
  theme_minimal()


Padul_NPP_brGDGTs


### e) CLIMATE-POTENTIAL NPP VS GDGT-DERIVED NPP ####


# --- Load data ---
Padul_Miami_climate <- read.xlsx("Dataset.xlsx", sheet = "Padul_climate")
head(Padul_Miami_climate)



# --- Propagate uncertainty from MAT and MAP 95% CI to NPP 95% CI --- 
# Parameters
n_iter <- 1000
n_samples <- nrow(Padul_df)

# Sotre NPP outputs
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
  geom_line(aes(y = NPP_mean), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = NPP_low, ymax = NPP_high), fill = "blue", alpha = 0.3) +
  scale_x_reverse(limits = c(35000, 8000), breaks = seq(35000, 8000, by = -5000)) +
  labs(
    x = "Age (yr cal BP)",
    y = "NPP (Miami)",
    title = "NPP de Miami with 95% CI"
  ) +
  theme_minimal()




# --- Human Appropriation of NPP (HANPP) --- 

# Interpolate NPP_Miami for the Age values in rf_long
Padul_Miami_interp <- approx(x = Padul_Miami_climate$Age,
                             y = Padul_Miami_climate$NPP_Miami,
                             xout = rf_long$Age)$y
rf_with_miami <- rf_long %>%
  mutate(NPP_Miami = Padul_Miami_interp,
         NPP_diff = NPP_Miami - value)


# Plot difference
HANPP <- ggplot() +
  # Líneas de cada árbol
  geom_line(data = rf_with_miami, aes(x = Age, y = NPP_diff, group = variable),
            alpha = 0.2, color = "gray") +
  # Línea de la media de la diferencia
  stat_summary(data = rf_with_miami, aes(x = Age, y = NPP_diff, group = 1),
               fun = mean, geom = "line", color = "red", linewidth = 1) +
  labs(x = "Age (yr cal BP)", y = "NPP Difference (Miami - RF)",
       title = "") +
  scale_x_reverse(limits = c(35000, 8000), 
                  breaks = seq(35000, 8000, by = -5000)) +
  theme_minimal()


HANPP




Padul_NPP_brGDGTs<-ggplot() +
  # Individual tree lines
  geom_line(data = rf_long, aes(x = Age, y = value, group = variable),
            alpha = 0.2, color = "gray") +
  # Mean prediction line
  geom_line(data = Padul_df, aes(x = `Age.(yr.cal.BP)`, y = NPP_predicted),
            color = "black", linewidth = 1) +
  # 95% confidence interval
  geom_line(data = Padul_df, aes(x = `Age.(yr.cal.BP)`, y = NPP_low),
            color = "orange", linewidth = 0.8) +
  geom_line(data = Padul_df, aes(x = `Age.(yr.cal.BP)`, y = NPP_high),
            color = "orange", linewidth = 0.8) +
  geom_line(data = Padul_Miami_climate, aes(x = Age, y = NPP_Miami),
            color = "blue", linewidth = 0.8, linetype = "dashed" ) +
  scale_x_reverse(limits = c(35000, 8000), 
                  breaks = seq(35000, 8000, by = -5000)) +
  # Miami climate NPP line

  labs(x = "Age (yr cal BP)", y = "NPP",
       title = " ") +
  theme_minimal()




grid.arrange(Padul_NPP_brGDGTs, HANPP, ncol=1)




### f) SUMMED PROBABILITY DISTRIBUTIONS ####

# ---- SPD by culture ----

# Delta values for marine samples
deltaR_val  <- 94
deltaR_sd   <- 61

# List of cultures
cultures <- c("Aurignacian", "Gravettian", "Solutrean", "Magdalenian", "Epipalaeolithic")

# Load data
SPD_df <- read.xlsx("Dataset.xlsx", sheet = "SPD", detectDates = FALSE)

# Dataframe to store outputs
Regions_SPD <- data.frame()

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
    h     = 300 # This was checked with bins of 200 and 500, with nearly identical results.
  )
  
  # Compute SPD
  spd_cult <- spd(caldates_ind, timeRange = c(35000, 6000), bins = bins)
  
  # Save SPD for this culture
  PrDens_df <- data.frame(
    PrDens  = spd_cult$grid$PrDens,
    Age     = spd_cult$grid$calBP,
    Culture = cult
  )
  
  Regions_SPD <- dplyr::bind_rows(Regions_SPD, PrDens_df)
}

# Save CSV / Optional
 # write.csv(Regions_SPD, "Regions_SPD_by_culture_combined_by_layer.csv", row.names = FALSE)

# Plot SPD
SPD_l <- ggplot(Regions_SPD, aes(x = Age, y = PrDens, fill = Culture)) +
  geom_area(alpha = 0.5, position = "identity", color = "black") +
  scale_x_reverse(breaks = seq(35000, 8000, by = -5000)) +
  theme_classic() +
  labs(x = "cal BP", y = "Probability Density", title = "Summed Probability Distribution by Culture") +
  scale_fill_brewer(palette = "Set2")

SPD_l


#NGRIP data for plot
NGRIP_df <- read.xlsx("Dataset.xlsx", sheet = "NGRIP", detectDates = FALSE)
head(NGRIP_df)

# Create df with events (for visual purposes)
events <- data.frame(
  Event = c("H3", "H2", "LGM", "H1",  "Younger Dryas"),
  start = c(31000, 24500, 21500, 17000,  12900),
  end   = c(29000, 23000, 19000, 15000, 11700)
)

# Plot
NGRIP_plot <- ggplot() +
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4) +
  geom_line(data = NGRIP_df, aes(x = Age * 1000, y = d18O),
            color = "blue", size = 1) +
  scale_x_reverse(breaks = seq(35000, 6000, by = -5000)) +
  labs(
    x = "Age (yr BP)",
    y = expression(delta^{18}*O~("\u2030")),
    title = " "
  ) +
  theme_minimal(base_size = 14)

NGRIP_plot


plot_grid(NGRIP_plot, Padul_NPP_brGDGTs, HANPP,SPD , ncol = 1, align = "v")






# ---- Pollen plot ----
# Load data
pollen_df <- read.xlsx("Dataset.xlsx", sheet = "Pollen", detectDates = FALSE)
head(pollen_df)

# Long formant 
pollen_long <- pollen_df %>%
  pivot_longer(
    cols = c(Medit_1505, Medit_III),
    names_to = "Record",
    values_to = "Medit_pollen"
  ) %>%
  mutate(Age = ifelse(Record == "Medit_1505", Age_1505, Age_III))

range_1505 <- range(pollen_df$Medit_1505, na.rm = TRUE)
range_III   <- range(pollen_df$Medit_III, na.rm = TRUE)
scale_factor <- diff(range_1505) / diff(range_III)

Medit_polllen <- ggplot() +
  geom_line(data = pollen_df, aes(x = Age_1505, y = Medit_1505, color = "1505")) +
  geom_line(data = pollen_df, aes(x = Age_III, y = Medit_III, color = "III")) +
  scale_y_continuous(
    name = "Mediterranean pollen (%)",
    sec.axis = sec_axis(~ ., name = "Mediterranean pollen (%)")  # misma escala
  ) +
  scale_x_reverse() +
  labs(x = "Age (years BP)", color = "Core") +
  theme_minimal()



# The same for the Precipitation Index
ip_long <- pollen_df %>%
  pivot_longer(
    cols = c(Ip_1505, Ip_III),
    names_to = "Record",
    values_to = "Ip"
  ) %>%
  mutate(Age = ifelse(Record == "Ip_1505", Age_1505, Age_III))

range_1505 <- range(pollen_df$Ip_1505, na.rm = TRUE)
range_III <- range(pollen_df$Ip_III, na.rm = TRUE)
scale_factor <- diff(range_1505) / diff(range_III)

lp<- ggplot() +
  geom_line(data = pollen_df, aes(x = Age_1505, y = Ip_1505, color = "1505")) +
  geom_line(data = pollen_df, aes(x = Age_III, y = Ip_III * scale_factor, color = "III")) +
  scale_y_continuous(
    name = "Ip 1505",
    sec.axis = sec_axis(~ . / scale_factor, name = "Ip III")
  ) +
  scale_x_reverse() +
  labs(x = "Age (years BP)", color = "Core") +
  theme_minimal()


# Final plot:

plot_grid(NGRIP_plot, Medit_polllen, lp, Padul_NPP_brGDGTs, HANPP,SPD , ncol = 1, align = "v")









