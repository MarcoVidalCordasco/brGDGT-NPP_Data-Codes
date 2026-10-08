# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# SCRIPT 1 — brGDGT-based NPP analyses and reconstructions
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#
# This script reproduces the analyses and figures associated
# with the brGDGT–NPP calibration, model validation, and NPP
# reconstruction presented in the manuscript.
#
# CONTENTS:
#
# 1. SUMMARY STATISTICS and EXPLORATORY PLOTS
#    Main Figures 1, 2, and 3A–C
#    Supplementary Tables 1 and 2
#    Supplementary Figure 3
#    Note: Supplementary Figures 1 and 4 were not generated in R.
#          Supplementary Figure 2 is included in the "Random Forest NPP model" section.
#          
# 2. RANDOM FOREST NPP MODEL 
#    Main Figure 4
#    Main Table 1
#    Supplementary Figures 2 and 5
#
# 3. INDEPENDENT MODEL VALIDATION
#    Main Figure 5
#    Supplementary Figure 6
#    
# 4. ENVIRONMENTAL NOVELTY
#    Supplementary Figure 7
#    
# 5. PADUL PALEOCLIMATE RECONSTRUCTION
#    Supplementary Figure 8
#
# 6. PADUL NPP RECONSTRUCTION
#    Main Figure 6
#    Supplementary Table 5
#    Supplementary Figures 9 and 13
#    12?
#
# 7. ARCHAEOLOGICAL ANALYSES
#    Main Figures 7–8
#    Supplementary Table 3 and 4
#    Supplementary Figures 10-11
#
# Additional herbivore biomass and consumer carrying-capacity
# analyses are provided in Script_2.R.
#

rm(list = ls()) # Clear all
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# REQUIRED PACKAGES
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 
library(terra)
library(dplyr)
library(tidyr)
library(RANN)
library(ggplot2)
library(maps)
library(scales)
library(openxlsx)
library(truncnorm)
library(vegan)
library(ggrepel)
library(grid)
library(compositions)
library(splines)
library(purrr)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(analogue)
library(corrplot)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(caret)
library(ggExtra)
library(randomForest)
library(reshape2)
library(rcarbon)
library(mgcv)
library(sf)
library(blockCV)
library(spdep)
library(ppcor)

# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##  1) SUMMARY STATISTICS & EXPLORATORY PLOTS ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# --- Read the data ---
data <- read.xlsx("Dataset_1.xlsx", sheet = "Dataset_R")


# --- Define compound groups ---
tetramethylated <- c("fIa", "fIb", "fIc")
pentamethylated <- c("fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_")
hexamethylated  <- c("fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_")

# --- Sum groups ---
data$Tetramethylated <- rowSums(data[, tetramethylated], na.rm = TRUE)
data$Pentamethylated <- rowSums(data[, pentamethylated], na.rm = TRUE)
data$Hexamethylated  <- rowSums(data[, hexamethylated], na.rm = TRUE)

# --- Compute CBT ---
 
data$CBT <- log10((data$fIc + data$fIIa_+data$fIIb_+data$fIIc_+data$fIIIa_+data$fIIIb_+data$fIIIc_) / 
                    (data$fIa + data$fIIa + data$fIIIa))

# --- Compute MBT ---
all_compounds <- c("fIa", "fIb", "fIc",
                   "fIIa", "fIIb", "fIIc",
                   "fIIIa")

data$MBT <- (data$fIa + data$fIb + data$fIc) / rowSums(data[, all_compounds], na.rm = TRUE)


# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
############### Main Figure 1 #############

# Plot specific individual compounds and NPP

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
  }, silent = TRUE)
  
  # 2. Linear model
  tryCatch({
    lm_fit <- lm(Value ~ NPP, data = datf)
    lm_summary <- summary(lm_fit)
    results$linear_r2 <- lm_summary$r.squared
    results$linear_p <- ifelse(nrow(coef(lm_summary)) >= 2, 
                               coef(lm_summary)[2, 4], NA)
  }, silent = TRUE)
  
  # 3. Quadratic model (U shape)
  tryCatch({
    quad_fit <- lm(Value ~ NPP + I(NPP^2), data = datf)
    quad_summary <- summary(quad_fit)
    results$quad_r2 <- quad_summary$r.squared
    results$quad_p <- ifelse(nrow(coef(quad_summary)) >= 3, 
                             coef(quad_summary)[3, 4], NA)
  }, silent = TRUE)
  
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
  }, silent = TRUE)
  
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
    }, silent = TRUE)
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
    }, silent = TRUE)
  }
  
  # 3. GAM model 
  # 3. GAM model
  if (!is.na(stats$gam_dev_expl) && stats$n >= 8) {
    try({
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
          linewidth = 1.5,
          linetype = "solid",
          alpha = 0.6
        )
      }
    }, silent = TRUE)
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
plots_tetra
plots_penta
plots_hexa

# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#### Supplementary Table 1 ####

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




#--- Shapiro Wilk test for normality of the variables ---

vars <- c(
  "NPP",
  "Tetramethylated", "Pentamethylated", "Hexamethylated",
  "MBT", "CBT", "PET",
  "BIO01", "BIO12",
  "pH", "Bulk", "Clay", "CEC",
  "Nitrogen", "Phosphorus",
  "SOC_0-5", "SOC_5-15",
  "Elevation"
)

# Names for the figure/table
labels <- c(
  BIO01 = "MAT",
  BIO12 = "MAP",
  MBT = "MBT'5ME",
  CBT = "CBT'"
)

# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#### Supplementary Table 2 ####

shapiro_table <- purrr::map_dfr(vars, function(v) {
  
  x <- data[[v]]
  x <- x[is.finite(x)]
  
  test <- shapiro.test(x)
  
  tibble(
    Variable = ifelse(v %in% names(labels), labels[v], v),
    W = round(unname(test$statistic), 3),
    `p-value` = ifelse(
      test$p.value < 0.001,
      "<0.001",
      sprintf("%.3f", test$p.value)
    )
  )
})

shapiro_table

# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#### Supplementary Figure 3 ####
# --- QQ plots ---

qq_data <- data %>%
  dplyr::select(all_of(vars)) %>%
  rename(
    MAT = BIO01,
    MAP = BIO12,
    `MBT'5ME` = MBT,
    `CBT'` = CBT
  ) %>%
  pivot_longer(
    everything(),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  drop_na()

ggplot(qq_data, aes(sample = Value)) +
  stat_qq(size = 1.3) +
  stat_qq_line(linewidth = 0.7) +
  facet_wrap(~ Variable, scales = "free", ncol = 4) +
  labs(
    x = "Theoretical quantiles",
    y = "Observed quantiles"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 11)
  )

# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
### Main Figure 2 ####

# --- Select variables for Spearman correlations, including CBT and MBT ---
selected_vars <- c("NPP","Tetramethylated", "Pentamethylated", "Hexamethylated", "MBT" , "CBT",
                   "PET","BIO01", "BIO12", "pH", 
                    "Bulk", "Clay", "CEC", "Nitrogen", "Phosphorus", "SOC_0-5", "SOC_5-15",
                   "Elevation" )


# --- Subset numeric ---
cor_data <- data[, selected_vars]
cor_data <- as.data.frame(lapply(cor_data, as.numeric))

# --- Compute correlation matrix ---
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs", method = "spearman")

# Change variable names for better visualization
colnames(cor_matrix)[colnames(cor_matrix) == "NPP"] <- "NPP"
colnames(cor_matrix)[colnames(cor_matrix) == "BIO01"] <- "MAT"
colnames(cor_matrix)[colnames(cor_matrix) == "BIO12"] <- "MAP"

rownames(cor_matrix) <- colnames(cor_matrix)  

# --- Plot ---

corrplot(
  cor_matrix,
  method = "square",
  type = "upper",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 1.1,
  number.cex = 0.85,
  cl.pos = "r",
  
  col = colorRampPalette(c(
    "firebrick3",
    "mistyrose",
    "white",
    "lightblue",
    "steelblue4"
  ))(200)
)


# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
### Main Figure 3A-C ####


# --- Variables used in partial correlations ---

br_vars <- c(
  "Tetramethylated",
  "Pentamethylated",
  "Hexamethylated"
)


# --- Function for partial correlation ---


partial_spearman <- function(var, controls) {
  
  df <- data %>%
    dplyr::select(all_of(c("NPP", var, controls))) %>%
    na.omit()
  
  partial <- ppcor::pcor.test(
    x = df$NPP,
    y = df[[var]],
    z = df[, controls, drop = FALSE],
    method = "spearman"
  )
  
  data.frame(
    Variable = var,
    rho = partial$estimate,
    p = partial$p.value
  )
}


# --- Zero order correlation ---


zero_order <- bind_rows(
  lapply(br_vars, function(var) {
    
    df <- data %>%
      dplyr::select(NPP, all_of(var)) %>%
      na.omit()
    
    x <- cor.test(
      df$NPP,
      df[[var]],
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      Variable = var,
      rho = unname(x$estimate),
      p = x$p.value
    )
  })
)


# --- Partial correlations ---


partial_SOC <- bind_rows(
  lapply(
    br_vars,
    partial_spearman,
    controls = c("SOC_0-5", "SOC_5-15")
  )
)

partial_pH <- bind_rows(
  lapply(
    br_vars,
    partial_spearman,
    controls = "pH"
  )
)

partial_Clay <- bind_rows(
  lapply(
    br_vars,
    partial_spearman,
    controls = "Clay"
  )
)

partial_Elevation <- bind_rows(
  lapply(
    br_vars,
    partial_spearman,
    controls = "Elevation"
  )
)

partial_MAT <- bind_rows(
  lapply(
    br_vars,
    partial_spearman,
    controls = "BIO01"
  )
)

partial_MAP <- bind_rows(
  lapply(
    br_vars,
    partial_spearman,
    controls = "BIO12"
  )
)


# --- Plot ---


plot_df <- bind_rows(
  
  zero_order %>%
    mutate(Adjustment = "Zero-order"),
  
  partial_SOC %>%
    mutate(Adjustment = "SOC"),
  
  partial_pH %>%
    mutate(Adjustment = "pH"),
  
  partial_Clay %>%
    mutate(Adjustment = "Clay"),
  
  partial_Elevation %>%
    mutate(Adjustment = "Elevation"),
  
  partial_MAT %>%
    mutate(Adjustment = "MAT"),
  
  partial_MAP %>%
    mutate(Adjustment = "MAP")
)



plot_df <- plot_df %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

# --- Order ---

plot_df$Adjustment <- factor(
  plot_df$Adjustment,
  levels = c(
    "Zero-order",
    "SOC",
    "pH",
    "Clay",
    "Elevation",
    "MAT",
    "MAP"
  )
)

plot_df$Variable <- factor(
  plot_df$Variable,
  levels = c(
    "Tetramethylated",
    "Hexamethylated",
    "Pentamethylated"
  )
)



partial <- ggplot(
  plot_df,
  aes(
    x = Adjustment,
    y = rho
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_segment(
    aes(
      x = Adjustment,
      xend = Adjustment,
      y = 0,
      yend = rho
    ),
    linewidth = 0.5
  ) +
  
  geom_point(
    size = 3.2
  ) +
  
  geom_text(
    aes(
      label = sig,
      y = rho + ifelse(rho >= 0, 0.045, -0.045)
    ),
    size = 4.2
  ) +
  
  facet_wrap(
    ~ Variable,
    nrow = 1
  ) +
  
  scale_y_continuous(
    limits = c(-0.82, 0.82),
    breaks = seq(-0.8, 0.8, 0.4)
  ) +
  
  labs(
    x = NULL,
    y = expression("Spearman " * rho)
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    panel.spacing = unit(1, "lines")
  )

partial




# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
### Main Figure 3B 

# --- RDA ---

# --- Define brGDGT compounds ---
brgdgt_vars <- c(
  "fIa", "fIb", "fIc",
  "fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_",
  "fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_"
)
# --- Select and rename explanatory variables ---
rda_data <- data %>%
  dplyr::select(
    all_of(brgdgt_vars),
    NPP,
    BIO01,
    BIO12,
    `SOC_0-5`,
    `SOC_5-15`,
    pH,
    Clay,
    Elevation
  ) %>%
  rename(
    MAT = BIO01,
    MAP = BIO12,
    SOC05 = `SOC_0-5`,
    SOC515 = `SOC_5-15`
  ) %>%
  filter(complete.cases(.))

br_clr <- as.matrix(
  rda_data[, brgdgt_vars]
)

# --- RDA with all 15 GDGT compounds ---

m_full <- vegan::rda(
  br_clr ~
    NPP +
    MAT +
    MAP +
    pH +
    SOC05 +
    SOC515 +
    Clay,
  data = rda_data
)

# --- Correlation coefficient -R2- --- 
vegan::RsquareAdj(m_full)

# --- Global significance ---
set.seed(123)
anova(
  m_full,
  permutations = 999
)

# --- Marginal significance of each predictor---
set.seed(123)
anova(
  m_full,
  by = "margin",
  permutations = 999
)

# Plot Figure 3B:


mod <- m_full 

# --- Scaling for visual purposes ---
scaling_use <- 2

# --- Get scores ---

site_scores <- as.data.frame(
  scores(mod, display = "sites", choices = 1:2, scaling = scaling_use)
)

s_scores <- as.data.frame(
  scores(mod, display = "species", choices = 1:2, scaling = scaling_use)
)

env_scores <- as.data.frame(
  scores(mod, display = "bp", choices = 1:2, scaling = scaling_use)
)

# --- Rename columns ---
colnames(site_scores)[1:2] <- c("Axis1", "Axis2")
colnames(s_scores)[1:2] <- c("Axis1", "Axis2")
colnames(env_scores)[1:2] <- c("Axis1", "Axis2")

# --- Labels ---
s_scores$Compound <- rownames(s_scores)
env_scores$label <- rownames(env_scores)

site_range_x <- diff(range(site_scores$Axis1, na.rm = TRUE))
site_range_y <- diff(range(site_scores$Axis2, na.rm = TRUE))

env_range_x <- diff(range(env_scores$Axis1, na.rm = TRUE))
env_range_y <- diff(range(env_scores$Axis2, na.rm = TRUE))

sp_range_x <- diff(range(s_scores$Axis1, na.rm = TRUE))
sp_range_y <- diff(range(s_scores$Axis2, na.rm = TRUE))

env_mult <- 0.85 * min(site_range_x / env_range_x,
                       site_range_y / env_range_y)

species_mult <- 0.35 * min(site_range_x / sp_range_x,
                           site_range_y / sp_range_y)

env_scores_plot <- env_scores %>%
  mutate(
    Axis1_plot = Axis1 * env_mult,
    Axis2_plot = Axis2 * env_mult
  )

s_scores_plot <- s_scores %>%
  mutate(
    Axis1_plot = Axis1 * species_mult,
    Axis2_plot = Axis2 * species_mult
  )


if (!is.null(mod$CCA) && length(mod$CCA$eig) >= 2) {
  axis_percent <- mod$CCA$eig[1:2] / sum(mod$CCA$eig) * 100
  
  x_lab <- paste0(
    "RDA1 (", round(axis_percent[1], 1), "% of constrained variance)"
  )
  
  y_lab <- paste0(
    "RDA2 (", round(axis_percent[2], 1), "% of constrained variance)"
  )
} else {
  x_lab <- "RDA1"
  y_lab <- "RDA2"
}


p_rda <- ggplot() +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    linetype = "dashed",
    colour = "grey60"
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.3,
    linetype = "dashed",
    colour = "grey60"
  ) +
  
  geom_point(
    data = site_scores,
    aes(x = Axis1, y = Axis2),
    size = 1.0,
    alpha = 0.22,
    shape = 21,
    fill = "grey75",
    colour = "grey40",
    stroke = 0.2
  ) +
  
  geom_segment(
    data = env_scores_plot,
    aes(
      x = 0, y = 0,
      xend = Axis1_plot,
      yend = Axis2_plot
    ),
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    linewidth = 0.8,
    colour = "black"
  ) +
  
  geom_text_repel(
    data = env_scores_plot,
    aes(
      x = Axis1_plot,
      y = Axis2_plot,
      label = label
    ),
    size = 3.5,
    colour = "black",
    segment.color = NA,
    max.overlaps = Inf
  ) +
  
  geom_segment(
    data = s_scores_plot,
    aes(
      x = 0, y = 0,
      xend = Axis1_plot,
      yend = Axis2_plot
    ),
    arrow = arrow(length = unit(0.08, "cm")),
    linewidth = 0.3,
    alpha = 0.55,
    colour = "red"
  ) +
  
  geom_text_repel(
    data = s_scores_plot,
    aes(
      x = Axis1_plot,
      y = Axis2_plot,
      label = Compound
    ),
    size = 2.4,
    colour = "blue",
    segment.color = NA,
    max.overlaps = Inf
  ) +
  
  labs(
    x = x_lab,
    y = y_lab
  ) +
  
  coord_equal() +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.margin = ggplot2::margin(
      t = 10, r = 20, b = 10, l = 10
    )
  )


p_rda

# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
### Main Figure 3C 

# Lambda 1 / Lambda 2 with all the 15 compounds
br_comp <- as.matrix(
  rda_data[, brgdgt_vars]
)

env_vars <- c("MAT", "NPP", "MAP", "pH",
              "Clay", "SOC515", "SOC05")

groups <- list(
  Tetra = c("fIa", "fIb", "fIc"),
  Penta = c("fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_"),
  Hexa  = c("fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_")
)

# ---Response matrices---
br_raw <- as.matrix(rda_data[, brgdgt_vars])

br_groups <- sapply(groups, function(x) {
  rowSums(rda_data[, x, drop = FALSE])
})

# ---Function to calculate lambda1/lambda2---
fit_rda <- function(Y, var) {
  
  mod <- rda(Y, rda_data[, var])
  
  lambda1 <- eigenvals(mod, model = "constrained")[1]
  lambda2 <- eigenvals(mod, model = "unconstrained")[1]
  
  tibble(
    Variable = var,
    Ratio = unname(lambda1 / lambda2),
    R2_adj = RsquareAdj(mod)$adj.r.squared
  )
}

# ---Run RDA for both compositions---
results <- imap_dfr(
  list("Individual brGDGTs" = br_raw,
       "Grouped brGDGTs" = br_groups),
  function(Y, name) {
    map_dfr(env_vars, ~ fit_rda(Y, .x)) %>%
      mutate(Composition = name)
  }
)

# ---Format variable labels and order---
results <- results %>%
  mutate(
    Variable = recode(
      Variable,
      SOC515 = "SOC 5–15 cm",
      SOC05 = "SOC 0–5 cm"
    ),
    Variable = factor(
      Variable,
      levels = rev(c(
        "MAT", "NPP", "MAP", "pH",
        "Clay", "SOC 5–15 cm",
        "SOC 0–5 cm"
      ))
    )
  )

# --- Plot Figure 3C ---

ggplot(results, aes(x = Ratio, y = Variable, fill = Composition)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.28,
    color = NA
  ) +
  geom_text(
    aes(label = ifelse(Ratio == 0, "0", sprintf("%.2f", Ratio))),
    position = position_dodge(width = 0.72),
    hjust = -0.15,
    size = 3
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    color = "#2E4DA7"
  ) +
  scale_fill_manual(
    values = c(
      "Grouped brGDGTs" = "#7A0000",   # dark red
      "Individual brGDGTs" = "grey55"  # grey
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1.35),
    breaks = c(0.0, 0.5, 1.0),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    x = expression(lambda[1] / lambda[2]),
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    panel.background = element_rect(fill = "grey92", color = NA),
    plot.background = element_rect(fill = "grey92", color = NA)
  )



# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 2) RANDOM FOREST NPP MODEL ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# --- Load data ---
data <- read.xlsx("Dataset_1.xlsx", sheet = "Dataset_R")
#Check
head(data)

# --- Select compounds ---
vars <- colnames(data)[12:26]

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


## RANDOM-FOREST ACROSS DIFERENT SAMPLE TYPES

# --- . Spearman’s rank correlation coefficients between the CRL-transformed relative abundance of individual brGDGT 
# compounds and NPP are shown for lacustrine sediments, peat, and soils---
unique(model_data$Sampletype)

gdgt <- colnames(clr_data)

cors_long <- model_data %>%
  filter(Sampletype_fixed %in% c("Soil", "Peat", "Lacustrine Sediment")) %>%
  group_by(Sampletype_fixed) %>%
  summarise(across(all_of(gdgt),
                   ~cor(.x, NPP, method = "spearman",
                        use = "complete.obs"))) %>%
  pivot_longer(-Sampletype_fixed,
               names_to = "brGDGT",
               values_to = "rho")
### Supplementary Figure 2 ####

ggplot(cors_long,
       aes(x = brGDGT, y = rho,
           group = Sampletype_fixed,
           colour = Sampletype_fixed)) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic() +
  labs(x = NULL,
       y = "Spearman rho with NPP",
       colour = "Sample type")


## --- Cross-environment Random Forest ---

#Cross-sample-type performance of the brGDGT-based NPP calibration. Model performance when trained on one sample type 
# and tested on another is summarized using the correlation coefficient (r) between observed and predicted NPP and the 
# root mean square error (RMSE):

types <- c("Soil",
           "Peat",
           "Lacustrine Sediment",
           "Aquatic_SPM_combined")

results <- expand.grid(
  Train = types,
  Test = types,
  stringsAsFactors = FALSE
) %>%
  filter(Train != Test)

results <- results %>%
  rowwise() %>%
  mutate(
    metrics = list({
      
      train <- model_data %>% filter(Sampletype_fixed == Train)
      test  <- model_data %>% filter(Sampletype_fixed == Test)
      
      rf <- randomForest(
        x = train[, gdgt],
        y = train$NPP
      )
      
      pred <- predict(rf, test[, gdgt])
      
      data.frame(
        r = cor(pred, test$NPP),
        RMSE = sqrt(mean((pred - test$NPP)^2))
      )
    })
  ) %>%
  unnest(metrics)

results



results$label <- paste0(
  "r = ", sprintf("%.2f", results$r),
  "\nRMSE = ", round(results$RMSE)
)

results$Test <- gsub("Aquatic_SPM_combined", "Aquatic SPM", results$Test)
results$Train <- gsub("Aquatic_SPM_combined", "Aquatic SPM", results$Train)

#### Supplementary Figure 5 ####
# --- Plot cross-environment performance ---

ggplot(results, aes(x = Test, y = Train, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient2(
    low = "white",
    mid = "lightblue",
    high = "steelblue4",
    midpoint = 0.5,
    limits = c(0, 0.8)
  ) +
  theme_classic() +
  labs(
    x = "Test sample type",
    y = "Training sample type",
    fill = "Correlation (r)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# --- RF model with all data --- ####
# All data is used for training and for testing

model_all_data <- train(
  formula,
  data = model_data,
  method = "rf",
  trControl = trainControl(method = "none"),  
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


#### Main Table 1 ####
#### 
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


#### Main Figure 4 ####
grid.arrange(p1,auto_r1$plot, p2,auto_r2$plot, p50,auto_r50$plot, p100,auto_r100$plot,
             p500,auto_r500$plot, p1000, auto_r1000$plot, ncol=4)





# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#   3) INDEPENDENT MODEL VALIDATION ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# --- Load data ---
m_NPP_df <- read.xlsx("Dataset_1.xlsx", sheet = "Measured_NPP_brGDGTs")

# --- Define fractions ---
tetramethylated <- c("fIa", "fIb", "fIc")
pentamethylated <- c("fIIa", "fIIa_", "fIIb", "fIIb_", "fIIc", "fIIc_")
hexamethylated  <- c("fIIIa", "fIIIa_", "fIIIb", "fIIIb_", "fIIIc", "fIIIc_")
all_fracs <- c(tetramethylated, pentamethylated, hexamethylated)

# --- Transform proportions and avoid 0 ---
frac_data <- m_NPP_df %>%
  dplyr::select(all_of(all_fracs)) %>%
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

####  Main Figure 5 #### 
# --- Plot results --- 

plot_data <- data %>%
  dplyr::select(
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


#### Supplementary Figure 6 ####
###
# --- Formula ---
formula <- as.formula(paste("NPP ~", paste(c(colnames(clr_data)), collapse = " + ")))
# Check
formula
# Now, for Supplementary Figure 6:
# 1. Load the same models trained without Sample type as predictor (or re-run the models without Sample type):
model_cv <- readRDS(paste0("WithoutSampleType/", "model_brGDGT-NPP_cv1.rds"))
model_spatial_cv_50<- readRDS(paste0("WithoutSampleType/", "spatial_cv_50km_models1.rds"))
model_spatial_cv_100<- readRDS(paste0("WithoutSampleType/", "spatial_cv_100km_models1.rds"))
model_spatial_cv_500<- readRDS(paste0("WithoutSampleType/", "spatial_cv_500km_models1.rds"))
model_spatial_cv_1000<- readRDS(paste0("WithoutSampleType/", "spatial_cv_1000km_models1.rds"))
# 2. Now re-run the independent validation chunk (Section 3)

# To avoid any errors, re-load the model that includes sample type as predictor:
model_cv <- readRDS(paste0("model_brGDGT-NPP_cv1.rds"))

# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 4) ENVIRONMENTAL NOVELTY ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# ---Load global environmental rasters ---
data <- read.xlsx("Dataset_1.xlsx", sheet = "Dataset_R")

# --- Raster files with environmental variables ---
r_elev <- rast("Elevation.tif")
r_bulk <- rast("rbulk_reproj.tif")
r_clay <- rast("rclay_reproj.tif")
r_ph   <- rast("rph_reproj.tif")
r_mat  <- rast("wc2.1_5m_bio_1.tif")
r_map  <- rast("wc2.1_5m_bio_12.tif")


# ---Put all rasters on the same CRS/grid ---

template <- r_mat

r_map  <- project(r_map,  template)
r_ph   <- project(r_ph,   template)
r_clay <- project(r_clay, template)
r_bulk <- project(r_bulk, template)
r_elev <- project(r_elev, template)

env_rast <- c(
  r_mat,
  r_map,
  r_ph,
  r_clay,
  r_bulk,
  r_elev
)

names(env_rast) <- c(
  "BIO01",
  "BIO12",
  "pH",
  "Clay",
  "Bulk",
  "Elevation"
)


# ---Calibration environmental space ---

vars <- names(env_rast)

train <- data %>%
  dplyr::select(all_of(vars)) %>%
  drop_na()

# Compare calibration and global raster ranges
rbind(
  Calibration_min  = sapply(train, min),
  Calibration_max  = sapply(train, max),
  Calibration_mean = sapply(train, mean),
  Raster_min = as.numeric(
    global(env_rast, "min", na.rm = TRUE)[,1]
  ),
  Raster_max = as.numeric(
    global(env_rast, "max", na.rm = TRUE)[,1]
  )
)

# ---Standardise all variables using calibration mean and SD---

mu  <- sapply(train, mean)
sdv <- sapply(train, sd)

train_z <- scale(
  train,
  center = mu,
  scale = sdv
)

env_z <- env_rast

for(i in seq_along(vars)){
  env_z[[i]] <- (env_rast[[i]] - mu[i]) / sdv[i]
}

# ---Convert raster cells to dataframe---

env_df <- as.data.frame(
  env_z,
  xy = TRUE,
  na.rm = TRUE
)

# ---Environmental novelty
#
# Euclidean distance in standardized multivariate environmental
# space to the nearest environmental analogue in the calibration
# dataset.


nn <- nn2(
  data = train_z,
  query = as.matrix(env_df[, vars]),
  k = 1
)

env_df$Novelty <- nn$nn.dists[,1]

# ---Cross-validated calibration-support threshold
#
# Each calibration sample is treated as withheld and compared
# with samples belonging to the remaining folds.

set.seed(123)

K <- 5

fold <- sample(
  rep(1:K, length.out = nrow(train_z))
)

cv_dist <- rep(NA_real_, nrow(train_z))

for(k in 1:K){
  
  ref <- train_z[
    fold != k,
    ,
    drop = FALSE
  ]
  
  test <- train_z[
    fold == k,
    ,
    drop = FALSE
  ]
  
  nn_cv <- nn2(
    data = ref,
    query = test,
    k = 1
  )
  
  cv_dist[fold == k] <- nn_cv$nn.dists[,1]
}

# Environmental-distance distribution under CV
cv_quantiles <- quantile(
  cv_dist,
  probs = c(0.50, 0.75, 0.90, 0.95, 0.99),
  na.rm = TRUE
)

cv_quantiles

# 95% calibration-support threshold
q95 <- quantile(
  cv_dist,
  0.95,
  na.rm = TRUE
)


# --- Novelty relative to calibration support ---
#
# <1 = within the 95% CV environmental-distance envelope
# >1 = environmental novelty exceeding that envelope


env_df$Novelty_relative <- env_df$Novelty / q95

env_df$Extrapolative <- env_df$Novelty_relative > 1

# Percentage of terrestrial/environmentally valid cells
# exceeding the threshold
mean(env_df$Extrapolative, na.rm = TRUE) * 100

# --- Reconstruct raster ---

novelty_map <- rast(
  env_df[, c("x", "y", "Novelty_relative")],
  type = "xyz",
  crs = crs(env_rast)
)

# Supplementary Figure 7 ####
#
# Cap only the DISPLAY at the 99th percentile so a few extreme
# cells do not dominate the colour gradient.
# Original values remain unchanged.

map_df <- env_df %>%
  dplyr::select(
    x,
    y,
    Novelty_relative
  )

display_cap <- quantile(
  map_df$Novelty_relative,
  0.99,
  na.rm = TRUE
)

map_df$Novelty_plot <- pmin(
  map_df$Novelty_relative,
  display_cap
)

world <- map_data("world")

ggplot() +
  
  geom_raster(
    data = map_df,
    aes(
      x = x,
      y = y,
      fill = Novelty_plot
    )
  ) +
  
  geom_polygon(
    data = world,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = NA,
    colour = "grey25",
    linewidth = 0.15
  ) +
  
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, display_cap),
    oob = squish,
    name = "Environmental novelty\n(relative to 95% CV threshold)"
  ) +
  
  coord_quickmap(
    xlim = c(-180, 180),
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_void(base_size = 11) +
  
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 5) PADUL PALEOCLIMATE RECONSTRUCTION ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

Padul_Miami_climate <- read.xlsx("Dataset_1.xlsx", sheet = "Regional_Paleoclimate")

colnames(Padul_Miami_climate)


NGRIP_df <- read.xlsx("Dataset_1.xlsx", sheet = "NGRIP", detectDates = FALSE)
head(NGRIP_df)

# Create df with events (for visual purposes)
events <- data.frame(
  Event = c( "H2", "LGM", "H1",  "Younger Dryas", "8.2"),
  start = c( 24500, 21500, 17000,  12900,  8400),
  end   = c( 23000, 19000, 15000, 11700,  8100)
)

#### Supplementary Figure 8 ####

# A. --- MAT. Pollen (This study) ---

ggplot(Padul_Miami_climate, aes(x = Age_P1, y = MAT_P1)) +
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4, inherit.aes = F) +
  geom_line(color = "black", linewidth = 1) +
  labs(title = "Mean Annual Temperature (MAT) from Pollen (This study)",
       x = "Age (cal BP)",
       y = "MAT (°C)") +
  geom_point(color = "black", size = 2) +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 15)) +
  scale_x_reverse()

# --- B. MAP. Pollen (This study) ---

ggplot(Padul_Miami_climate, aes(x = Age_P1, y = MAP_P1)) +
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4, inherit.aes = F) +
  geom_line(color = "black", linewidth = 1) +
  labs(title = "Mean Annual Precipitation (MAP) from Pollen (This study)",
       x = "Age (cal BP)",
       y = "MAP (mm/yr)") +
  geom_point(color = "black", size = 2) +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 800)) +
  scale_x_reverse()


# --- C. MAT. HadCM3B-M2.1 --- 

ggplot(Padul_Miami_climate, aes(x = Age_GCM, y = MAT_GCM)) +
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4, inherit.aes = F) +
  geom_line(color = "black", linewidth = 1) +
  labs(title = "Mean Annual Temperature (MAT) from HadCM3B-M2.1 ",
       x = "Age (cal BP)",
       y = "MAT (°C)") +
  geom_point(color = "black", size = 2) +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 15)) +
  scale_x_reverse()

# --- D. MAP. HadCM3B-M2.1 ---

ggplot(Padul_Miami_climate, aes(x = Age_GCM, y = MAP_GCM)) +
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4, inherit.aes = F) +
  geom_line(color = "black", linewidth = 1) +
  labs(title = "Mean Annual Precipitation (MAP) from HadCM3B-M2.1",
       x = "Age (cal BP)",
       y = "MAP (mm/yr)") +
  geom_point(color = "black", size = 2) +
  theme_minimal()+
  scale_y_continuous(limits = c(200, 800)) +
  scale_x_reverse()


# --- E. MAT CHELSA-TraCE21k ---

ggplot(Padul_Miami_climate, aes(x = Age_CHELSA, y = MAT_Chelsa)) +
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4, inherit.aes = F) +
  geom_line(color = "black", linewidth = 1) +
  labs(title = "Mean Annual Temperature (MAT) from CHELSA-TraCE21k ",
       x = "Age (cal BP)",
       y = "MAT (°C)") +
  geom_point(color = "black", size = 2) +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 10)) +
  scale_x_reverse()

# --- F. MAP CHELSA-TraCE21k ---

ggplot(Padul_Miami_climate, aes(x = Age_CHELSA, y = MAP_Chelsa)) +
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4, inherit.aes = F) +
  geom_line(color = "black", linewidth = 1) +
  labs(title = "Mean Annual Precipitation (MAP) from CHELSA-TraCE21k",
       x = "Age (cal BP)",
       y = "MAP (mm/yr)") +
  geom_point(color = "black", size = 2) +
  theme_minimal()+
  scale_y_continuous(limits = c(400, 900)) +
  scale_x_reverse()

# --- G. MAT Alboran Sea Temperature---

ggplot(Padul_Miami_climate, aes(x = Age_SST, y = SST)) +
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4, inherit.aes = F) +
  geom_line(color = "black", linewidth = 1) +
  labs(title = "Mean Temperature from Alboran Sea Temperature ",
       x = "Age (cal BP)",
       y = "MAT (°C)") +
  geom_point(color = "black", size = 2) +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 22)) +
  scale_x_reverse()

# --- H. MAP Camuera et al. 2022 ---

ggplot(Padul_Miami_climate, aes(x = Age_P2, y = Camuerta.et.al..2022)) +
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4, inherit.aes = F) +
  geom_line(color = "black", linewidth = 1) +
  labs(title = "Mean Annual Precipitation (MAP) from Camuera et al. 2022",
       x = "Age (cal BP)",
       y = "MAP (mm/yr)") +
  geom_point(color = "black", size = 2) +
  theme_minimal()+
  scale_y_continuous(limits = c(200, 800)) +
  scale_x_reverse()



# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 6) PADUL NPP RECONSTRUCTION ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

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


# --- CALCULATE 95% CONFIDENCE INTERVALS ---

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
  scale_x_reverse(limits = c(29000, 5000),
                  breaks = seq(29000, 5000, by = -2000)) +
  labs(x = "Age (yr cal BP)", y = "NPP",
       title = "brGDGT-NPP") +
  theme_minimal()



print(Padul_NPP_plot_simple)


# Check whether predictions change substantially when Sample type changes:
#  Supplementary Figure 9 ####
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


### Climate-derived (i.e., potential) vs brGDGT (i.e., actual) NPP ####

# --- Load data ---
# To use the pollen-derived data:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
Padul_Miami_climate <- read.xlsx("Dataset_1.xlsx", sheet = "Padul_climate")
head(Padul_Miami_climate)

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

Padul_Miami_climate <- read.xlsx("Dataset_1.xlsx", sheet = "Padul_climate")
head(Padul_Miami_climate)
colnames(Padul_Miami_climate)

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




# ---Propagated uncertainty for DeltaNPP / NPP Deficit ---
# ---Load climate-derived potential NPP data ---

Padul_Miami_climate <- read.xlsx(
  "Dataset_1.xlsx",
  sheet = "Padul_climate"
)

head(Padul_Miami_climate)
colnames(Padul_Miami_climate)


# --- Check variables ---

summary(Padul_Miami_climate$Predicted_MAT)
summary(Padul_Miami_climate$MAT_RMSE)

summary(Padul_Miami_climate$Predicted_MAP)
summary(Padul_Miami_climate$RMSE_MAP)

anyNA(Padul_Miami_climate$Predicted_MAT)
anyNA(Padul_Miami_climate$MAT_RMSE)

anyNA(Padul_Miami_climate$Predicted_MAP)
anyNA(Padul_Miami_climate$RMSE_MAP)


# --- mean brGDGT-derived NPP values---

rf_mean <- tree_long %>%
  dplyr::group_by(Age_calBP) %>%
  dplyr::summarise(
    NPP_rf = mean(NPP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Age_calBP)


# --- Mean Miami NPP interpolated to brGDGT ages ---

miami_interp <- approx(
  x = Padul_Miami_climate$Age,
  y = Padul_Miami_climate$NPP_Miami,
  xout = rf_mean$Age_calBP,
  rule = 2
)$y


# --- DeltaNPP---

DeltaNPP_central <- miami_interp - rf_mean$NPP_rf

delta_central_df <- data.frame(
  Age_calBP = rf_mean$Age_calBP,
  DeltaNPP = DeltaNPP_central
)


# --- Monte Carlo for potential NPP: sample-specific +/- RMSE for MAT and MAT ---

set.seed(123)

n_iter <- 5000
n_clim <- nrow(Padul_Miami_climate)

NPPpot_mc <- matrix(
  NA_real_,
  nrow = n_clim,
  ncol = n_iter
)


for(i in seq_len(n_clim)) {
  
  # ---MAT uncertainty: +/- specific RMSE---
  
  mat_low <-
    Padul_Miami_climate$Predicted_MAT[i] -
    Padul_Miami_climate$MAT_RMSE[i]
  
  mat_high <-
    Padul_Miami_climate$Predicted_MAT[i] +
    Padul_Miami_climate$MAT_RMSE[i]
  
  mat_sim <- rtruncnorm(
    n = n_iter,
    a = mat_low,
    b = mat_high,
    mean = Padul_Miami_climate$Predicted_MAT[i],
    sd = Padul_Miami_climate$MAT_RMSE[i]
  )
  
  
  # ---MAP uncertainty: +/- specific RMSE---
  
  map_low <- max(
    0,
    Padul_Miami_climate$Predicted_MAP[i] -
      Padul_Miami_climate$RMSE_MAP[i]
  )
  
  map_high <-
    Padul_Miami_climate$Predicted_MAP[i] +
    Padul_Miami_climate$RMSE_MAP[i]
  
  map_sim <- rtruncnorm(
    n = n_iter,
    a = map_low,
    b = map_high,
    mean = Padul_Miami_climate$Predicted_MAP[i],
    sd = Padul_Miami_climate$RMSE_MAP[i]
  )
  
  
  # --- Miami model ---
  
  npp_miami_mat <- 3000 /
    (1 + exp(1.315 - 0.119 * mat_sim))
  
  npp_miami_map <- 3000 *
    (1 - exp(-0.000664 * map_sim))
  
  NPPpot_mc[i, ] <- pmin(
    npp_miami_mat,
    npp_miami_map
  )
}


# --- Uncertainty in potential-NPP ---

NPPpot_low <- apply(
  NPPpot_mc,
  1,
  quantile,
  probs = 0.025,
  na.rm = TRUE
)

NPPpot_high <- apply(
  NPPpot_mc,
  1,
  quantile,
  probs = 0.975,
  na.rm = TRUE
)

NPPpot_median <- apply(
  NPPpot_mc,
  1,
  median,
  na.rm = TRUE
)

summary(
  NPPpot_high - NPPpot_low
)


# --- Dataframe for potential NPP ---

df_NPPpot <- data.frame(
  Age = Padul_Miami_climate$Age,
  
  # Keep original Miami estimate as central line
  NPP_mean = Padul_Miami_climate$NPP_Miami,
  
  NPP_low = NPPpot_low,
  NPP_high = NPPpot_high,
  
  # diagnostic only
  NPP_MC_median = NPPpot_median
)


# --- Prepare for individual tree predictions---

rf_wide <- tree_long %>%
  dplyr::select(
    Age_calBP,
    Tree,
    NPP
  ) %>%
  tidyr::pivot_wider(
    names_from = Tree,
    values_from = NPP
  ) %>%
  dplyr::arrange(Age_calBP)


rf_age <- rf_wide$Age_calBP

RF_matrix <- as.matrix(
  rf_wide[, -1]
)

NPPact_central <- rowMeans(
  RF_matrix,
  na.rm = TRUE
)


# --- Mean potential NPP

NPPpot_central_rf <- approx(
  x = Padul_Miami_climate$Age,
  y = Padul_Miami_climate$NPP_Miami,
  xout = rf_age,
  rule = 2
)$y


# --- Mean DeltaNPP (difference between potential/climate-derived and actual/brGDGT-derived NPP)---

DeltaNPP_central <-
  NPPpot_central_rf -
  NPPact_central


# --- Propagate potential + actual NPP uncertainty ---

set.seed(123)

n_delta <- 5000

DeltaNPP_mc <- matrix(
  NA_real_,
  nrow = length(rf_age),
  ncol = n_delta
)


for(m in seq_len(n_delta)) {
  
  # ---Select one complete potential-NPP Monte Carlo ---
  
  pot_id <- sample(
    seq_len(ncol(NPPpot_mc)),
    size = 1
  )
  
  NPPpot_draw <- approx(
    x = Padul_Miami_climate$Age,
    y = NPPpot_mc[, pot_id],
    xout = rf_age,
    rule = 2
  )$y
  
  
  # ----Select one RF tree ---
  
  rf_id <- sample(
    seq_len(ncol(RF_matrix)),
    size = 1
  )
  
  NPPact_draw <- RF_matrix[, rf_id]
  
  
  # ---Delta NPP ---
  
  DeltaNPP_mc[, m] <-
    NPPpot_draw -
    NPPact_draw
}


# ---Propagated DeltaNPP uncertainty---

DeltaNPP_low <- apply(
  DeltaNPP_mc,
  1,
  quantile,
  probs = 0.025,
  na.rm = TRUE
)

DeltaNPP_high <- apply(
  DeltaNPP_mc,
  1,
  quantile,
  probs = 0.975,
  na.rm = TRUE
)

DeltaNPP_median <- apply(
  DeltaNPP_mc,
  1,
  median,
  na.rm = TRUE
)


# ---Final DeltaNPP dataframe---

hanpp_mean <- data.frame(
  Age_calBP = rf_age,
  
  # Central estimate
  HANPP_mean = DeltaNPP_central,
  
  # Propagated uncertainty envelope
  HANPP_low = DeltaNPP_low,
  HANPP_high = DeltaNPP_high,
  
  # Diagnostic only
  HANPP_MC_median = DeltaNPP_median
)

head(hanpp_mean)


#### Main Figure 6B and 6C ####
#### 
# --- NGRIP ---

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


# --- Pollen plot ---
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


#### Main Figure 6D and 6E ####
#### 
# --- Plot NPP direved from climate and derived from brGDGTs 

NPP_plot_uncertainty <- ggplot() +
  
  # Potential NPP uncertainty
  geom_ribbon(
    data = df_NPPpot,
    aes(
      x = Age,
      ymin = NPP_low,
      ymax = NPP_high
    ),
    fill = "black",
    alpha = 0.10
  ) +
  
  # Potential NPP central estimate
  geom_line(
    data = df_NPPpot,
    aes(
      x = Age,
      y = NPP_mean
    ),
    color = "black",
    linewidth = 1
  ) +
  
  # Actual NPP uncertainty
  geom_ribbon(
    data = results_complete,
    aes(
      x = Age_calBP,
      ymin = NPP_CI95_lower,
      ymax = NPP_CI95_upper
    ),
    fill = "darkgreen",
    alpha = 0.10
  ) +
  
  # Actual NPP central estimate
  geom_line(
    data = results_complete,
    aes(
      x = Age_calBP,
      y = NPP_pred_mean
    ),
    color = "darkgreen",
    linewidth = 1
  ) +
  
  scale_x_reverse(
    limits = c(27000, 5000),
    breaks = seq(
      25000,
      5000,
      by = -5000
    )
  ) +
  
  labs(
    x = "Age (yr cal BP)",
    y = expression(
      NPP~(g~C~m^{-2}~yr^{-1})
    )
  ) +
  
  theme_minimal()


NPP_plot_uncertainty



DNPP_plot2 <- ggplot(
  hanpp_mean,
  aes(x = Age_calBP)
) +
  
  # Propagated uncertainty
  geom_ribbon(
    aes(
      ymin = HANPP_low,
      ymax = HANPP_high
    ),
    fill = "blue",
    alpha = 0.15
  ) +
  
  # Central DeltaNPP
  geom_line(
    aes(y = HANPP_mean),
    color = "blue",
    linewidth = 0.7
  ) +
  
  # Zero reference
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  
  scale_x_reverse(
    limits = c(27000, 5000),
    breaks = seq(
      25000,
      5000,
      by = -5000
    )
  ) +
  ylim(c(-400, 900))+
  
  labs(
    x = "Age (cal yr BP)",
    y = expression(
      Delta*NPP~
        (g~C~m^{-2}~yr^{-1})
    )
  ) +
  
  theme_minimal()


DNPP_plot2



#### Supplementary Figure 13 ####

# ---- Sensitivity of DNPP to temporal lags ---
# 
# Sensitivity of reconstructed NPP deficit (ΔNPP) to temporal offsets between pollen/climate-derived potential NPP 
# and brGDGT-derived actual NPP. 

lags <- c(-1000, -500, -250, -100, 0, 100, 250, 500, 1000)

lag_test <- bind_rows(lapply(lags, function(L){
  
  tmp <- rf_mean %>%
    mutate(
      Lag = L,
      NPP_Miami = approx(
        Padul_Miami_climate$Age + L,
        Padul_Miami_climate$NPP_Miami,
        xout = Age_calBP,
        rule = 1
      )$y,
      HANPP = NPP_Miami - NPP_rf
    )
  
  tmp
}))

# Compare each lag scenario with the original (0-year lag)
base <- lag_test %>%
  filter(Lag == 0) %>%
  dplyr::select(Age_calBP, HANPP_0 = HANPP)

lag_summary <- lag_test %>%
  left_join(base, by = "Age_calBP") %>%
  group_by(Lag) %>%
  summarise(
    rho = cor(HANPP, HANPP_0,
              method = "spearman",
              use = "complete.obs"),
    sign_agreement = mean(
      sign(HANPP) == sign(HANPP_0),
      na.rm = TRUE
    ) * 100
  )

lag_summary

ggplot(
  lag_test,
  aes(
    x = Age_calBP,
    y = HANPP,
    color = Lag,
    group = Lag
  )
) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "grey30",
    high = "red",
    midpoint = 0,
    breaks = c(-1000, -500, 0, 500, 1000),
    name = "Lag (yr)"
  ) +
  scale_x_reverse(
    limits = c(27000, 5000),
    breaks = seq(25000, 5000, by = -5000)
  ) +
  labs(
    x = "Age (cal yr BP)",
    y = expression("NPP Deficit (g C m"^{-2}*" yr"^{-1}*")")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right"
  )


ggplot(
  lag_test,
  aes(
    x = Age_calBP,
    y = HANPP,
    color = Lag,
    group = Lag
  )
) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  
  geom_line(
    data = subset(lag_test, Lag == 0),
    aes(x = Age_calBP, y = HANPP),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.3
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  
  scale_color_gradient2(
    low = "blue",
    mid = "grey70",
    high = "red",
    midpoint = 0,
    breaks = c(-1000, -500, 0, 500, 1000),
    name = "Lag (yr)"
  ) +
  
  scale_x_reverse(
    limits = c(27000, 5000),
    breaks = seq(25000, 5000, by = -5000)
  ) +
  
  labs(
    x = "Age (cal yr BP)",
    y = expression("NPP Deficit (g C m"^{-2}*" yr"^{-1}*")")
  ) +
  
  theme_minimal(base_size = 12)




#### Supplementary Table 5 ####

# --- Extract the original (zero-lag) HANPP reconstruction ---
base <- lag_test %>%
  dplyr::filter(Lag == 0) %>%
  dplyr::select(
    Age_calBP,
    HANPP_0 = HANPP
  )

# --- Compare every lag scenario against the zero-lag reconstruction ---
lag_sensitivity <- lag_test %>%
  dplyr::left_join(base, by = "Age_calBP") %>%
  dplyr::group_by(Lag) %>%
  dplyr::summarise(
    
    # Number of paired observations
    n = sum(complete.cases(HANPP, HANPP_0)),
    
    # Spearman correlation
    Spearman_rho = cor(
      HANPP,
      HANPP_0,
      method = "spearman",
      use = "complete.obs"
    ),
    
    # Percentage retaining the same sign
    Sign_agreement = mean(
      sign(HANPP) == sign(HANPP_0),
      na.rm = TRUE
    ) * 100,
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(Lag)

lag_sensitivity



### Supplementary Figure 12 ####

# --- Additional sensitivity test: effect of using different paleoclimate datasets for potential NPP reconstruction ---

Padul_Miami_climate <- read.xlsx("Dataset_1.xlsx", sheet = "Regional_Paleoclimate")
head(Padul_Miami_climate)


# --- A. Pollen-based paleoclimate ---
# 1. Calculate mean brGDGT-NPP for each age (from all 1000 trees)
rf_mean <- tree_long %>%
  group_by(Age_calBP) %>%
  summarise(NPP_rf = mean(NPP, na.rm = TRUE)) %>%
  arrange(desc(Age_calBP))


# 2. Interpolate Miami NPP to match brGDGT ages
miami_interp <- approx(x = Padul_Miami_climate$Age_P1,
                       y = Padul_Miami_climate$NPP_Miami,
                       xout = rf_mean$Age_calBP)$y

# 3. Calculate DNPP (HANPP) (Miami - brGDGT) for each tree
hanpp_data <- tree_long %>%
  mutate(NPP_Miami = approx(Padul_Miami_climate$Age_P1,
                            Padul_Miami_climate$NPP_Miami,
                            xout = Age_calBP)$y) %>%
  mutate(HANPP =   NPP_Miami -NPP  ) %>%  
  arrange(desc(Age_calBP))

# 4. Calculate mean HANPP across trees
hanpp_mean <- hanpp_data %>%
  group_by(Age_calBP) %>%
  summarise(HANPP_mean = mean(HANPP, na.rm = TRUE),
            HANPP_sd = sd(HANPP, na.rm = TRUE))

# 5. Plot DNPP/HANPP with individual trees and mean
DNPP_plot1 <- ggplot() +
  geom_line(data = hanpp_data,
            aes(x = Age_calBP, y = HANPP, group = Tree),
            alpha = 0.1, color = "darkgrey", linewidth = 0.3) +
  geom_line(data = hanpp_mean,
            aes(x = Age_calBP, y = HANPP_mean),
            color = "darkgreen", linewidth = 0.5) +
  scale_x_reverse(limits = c(27000, 5000),
                  breaks = seq(35000, 5000, by = -5000)) +
  labs(x = "Age (cal yr BP)",
       y = expression(paste("NPP Deficit (g C m"^{-2}, " yr"^{-1}, ")")),
       title = " ") +
  theme_minimal()
#
DNPP_plot1


# --- B. HadCM3B-M2.1 paleoclimate ---
# 
# 1. Calculate mean brGDGT-NPP for each age (from all 1000 trees)
rf_mean <- tree_long %>%
  group_by(Age_calBP) %>%
  summarise(NPP_rf = mean(NPP, na.rm = TRUE)) %>%
  arrange(desc(Age_calBP))


# 2. Interpolate Miami NPP to match brGDGT ages
miami_interp <- approx(x = Padul_Miami_climate$Age_GCM,
                       y = Padul_Miami_climate$NPP_Miami_GCM ,
                       xout = rf_mean$Age_calBP)$y

# 3. Calculate DNPP (HANPP) (Miami - brGDGT) for each tree
hanpp_data <- tree_long %>%
  mutate(NPP_Miami = approx(Padul_Miami_climate$Age_GCM,
                            Padul_Miami_climate$NPP_Miami_GCM ,
                            xout = Age_calBP)$y) %>%
  mutate(HANPP =   NPP_Miami -NPP  ) %>%  
  arrange(desc(Age_calBP))

# 4. Calculate mean HANPP across trees
hanpp_mean <- hanpp_data %>%
  group_by(Age_calBP) %>%
  summarise(HANPP_mean = mean(HANPP, na.rm = TRUE),
            HANPP_sd = sd(HANPP, na.rm = TRUE))

# 5. Plot DNPP/HANPP with individual trees and mean
DNPP_plot2 <- ggplot() +
  geom_line(data = hanpp_data,
            aes(x = Age_calBP, y = HANPP, group = Tree),
            alpha = 0.1, color = "darkgrey", linewidth = 0.3) +
  geom_line(data = hanpp_mean,
            aes(x = Age_calBP, y = HANPP_mean),
            color = "darkgreen", linewidth = 0.5) +
  scale_x_reverse(limits = c(27000, 5000),
                  breaks = seq(35000, 5000, by = -5000)) +
  labs(x = "Age (cal yr BP)",
       y = expression(paste("NPP Deficit (g C m"^{-2}, " yr"^{-1}, ")")),
       title = " ") +
  theme_minimal()
#
DNPP_plot2

# --- C. CHELSA-TraCE21k paleoclimate ---
miami_interp <- approx(x = Padul_Miami_climate$Age_CHELSA ,
                       y = Padul_Miami_climate$NPP_Miami_CHELSA ,
                       xout = rf_mean$Age_calBP)$y

# 3. Calculate DNPP (HANPP) (Miami - brGDGT) for each tree
hanpp_data <- tree_long %>%
  mutate(NPP_Miami = approx(Padul_Miami_climate$Age_CHELSA ,
                            Padul_Miami_climate$NPP_Miami_CHELSA ,
                            xout = Age_calBP)$y) %>%
  mutate(HANPP =   NPP_Miami -NPP  ) %>%  
  arrange(desc(Age_calBP))

# 4. Calculate mean HANPP across trees
hanpp_mean <- hanpp_data %>%
  group_by(Age_calBP) %>%
  summarise(HANPP_mean = mean(HANPP, na.rm = TRUE),
            HANPP_sd = sd(HANPP, na.rm = TRUE))

# 5. Plot DNPP/HANPP with individual trees and mean
DNPP_plot3 <- ggplot() +
  geom_line(data = hanpp_data,
            aes(x = Age_calBP, y = HANPP, group = Tree),
            alpha = 0.1, color = "darkgrey", linewidth = 0.3) +
  geom_line(data = hanpp_mean,
            aes(x = Age_calBP, y = HANPP_mean),
            color = "darkgreen", linewidth = 0.5) +
  scale_x_reverse(limits = c(27000, 5000),
                  breaks = seq(35000, 5000, by = -5000)) +
  labs(x = "Age (cal yr BP)",
       y = expression(paste("NPP Deficit (g C m"^{-2}, " yr"^{-1}, ")")),
       title = " ") +
  theme_minimal()
#
DNPP_plot3



grid.arrange(DNPP_plot1,
             DNPP_plot2,
             DNPP_plot3)

# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 7)  ARCHAEOLOGICAL ANALYSES ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# --- SPD by culture ---

# Delta values for marine samples
deltaR_val  <- 94
deltaR_sd   <- 61

# List of cultures
cultures <- c("Gravettian", "Solutrean", "Magdalenian", "Epipalaeolithic", "Neolithic")


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


#### Main Figure 6F ####
#### 
# --- SPD ---
# Charcoal concentration (n charcoals/cm3) throughout the sequence:
# Scaling factor to put charcoal on the same plotting range as PrDens
scale_factor <- max(Regions_SPD$PrDens, na.rm = TRUE) /
  max(pollen_df$Charcoals_cm3, na.rm = TRUE)

SPD_l3_charcoal <- ggplot() +
  
  # SPD by culture
  geom_area(
    data = Regions_SPD,
    aes(x = Age, y = PrDens, fill = Culture),
    alpha = 0.5,
    position = "identity",
    color = "black"
  ) +
  
  # Charcoal record
  geom_point(
    data = pollen_df,
    aes(x = Age, y = Charcoals_cm3 * scale_factor),
    color = "black",
    size = 1.5
  ) +
  
  scale_x_reverse(
    limits = c(27000, 5000),
    breaks = seq(25000, 5000, by = -5000)
  ) +
  
  scale_y_continuous(
    name = "Probability Density",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = expression("Charcoal (particles cm"^{-3}*")")
    )
  ) +
  
  scale_fill_brewer(palette = "Set2") +
  
  labs(
    x = "cal BP",
    title = "Summed Probability Distribution and Charcoal"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    axis.title.y.right = element_text(),
    axis.text.y.right = element_text()
  )

SPD_l3_charcoal



# --- Create a dataframe that combines PrDens, NPP_brGDGT and DNPP ---

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



# --- Correlations with mean values of each culture ---

spd_ages <- Regions_SPD %>%
  filter(PrDens >= 0.001) %>%
  dplyr::select(Age, Culture, PrDens) %>%
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

#### Main Figure 7A-C ####
#--- Summary statistics and correlations ---

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




# --- Boxplot NPP_brGDGT / Culture ---
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


p3 <- ggplot(all_data_complete, aes(x = reorder(Culture, -Age), y = HANPP, fill = Culture)) +
  geom_boxplot(alpha = 0.7,   outlier.color = "gray", 
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

#### Supplementary Table 4 ####
print(stats_by_culture)



# --- Filter PrDens >= 0.001 ---
all_data_fixed_thresh <- all_data %>%
  filter(PrDens >= 0.001)

#--- Correlation between PrDens and Defcit of NPP (DNPP/HANPP) (Figure 8D) ---
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

# --- Same ages for DNPP and NPP_rf ---
hanpp_interp <- approx(
  x = hanpp_data$Age_calBP,
  y = hanpp_data$HANPP,
  xout = rf_mean$Age_calBP,
  rule = 2 
)$y

# --- The same for Pollen_df (Pollen_cm3, Carbones_cm3, PAR) ---
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

#### Supplementary Figure 10A-C ####
# --- Correlation potential vs actual NPP ---
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






#### Supplementary Figure 11A-D ####
# Now check association between charcoal concentration (charcoals/cm3) and
# deficit of NPP for the whole sequence, and since the Magdalenian:

combined_df_P<- subset(combined_df, combined_df$Age_calBP<16000)
colnames(combined_df)

Pot_Carchoal_All<-ggplot(combined_df, aes(x = combined_df$Carbones_cm3, y = NPP_climate)) +
  # Points colored by Age
  geom_point(aes(color = Age_calBP ), size = 3, alpha = 0.8) +
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

Act_Carchoal_All<-ggplot(combined_df, aes(x = combined_df$Carbones_cm3, y = NPP_rf)) +
  # Points colored by Age
  geom_point(aes(color = Age_calBP ), size = 3, alpha = 0.8) +
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

Act_Carchoal_All




Act_Carchoal_Magd<-ggplot(combined_df_P, aes(x = combined_df_P$Carbones_cm3, y = NPP_rf)) +
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

Act_Carchoal_Magd

Pot_Carchoal_Magd<-ggplot(combined_df_P, aes(x = combined_df_P$Carbones_cm3, y = NPP_climate)) +
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



grid.arrange(Pot_Carchoal_All, Pot_Carchoal_Magd, Act_Carchoal_All, Act_Carchoal_Magd, ncol=2)


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




#### Supplementary Table 3 ####

# --- ANOVA ---
anova_brgdgt <- aov(
  NPP_brGDGT ~ Culture,
  data = all_data_complete
)

anova_miami <- aov(
  NPP_Miami ~ Culture,
  data = all_data_complete
)

anova_hanpp <- aov(
  HANPP ~ Culture,
  data = all_data_complete
)


# --- Tukey hsd ---

tukey_brgdgt <- TukeyHSD(anova_brgdgt)
tukey_miami  <- TukeyHSD(anova_miami)
tukey_hanpp  <- TukeyHSD(anova_hanpp)


# Function to convert tukey output to dataframe ---

tidy_tukey <- function(tukey_object, variable_name) {
  
  out <- as.data.frame(tukey_object$Culture)
  
  out$Comparison <- rownames(out)
  rownames(out) <- NULL
  
  out %>%
    transmute(
      Variable = variable_name,
      Comparison = Comparison,
      Mean_difference = diff,
      CI_95_lower = lwr,
      CI_95_upper = upr,
      p_adjusted = `p adj`
    )
}


# --- Create tables ---

table_brgdgt <- tidy_tukey(
  tukey_brgdgt,
  "brGDGT-derived NPP"
)

table_miami <- tidy_tukey(
  tukey_miami,
  "Climate-derived NPP"
)

table_hanpp <- tidy_tukey(
  tukey_hanpp,
  "NPP deficit"
)

# --- Combine ---
tukey_table <- bind_rows(
  table_brgdgt,
  table_miami,
  table_hanpp
)


tukey_table


# --- Correlations between prDens and NPP ---

# --- SPD ---
spd_ages <- Regions_SPD %>%
  filter(PrDens >= 0.001) %>%
  dplyr::select(Age, Culture, PrDens) %>%
  distinct()

# --- Same ages ---
hanpp_by_age <- hanpp_data %>%
  group_by(Age_calBP) %>%
  summarise(HANPP_mean = mean(HANPP, na.rm = TRUE))

hanpp_interp <- approx(
  x = hanpp_by_age$Age_calBP,
  y = hanpp_by_age$HANPP_mean,
  xout = spd_ages$Age,
  rule = 2
)

# --- Create DF ---
all_data_complete <- spd_ages %>%
  mutate(HANPP = hanpp_interp$y) %>%
  filter(!is.na(HANPP) & !is.na(PrDens))

# --- Summary ---

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


# --- Correlation tests and plots ---

cultures_corr <- c(
  "Aurignacian",
  "Gravettian",
  "Solutrean",
  "Magdalenian",
  "Epipalaeolithic",
  "Neolithic"
)

Regions_SPD_corr <- data.frame()

for(cult in cultures_corr){
  
  # Filter by culture and non-missing data
  cult_df <- SPD_df %>%
    dplyr::filter(
      Included == "y",
      Culture == cult
    ) %>%
    dplyr::filter(
      !is.na(Age),
      !is.na(s.dev.),
      !is.na(Curve),
      !is.na(Labref),
      !is.na(layer)
    )
  
  if(nrow(cult_df) == 0){
    message("No dates for culture: ", cult)
    next
  }
  
  # Calibration curves and deltaR
  curves <- ifelse(
    cult_df$Curve == "Terrestrial",
    "intcal20",
    "marine20"
  )
  
  deltaR_vals <- ifelse(
    cult_df$Curve == "Marine",
    deltaR_val,
    0
  )
  
  deltaR_sds <- ifelse(
    cult_df$Curve == "Marine",
    deltaR_sd,
    0
  )
  
  # Calibrate individual dates
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
    dplyr::group_by(layer) %>%
    dplyr::summarise(
      combined_rads = list(
        combine(
          caldates_ind[
            which(cult_df$layer == unique(layer))
          ],
          fixIDs = TRUE
        )$date
      ),
      .groups = "drop"
    )
  
  # Prepare vectors for SPD
  ages <- unlist(
    lapply(combined$combined_rads, function(x) x$age)
  )
  
  sdev <- unlist(
    lapply(combined$combined_rads, function(x) x$error)
  )
  
  # 300-year bins
  bins <- binPrep(
    sites = rep(
      combined$layer,
      times = sapply(combined$combined_rads, length)
    ),
    ages = ages,
    h = 300
  )
  
  # Compute SPD
  spd_cult <- spd(
    caldates_ind,
    timeRange = c(35000, 2000),
    bins = bins
  )
  
  # Add to correlation-specific SPD dataframe
  Regions_SPD_corr <- rbind(
    Regions_SPD_corr,
    data.frame(
      PrDens = spd_cult$grid$PrDens,
      Age = spd_cult$grid$calBP,
      Culture = cult
    )
  )
}


spd_ages <- Regions_SPD_corr %>%
  filter(PrDens >= 0.001) %>%
  dplyr::select(Age, Culture, PrDens) %>%
  distinct()

hanpp_by_age <- hanpp_data %>%
  group_by(Age_calBP) %>%
  summarise(
    HANPP_mean = mean(HANPP, na.rm = TRUE),
    .groups = "drop"
  )

hanpp_interp <- approx(
  x = hanpp_by_age$Age_calBP,
  y = hanpp_by_age$HANPP_mean,
  xout = spd_ages$Age,
  rule = 2
)

all_data_complete <- spd_ages %>%
  mutate(
    HANPP = hanpp_interp$y
  ) %>%
  filter(
    !is.na(HANPP),
    !is.na(PrDens)
  )

culture_stats <- all_data_complete %>%
  group_by(Culture) %>%
  summarise(
    mean_PrDens = mean(PrDens, na.rm = TRUE),
    sd_PrDens   = sd(PrDens, na.rm = TRUE),
    mean_HANPP  = mean(HANPP, na.rm = TRUE),
    sd_HANPP    = sd(HANPP, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )


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

#### Main Figure 8A-E ####
#### 
grid.arrange(p, p_l, p_d,
             ncol = 3, nrow = 1)
grid.arrange(DNPP_Dens_plot,DNPP_Charcoal_cm3,
             ncol = 2, nrow = 1)


############# ADDITIONAL INFORMATION ------------------------------

# The rest of the code is available in the Script_2.R file. It includes all the modelling of
# herbivore biomass and carrying capapcity of secondary consumers.
#


