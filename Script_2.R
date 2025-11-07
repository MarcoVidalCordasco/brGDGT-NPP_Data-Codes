
rm(list = ls()) # Clear all
# SETUP ####

# the following command sets the working directory to the folder where this
# script is located (similar to RStudio menu "Session - Set Working Directory -
# To Source File Location"):
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####
library(openxlsx) 
library(gridExtra)
library(beanplot)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(magick)
library(robustbase)
library(RColorBrewer)
library(tidyr)
library(igraph)

df<-read.xlsx("Dataset_CC.xlsx", rowNames=FALSE, 
              colNames=TRUE, sheet="Fauna_Padul_Pleistocene")
head(df)

### a) HERBIVORE DENSITIES ####

# ---- 1. Load data from Excel file ----

# Read data from the "Fauna_Padul_Pleistocene" sheet
df <- read.xlsx("Dataset_CC.xlsx",
                rowNames = FALSE,
                colNames = TRUE,
                sheet = "Fauna_Padul_Pleistocene")

# Check 
head(df)


# ---- 2. Compute THB values (actual and potential) ----

# THB_act = Total Herbivore biomass based on actual (brGDGT-derived) NPP
df$THB_act <- 10^(1.401 * log10(df$NPP_act) - 0.642)

# THB_pot = Total Herbivore biomass based on potential (brGDGT-derived) NPP
df$THB_pot <- 10^(1.401 * log10(df$NPP_pot) - 0.642)

# Verify that new columns are created
head(df)


# ---- 3. Create an empty data frame to store results ----

# Start with the species and their body mass
outputs <- data.frame(
  Primary_Consumers = df$Primary_Consumers,
  BM_pc = df$BM_pc
)



# ---- 4. Get list of unique cultures (excluding NAs) ----

cultures <- unique(df$Culture)
cultures <- cultures[!is.na(cultures)]
cultures


# ---- 5. Loop through each culture ----

# For every culture, calculate densities (actual and potential)
# and add them as new columns in 'outputs'

for (i in seq_along(cultures)) {
  
  # Get the name of the culture
  culture_name <- cultures[i]
  cat("Processing culture:", culture_name, "\n")
  
  # Select the rows corresponding to that culture
  df_cult <- df[df$Culture == culture_name, ]
  

  # First, densities based on actual THB
  # Step 1: get THB_act for that culture
  THB_act_value <- df_cult$THB_act[1]
  
  # Step 2: compute the sum of BM_pc^0.25
  sum_BM25 <- sum(df$BM_pc^0.25, na.rm = TRUE)
  
  # Step 3: compute c
  c_act <- THB_act_value / sum_BM25
  
  # Step 4: compute density for each species (actual)
  Density_act <- c_act * (df$BM_pc^-0.75)
  
  # Add this new column to 'outputs'
  colname_act <- paste0(culture_name, "_act")
  outputs[[colname_act]] <- Density_act
  
  
  # Now, densities based on potential THB
  # Step 1: get THB_pot for that culture
  THB_pot_value <- df_cult$THB_pot[1]
  
  # Step 2: compute c
  c_pot <- THB_pot_value / sum_BM25
  
  # Step 3: compute density for each species (potential)
  Density_pot <- c_pot * (df$BM_pc^-0.75)
  
  # Add this new column to 'outputs'
  colname_pot <- paste0(culture_name, "_pot")
  outputs[[colname_pot]] <- Density_pot
}


# ---- 6. Check the final output table ----

head(outputs)

# ---- 7. Save results to a CSV file ----

write.csv(outputs, "Outputs_Density_Pleistocene.csv", row.names = FALSE)


## Repeat for the Holocene mammal composition following these steps:

# 1) In line 38, replace "Fauna_Padul_Pleistocene" by "Padul_Holocene"
# 2) Run again from line 35 to line 120
# 3) Replace in line 131 "Outputs_Density_Pleistocene.csv" by "Outputs_Density_Holocene.csv"




# b) HERBIVORE POPULATION STRUCTURE & BIOMASS BY WEIGHT CLASS ####
### GET POPULATION STRUCTURE & BIOMASS BY WEIGHT CLASSES ####

# ---- 1. Function to reconstruct stable and stationary population ----
stability_function <- function(params, a, beta) {
  alpha <- params[1]
  sum(a * exp(-((1:length(a) - 1) / alpha)^beta)) - 1
}

# ---- 2. List of species ----
 # Modify the list for Holocene scenario 

species_list<- c("Mammuthus.primigenius",
                 "Equus.hydruntinus",
                 "Equus.ferus",
                 "Bison.priscus",
                 "Bos.primigenius",
                 "Cervus.elaphus",
                 "Capreolus.capreolus",
                 "Capra.pyrenaica",
                 "Rupicapra.rupicapra",
                 "Sus.scrofa", 
                 "Oryctolagus")

# ---- 3. Density data ----
# These are the outputs obtained in a) HERBIVORE DENSITIES and saved as 'Outputs_Density_Pleistocene.csv' or 'Outputs_Density_Holocene.csv'
df_PC <- read.xlsx("Dataset_CC.xlsx",
                rowNames = FALSE,
                colNames = TRUE,
                sheet = "Actual_Pleistocene")
head(df_PC)

chrons <- c("Aurignacian_act",    "Gravettian_act" ,   
             "Solutrean_act",     "Magdalenian_act",  "Epipaleolithic_act")

# ---- 4. Dataframe where outputs will be stored ----
outputs_bm <- data.frame(
  Species = character(),
  Chron = character(),
  Small_Biomass = numeric(),
  Medium_Biomass = numeric(),
  MediumLarge_Biomass = numeric(),
  Large_Biomass = numeric(),
  stringsAsFactors = FALSE
)


# Example with one species:
df<-read.xlsx("Dataset_CC.xlsx", rowNames=FALSE, 
              colNames=TRUE, sheet="Mammuthus.primigenius") # Replace with any species

Species <- na.omit(df$Species[1])
a<- na.omit(df$Fecundity)
# Generate 100 radom values with beta between 0 and 1
beta_values <- runif(100, min = 0, max = 1)
survival_curves <- list()

# Loop 100 survival curves
for (beta in beta_values) {
  # Find valid values of alpha
  test_interval <- function(lower, upper, stability_function, a, beta) {
    tryCatch({
      alpha_estimate <- uniroot(stability_function, interval = c(lower, upper), a = a, beta = beta)$root
      return(alpha_estimate)
    }, error = function(e) {
      return(NA)
    })
  }
  
  search_range <- seq(0, 10, by = 0.1)
  valid_intervals <- lapply(1:(length(search_range) - 1), function(i) {
    lower <- search_range[i]
    upper <- search_range[i + 1]
    alpha_estimate <- test_interval(lower, upper, stability_function, a, beta)
    if (!is.na(alpha_estimate)) {
      return(list(lower = lower, upper = upper, alpha_estimate = alpha_estimate))
    } else {
      return(NULL)
    }
  })
  
  valid_intervals <- Filter(Negate(is.null), valid_intervals)
  valid_intervals
  
  # If there are not valid values of alpha, try with another beta value
  if (length(valid_intervals) == 0) next
  
  # Get the first alpha value
  alpha_estimate <- valid_intervals[[1]]$alpha_estimate
  
  # Compute the survival curve for that beta and alpha
  proportion_alive <- exp(-(((1:length(a) - 1) / alpha_estimate) ^ beta))
  survival_curves[[length(survival_curves) + 1]] <- list(proportion_alive = proportion_alive, beta = beta, alpha = alpha_estimate)
}

# Compute the mean survival curve 
proportion_matrix <- do.call(cbind, lapply(survival_curves, function(curve) curve$proportion_alive))

# Survival curve by age
mean_proportion_alive <- rowMeans(proportion_matrix)

# Plot the 100 population survival curves 
plot(1, type = "n", xlab = "Age", ylab = "Alive (%)", xlim = c(1, length(a)), ylim = c(0, 1), main = paste("", Species))
for (curve in survival_curves) {
  lines(curve$proportion_alive, col = "grey")
  lines(mean_proportion_alive, col = "red", lwd = 3) # Average
}






# ---- 5. Compute survival curves for all herbivore species ----
for (j in species_list){

  df<-read.xlsx("Dataset_CC.xlsx", rowNames=FALSE, 
                colNames=TRUE, sheet=j) # j
  head(df)
  Species <- na.omit(df$Species[1])
  a<- na.omit(df$Fecundity)
  # Generate 100 radom values with beta between 0 and 1
  beta_values <- runif(100, min = 0, max = 1)
  survival_curves <- list()
  
  for (beta in beta_values) {
    # Find valid values of alpha
    test_interval <- function(lower, upper, stability_function, a, beta) {
      tryCatch({
        alpha_estimate <- uniroot(stability_function, interval = c(lower, upper), a = a, beta = beta)$root
        return(alpha_estimate)
      }, error = function(e) {
        return(NA)
      })
    }
    
    search_range <- seq(0, 10, by = 0.1)
    valid_intervals <- lapply(1:(length(search_range) - 1), function(i) {
      lower <- search_range[i]
      upper <- search_range[i + 1]
      alpha_estimate <- test_interval(lower, upper, stability_function, a, beta)
      if (!is.na(alpha_estimate)) {
        return(list(lower = lower, upper = upper, alpha_estimate = alpha_estimate))
      } else {
        return(NULL)
      }
    })
    
    valid_intervals <- Filter(Negate(is.null), valid_intervals)
    valid_intervals
    
    # If there are not valid values of alpha, try with another beta value
    if (length(valid_intervals) == 0) next
    
    # Get the first alpha value
    alpha_estimate <- valid_intervals[[1]]$alpha_estimate
    
    # Compute the survival curve for that beta and alpha
    proportion_alive <- exp(-(((1:length(a) - 1) / alpha_estimate) ^ beta))
    survival_curves[[length(survival_curves) + 1]] <- list(proportion_alive = proportion_alive, beta = beta, alpha = alpha_estimate)
  }
  
  # Compute the mean survival curve 
  proportion_matrix <- do.call(cbind, lapply(survival_curves, function(curve) curve$proportion_alive))
  
  # Survival curve by age
  mean_proportion_alive <- rowMeans(proportion_matrix)
  
  
  # Plot the 100 population survival curves 
  plot(1, type = "n", xlab = "Age", ylab = "Alive (%)", xlim = c(1, length(a)), ylim = c(0, 1), main = paste("", Species))
  for (curve in survival_curves) {
    lines(curve$proportion_alive, col = "grey")
    lines(mean_proportion_alive, col = "red", lwd = 3) # Average
  }
  
  # Use the mean survival curves in the following chunks
  proportion_alive <- mean_proportion_alive
  
  # Get alpha and beta values of each curve
  alpha_values <- sapply(survival_curves, function(curve) curve$alpha)
  beta_values <- sapply(survival_curves, function(curve) curve$beta)
  
  # Compute mean
  mean_alpha <- mean(alpha_values, na.rm = TRUE)
  mean_beta <- mean(beta_values, na.rm = TRUE)
  
  alpha_estimate <- mean_alpha
  beta <- mean_beta
  
  # Verify that population is stable and stationary (sum(e) = 1)
  Fec <- a
  e <- proportion_alive * a
  print(round(sum(e)))
  
  
  # ---- 5. Compute biomass of small, medium, mediumlarge and large herbivores ----
  for (chron_value in chrons) {

    # Compute individuals (ind/km2) that survive each year and biomass
    
    # Find the row corresponding to that species
    row_index <- which(df_PC$Primary_Consumers == j)
    
    # Skip if species not found
    if (length(row_index) == 0) {
      cat("Check species' name: species not found in df_PC:", j, "\n")
      next
    }
    
    # Extract density
    Density <- df_PC[[chron_value]][row_index]
    
    # Skip if density is NA
    if (is.na(Density)) {
      cat("Check species' name: no density for", j, "in", chron_value, "\n")
      next
    }
    
  #Create file where biomass outputs will be stored
  
    df$factor_death <- abs(c(proportion_alive[1], diff(proportion_alive)))
    df$factor_death<- proportion_alive / sum(proportion_alive) 
    sum(df$factor_death)
    plot(df$factor_death)
    
    head(df)
    
    
    # Create weight size categories in df
    df$Category[df$Weight < 10] <- "Small"
    df$Category[df$Weight >= 10 & df$Weight < 100] <- "Medium"
    df$Category[df$Weight >= 100 & df$Weight < 500] <- "MediumLarge"
    df$Category[df$Weight >= 500] <- "Large"
    
head(df)
Small_f <- sum(subset(df, Category=="Small")$factor_death * subset(df, Category=="Small")$Weight)
Medium_f <- sum(subset(df, Category=="Medium")$factor_death * subset(df, Category=="Medium")$Weight)
MediumLarge_f <- sum(subset(df, Category=="MediumLarge")$factor_death * subset(df, Category=="MediumLarge")$Weight)
Large_f <- sum(subset(df, Category=="Large")$factor_death * subset(df, Category=="Large")$Weight)

Small_Biomass <- Small_f * Density
Medium_Biomass <- Medium_f * Density
MediumLarge_Biomass <- MediumLarge_f * Density
Large_Biomass <- Large_f * Density

# Add a row to the outputs table
outputs_bm <- rbind(outputs_bm, data.frame(
  Species = j,
  Chron = chron_value,
  Small_Biomass = Small_Biomass,
  Medium_Biomass = Medium_Biomass,
  MediumLarge_Biomass = MediumLarge_Biomass,
  Large_Biomass = Large_Biomass,
  stringsAsFactors = FALSE
))
    
    
  }
}

head(outputs_bm)
# ---- 6. Save outputs ----
write.csv(outputs_bm, "outputs__Biomass.csv")


# ---- 7. Summary statistics ----
# Convert to data.frame
outputs_bm <- as.data.frame(outputs_bm)
outputs_bm$Chron <- as.character(outputs_bm$Chron)

# Compute mean and SD
summary_biomass <- aggregate(
  cbind(Small_Biomass, Medium_Biomass, MediumLarge_Biomass, Large_Biomass) ~ Chron,
  data = outputs_bm,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
)
summary_biomass <- do.call(data.frame, summary_biomass)

# Rename columns
colnames(summary_biomass) <- c(
  "Chron",
  "Small_mean", "Small_sd",
  "Medium_mean", "Medium_sd",
  "MediumLarge_mean", "MediumLarge_sd",
  "Large_mean", "Large_sd"
)

# Order
summary_biomass <- summary_biomass[order(as.character(summary_biomass$Chron)), , drop = FALSE]

# Results
print(summary_biomass)

#Save
write.csv(summary_biomass, "summary_biomass.csv") ## These outputs are available in sheet Biomass in Dataset_CC.xlsx






# ---- 8. Load data ----
df_biomass <- read.xlsx("Dataset_CC.xlsx", sheet= "Biomass")

head(df_biomass)

# Check the unique combinations of Composition and Chron
unique_combinations <- df_biomass %>%
  distinct(Composition, Chron)

print(unique_combinations)

# Create grouping variable
group_vars <- with(df_biomass, interaction(Composition, Chron, drop = TRUE))
group_levels <- levels(group_vars)

# ---- 9. Calculate statistics by group ----

result <- by(df_biomass, group_vars, function(x) {
  c(
    Small_mean = mean(x$Small_Biomass, na.rm = TRUE),
    Small_sd = sd(x$Small_Biomass, na.rm = TRUE),
    Medium_mean = mean(x$Medium_Biomass, na.rm = TRUE),
    Medium_sd = sd(x$Medium_Biomass, na.rm = TRUE),
    MediumLarge_mean = mean(x$MediumLarge_Biomass, na.rm = TRUE),
    MediumLarge_sd = sd(x$MediumLarge_Biomass, na.rm = TRUE),
    Large_mean = mean(x$Large_Biomass, na.rm = TRUE),
    Large_sd = sd(x$Large_Biomass, na.rm = TRUE)
  )
})

# Extract group names and create proper data frame
group_names <- strsplit(as.character(group_levels), "\\.")
Composition <- sapply(group_names, `[`, 1)
Chron <- sapply(group_names, `[`, 2)

# Create final summary statistics data frame
summary_stats_correct <- data.frame(
  Composition = Composition,
  Chron = Chron,
  do.call(rbind, result),
  row.names = NULL
)

# ---- 10. Print, save and plot results ----

print(summary_stats_correct)

write.csv(summary_stats_correct, "summary_stats_HerbivoreBiomass.csv") # This is available in the sheed Summary_Biomass of Dataset_CC.xlsx


# Define correct chronological order
chronological_order <- c("Aurignacian", "Gravettian", "Solutrean", "Magdalenian", "Epipaleolithic")

transition_data <- summary_stats_correct %>%
  mutate(
    Chron = factor(Chron, levels = chronological_order),
    Period_Type = case_when(
      Chron %in% c("Aurignacian", "Gravettian", "Solutrean") ~ "Definite Pleistocene",
      Chron == "Magdalenian" ~ "Transition", 
      Chron == "Epipaleolithic" ~ "Definite Holocene"
    ),
    Period_Type = factor(Period_Type, 
                        levels = c("Definite Pleistocene", "Transition", "Definite Holocene"))
  ) %>%
  pivot_longer(
    cols = c(Small_mean, Medium_mean, MediumLarge_mean, Large_mean),
    names_to = "Biomass_Type",
    values_to = "Mean"
  ) %>%
  mutate(
    Biomass_Type = gsub("_mean", "", Biomass_Type),
    SD = case_when(
      Biomass_Type == "Small" ~ Small_sd,
      Biomass_Type == "Medium" ~ Medium_sd,
      Biomass_Type == "MediumLarge" ~ MediumLarge_sd,
      Biomass_Type == "Large" ~ Large_sd
    ),
    Biomass_Type = factor(Biomass_Type,
                         levels = c("Small", "Medium", "MediumLarge", "Large"),
                         labels = c("Small", "Medium", "Medium-Large", "Large")),
    lower_error = pmax(Mean - SD, 0),  
    upper_error = Mean + SD
  ) %>%
  select(Composition, Chron, Biomass_Type, Mean, SD, lower_error, upper_error, Period_Type)

# Plot
ggplot(transition_data, aes(x = Chron, y = Mean, fill = Biomass_Type)) +
  geom_col(position = "dodge", alpha = 0.8) +
  facet_wrap(~ Composition, ncol = 1, scales = "free_y", 
             labeller = as_labeller(c("Actual_Pleistocene" = "Pleistocene Composition", 
                                      "Actual_Holocene" = "Holocene Composition"))) +
  scale_fill_manual(values = c("Small" = "#66C2A5", "Medium" = "#FC8D62", 
                               "Medium-Large" = "#8DA0CB", "Large" = "#E78AC3")) +
  labs(title = "Herbivore Biomass by Size Class and Faunal Composition",
       subtitle = "Comparison of Pleistocene (top) vs Holocene (bottom) scenarios across chronological periods",
       x = "Chronological Period", 
       y = "Mean Biomass",
       fill = "Size Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing = unit(1, "lines"))



# Compute % change between Pleistocene and Holocene
change_data <- summary_stats_correct %>%
  select(Composition, Chron, Small_mean, Medium_mean, MediumLarge_mean, Large_mean) %>%
  pivot_wider(names_from = Composition, values_from = c(Small_mean, Medium_mean, MediumLarge_mean, Large_mean)) %>%
  mutate(
    Small_change = (Small_mean_Actual_Holocene - Small_mean_Actual_Pleistocene) / Small_mean_Actual_Pleistocene * 100,
    Medium_change = (Medium_mean_Actual_Holocene - Medium_mean_Actual_Pleistocene) / Medium_mean_Actual_Pleistocene * 100,
    MediumLarge_change = (MediumLarge_mean_Actual_Holocene - MediumLarge_mean_Actual_Pleistocene) / MediumLarge_mean_Actual_Pleistocene * 100,
    Large_change = (Large_mean_Actual_Holocene - Large_mean_Actual_Pleistocene) / Large_mean_Actual_Pleistocene * 100
  ) %>%
  select(Chron, ends_with("change")) %>%
  pivot_longer(cols = -Chron, names_to = "Size_Class", values_to = "Percent_Change") %>%
  mutate(
    Size_Class = gsub("_change", "", Size_Class),
    Size_Class = factor(Size_Class, levels = c("Small", "Medium", "MediumLarge", "Large")),
    Chron = factor(Chron, levels = chronological_order)
  )

ggplot(change_data, aes(x = Chron, y = Percent_Change, fill = Size_Class)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Percentage Change in Biomass: Holocene vs Pleistocene",
       subtitle = "Positive values = Increase in Holocene | Negative values = Decrease in Holocene",
       x = "Chronological Period", 
       y = "Percentage Change (%)",
       fill = "Size Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#### c) CARRYING CAPACITY FOR SECONDARY CONSUMERS ####

# ---- 1. Load data from Excel file ----
df_biomass<- read.xlsx("Dataset_CC.xlsx", sheet = "Summary_Biomass")

head(df_biomass)

df_biomass<- subset(df_biomass, df_biomass$Composition=="Actual_Holocene") # Replace with Pleistocene and the % of human meat intake [30-60, by 10]

Prey_preferences<- read.xlsx("Prey_preferences.xlsx", rowNames=FALSE, 
                             colNames=TRUE, sheet="Diet_30_Holocene")  # Replace with Pleistocene and the % of human meat intake [30-60, by 10]
head(Prey_preferences)


db_values <- data.frame(
  Predator = colnames(Prey_preferences)[-1],  # Exclude Weight_Category column 
  AnnualBiomass = as.numeric(Prey_preferences[2, -1])  # Remove first column
)

# Verify
print(db_values)


# ---- 2. Wastage Factor ----
WastageFactor <- c(Small = 0.8, Medium = 0.75, MedLarge = 0.65, Large = 0.6)

# ---- 3. Apply Wastage Factors to mean values of biomass  ----
df_biomass_corrected <- df_biomass
df_biomass_corrected$Small_mean      <- df_biomass$Small_mean      * WastageFactor["Small"]
df_biomass_corrected$Medium_mean     <- df_biomass$Medium_mean     * WastageFactor["Medium"]
df_biomass_corrected$MediumLarge_mean <- df_biomass$MediumLarge_mean* WastageFactor["MedLarge"]
df_biomass_corrected$Large_mean      <- df_biomass$Large_mean      * WastageFactor["Large"]

# Filter
prey_subset <- Prey_preferences[Prey_preferences$Weight_Category %in% c("Small","Medium","MedLarge","Large"), ]

# Convert into numeric
prey_numeric <- prey_subset
for (col in colnames(prey_numeric)[-1]) {
  prey_numeric[[col]] <- as.numeric(as.character(prey_numeric[[col]]))
}

predators <- colnames(prey_numeric)[-1]
prey_long <- data.frame(
  SizeCategory = rep(prey_numeric$Weight_Category, times = length(predators)),
  Predator = rep(predators, each = nrow(prey_numeric)),
  DietFraction = unlist(prey_numeric[, predators]),
  stringsAsFactors = FALSE
)

# --- 4. Compute B for each cultural period ----

Dataset_CC_list <- list()

for (culture in df_biomass_corrected$Culture) {
  
  # Biomass in each cultural period
  EB_row <- df_biomass_corrected[df_biomass_corrected$Culture == culture,
                                 c("Small_mean","Medium_mean","MediumLarge_mean","Large_mean")]
  EB <- data.frame(
    SizeCategory = c("Small","Medium","MedLarge","Large"),
    EdibleBiomass = as.numeric(EB_row[1,])
  )
  
  # Inicialise 
  B_table <- data.frame()
  
  # For each herbivore weight class
  for (size in EB$SizeCategory) {
    edible <- EB$EdibleBiomass[EB$SizeCategory == size]
    
    # Filter diet
    diet_row <- prey_long[prey_long$SizeCategory == size, ]
    
    # Total deman to normalise
    total_demand <- sum(diet_row$DietFraction, na.rm = TRUE)
    
    # Available biomass for predator
    diet_row$BiomassAvailable <- edible * diet_row$DietFraction / total_demand
    
    B_table <- rbind(B_table, diet_row)
  }
  
  # Compute for each predator
  Dataset_CC_culture <- aggregate(BiomassAvailable ~ Predator, data = B_table, sum)
  Dataset_CC_culture$AnnualBiomass <- db_values$AnnualBiomass[match(Dataset_CC_culture$Predator, db_values$Predator)]
  Dataset_CC_culture$Dataset_CC <- Dataset_CC_culture$BiomassAvailable / Dataset_CC_culture$AnnualBiomass
  Dataset_CC_culture$Culture <- culture
  
  # Save
  Dataset_CC_list[[culture]] <- Dataset_CC_culture[, c("Culture","Predator","BiomassAvailable","AnnualBiomass","Dataset_CC")]
}

# Bind cultural periods
Dataset_CC_final <- do.call(rbind, Dataset_CC_list)

# --- 5. Print, save and plot results ----
print(Dataset_CC_final)

write.csv(Dataset_CC_final, "Output_CC.csv")



# Load data
Dataset_CC_human<-read.xlsx("Dataset_CC.xlsx", sheet= "Human_CC") # This includes the Output_CC.csv file (line 678) for Pleistocene and Holocene
# mammal compositions and for a 30, 40, 50 and 60% of human mean intake 
head(Dataset_CC_human)

Dataset_CC_human_plot <- Dataset_CC_human %>%
  mutate(
    Composition_Label = ifelse(Scenario == "Pleistocene_Actual", "Pre-extinctions", "Post-extinctions"),
    Chron = factor(Culture, levels = chronological_order),
    Diet = as.factor(Diet)  # As factor
  )
y_range_Dataset_CC<- c(0.001, 0.2)
# Plot 
plot_Dataset_CC <- Dataset_CC_human_plot %>%
  ggplot(aes(x = Chron, y = CC, color = as.factor(Diet), group = as.factor(Diet))) +
  geom_point(size = 4, alpha = 0.8) +
  geom_line(size = 1, alpha = 0.6) +
  facet_wrap(~ Composition_Label, ncol = 1) +
  scale_color_manual(values = c("30" = "#1b9e77", "40" = "#d95f02", "50" = "#7570b3", "60" = "#e7298a"),
                     name = "Diet (% meat)") +
  labs(title = "Human Carrying Capacity by Diet Percentage",
       subtitle = "Top: Pre-extinction | Bottom: Post-extinction\nColors indicate meat percentage in diet",
       x = "Chronological Period", 
       y = "Carrying Capacity") +
  scale_y_continuous(limits = y_range_Dataset_CC) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        strip.background = element_rect(fill = "#f0f0f0", color = "grey50"),
        strip.text = element_text(face = "bold", size = 10, color = "black"),
        panel.spacing = unit(1.2, "lines"),
        legend.position = "bottom",
        plot.subtitle = element_text(color = "gray40", size = 10))

# Combine
combined_plot <- plot_biomass + plot_Dataset_CC +
  plot_annotation(
    title = "Herbivore Biomass and Human Carrying Capacity Through Time",
    subtitle = "Left: Herbivore biomass by size class | Right: Human carrying capacity by diet percentage",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

print(combined_plot)






# Order of cultural periods
chronological_order <- c("Aurignacian", "Gravettian", "Solutrean", "Magdalenian", "Epipaleolithic")

# Combine biomass data to assess the scenario in which extinctions occur between the Magdalenian and Epipaleolithic
transition_biomass <- transition_data_labeled %>%
  mutate(Chron = factor(Chron, levels = chronological_order)) %>%
  # Filtrar: Pleistoceno para primeros 3, Holoceno para últimos 2
  filter(
    (Chron %in% c("Aurignacian", "Gravettian", "Solutrean") & Composition_Label == "Pre-extinctions") |
      (Chron %in% c("Magdalenian", "Epipaleolithic") & Composition_Label == "Post-extinctions")
  ) %>%
  mutate(Scenario = ifelse(Composition_Label == "Pre-extinctions", "Pleistocene", "Holocene"))

# Plot 
plot_biomass_realistic <- ggplot(transition_biomass, aes(x = Chron, y = Mean, fill = Biomass_Type)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  scale_fill_manual(values = c("Large" = "#E78AC3", "MediumLarge" = "#8DA0CB", 
                               "Medium" = "#FC8D62", "Small" = "#66C2A5"),
                    name = "Herbivore Size Class") +
  labs(title = "Realistic Herbivore Biomass Transition",
       subtitle = "Pleistocene fauna → Holocene fauna transition at Magdalenian",
       x = "Chronological Period", 
       y = "Mean Biomass") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom",
        plot.subtitle = element_text(color = "gray40", size = 10))

# Create combined data
transition_Dataset_CC <- Dataset_CC_human %>%
  mutate(
    Chron = factor(Culture, levels = chronological_order),
    Scenario = ifelse(Scenario == "Pleistocene_Actual", "Pleistocene", "Holocene")
  ) %>%
  filter(
    (Chron %in% c("Aurignacian", "Gravettian", "Solutrean") & Scenario == "Pleistocene") |
      (Chron %in% c("Magdalenian", "Epipaleolithic") & Scenario == "Holocene")
  )

y_range_Dataset_CC_realistic<- c(0, 0.18) # Range of Y axis labels
culturas <- unique(transition_Dataset_CC$Chron)
resultados <- data.frame()

for(cultura in culturas) {
  datos_cultura <- transition_Dataset_CC %>% filter(Chron == cultura)
  
  media_Dataset_CC <- mean(datos_cultura$CC, na.rm = TRUE)
  sd_Dataset_CC <- sd(datos_cultura$CC, na.rm = TRUE)
  
  nueva_fila <- data.frame(
    Chron = cultura,
    Dataset_CC_mean = media_Dataset_CC,
    Dataset_CC_sd = sd_Dataset_CC
  )
  
  resultados <- rbind(resultados, nueva_fila)
}

resultados <- resultados %>%
  mutate(Chron = factor(Chron, levels = chronological_order))

# See results
print("Mean and SD by culture:")
print(resultados)

# Plot
plot_Dataset_CC_realistic <- ggplot(resultados, aes(x = Chron, y = Dataset_CC_mean)) +
  geom_point(size = 4, color = "#377EB8", alpha = 0.8) +
  geom_line(aes(group = 1), size = 1, color = "#377EB8", alpha = 0.6) +
  geom_errorbar(aes(ymin = Dataset_CC_mean - Dataset_CC_sd, ymax = Dataset_CC_mean + Dataset_CC_sd), 
                width = 0.2, color = "#377EB8", size = 0.8, alpha = 0.7) +
  labs(title = "Human Carrying Capacity (Mean ± SD)",
       subtitle = "Pleistocene Dataset_CC → Holocene Dataset_CC transition at Magdalenian\nShowing mean and standard deviation across all diets", 
       x = "Chronological Period", 
       y = "Carrying Capacity") +
  scale_y_continuous(limits = c(0, 0.17)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        legend.position = "none",
        plot.subtitle = element_text(color = "gray40", size = 10))

# Combine
combined_realistic_plot <- plot_biomass_realistic + plot_Dataset_CC_realistic +
  plot_annotation(
    title = "Transitional scenario: Pleistocene to Holocene Faunal Composition",
    subtitle = "Vertical red line marks the transition point at Magdalenian period",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

print(combined_realistic_plot)



df<- read.xlsx("Prey_preferences.xlsx", sheet= "dB")
head(df)

# Load data
df <- read.xlsx("Prey_preferences.xlsx", sheet = "dB")
head(df)

# Define species order by body size
species_order <- c("Ursus.arctos", "Crocuta.crocuta", "Panthera.pardus", 
                   "Homo.sapiens", "Canis.lupus", "Cuon.alpinus", 
                   "Lynx.pardinus", "Felis.silvestris", "Meles.meles", 
                   "Vulpes.vulpes", "Martes.martes")

# Convert data to long format
df_long <- df %>%
  pivot_longer(cols = -X1, names_to = "Species", values_to = "Value") %>%
  filter(Value > 0)  # Keep only connections with Value > 0

# Create edge dataframe
edges <- data.frame(
  from = df_long$X1,
  to = df_long$Species,
  weight = df_long$Value
)

# Create node dataframe
size_categories <- unique(df$X1)
species <- colnames(df)[-1]  # All species except the first column (X1)

nodes <- data.frame(
  name = c(size_categories, species),
  type = c(rep("Category", length(size_categories)), 
           rep("Species", length(species)))
)

# Create graph
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# Define manual layout to order species
layout_matrix <- matrix(0, nrow = length(nodes$name), ncol = 2)

# Positions for size categories (left column)
category_positions <- seq(1, 0, length.out = length(size_categories))
for(i in 1:length(size_categories)) {
  idx <- which(nodes$name == size_categories[i])
  layout_matrix[idx, ] <- c(0, category_positions[i])
}

# Positions for species (right column) in the specified order
species_positions <- seq(1, 0, length.out = length(species_order))
for(i in 1:length(species_order)) {
  idx <- which(nodes$name == species_order[i])
  if(length(idx) > 0) {
    layout_matrix[idx, ] <- c(1, species_positions[i])
  }
}

# Set node colors and shapes
V(g)$color <- ifelse(V(g)$type == "Category", "lightblue", "lightcoral")
V(g)$shape <- ifelse(V(g)$type == "Category", "square", "circle")
V(g)$size <- ifelse(V(g)$type == "Category", 15, 12)

# Set edge thickness according to weight
E(g)$width <- E(g)$weight / 100  # Adjusted for better visualization

# Create the plot
plot(g, 
     layout = layout_matrix,
     vertex.label = V(g)$name,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.frame.color = "gray",
     edge.color = "gray50",
     edge.curved = 0.2,
     main = "Demanded biomass",
     margin = c(0, 0, 0, 0))

# Add legend
legend("topleft",
       legend = c("Body Size Categories", "Species"),
       pch = c(15, 19),
       col = c("lightblue", "lightcoral"),
       pt.cex = 1.5,
       cex = 0.8,
       bty = "n")
