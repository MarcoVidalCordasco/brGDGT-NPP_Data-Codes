# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# SCRIPT 2 — Herbivore density estimations
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#
# This script reproduces the analyses and figures associated
# with the herbivore biomass and human carrying capacity estimations
# based on 1) the composition primary and secondary consumer guilds, 2) Their
# body masses and energy requirements, 3) prey-size preferences, and 4) NPP
# estimated from paleoclimate reconstructions and brGDGTs recovered from sediments.
# 
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


### a) HERBIVORE DENSITIES ####

# ---- 1. Load data from Excel file ----

# Read data from the "Fauna_Padul_Pleistocene" sheet
df <- read.xlsx("Dataset_2.xlsx",
                rowNames = FALSE,
                colNames = TRUE,
                sheet = "Fauna_Padul_Holocene")

# Check 
head(df)


# ---- 2. Compute THB values (actual and potential) ----

# THB_act = Total Herbivore biomass based on actual (brGDGT-derived) NPP 
# The formula is derived from the regression in the paper Vidal-Cordasco et al., 2022

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

species_list<- c("Mammuthus.primigenius",  # For Holocene scenario, replace this species from this list
                 "Equus.hydruntinus", # For Holocene scenario, replace this species from this list
                 "Equus.ferus",
                 "Bison.priscus", # For Holocene scenario, replace this species from this list
                 "Bos.primigenius",
                 "Cervus.elaphus",
                 "Capreolus.capreolus",
                 "Capra.pyrenaica",
                 "Rupicapra.rupicapra",
                 "Sus.scrofa", 
                 "Oryctolagus",
                 "Lepus.granatensis")

# ---- 3. Density data ----
# These are the outputs obtained in a) HERBIVORE DENSITIES and saved as 'Outputs_Density_Pleistocene.csv' or 'Outputs_Density_Holocene.csv'
df_PC <- read.xlsx("Dataset_2.xlsx",
                rowNames = FALSE,
                colNames = TRUE,
                sheet = "Potential_Pleistocene")
head(df_PC)

# For actual, brGDGT-derived NPP:
chrons <- c(   "Gravettian_act" ,   
             "Solutrean_act",     "Magdalenian_act",  "Epipaleolithic_act")

# For potential, climate-derived NPP:
chrons <- c(   "Gravettian_pot" ,   
               "Solutrean_pot",     "Magdalenian_pot",  "Epipaleolithic_pot")

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


#### Supplementary Figure 15 ####
# Example with one species:
df<-read.xlsx("Dataset_2.xlsx", rowNames=FALSE, 
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

  df<-read.xlsx("Dataset_2.xlsx", rowNames=FALSE, 
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
  
  
  # ---- 6. Compute biomass of small, medium, mediumlarge and large herbivores ----
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
# ---- 7. Save outputs ----
write.csv(outputs_bm, "outputs__Biomass_Potential_Pleisto.csv")


#### Supplementary Tables 6, 7 and 8 ####
# ---- 8. Summary statistics ----
# Convert to data.frame

Biomass_df<- read.xlsx("Dataset_2.xlsx", sheet = "Biomass")
head(Biomass_df)

Biomass_df<- subset(Biomass_df, Biomass_df$Composition=="Holocene")
Biomass_df<- subset(Biomass_df, Biomass_df$NPP=="ACT")
outputs_bm <- as.data.frame(Biomass_df)
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
write.csv(summary_biomass, "summary_biomass.csv") ## These outputs are available in sheet Biomass in Dataset_2.xlsx





#### c) CARRYING CAPACITY FOR SECONDARY CONSUMERS ####

# ---- 1. Load data from Excel file ----
df_biomass<- read.xlsx("Dataset_2.xlsx", sheet = "Summary_Biomass")

head(df_biomass)

df_biomass<- subset(df_biomass, df_biomass$Composition=="Pleistocene") # Replace with Holocene if necessary 

Prey_preferences<- read.xlsx("Prey_preferences.xlsx", rowNames=FALSE, 
                             colNames=TRUE, sheet="Diet_60_Pleistocene")  # Replace with % of human meat intake [30-60, by 10]
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

write.csv(Dataset_CC_final, "Output_CC_Pleistocene_60.csv")



# Load and summarise data 


CC_df<-read.xlsx("Dataset_2.xlsx", sheet= "Carnivores_CC")

head(CC_df)

summary_table_df <- CC_df %>%
  group_by(Paleocommunity, NPP, Culture, Species) %>%
  summarise(mean_CC = mean(CC, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = Species,
    values_from = mean_CC
  )

write.xlsx(summary_table_df, "summary_stats_cc.xlsx")


#### Main Figure 9 ####
# BIOMASS PLOT

getwd()
Biomass_herbivores_df<- read.xlsx("Dataset_2.xlsx", sheet="Summary_Biomass")
head(Biomass_herbivores_df)


# 1. Filter NPP = brGDGT
df_plot <- Biomass_herbivores_df %>%
  filter(NPP == "brGDGT") %>%
  
  # 2. Fix culture order
  mutate(Culture = factor(Culture,
                          levels = c("Gravettian",
                                     "Solutrean",
                                     "Magdalenian",
                                     "Epipaleolithic")))

# 3. Reshape to long format (only the 3 categories you want)
df_long <- df_plot %>%
  dplyr::select(Composition, Culture,
         Small_mean, Medium_mean, MediumLarge_mean, Large_mean) %>%
  
  pivot_longer(cols = c(Small_mean, Medium_mean, MediumLarge_mean, Large_mean),
               names_to = "Size",
               values_to = "Biomass") %>%
  
  mutate(Size = recode(Size,
                       Small_mean = "Small",
                       Medium_mean = "Medium",
                       MediumLarge_mean = "Medium-large",
                       Large_mean="Large"))

# 4. Plot
ggplot(df_long, aes(x = Culture, y = Biomass, fill = Size)) +
  
  geom_col(position = position_dodge(width = 0.8)) +
  
  facet_wrap(~ Composition) +
  
  labs(x = "Culture",
       y = "Herbivore biomass (kg km⁻² yr⁻¹)",
       fill = "Body size class") +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )




# Filter Homo sapiens: 
df_hs <- CC_df %>%
  filter(Species == "Homo.sapiens")

# Compute mean and SD:
summary_df <- df_hs %>%
  group_by(Culture, NPP, Paleocommunity) %>%
  summarise(
    mean_CC = mean(CC, na.rm = TRUE),
    sd_CC   = sd(CC, na.rm = TRUE),
    .groups = "drop"
  )

summary_df<- subset(summary_df, summary_df$NPP=="brGDGT")

# Plot
ggplot(summary_df, aes(x = Culture, y = mean_CC, fill = Paleocommunity)) +
  geom_point()+
  geom_errorbar(
    aes(ymin = mean_CC - sd_CC, ymax = mean_CC + sd_CC),
    width = 0.2
  ) +
  facet_wrap(~ NPP) +
  labs(
    x = "Culture",
    y = "Carrying Capacity (CC)",
    fill = "Paleocommunity"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )





# Filter for Homo sapiens only
homo_data <- CC_df %>%
  filter(Species == "Homo.sapiens") %>%
  # Clean up culture names if needed (remove "_act" suffix)
  mutate(Culture = gsub("_act", "", Culture)) %>%
  # Ensure the culture order
  mutate(Culture = factor(Culture, 
                          levels = c("Gravettian", "Solutrean", 
                                     "Magdalenian", "Epipaleolithic")))

# Calculate summary statistics
summary_homo <- homo_data %>%
  group_by(Paleocommunity, Culture) %>%
  summarise(
    Mean_CC = mean(CC, na.rm = TRUE),
    SD_CC = sd(CC, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot with error bars
ggplot(summary_homo, aes(x = Culture, y = Mean_CC, color = Paleocommunity, group = Paleocommunity)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Mean_CC - SD_CC, ymax = Mean_CC + SD_CC), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  labs(
    title = "Homo sapiens CC across Cultures",
    subtitle = "Comparison between Paleocommunity A and B",
    x = "Culture",
    y = "Mean CC (± SD)",
    color = "Paleocommunity"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )




# Load data for food-web plot
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



# BIOMASS
Biomass_df<- read.xlsx("Dataset_2.xlsx", sheet= "Biomass")
head(Biomass_df)

# Reshape the data to long format for easier manipulation
biomass_long <- Biomass_df %>%
  # Select relevant columns
  dplyr::select(Composition, Species, Chron, Small_Biomass, Medium_Biomass, 
         MediumLarge_Biomass, Large_Biomass) %>%
  # Pivot to long format
  pivot_longer(cols = c(Small_Biomass, Medium_Biomass, MediumLarge_Biomass, Large_Biomass),
               names_to = "Size_Class",
               values_to = "Biomass") %>%
  # Clean size class names
  mutate(Size_Class = case_when(
    Size_Class == "Small_Biomass" ~ "Small",
    Size_Class == "Medium_Biomass" ~ "Medium",
    Size_Class == "MediumLarge_Biomass" ~ "MediumLarge",
    Size_Class == "Large_Biomass" ~ "Large"
  )) %>%
  # Filter out zero biomass values to clean up pies
  filter(Biomass > 0)

# Aggregate across all periods and compositions
biomass_aggregated <- biomass_long %>%
  group_by(Size_Class, Species) %>%
  summarise(Total_Biomass = sum(Biomass, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Size_Class) %>%
  mutate(Percentage = Total_Biomass / sum(Total_Biomass) * 100) %>%
  arrange(Size_Class, desc(Percentage))

# Separate pies for each Composition (Actual_Pleistocene vs Actual_Holocene)
biomass_by_composition <- biomass_long %>%
  group_by(Composition, Size_Class, Species) %>%
  summarise(Total_Biomass = sum(Biomass, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Composition, Size_Class) %>%
  mutate(Percentage = Total_Biomass / sum(Total_Biomass) * 100) %>%
  arrange(Composition, Size_Class, desc(Percentage))

# Separate pies for each Chron (time period)
biomass_by_chron <- biomass_long %>%
  group_by(Chron, Size_Class, Species) %>%
  summarise(Total_Biomass = sum(Biomass, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Chron, Size_Class) %>%
  mutate(Percentage = Total_Biomass / sum(Total_Biomass) * 100) %>%
  arrange(Chron, Size_Class, desc(Percentage))

# Create pie charts 
create_pie <- function(data, size_class) {
  data %>%
    filter(Size_Class == size_class) %>%
    ggplot(aes(x = "", y = Percentage, fill = Species)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    labs(title = size_class) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) 
}

# Create all four pies
pie_small <- create_pie(biomass_aggregated, "Small")
pie_medium <- create_pie(biomass_aggregated, "Medium")
pie_mediumlarge <- create_pie(biomass_aggregated, "MediumLarge")
pie_large <- create_pie(biomass_aggregated, "Large")

grid.arrange(pie_small, pie_medium, pie_mediumlarge, pie_large)
