# brGDGT-based NPP Modelling and Analysis

## 📌 Overview

This repository contains an R workflow for analysing branched GDGT (brGDGT) compounds and modelling **Net Primary Productivity (NPP)** using statistical and machine learning approaches.
It also includes the full ecological modelling workflow to recosntruct herbivore biomass and carrying capacity of secondary consumers.

There are two scripts:

* Script_1
* Script_2

The Script_1 performs:

* Exploratory data analysis
* Statistical modelling (linear, quadratic, GAM)
* Random Forest modelling with three validation strategies
* Independent dataset validation
* NPP reconstruction from Padul-15-05 core
* Comparison with climate-derived NPP (Miami model)
* Summed Probability Distribution (SPD) analysis for archaeological data
* Correlation analysis
* Summary statistics and plots


The Script_2 performs:

* Herbivore density estimations from primary productivity
* Population structure modelling using survival curves
* Biomass partitioning across body size classes
* Secondary consumers carrying capacity calculations
* Summary statistics and plots

---

## 📂 Data Requirements

The Script_1 expects an Excel file:

```
Dataset_1.xlsx

This contains:

* A global dataset of present-day measured brGDGTs and associated climatic, edaphic and environmental variables
* Areas with in-field NPP measurements
* Climate conditions reconstructed in Padul
* Archaeological sites, levels and dates for the summed probability distributions

```

The Script_2 expects two additional Excel files:

```
Dataset_2.xlsx

This contains:

* Herbivore traits (species, body mass, life history traits, etc.),
* The wild mammal species list recovered in the region during the Late Pleistocene and Early Holocene (including sites and references)
* NPP values (actual and potential) obtained from Script_1
* Primary and secondary consumer biomass summaries

```
```
Prey_preferences.xlsx

This contains:

* Predator diet composition and annual biomass requirements

```

## 📦 R Packages

A package to apply the random-forest NPP predictions from brGDGTs is under construction

---

## ⚠️ Notes

* Ensure working directory is set correctly:

```r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
```

* All files available in this repository should be saved into the same folder/directory

---

## 📬 Contact 

# mav54#cam.ac.uk
