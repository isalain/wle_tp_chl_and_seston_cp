################################################################################
# Hierarchical Bayesian models for TP-Chlorophyll relationships and seston C:P
################################################################################

# SNIPPET 1: DATA PREPROCESSING
####################################

# list of R packages used
library(dplyr)
library(lubridate)
library(brms)
library(tidyverse)
library(posterior)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(viridis)
library(ellipse)

# ==============================================================================
# WORKING DIRECTORY AND DOWNLOAD DATA
# ==============================================================================

# the local working directory that needs to be adjusted 
# at a different computer
work_dir <- "C:/Users/aisab/Desktop/tp_chl"
if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}
setwd(work_dir)

# permanent NCEI links where data are archived
# except "Lake_Erie_Merged_Data_2023.csv" the 5th file which
# was circulated internally

urls <- c(
  "https://www.nodc.noaa.gov/archive/arc0204/0254720/1.1/data/0-data/noaa-glerl-erie-habs-field-sampling-results-2020-2021.csv",
  "https://www.nodc.noaa.gov/archive/arc0152/0209116/1.1/data/0-data/lake_erie_habs_field_sampling_results_2019.csv",
  "https://www.nodc.noaa.gov/archive/arc0135/0187718/2.2/data/0-data/lake_erie_habs_field_sampling_results_2012_2018_v2.csv",
  "https://www.ncei.noaa.gov/data/oceans/archive/arc0225/0292222/1.1/data/0-data/noaa-glerl-erie-habs-field-sampling-results-2022.csv"
)

for (i in seq_along(urls)) {
  filename <- paste0("file_", i, ".csv")
  destination <- file.path(work_dir, filename)
  if (!file.exists(destination)) {
    download.file(urls[i], destination, mode = "wb")
  }
}

# ==============================================================================
# DATA LOADING AND PREPROCESSING
# ==============================================================================

process_dataset <- function(data, select_cols, simp_names) {
  for (i in seq_along(select_cols)) {
    colnames(data)[select_cols[i]] <- simp_names[i]
  }
  return(subset(data, select = select_cols))
}

# the numeric values depict a column within which each variable is located
# for example 1 represents 'Date' for all the five files and 
# 26 represent SRP for the 1st file, while SRP is in column 28 for the 5th file
SelectCol1 <- c(1, 2, 5, 23, 24, 25, 26, 29)
SelectCol2 <- c(1, 2, 5, 23, 24, 25, 26, 29)
SelectCol3 <- c(1, 2, 5, 24, 25, 26, 27, 31)
SelectCol4 <- c(1, 2, 5, 22, 25, 26, 27, 30)
SelectCol5 <- c(1, 2, 3, 23, 26, 27, 28, 31)

SimpNames <- c("Date", "Site", "Category", "Chla", "TP", "TDP", "SRP", "POC")

RawData1 <- read.csv(file.path(work_dir, "file_1.csv"), check.names = FALSE)
RawData1 <- process_dataset(RawData1, SelectCol1, SimpNames)

RawData2 <- read.csv(file.path(work_dir, "file_2.csv"), check.names = FALSE)
RawData2 <- process_dataset(RawData2, SelectCol2, SimpNames)

RawData3 <- read.csv(file.path(work_dir, "file_3.csv"), check.names = FALSE)
RawData3 <- process_dataset(RawData3, SelectCol3, SimpNames)

RawData4 <- read.csv(file.path(work_dir, "file_4.csv"), check.names = FALSE)
RawData4 <- process_dataset(RawData4, SelectCol4, SimpNames)

if (file.exists(file.path(work_dir, "Lake_Erie_Merged_Data_2023.csv"))) {
  RawData5 <- read.csv(file.path(work_dir, "Lake_Erie_Merged_Data_2023.csv"), check.names = FALSE)
  RawData5 <- process_dataset(RawData5, SelectCol5, SimpNames)
  RawData <- rbind(RawData1, RawData2, RawData3, RawData4, RawData5)
} else {
  RawData <- rbind(RawData1, RawData2, RawData3, RawData4)
}

# ==============================================================================
# DATA FILTERING
# ==============================================================================

RawData <- RawData[RawData$Category == "Surface", ]

known_sites <- c("WE2", "WE4", "WE6", "WE8", "WE9", 
                 "WE12", "WE13", "WE14", "WE15", "WE16")

RawData <- RawData[RawData$Site %in% known_sites, ]

# ==============================================================================
# DATA TYPE CONVERSIONS
# ==============================================================================

# creating new colums for dates

RawData$Date <- as.Date(RawData$Date, format = "%m/%d/%Y")
RawData$Year <- year(RawData$Date)  # kept for reference only
RawData$Month <- factor(format(RawData$Date, "%b"), 
                        levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))

convert_to_numeric <- function(x) {
  x[grep("<", x)] <- NA
  as.numeric(x)
}

numeric_cols <- c("Chla", "TP", "TDP", "SRP", "POC")
for (col in numeric_cols) {
  if (col %in% names(RawData)) {
    RawData[[col]] <- convert_to_numeric(RawData[[col]])
  }
}

# ==============================================================================
# DERIVED VARIABLES
# ==============================================================================
# e.g., log transformations and Seston C:P

RawData <- RawData %>%
  mutate(
    PP = ifelse(TP >= TDP, TP - TDP, NA),
    Seston_C_to_P = ((POC / 12.0107) * 1000) / ((PP) / 30.973761),
    log_Chla = log(Chla),
    log_TP = log(TP),
    log_CP = log10(Seston_C_to_P),
    Site_factor = factor(Site, levels = known_sites),
    Month_factor = factor(Month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  )

# ==============================================================================
# OUTLIER REMOVAL
# ==============================================================================

RawData <- RawData %>% 
  filter(TP < 1000 | is.na(TP)) %>% 
  filter(Chla < 1000 | is.na(Chla))

# ==============================================================================
# CREATE ANALYSIS DATASETS
# ==============================================================================

RawData_clean <- RawData %>%
  filter(!is.na(Chla) & !is.na(TP) & !is.na(POC) & !is.na(PP) & 
           !is.na(Month_factor) & !is.na(Site_factor)) %>%
  select(Date, Site, Year, Chla, TP, TDP, PP, POC,
         Month, Seston_C_to_P,
         log_Chla, log_TP, log_CP,
         Site_factor, Month_factor)

model_data_tpchl <- RawData_clean %>%
  filter(!is.na(log_Chla) & !is.na(log_TP))

seston_data <- RawData_clean %>%
  filter(!is.na(Seston_C_to_P) & !is.na(log_CP))

# ==============================================================================
# SETTING UP COLOR PALETTES
# ==============================================================================

site_colors <- c(
  "WE2" = "#440154", "WE4" = "#3B528B", "WE6" = "#21908C",
  "WE8" = "#5DC863", "WE9" = "#FDE725", "WE12" = "#E8601C",
  "WE13" = "#7A0403", "WE14" = "#C73E4C", "WE15" = "#8B4789", 
  "WE16" = "#2C7FB8"
)

month_colors <- c(
  "Apr" = "#440154", "May" = "#3B528B", "Jun" = "#21908C", 
  "Jul" = "#5DC863", "Aug" = "#FDE725", "Sep" = "#E8601C", 
  "Oct" = "#7A0403"
)
