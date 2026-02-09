rm(list = ls())

# Load packages -----------------------------------------------------------

library(terra)
library(dtplyr)
library(dplyr)
library(usdm)
library(tidyr)
library(tidyterra)
library(purrr)
library(data.table)
library(ggplot2)
library(ggfortify)
library(viridis)
library(visdat)
library(foreach)
library(doFuture)
library(exactextractr)
library(sf)
library(stringr)
library(mclust)
library(RJDemetra)
library(mgcv)
library(mgcViz)
library(corrplot)
library(Boruta)
library(xgboost)
library(skimr)
library(missRanger)
library(fitdistrplus)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
library(ranger)
library(pdp)
library(gridExtra)
library(ggalluvial)
library(DescTools)
library(fastglm)
library(blockCV)
library(pROC)
library(patchwork)
library(geobr)
library(tidyverse)
library(cowplot)
library(rnaturalearth)
library(ggspatial) # For North Arrow & Scale Bar
library(scales) # For axis label formatting

# Load the Cerrado geometry
Cerrado <- readRDS("Data/Cerrado.rds")
print(Cerrado)

# --- Define a metric projection for Brazil ---
# SIRGAS 2000 / Brazil Polyconic (EPSG:5880) is a good choice for the whole country.
metric_crs <- "EPSG:5880"

# Select only the geometry object from "Cerrado"
Cerrado <- Cerrado %>% dplyr::select(geometry)
Cerrado <- sf::st_transform(Cerrado, metric_crs)

# Organizing predictors ---------------------------------------------------

## Land use and land cover -------------------------------------------------

legend_dt <- read.delim2(
  "Data/Codigos-da-legenda-colecao-10.csv"
)

lulc_df <- readRDS(
  "Data/lulc_predictors.rds"
)

grid_polygons <- terra::vect("grid_polygons.shp")
# grid_polygons <- terra::project(grid_polygons, "EPSG:5880")
centroids_grid <- centroids(grid_polygons)
grid_df <- terra::as.data.frame(centroids_grid, geom = "xy")
names(grid_df)[1] <- "plot_id"
lulc_df <- left_join(lulc_df, grid_df, by = "plot_id")

lulc_df %>%
  group_by(x, y, plot_id) %>%
  summarise(savanna = mean(savanna, na.rm = T)) %>%
  ggplot(aes(x = x, y = y, z = savanna, fill = savanna)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")

# Load one of the reprojected binary rasters to use as a template
template_raster <- rast("~/Documents/ERA5/precip.grib")[[1]]
print(template_raster)

# SIRGAS 2000 / Brazil Polyconic (EPSG:5880) is a good choice for the whole country.
metric_crs <- "EPSG:5880"
template_raster <- terra::project(
  template_raster,
  metric_crs,
  threads = T,
  method = "bilinear"
)

template_raster <- crop(template_raster, Cerrado) %>% mask(Cerrado)


# Get list of years and columns to process
years <- sort(unique(lulc_df$year))
# Identify metric columns (exclude ID/Coords/Year)
metric_cols <- setdiff(names(lulc_df), c("plot_id", "x", "y", "year"))

# THE RE-ALIGNMENT FUNCTION

process_year_lulc <- function(curr_year) {
  message(paste("Processing Year:", curr_year, "..."))

  # A. Filter data for this year
  df_year <- lulc_df %>% filter(year == curr_year)

  # B. Rasterize the "Dirty" Data
  # We convert the points to a raster.
  # terra::rasterize allows creating a multi-layer raster (one for each metric) at once.
  v_year <- vect(df_year, geom = c("x", "y"), crs = "EPSG:4326")

  # 3. Define the Correct Extent (Padding for Centers)
  # If x is the center, the left edge is x - res/2
  min_x <- min(df_year$x) - (res_x / 2)
  max_x <- max(df_year$x) + (res_x / 2)
  min_y <- min(df_year$y) - (res_y / 2)
  max_y <- max(df_year$y) + (res_y / 2)

  # 4. Create the Template Raster
  # We specify the extent manually to ensure perfect alignment
  r_temp_dirty <- rast(
    xmin = min_x,
    xmax = max_x,
    ymin = min_y,
    ymax = max_y,
    resolution = c(0.1, 0.1),
    crs = "EPSG:4326"
  )

  # Rasterize all metrics at once
  r_year_dirty <- rasterize(
    v_year,
    r_temp_dirty,
    field = metric_cols,
    fun = mean
  )

  # C. Resample to the MASTER Grid
  # This fixes the misalignment.
  # Use 'bilinear' because land use percentages are continuous.
  # Use 'near' if you had categorical classes.
  r_year_projected <- terra::project(r_year_dirty, metric_crs, threads = T)
  r_year_aligned <- resample(
    r_year_projected,
    template_raster,
    method = "bilinear"
  )

  v_target <- terra::as.data.frame(template_raster, xy = TRUE, cells = TRUE)
  # D. Extract at Correct Locations
  # This pulls the values from the aligned raster into our Target Points
  extracted_values <- terra::extract(
    r_1985_aligned,
    v_target[, c("x", "y")],
    ID = FALSE
  )

  # E. Re-assemble
  # Combine with the correct plot_ids and add the year
  result_df <- cbind(
    plot_id = v_target$cell,
    year = curr_year,
    extracted_values,
    x = v_target$x,
    y = v_target$y
  )

  return(as_tibble(result_df))
}

# Use map_dfr to run the function for every year and bind rows automatically
corrected_lulc_df <- map_dfr(years, process_year_lulc)

print(head(corrected_lulc_df))
skim(corrected_lulc_df)

corrected_lulc_df %>%
  group_by(x, y, plot_id) %>%
  summarise(savanna = mean(savanna_mean, na.rm = T)) %>%
  ggplot(aes(x = x, y = y, z = savanna, fill = savanna)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")

lulc_df <- corrected_lulc_df


skim(lulc_df)
usdm::vif(lulc_df[, c(3:18)])
vif_lulc <- vifstep(lulc_df[, c(3:18)], th = 2)
print(vif_lulc)

# Correlation plot
correlation_matrix <- cor(
  lulc_df[, c(3:18)],
  use = "pairwise.complete.obs",
  method = "spearman"
)

# Create the correlation plot using ellipses
# A combined plot: ellipses in the upper triangle, numbers in the lower
pdf("Figs/corrplot_lulc.pdf", paper = "a4r", width = 0, height = 0)
corrplot(
  correlation_matrix,
  method = "ellipse",
  type = "upper",
  order = "hclust", # Reorders variables to group similar ones
  addCoef.col = "black", # Adds correlation coefficients to the lower half
  tl.col = "black",
  tl.srt = 45, # Tweak text labels
  diag = FALSE
) # Hides the diagonal
dev.off()

## NDVI --------------------------------------------------------------------

# Read data
ndvi_df <- readRDS("Data/ndvi_predictor.rds")
skim(ndvi_df)

ndvi_df %>%
  group_by(longitude, latitude, plot_id) %>%
  summarise(NDVI = mean(NDVI, na.rm = T)) %>%
  ggplot(aes(x = longitude, y = latitude, z = NDVI, fill = NDVI)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")

## Population density ------------------------------------------------------

# Read data
pop_df <- readRDS("Data/pop_predictor.rds")
skim(pop_df)

pop_df %>%
  group_by(longitude, latitude, plot_id) %>%
  summarise(pop = mean(pop, na.rm = T)) %>%
  ggplot(aes(x = longitude, y = latitude, z = pop, fill = pop)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")

## Climate -----------------------------------------------------------------

# Read results
climate_df <- readRDS(
  "Data/climate_predictors.rds"
)
summary(climate_df[, 27])
climate_df <- climate_df[, -27]

usdm::vif(climate_df[, c(6:44)])
vif_climate <- vifstep(climate_df[, c(6:44)], th = 2)
print(vif_climate)

climate_df <- usdm::exclude(climate_df, vif_climate)
skim(climate_df)

climate_df %>%
  group_by(x, y, plot_id) %>%
  summarise(precip_min = mean(precip_min, na.rm = T)) %>%
  ggplot(aes(x = x, y = y, z = precip_min, fill = precip_min)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")

# Correlation plot
correlation_matrix <- cor(
  climate_df[sample(1:nrow(climate_df), 100000), c(6:15)],
  use = "pairwise.complete.obs",
  method = "spearman"
)

# Create the correlation plot using ellipses
# A combined plot: ellipses in the upper triangle, numbers in the lower
pdf("Figs/corrplot_climate.pdf", paper = "a4r", width = 0, height = 0)
corrplot(
  correlation_matrix,
  method = "ellipse",
  type = "upper",
  order = "hclust", # Reorders variables to group similar ones
  addCoef.col = "black", # Adds correlation coefficients to the lower half
  tl.col = "black",
  tl.srt = 45, # Tweak text labels
  diag = FALSE
) # Hides the diagonal
dev.off()

## Terrain variables -------------------------------------------------------
terrain_df <- readRDS("Data/terrain_predictors.rds")

usdm::vif(terrain_df[, c(4:8)])
vif_terrain <- vifstep(terrain_df[, c(4:8)], th = 2)
print(vif_terrain)

terrain_df <- usdm::exclude(terrain_df, vif_terrain)
skim(terrain_df)

ggplot(terrain_df, aes(x = x, y = y, z = elevation, fill = elevation)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")

# Correlation plot
correlation_matrix <- cor(
  terrain_df[, c(4:7)],
  use = "pairwise.complete.obs",
  method = "spearman"
)

# Create the correlation plot using ellipses
# A combined plot: ellipses in the upper triangle, numbers in the lower
pdf("Figs/corrplot_terrain.pdf", paper = "a4r", width = 0, height = 0)
corrplot(
  correlation_matrix,
  method = "ellipse",
  type = "upper",
  order = "hclust", # Reorders variables to group similar ones
  addCoef.col = "black", # Adds correlation coefficients to the lower half
  tl.col = "black",
  tl.srt = 45, # Tweak text labels
  diag = FALSE
) # Hides the diagonal
dev.off()


# Merge fire variables with predictors ----------------------------------

## Fire occurrence ---------------------------------------------------------

fire_occ_df <- readRDS(
  "Data/final_fire_occ_data.rds"
)

fire_occ_df <- dtplyr::lazy_dt(fire_occ_df, immutable = FALSE)

# Standardize names and variables
fire_occ_df <- dplyr::rename(fire_occ_df, x = longitude)
fire_occ_df <- dplyr::rename(fire_occ_df, y = latitude)
print(fire_occ_df)

climate_df <- climate_df %>%
  mutate(year = as.character(year))

lulc_df <- lulc_df %>%
  mutate(year = as.character(year))

# Join datasets
fire_occ_df <- fire_occ_df %>%
  left_join(climate_df, by = c("plot_id", "year", "month")) %>%
  left_join(lulc_df, by = c("plot_id", "year")) %>%
  left_join(ndvi_df, by = c("plot_id", "year")) %>%
  left_join(pop_df, by = c("plot_id", "year")) %>%
  left_join(terrain_df, by = c("plot_id")) %>%
  as_tibble()
fire_occ_df <- as_tibble(fire_occ_df)


fire_occ_df %>%
  group_by(x.x, y.x, plot_id) %>%
  summarise(savanna = mean(savanna_mean, na.rm = T)) %>%
  ggplot(aes(x = x.x, y = y.x, z = savanna, fill = savanna)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")

skim(fire_occ_df)

fire_occ_df <- dtplyr::lazy_dt(fire_occ_df)
fire_occ_df <- fire_occ_df %>%
  dplyr::rename(x = x.x) %>%
  dplyr::rename(y = y.x) %>%
  dplyr::select(
    !c(
      "x.y",
      "y.y",
      "x.x.x",
      "y.x.x",
      "longitude.x",
      "latitude.x",
      "longitude.y",
      "latitude.y",
      "x.y.y",
      "y.y.y"
    )
  ) %>%
  as_tibble()


fire_occ_df <- fire_occ_df[complete.cases(fire_occ_df[, 9:33]), ]
skim(fire_occ_df)
fire_occ_df$year <- as.integer(fire_occ_df$year)
fire_occ_df <- dplyr::select(fire_occ_df, !"decade")

### Impute missing values ---------------------------------------------------
missRanger_fire_occ <- missRanger(
  as.data.frame(fire_occ_df),
  num.trees = 200,
  max.depth = 15,
  sample.fraction = 0.2,
  min.node.size = 10,
  data_only = FALSE,
  returnOOB = TRUE,
  keep_forests = FALSE
)
print(missRanger_fire_occ)
# - best average OOB imputation error: 0.1058643

# Rename land use variables
fire_occ_df <- missRanger_fire_occ$data

fire_occ_df <- fire_occ_df %>%
  dplyr::rename(savanna = savanna_mean) %>%
  dplyr::rename(forest_plantation = forest_plantation_mean) %>%
  dplyr::rename(grassland = grassland_mean) %>%
  dplyr::rename(pasture = pasture_mean) %>%
  dplyr::rename(sugar_cane = sugar_cane_mean) %>%
  dplyr::rename(mosaic_uses = mosaic_uses_mean) %>%
  dplyr::rename(urban = urban_mean) %>%
  dplyr::rename(soybean = soybean_mean) %>%
  dplyr::rename(non_vegetated = non_vegetated_mean) %>%
  dplyr::rename(forest = forest_mean) %>%
  dplyr::rename(temporary_crop = temporary_crop_mean) %>%
  dplyr::rename(perennial_crop = perennial_crop_mean) %>%
  dplyr::rename(water = water_mean) %>%
  dplyr::rename(wet_herbaceous = wet_herbaceous_mean) %>%
  dplyr::rename(native_area = native_area_mean) %>%
  dplyr::rename(anthropic_area = anthropic_area_mean) %>%
  as_tibble()

saveRDS(fire_occ_df, "Output/fire_occ_df_modeling.rds")

## Fire regime metrics -----------------------------------------------------

fire_metrics_df <- readRDS(
  "Data/final_clean_data.rds"
)

fire_metrics_df <- dtplyr::lazy_dt(fire_metrics_df, immutable = FALSE)

# Standardize names and variables
fire_metrics_df <- dplyr::rename(fire_metrics_df, x = longitude)
fire_metrics_df <- dplyr::rename(fire_metrics_df, y = latitude)
fire_metrics_df <- fire_metrics_df %>%
  mutate(year = as.integer(year))

# Join datasets
fire_metrics_df <- left_join(
  fire_metrics_df,
  missRanger_fire_occ$data,
  by = c("plot_id", "year", "month")
) %>%
  as_tibble()

skim(fire_metrics_df)

# Check if join preserved spatial properties
fire_metrics_df %>%
  group_by(x.x, y.x, plot_id) %>%
  summarise(lsm_c_cai_sd = mean(lsm_c_cai_sd)) %>%
  ggplot(aes(x = x.x, y = y.x, z = lsm_c_cai_sd, fill = lsm_c_cai_sd)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")

# Remove NAs
fire_metrics_df <- fire_metrics_df[complete.cases(fire_metrics_df[, 13:47]), ]
skim(fire_metrics_df)

# Remove duplicated variables
fire_metrics_df <- dtplyr::lazy_dt(fire_metrics_df)
fire_metrics_df <- fire_metrics_df %>%
  dplyr::rename(x = x.x) %>%
  dplyr::rename(y = y.x) %>%
  dplyr::select(!c("x.y", "y.y")) %>%
  as_tibble()

# Rename land use variables
fire_metrics_df <- fire_metrics_df %>%
  dplyr::rename(savanna = savanna_mean) %>%
  dplyr::rename(forest_plantation = forest_plantation_mean) %>%
  dplyr::rename(grassland = grassland_mean) %>%
  dplyr::rename(pasture = pasture_mean) %>%
  dplyr::rename(sugar_cane = sugar_cane_mean) %>%
  dplyr::rename(mosaic_uses = mosaic_uses_mean) %>%
  dplyr::rename(urban = urban_mean) %>%
  dplyr::rename(soybean = soybean_mean) %>%
  dplyr::rename(non_vegetated = non_vegetated_mean) %>%
  dplyr::rename(forest = forest_mean) %>%
  dplyr::rename(temporary_crop = temporary_crop_mean) %>%
  dplyr::rename(perennial_crop = perennial_crop_mean) %>%
  dplyr::rename(water = water_mean) %>%
  dplyr::rename(wet_herbaceous = wet_herbaceous_mean) %>%
  dplyr::rename(native_area = native_area_mean) %>%
  dplyr::rename(anthropic_area = anthropic_area_mean) %>%
  as_tibble()

saveRDS(fire_metrics_df, "Data/fire_metrics_df_modeling.rds")

# VIF analysis -----------------------------------------------
fire_occ_df <- readRDS("Data/fire_occ_df_modeling.rds")
vif_analysis <- vifstep(as.data.frame(fire_occ_df[, c(8:16, 33:38)]), th = 2)
print(vif_analysis)

# Time series decomposition -----------------------------------------------

## Fire occurrence ---------------------------------------------------------

fire_occ_df <- readRDS("Data/fire_occ_df_modeling.rds")
fire_occ_df <- lazy_dt(fire_occ_df)

fire_occ_ts <- fire_occ_df %>%
  group_by(year, month) %>%
  summarise(fire_occurrence = mean(fire_occurrence, na.rm = T)) %>%
  as_tibble() %>%
  pull(fire_occurrence) %>%
  ts(start = c(1985, 01), end = c(2024, 12), frequency = 12)

print(fire_occ_ts)

x13_fire_occ <- x13(fire_occ_ts)
print(x13_fire_occ)

pdf("Figs/decomp_x13_fire_occ.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_fire_occ)
dev.off()

pdf("Figs/residuals_x13_fire_occ.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(3, 2))
plot(x13_fire_occ$regarima)
par(mfrow = c(1, 1))
dev.off()

pdf("Figs/SIratio_x13_fire_occ.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_fire_occ$decomposition)
dev.off()

## Fire metrics ---------------------------------------------------------

fire_metrics_df <- readRDS("Data/fire_metrics_df_modeling.rds")

fire_metrics_df <- lazy_dt(fire_metrics_df)

# Aggregate data
fire_metrics_ts <- fire_metrics_df %>%
  group_by(year, month) %>%
  summarise(
    lsm_c_cai_sd = mean(lsm_c_cai_sd, na.rm = T),
    lsm_c_circle_mn = mean(lsm_c_circle_mn, na.rm = T),
    lsm_c_core_mn = mean(lsm_c_core_mn, na.rm = T),
    lsm_c_dcore_sd = mean(lsm_c_dcore_sd, na.rm = T),
    lsm_c_enn_sd = mean(lsm_c_enn_sd, na.rm = T),
    lsm_c_pd = mean(lsm_c_pd, na.rm = T),
    lsm_c_split = mean(lsm_c_split, na.rm = T)
  ) %>%
  as_tibble()

print(fire_metrics_ts)


### Fire Scar Integrity Variation (lsm_c_cai_sd) ----------------------------

cai_sd_ts <- fire_metrics_ts %>%
  pull(lsm_c_cai_sd) %>%
  ts(start = c(1985, 01), end = c(2024, 12), frequency = 12)

print(cai_sd_ts)

x13_cai_sd <- x13(cai_sd_ts)
print(x13_cai_sd)

pdf("Figs/decomp_x13_cai_sd.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_cai_sd)
dev.off()

pdf("Figs/residuals_x13_cai_sd.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(3, 2))
plot(x13_cai_sd$regarima)
par(mfrow = c(1, 1))
dev.off()

pdf("Figs/SIratio_x13_cai_sd.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_cai_sd$decomposition)
dev.off()

### Fire Scar Shape Regularity (lsm_c_circle_mn) ----------------------------

circle_mn_ts <- fire_metrics_ts %>%
  pull(lsm_c_circle_mn) %>%
  ts(start = c(1985, 01), end = c(2024, 12), frequency = 12)

print(circle_mn_ts)

x13_circle_mn <- x13(circle_mn_ts)
print(x13_circle_mn)

pdf("Figs/decomp_x13_circle_mn.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_circle_mn)
dev.off()

pdf("Figs/residuals_x13_circle_mn.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(3, 2))
plot(x13_circle_mn$regarima)
par(mfrow = c(1, 1))
dev.off()

pdf("Figs/SIratio_x13_circle_mn.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_circle_mn$decomposition)
dev.off()

### Mean Fire Scar Core Size (lsm_c_core_mn) ----------------------------

core_mn_ts <- fire_metrics_ts %>%
  pull(lsm_c_core_mn) %>%
  ts(start = c(1985, 01), end = c(2024, 12), frequency = 12)

print(core_mn_ts)

x13_core_mn <- x13(core_mn_ts)
print(x13_core_mn)

pdf("Figs/decomp_x13_core_mn.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_core_mn)
dev.off()

pdf("Figs/residuals_x13_core_mn.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(3, 2))
plot(x13_core_mn$regarima)
par(mfrow = c(1, 1))
dev.off()

pdf("Figs/SIratio_x13_core_mn.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_core_mn$decomposition)
dev.off()

### Fire Scar Cohesion Variation (lsm_c_dcore_sd) ----------------------------

dcore_sd_ts <- fire_metrics_ts %>%
  pull(lsm_c_dcore_sd) %>%
  ts(start = c(1985, 01), end = c(2024, 12), frequency = 12)

print(dcore_sd_ts)

x13_dcore_sd <- x13(dcore_sd_ts)
print(x13_dcore_sd)

pdf("Figs/decomp_x13_dcore_sd.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_dcore_sd)
dev.off()

pdf("Figs/residuals_x13_dcore_sd.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(3, 2))
plot(x13_dcore_sd$regarima)
par(mfrow = c(1, 1))
dev.off()

pdf("Figs/SIratio_x13_dcore_sd.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_dcore_sd$decomposition)
dev.off()

### Fire Scar Isolation Variation (lsm_c_enn_sd) ----------------------------

enn_sd_ts <- fire_metrics_ts %>%
  pull(lsm_c_enn_sd) %>%
  ts(start = c(1985, 01), end = c(2024, 12), frequency = 12)

print(enn_sd_ts)

x13_enn_sd <- x13(enn_sd_ts)
print(x13_enn_sd)

pdf("Figs/decomp_x13_enn_sd.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_enn_sd)
dev.off()

pdf("Figs/residuals_x13_enn_sd.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(3, 2))
plot(x13_enn_sd$regarima)
par(mfrow = c(1, 1))
dev.off()

pdf("Figs/SIratio_x13_enn_sd.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_enn_sd$decomposition)
dev.off()

### Fire Scar Density (lsm_c_pd) ----------------------------

pd_ts <- fire_metrics_ts %>%
  pull(lsm_c_pd) %>%
  ts(start = c(1985, 01), end = c(2024, 12), frequency = 12)

print(pd_ts)

x13_pd <- x13(pd_ts)
print(x13_pd)

pdf("Figs/decomp_x13_pd.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_pd)
dev.off()

pdf("Figs/residuals_x13_pd.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(3, 2))
plot(x13_pd$regarima)
par(mfrow = c(1, 1))
dev.off()

pdf("Figs/SIratio_x13_pd.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_pd$decomposition)
dev.off()

### Fire Scar Fragmentation (lsm_c_split) ----------------------------

split_ts <- fire_metrics_ts %>%
  pull(lsm_c_split) %>%
  ts(start = c(1985, 01), end = c(2024, 12), frequency = 12)

print(split_ts)

x13_split <- x13(split_ts)
print(x13_split)

pdf("Figs/decomp_x13_split.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_split)
dev.off()

pdf("Figs/residuals_x13_split.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(3, 2))
plot(x13_split$regarima)
par(mfrow = c(1, 1))
dev.off()

pdf("Figs/SIratio_x13_split.pdf", paper = "a4r", width = 0, height = 0)
plot(x13_split$decomposition)
dev.off()

# Variable Selection (Boruta) ---------------------------------------------

## Fire occurrence ---------------------------------------------------------

fire_occ_boruta <- Boruta(
  fire_occurrence ~ dewpoint_temp_sd +
    pot_evap_sd +
    precip_min +
    solrad_mean +
    total_evap_max +
    wind_min +
    wind_sd +
    water_soil_mean +
    water_soil_sd +
    savanna +
    forest_plantation +
    grassland +
    pasture +
    sugar_cane +
    mosaic_uses +
    urban +
    soybean +
    non_vegetated +
    forest +
    temporary_crop +
    perennial_crop +
    water +
    wet_herbaceous +
    native_area +
    anthropic_area +
    NDVI +
    pop +
    elevation +
    aspect +
    TPI +
    TRI,
  data = as.data.frame(fire_occ_df),
  doTrace = 1,
  getImp = getImpXgboost
)

print(fire_occ_boruta)
# Boruta performed 12 iterations in 5.765394 mins.
# 29 attributes confirmed important: aspect, dewpoint_temp_sd, elevation, forest, forest_plantation and 24 more;
# 2 attributes confirmed unimportant: anthropic_area, sugar_cane;
attStats(fire_occ_boruta)

pdf("Figs/boruta_fire_occ.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_boruta, las = 2, bty = "n", xlab = "", xaxt = "n")

lz <- lapply(1:ncol(fire_occ_boruta$ImpHistory), function(i) {
  fire_occ_boruta$ImpHistory[is.finite(fire_occ_boruta$ImpHistory[, i]), i]
})
names(lz) <- colnames(fire_occ_boruta$ImpHistory)
Labels <- sort(sapply(lz, median))

axis(
  side = 1,
  las = 2,
  labels = names(Labels),
  at = 1:ncol(fire_occ_boruta$ImpHistory),
  cex.axis = 0.7
)
dev.off()

saveRDS(fire_occ_boruta, "Output/fire_occ_boruta.rds")

## Fire Scar Integrity Variation (lsm_c_cai_sd) ----------------------------
fire_cai_sd_boruta <- Boruta(
  lsm_c_cai_sd ~ dewpoint_temp_sd +
    pot_evap_sd +
    precip_min +
    solrad_mean +
    total_evap_max +
    wind_min +
    wind_sd +
    water_soil_mean +
    water_soil_sd +
    savanna +
    forest_plantation +
    grassland +
    pasture +
    sugar_cane +
    mosaic_uses +
    urban +
    soybean +
    non_vegetated +
    forest +
    temporary_crop +
    perennial_crop +
    water +
    wet_herbaceous +
    native_area +
    anthropic_area +
    NDVI +
    pop +
    elevation +
    aspect +
    TPI +
    TRI,
  data = as.data.frame(fire_metrics_df),
  doTrace = 1,
  getImp = getImpXgboost
)

print(fire_cai_sd_boruta)
# Boruta performed 12 iterations in 1.694677 mins.
# 30 attributes confirmed important: anthropic_area, aspect, dewpoint_temp_sd, elevation, forest and 25 more;
# 1 attributes confirmed unimportant: urban;
attStats(fire_cai_sd_boruta)

pdf("Figs/boruta_fire_cai_sd.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_cai_sd_boruta, las = 2, bty = "n", xlab = "", xaxt = "n")

lz <- lapply(1:ncol(fire_cai_sd_boruta$ImpHistory), function(i) {
  fire_cai_sd_boruta$ImpHistory[
    is.finite(fire_cai_sd_boruta$ImpHistory[, i]),
    i
  ]
})
names(lz) <- colnames(fire_cai_sd_boruta$ImpHistory)
Labels <- sort(sapply(lz, median))

axis(
  side = 1,
  las = 2,
  labels = names(Labels),
  at = 1:ncol(fire_cai_sd_boruta$ImpHistory),
  cex.axis = 0.7
)
dev.off()

saveRDS(fire_cai_sd_boruta, "Output/fire_cai_sd_boruta.rds")

## Fire Scar Shape Regularity (lsm_c_circle_mn) ----------------------------
fire_circle_mn_boruta <- Boruta(
  lsm_c_circle_mn ~ dewpoint_temp_sd +
    pot_evap_sd +
    precip_min +
    solrad_mean +
    total_evap_max +
    wind_min +
    wind_sd +
    water_soil_mean +
    water_soil_sd +
    savanna +
    forest_plantation +
    grassland +
    pasture +
    sugar_cane +
    mosaic_uses +
    urban +
    soybean +
    non_vegetated +
    forest +
    temporary_crop +
    perennial_crop +
    water +
    wet_herbaceous +
    native_area +
    anthropic_area +
    NDVI +
    pop +
    elevation +
    aspect +
    TPI +
    TRI,
  data = as.data.frame(fire_metrics_df),
  doTrace = 1,
  getImp = getImpXgboost
)

print(fire_circle_mn_boruta)
# Boruta performed 16 iterations in 2.093484 mins.
# 30 attributes confirmed important: anthropic_area, aspect, dewpoint_temp_sd, elevation, forest and 25 more;
# 1 attributes confirmed unimportant: sugar_cane;
attStats(fire_circle_mn_boruta)

pdf("Figs/boruta_fire_circle_mn.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_circle_mn_boruta, las = 2, bty = "n", xlab = "", xaxt = "n")

lz <- lapply(1:ncol(fire_circle_mn_boruta$ImpHistory), function(i) {
  fire_circle_mn_boruta$ImpHistory[
    is.finite(fire_circle_mn_boruta$ImpHistory[, i]),
    i
  ]
})
names(lz) <- colnames(fire_circle_mn_boruta$ImpHistory)
Labels <- sort(sapply(lz, median))

axis(
  side = 1,
  las = 2,
  labels = names(Labels),
  at = 1:ncol(fire_circle_mn_boruta$ImpHistory),
  cex.axis = 0.7
)
dev.off()

saveRDS(fire_circle_mn_boruta, "Output/fire_circle_mn_boruta.rds")

## Mean Fire Scar Core Size (lsm_c_core_mn) ----------------------------
fire_core_mn_boruta <- Boruta(
  lsm_c_core_mn ~ dewpoint_temp_sd +
    pot_evap_sd +
    precip_min +
    solrad_mean +
    total_evap_max +
    wind_min +
    wind_sd +
    water_soil_mean +
    water_soil_sd +
    savanna +
    forest_plantation +
    grassland +
    pasture +
    sugar_cane +
    mosaic_uses +
    urban +
    soybean +
    non_vegetated +
    forest +
    temporary_crop +
    perennial_crop +
    water +
    wet_herbaceous +
    native_area +
    anthropic_area +
    NDVI +
    pop +
    elevation +
    aspect +
    TPI +
    TRI,
  data = as.data.frame(fire_metrics_df),
  doTrace = 1,
  maxRuns = 200,
  getImp = getImpXgboost
)

print(fire_core_mn_boruta)
# Boruta performed 199 iterations in 7.881913 mins.
# No attributes deemed important.
# 29 attributes confirmed unimportant: anthropic_area, aspect, dewpoint_temp_sd, elevation, forest and 24 more;
# 2 tentative attributes left: solrad_mean, water_soil_sd;
TentativeRoughFix(fire_core_mn_boruta)
# Boruta performed 199 iterations in 7.881913 mins.
# Tentatives roughfixed over the last 199 iterations.
# 1 attributes confirmed important: solrad_mean;
# 30 attributes confirmed unimportant: anthropic_area, aspect, dewpoint_temp_sd, elevation, forest and 25 more;
attStats(fire_core_mn_boruta)

pdf("Figs/boruta_fire_core_mn.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_core_mn_boruta, las = 2, bty = "n", xlab = "", xaxt = "n")

lz <- lapply(1:ncol(fire_core_mn_boruta$ImpHistory), function(i) {
  fire_core_mn_boruta$ImpHistory[
    is.finite(fire_core_mn_boruta$ImpHistory[, i]),
    i
  ]
})
names(lz) <- colnames(fire_core_mn_boruta$ImpHistory)
Labels <- sort(sapply(lz, median))

axis(
  side = 1,
  las = 2,
  labels = names(Labels),
  at = 1:ncol(fire_core_mn_boruta$ImpHistory),
  cex.axis = 0.7
)
dev.off()

saveRDS(fire_core_mn_boruta, "Output/fire_core_mn_boruta.rds")

## Fire Scar Cohesion Variation (lsm_c_dcore_sd) ----------------------------
fire_dcore_sd_boruta <- Boruta(
  lsm_c_dcore_sd ~ dewpoint_temp_sd +
    pot_evap_sd +
    precip_min +
    solrad_mean +
    total_evap_max +
    wind_min +
    wind_sd +
    water_soil_mean +
    water_soil_sd +
    savanna +
    forest_plantation +
    grassland +
    pasture +
    sugar_cane +
    mosaic_uses +
    urban +
    soybean +
    non_vegetated +
    forest +
    temporary_crop +
    perennial_crop +
    water +
    wet_herbaceous +
    native_area +
    anthropic_area +
    NDVI +
    pop +
    elevation +
    aspect +
    TPI +
    TRI,
  data = as.data.frame(fire_metrics_df),
  doTrace = 1,
  getImp = getImpXgboost
)

print(fire_dcore_sd_boruta)
# Boruta performed 99 iterations in 10.61709 mins.
# 19 attributes confirmed important: anthropic_area, dewpoint_temp_sd, elevation, forest, grassland and 14 more;
# 8 attributes confirmed unimportant: aspect, forest_plantation, non_vegetated, perennial_crop, sugar_cane and 3 more;
# 4 tentative attributes left: mosaic_uses, soybean, TPI, wind_sd;
TentativeRoughFix(fire_dcore_sd_boruta)
# Boruta performed 99 iterations in 10.61709 mins.
# Tentatives roughfixed over the last 99 iterations.
# 22 attributes confirmed important: anthropic_area, dewpoint_temp_sd, elevation, forest, grassland and 17 more;
# 9 attributes confirmed unimportant: aspect, forest_plantation, mosaic_uses, non_vegetated, perennial_crop and 4 more;
attStats(fire_dcore_sd_boruta)

pdf("Figs/boruta_fire_dcore_sd.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_dcore_sd_boruta, las = 2, bty = "n", xlab = "", xaxt = "n")

lz <- lapply(1:ncol(fire_dcore_sd_boruta$ImpHistory), function(i) {
  fire_dcore_sd_boruta$ImpHistory[
    is.finite(fire_dcore_sd_boruta$ImpHistory[, i]),
    i
  ]
})
names(lz) <- colnames(fire_dcore_sd_boruta$ImpHistory)
Labels <- sort(sapply(lz, median))

axis(
  side = 1,
  las = 2,
  labels = names(Labels),
  at = 1:ncol(fire_dcore_sd_boruta$ImpHistory),
  cex.axis = 0.7
)
dev.off()

saveRDS(fire_core_mn_boruta, "Output/fire_dcore_sd_boruta.rds")

## Fire Scar Isolation Variation (lsm_c_enn_sd) ----------------------------
fire_enn_sd_boruta <- Boruta(
  lsm_c_enn_sd ~ dewpoint_temp_sd +
    pot_evap_sd +
    precip_min +
    solrad_mean +
    total_evap_max +
    wind_min +
    wind_sd +
    water_soil_mean +
    water_soil_sd +
    savanna +
    forest_plantation +
    grassland +
    pasture +
    sugar_cane +
    mosaic_uses +
    urban +
    soybean +
    non_vegetated +
    forest +
    temporary_crop +
    perennial_crop +
    water +
    wet_herbaceous +
    native_area +
    anthropic_area +
    NDVI +
    pop +
    elevation +
    aspect +
    TPI +
    TRI,
  data = as.data.frame(fire_metrics_df),
  doTrace = 1,
  getImp = getImpXgboost
)

print(fire_enn_sd_boruta)
# Boruta performed 99 iterations in 12.36234 mins.
# 29 attributes confirmed important: anthropic_area, aspect, dewpoint_temp_sd, elevation, forest and 24 more;
# 1 attributes confirmed unimportant: forest_plantation;
# 1 tentative attributes left: wind_min;
TentativeRoughFix(fire_enn_sd_boruta)
attStats(fire_enn_sd_boruta)

pdf("Figs/boruta_fire_enn_sd.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_enn_sd_boruta, las = 2, bty = "n", xlab = "", xaxt = "n")

lz <- lapply(1:ncol(fire_enn_sd_boruta$ImpHistory), function(i) {
  fire_enn_sd_boruta$ImpHistory[
    is.finite(fire_enn_sd_boruta$ImpHistory[, i]),
    i
  ]
})
names(lz) <- colnames(fire_enn_sd_boruta$ImpHistory)
Labels <- sort(sapply(lz, median))

axis(
  side = 1,
  las = 2,
  labels = names(Labels),
  at = 1:ncol(fire_enn_sd_boruta$ImpHistory),
  cex.axis = 0.7
)
dev.off()

saveRDS(fire_enn_sd_boruta, "Output/fire_enn_sd_boruta.rds")

## Fire Scar Density (lsm_c_pd) ----------------------------
fire_pd_boruta <- Boruta(
  lsm_c_pd ~ dewpoint_temp_sd +
    pot_evap_sd +
    precip_min +
    solrad_mean +
    total_evap_max +
    wind_min +
    wind_sd +
    water_soil_mean +
    water_soil_sd +
    savanna +
    forest_plantation +
    grassland +
    pasture +
    sugar_cane +
    mosaic_uses +
    urban +
    soybean +
    non_vegetated +
    forest +
    temporary_crop +
    perennial_crop +
    water +
    wet_herbaceous +
    native_area +
    anthropic_area +
    NDVI +
    pop +
    elevation +
    aspect +
    TPI +
    TRI,
  data = as.data.frame(fire_metrics_df),
  doTrace = 1,
  getImp = getImpXgboost
)

print(fire_pd_boruta)
# Boruta performed 99 iterations in 12.75684 mins.
# 28 attributes confirmed important: anthropic_area, aspect, dewpoint_temp_sd, elevation, forest and 23 more;
# 1 attributes confirmed unimportant: urban;
# 2 tentative attributes left: soybean, TPI;
TentativeRoughFix(fire_pd_boruta)
attStats(fire_pd_boruta)

pdf("Figs/boruta_fire_pd.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_pd_boruta, las = 2, bty = "n", xlab = "", xaxt = "n")

lz <- lapply(1:ncol(fire_pd_boruta$ImpHistory), function(i) {
  fire_pd_boruta$ImpHistory[is.finite(fire_pd_boruta$ImpHistory[, i]), i]
})
names(lz) <- colnames(fire_pd_boruta$ImpHistory)
Labels <- sort(sapply(lz, median))

axis(
  side = 1,
  las = 2,
  labels = names(Labels),
  at = 1:ncol(fire_pd_boruta$ImpHistory),
  cex.axis = 0.7
)
dev.off()

saveRDS(fire_pd_boruta, "Output/fire_pd_boruta.rds")

## Fire Scar Fragmentation (lsm_c_split) ----------------------------
fire_split_boruta <- Boruta(
  lsm_c_split ~ dewpoint_temp_sd +
    pot_evap_sd +
    precip_min +
    solrad_mean +
    total_evap_max +
    wind_min +
    wind_sd +
    water_soil_mean +
    water_soil_sd +
    savanna +
    forest_plantation +
    grassland +
    pasture +
    sugar_cane +
    mosaic_uses +
    urban +
    soybean +
    non_vegetated +
    forest +
    temporary_crop +
    perennial_crop +
    water +
    wet_herbaceous +
    native_area +
    anthropic_area +
    NDVI +
    pop +
    elevation +
    aspect +
    TPI +
    TRI,
  data = as.data.frame(fire_metrics_df),
  doTrace = 1,
  getImp = getImpXgboost
)

print(fire_split_boruta)
attStats(fire_split_boruta)

pdf("Figs/boruta_fire_split.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_split_boruta, las = 2, bty = "n", xlab = "", xaxt = "n")

lz <- lapply(1:ncol(fire_split_boruta$ImpHistory), function(i) {
  fire_split_boruta$ImpHistory[is.finite(fire_split_boruta$ImpHistory[, i]), i]
})
names(lz) <- colnames(fire_split_boruta$ImpHistory)
Labels <- sort(sapply(lz, median))

axis(
  side = 1,
  las = 2,
  labels = names(Labels),
  at = 1:ncol(fire_split_boruta$ImpHistory),
  cex.axis = 0.7
)
dev.off()

saveRDS(fire_split_boruta, "Output/fire_split_boruta.rds")


# GAM ---------------------------------------------------------------------

## Fire occurrence ---------------------------------------------------------

# 1. Get the starting year and month (requires collecting a small part)
start_info <- fire_occ_df %>%
  summarise(start_year = min(year), .by = NULL) %>% # Need .by=NULL or similar in older dplyr?
  as.data.table() # Collect the minimum year
start_year <- start_info$start_year
start_month <- fire_occ_df %>%
  filter(year == start_year) %>%
  summarise(start_month = min(month), .by = NULL) %>%
  as.data.table() # Collect the minimum month for the start year
start_month <- start_month$start_month

# 2. Add the 't' column using mutate
fire_occ_df <- fire_occ_df %>%
  mutate(t = (year - start_year) * 12 + (month - start_month) + 1)

# 1. Create the lazy table
fire_lazy <- lazy_dt(fire_occ_df)

# Define your saturation point (e.g., 4 years)
saturation_point <- 48

# 2. Perform the calculation pipeline
fire_dt_with_tslf <- fire_lazy %>%
  group_by(plot_id) %>%
  arrange(t) %>%
  # --- THE FIX ---
  # We define the interval based on the LAG of fire occurrence.
  # This means if a fire happens at t=5, the new group starts at t=6.
  mutate(
    lag_fire = lag(fire_occurrence, default = 0),
    fire_interval_id = cumsum(lag_fire)
  ) %>%
  group_by(plot_id, fire_interval_id) %>%
  mutate(
    # Now, if a fire happens at month 36, this value will be 36.
    # It resets to 1 in the next month.
    months_since_fire = row_number()
  ) %>%
  # --- Saturation Logic (Same as before) ---
  ungroup() %>%
  mutate(
    # Handle pre-history (assume saturated if we haven't seen a previous fire)
    months_since_fire = if_else(
      fire_interval_id == 0,
      as.numeric(saturation_point),
      as.numeric(months_since_fire)
    ),
    # Cap the maximum
    months_since_fire = if_else(
      months_since_fire > saturation_point,
      as.numeric(saturation_point),
      months_since_fire
    )
  ) %>%
  dplyr::select(-lag_fire, -fire_interval_id) %>% # Clean up helper columns
  as_tibble()

# Sanity Check (Very Important!)
# Look at a single plot that had a fire.
# You should see 'months' increasing up to the fire event.
View(
  fire_dt_with_tslf %>%
    filter(plot_id == 1528) %>%
    dplyr::select(t, fire_occurrence, months_since_fire)
)
# Check the result
glimpse(fire_dt_with_tslf)

# To see the actual result, collect it
fire_occ_df <- fire_dt_with_tslf %>% as_tibble()
head(fire_occ_df)
print(tail(fire_occ_df))

# Null model
fire_occ_bam_null <- bam(
  fire_occurrence ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15),
  data = fire_occ_df,
  family = "binomial",
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)

summary(fire_occ_bam_null)

check(getViz(fire_occ_bam_null))

viz_null <- getViz(fire_occ_bam_null)

# Plot the Fuel Accumulation Curve
pdf(
  "Figs/fire_occ_bam_null_months_since_fire.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)

# trans = plogis converts Log-Odds -> Probability
plot(sm(viz_null, 3)) + # 3 is the index of s(months_since_fire)
  l_fitLine() +
  l_ciPoly() +
  theme_classic() +
  labs(
    title = "Fuel Limitation Effect",
    x = "Months Since Last Fire",
    y = "Effect on Fire Probability"
  )
dev.off()

# Full model
print(fire_occ_boruta$finalDecision)

fire_occ_bam_full <- bam(
  fire_occurrence ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +

    # Environmental drivers
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +

    # Land use
    s(savanna, bs = "cr") +
    s(forest_plantation, bs = "cr") +
    s(grassland, bs = "cr") +
    s(pasture, bs = "cr") +
    s(mosaic_uses, bs = "cr") +
    s(urban, bs = "cr") +
    s(soybean, bs = "cr") +
    s(non_vegetated, bs = "cr") +
    s(forest, bs = "cr") +
    s(temporary_crop, bs = "cr") +
    s(perennial_crop, bs = "cr") +
    s(water, bs = "cr") +
    s(wet_herbaceous, bs = "cr") +
    s(pop, bs = "cr") +

    # Biophysical
    s(NDVI, bs = "cr") +
    s(elevation, bs = "cr") +
    s(aspect, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_occ_df,
  family = "binomial",
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_occ_bam_full)
saveRDS(fire_occ_bam_full, "Output/fire_occ_bam_full.rds")

# Simpler model (lulc only native and anthropic areas)

fire_occ_bam_simple <- bam(
  fire_occurrence ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(native_area, bs = "cr") +
    s(anthropic_area, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(aspect, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_occ_df,
  family = "binomial",
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_occ_bam_simple)

# Compare models
anova(fire_occ_bam_null, fire_occ_bam_full, test = "Chisq")

AIC(fire_occ_bam_null, fire_occ_bam_full, fire_occ_bam_simple)
bbmle::AICctab(fire_occ_bam_null, fire_occ_bam_full, fire_occ_bam_simple)


# Plots
# Full model
fire_occ_bam_full <- readRDS("Output/fire_occ_bam_full.rds")

fire_occ_bam_full_viz <- getViz(fire_occ_bam_full)

pdf("Figs/fire_occ_bam_full_all.pdf", paper = "a4r", width = 0, height = 0)
print(
  plot(fire_occ_bam_full_viz, allTerms = TRUE) +
    l_fitLine() +
    l_ciLine() +
    labs(y = "") +
    # Rotate x-axis text 45 degrees and align it to the right
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  pages = 1
)
dev.off()

pl <- plot(fire_occ_bam_full_viz, allTerms = T) +
  theme_get() +
  labs(title = NULL)

pdf(
  "Figs/fire_occ_bam_full_all_contours.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
print(pl, pages = 1)
dev.off()

pdf("Figs/fire_occ_bam_full_trend.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 2, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_TSLF.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 3, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_dewpoint_temp_sd.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 4, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_pot_evap_sd.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 5) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_precip_min.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 6) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_solrad_mean.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 7) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_total_evap_max.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 8) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_wind_min.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 9) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_wind_sd.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 10) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_water_soil_mean.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 11) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_water_soil_sd.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 12) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_savanna.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 13) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()


pdf(
  "Figs/fire_occ_bam_full_forest_plantation.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 14) + l_fitLine() + l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_grassland.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 15) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_pasture.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 16) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_mosaic_uses.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 17) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_urban.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 18) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_soybean.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 19) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_non_vegetated.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 20) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_forest.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 21) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_temporary_crop.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 22) +
  l_dens(type = "cond") +
  l_rug() +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_perennial_crop.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 23) + l_rug() + l_fitLine() + l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_water.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 24) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_occ_bam_full_wet_herbaceous.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 25) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_pop.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 26) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_NDVI.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 27) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()


pdf(
  "Figs/fire_occ_bam_full_elevation.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_occ_bam_full_viz, select = 28) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_aspect.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 29) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_TPI.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 30) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_occ_bam_full_TRI.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_occ_bam_full_viz, select = 31) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()


# January
pdf("Figs/fire_occ_bam_full_Jan.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# February
pdf("Figs/fire_occ_bam_full_Feb.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# March
pdf("Figs/fire_occ_bam_full_Mar.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# April
pdf("Figs/fire_occ_bam_full_Apr.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# May
pdf("Figs/fire_occ_bam_full_May.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# June
pdf("Figs/fire_occ_bam_full_Jun.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# July
pdf("Figs/fire_occ_bam_full_Jul.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# August
pdf("Figs/fire_occ_bam_full_Aug.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# September
pdf("Figs/fire_occ_bam_full_Sep.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# October
pdf("Figs/fire_occ_bam_full_Oct.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# November
pdf("Figs/fire_occ_bam_full_Nov.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# December
pdf("Figs/fire_occ_bam_full_Dec.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour()

dev.off()

pdf("Figs/fire_occ_bam_full_check.pdf", paper = "a4r", width = 0, height = 0)
check(fire_occ_bam_full_viz)
dev.off()

## Fire Scar Integrity Variation (lsm_c_cai_sd) ----------------------------

# 1. Get the starting year and month (requires collecting a small part)
start_info <- fire_metrics_df %>%
  summarise(start_year = min(year), .by = NULL) %>% # Need .by=NULL or similar in older dplyr?
  as.data.table() # Collect the minimum year
start_year <- start_info$start_year
start_month <- fire_metrics_df %>%
  filter(year == start_year) %>%
  summarise(start_month = min(month), .by = NULL) %>%
  as.data.table() # Collect the minimum month for the start year
start_month <- start_month$start_month

# 2. Add the 't' column using mutate
fire_metrics_df <- fire_metrics_df %>%
  mutate(t = (year - start_year) * 12 + (month - start_month) + 1)

# Select only keys and the variable we want to transfer to save memory
tslf_lookup <- fire_occ_df %>%
  dplyr::select(plot_id, t, months_since_fire)

fire_metrics_df <- fire_metrics_df %>%
  inner_join(tslf_lookup, by = c("plot_id", "t"))

# Check the (lazy) result structure
print(fire_metrics_df)

# To see the actual result, collect it
fire_metrics_df <- fire_metrics_df %>% as_tibble()
head(fire_metrics_df)
print(tail(fire_metrics_df))

# Check distribution

pdf("Figs/histogram_lsm_c_cai_sd.pdf", paper = "a4", width = 0, height = 0)
hist(fire_metrics_df$lsm_c_cai_sd)
dev.off()

fire_cai_sd_bam_tweedie <- bam(
  lsm_c_cai_sd ~ s(t),
  data = fire_metrics_df,
  family = Tweedie(p = 1.5, link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_cai_sd_bam_tweedie)


fire_cai_sd_bam_scat <- bam(
  lsm_c_cai_sd ~ s(t),
  data = fire_metrics_df,
  family = scat(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_cai_sd_bam_scat)

fire_cai_sd_bam_gaussian <- bam(
  lsm_c_cai_sd ~ s(t),
  data = fire_metrics_df,
  family = gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_cai_sd_bam_gaussian)

fire_cai_sd_bam_loggaussian <- bam(
  log1p(lsm_c_cai_sd) ~ s(t),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_cai_sd_bam_loggaussian)

AIC(
  fire_cai_sd_bam_tweedie,
  fire_cai_sd_bam_scat,
  fire_cai_sd_bam_gaussian,
  fire_cai_sd_bam_loggaussian
)

bbmle::AICtab(
  fire_cai_sd_bam_tweedie,
  fire_cai_sd_bam_scat,
  fire_cai_sd_bam_gaussian,
  fire_cai_sd_bam_loggaussian
)


pdf(
  "Figs/gamcheck_loggaussian_lsm_c_cai_sd.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
par(mfrow = c(2, 2))
check(getViz(fire_cai_sd_bam_loggaussian))
dev.off()

# Null model
fire_cai_sd_bam_null <- bam(
  log1p(lsm_c_cai_sd) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)


summary(fire_cai_sd_bam_null)

# Full model
sort(TentativeRoughFix(fire_cai_sd_boruta)$finalDecision)

fire_cai_sd_bam_full <- bam(
  log1p(lsm_c_cai_sd) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(savanna, bs = "cr") +
    s(forest_plantation, bs = "cr") +
    s(grassland, bs = "cr") +
    s(pasture, bs = "cr") +
    s(sugar_cane, bs = "cr") +
    s(mosaic_uses, bs = "cr") +
    s(soybean, bs = "cr") +
    s(non_vegetated, bs = "cr") +
    s(forest, bs = "cr") +
    s(temporary_crop, bs = "cr") +
    s(perennial_crop, bs = "cr") +
    s(water, bs = "cr") +
    s(wet_herbaceous, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(aspect, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_cai_sd_bam_full)
saveRDS(fire_cai_sd_bam_full, "Output/fire_cai_sd_bam_full.rds")

# Simpler model (lulc only native and anthropic areas)

fire_cai_sd_bam_simple <- bam(
  log1p(lsm_c_cai_sd) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(native_area, bs = "cr") +
    s(anthropic_area, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(aspect, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_cai_sd_bam_simple)

# Compare models
anova(fire_cai_sd_bam_null, fire_cai_sd_bam_full, test = "Chisq")

AIC(fire_cai_sd_bam_null, fire_cai_sd_bam_full, fire_cai_sd_bam_simple)
bbmle::AICctab(
  fire_cai_sd_bam_null,
  fire_cai_sd_bam_full,
  fire_cai_sd_bam_simple
)

# Plots

# Full model
fire_cai_sd_bam_full <- readRDS("Output/fire_cai_sd_bam_full.rds")
fire_cai_sd_bam_full_viz <- getViz(fire_cai_sd_bam_full)

pdf("Figs/fire_cai_sd_bam_full_all.pdf", paper = "a4r", width = 0, height = 0)
print(
  plot(fire_cai_sd_bam_full_viz, allTerms = TRUE) +
    l_fitLine() +
    l_ciLine() +
    labs(y = "") +
    # Rotate x-axis text 45 degrees and align it to the right
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  pages = 1
)
dev.off()

pl <- plot(fire_cai_sd_bam_full_viz, allTerms = T) +
  theme_get() +
  labs(title = NULL)

pdf(
  "Figs/fire_cai_sd_bam_full_all_contours.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
print(pl, pages = 1)
dev.off()

pdf("Figs/fire_cai_sd_bam_full_trend.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_cai_sd_bam_full_viz, select = 2, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_cai_sd_bam_full_TSLF.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_cai_sd_bam_full_viz, select = 3, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

# pdf("Figs/fire_cai_sd_bam_full_dewpoint_temp_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 4, type = "response") + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_pot_evap_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 5) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_precip_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 6) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_solrad_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 7) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_total_evap_max.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 8) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_wind_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 9) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_wind_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 10) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_water_soil_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 11) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_water_soil_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 12) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_savanna.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 14) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
#
# pdf("Figs/fire_cai_sd_bam_full_forest_plantation.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 15) + l_fitLine() + l_ciLine() # scale estimate is zero for input data
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_grassland.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 16) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_pasture.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 17) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_sugar_cane.pdf", paper = "a4r", width = 0, height = 0) #scale estimate is zero for input data
# plot(fire_cai_sd_bam_full_viz, select = 18) + l_rug() + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_soybean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 19) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_non_vegetated.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 20) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_forest.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 21) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_NDVI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 22) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_pop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 23) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_elevation.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 24) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_TPI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 25) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_cai_sd_bam_full_TRI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_cai_sd_bam_full_viz, select = 26) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()

# January
pdf("Figs/fire_cai_sd_bam_full_Jan.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# February
pdf("Figs/fire_cai_sd_bam_full_Feb.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# March
pdf("Figs/fire_cai_sd_bam_full_Mar.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# April
pdf("Figs/fire_cai_sd_bam_full_Apr.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# May
pdf("Figs/fire_cai_sd_bam_full_May.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# June
pdf("Figs/fire_cai_sd_bam_full_Jun.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# July
pdf("Figs/fire_cai_sd_bam_full_Jul.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# August
pdf("Figs/fire_cai_sd_bam_full_Aug.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# September
pdf("Figs/fire_cai_sd_bam_full_Sep.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# October
pdf("Figs/fire_cai_sd_bam_full_Oct.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# November
pdf("Figs/fire_cai_sd_bam_full_Nov.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# December
pdf("Figs/fire_cai_sd_bam_full_Dec.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# Diagnostics
pdf("Figs/fire_cai_sd_bam_full_check.pdf", paper = "a4r", width = 0, height = 0)
check(fire_cai_sd_bam_full_viz)
dev.off()

## Fire Scar Shape Regularity (lsm_c_circle_mn) ----------------------------

head(fire_metrics_df)
print(tail(fire_metrics_df))
skim(fire_metrics_df$lsm_c_circle_mn)

# Check distribution

pdf("Figs/histogram_lsm_c_circle_mn.pdf", paper = "a4", width = 0, height = 0)
hist(fire_metrics_df$lsm_c_circle_mn)
dev.off()

fire_circle_mn_bam_tweedie <- bam(
  lsm_c_circle_mn ~ s(t),
  data = fire_metrics_df,
  family = Tweedie(p = 1.5, link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_circle_mn_bam_tweedie)


fire_circle_mn_bam_scat <- bam(
  lsm_c_circle_mn ~ s(t),
  data = fire_metrics_df,
  family = scat(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_circle_mn_bam_scat)

fire_circle_mn_bam_gamma <- bam(
  lsm_c_circle_mn ~ s(t),
  data = fire_metrics_df,
  family = Gamma(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_circle_mn_bam_gamma)

fire_circle_mn_bam_gaussian <- bam(
  lsm_c_circle_mn ~ s(t),
  data = fire_metrics_df,
  family = gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_circle_mn_bam_gaussian)

fire_circle_mn_bam_invgaussian <- bam(
  lsm_c_circle_mn ~ s(t),
  data = fire_metrics_df,
  family = inverse.gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_circle_mn_bam_invgaussian)

fire_circle_mn_bam_loggaussian <- bam(
  log1p(lsm_c_circle_mn) ~ s(t),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_circle_mn_bam_loggaussian)


AIC(
  fire_circle_mn_bam_tweedie,
  fire_circle_mn_bam_scat,
  fire_circle_mn_bam_gamma,
  fire_circle_mn_bam_gaussian,
  fire_circle_mn_bam_invgaussian,
  fire_circle_mn_bam_loggaussian
)

bbmle::AICctab(
  fire_circle_mn_bam_tweedie,
  fire_circle_mn_bam_scat,
  fire_circle_mn_bam_gamma,
  fire_circle_mn_bam_gaussian,
  fire_circle_mn_bam_invgaussian,
  fire_circle_mn_bam_loggaussian
)


pdf(
  "Figs/gamcheck_loggaussian_lsm_c_circle_mn.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
par(mfrow = c(2, 2))
check(getViz(fire_circle_mn_bam_loggaussian))
dev.off()

# Null model
fire_circle_mn_bam_null <- bam(
  log1p(lsm_c_circle_mn) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)

summary(fire_circle_mn_bam_null)

# Full model
sort(fire_circle_mn_boruta$finalDecision)

fire_circle_mn_bam_full <- bam(
  log1p(lsm_c_circle_mn) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(savanna, bs = "cr") +
    s(forest_plantation, bs = "cr") +
    s(grassland, bs = "cr") +
    s(pasture, bs = "cr") +
    s(mosaic_uses, bs = "cr") +
    s(urban, bs = "cr") +
    s(soybean, bs = "cr") +
    s(non_vegetated, bs = "cr") +
    s(forest, bs = "cr") +
    s(temporary_crop, bs = "cr") +
    s(perennial_crop, bs = "cr") +
    s(water, bs = "cr") +
    s(wet_herbaceous, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(aspect, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_circle_mn_bam_full)
saveRDS(fire_circle_mn_bam_full, "Output/fire_circle_mn_bam_full.rds")

# Simpler model (lulc only native and anthropic areas)

fire_circle_mn_bam_simple <- bam(
  log1p(lsm_c_circle_mn) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(native_area, bs = "cr") +
    s(anthropic_area, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(aspect, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_circle_mn_bam_simple)

# Compare models
anova(fire_circle_mn_bam_null, fire_circle_mn_bam_full, test = "Chisq")

AIC(fire_circle_mn_bam_null, fire_circle_mn_bam_full, fire_circle_mn_bam_simple)
bbmle::AICctab(
  fire_circle_mn_bam_null,
  fire_circle_mn_bam_full,
  fire_circle_mn_bam_simple
)

# Plots

# Full model
fire_circle_mn_bam_full <- readRDS("Output/fire_circle_mn_bam_full.rds")
fire_circle_mn_bam_full_viz <- getViz(fire_circle_mn_bam_full)

pdf(
  "Figs/fire_circle_mn_bam_full_all.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
print(
  plot(fire_circle_mn_bam_full_viz, allTerms = TRUE) +
    l_fitLine() +
    l_ciLine() +
    labs(y = "") +
    # Rotate x-axis text 45 degrees and align it to the right
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  pages = 1
)
dev.off()

pl <- plot(fire_circle_mn_bam_full_viz, allTerms = T) +
  theme_get() +
  labs(title = NULL)

pdf(
  "Figs/fire_circle_mn_bam_full_all_contours.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
print(pl, pages = 1)
dev.off()

pdf(
  "Figs/fire_circle_mn_bam_full_trend.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_circle_mn_bam_full_viz, select = 2, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_circle_mn_bam_full_TSLF.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_circle_mn_bam_full_viz, select = 3, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

# pdf("Figs/fire_circle_mn_bam_full_dewpoint_temp_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 4, type = "response") + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_temp_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 5) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_pot_evap_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 6) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_precip_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 7) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_solrad_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 8) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_total_evap_max.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 9) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_wind_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 10) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_wind_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 11) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_water_soil_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 12) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_water_soil_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 13) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_savanna.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 14) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_forest_plantation.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 15) + l_fitLine() + l_ciLine() # scale estimate is zero for input data
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_grassland.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 16) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_pasture.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 17) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_sugar_cane.pdf", paper = "a4r", width = 0, height = 0) #scale estimate is zero for input data
# plot(fire_circle_mn_bam_full_viz, select = 18) + l_rug() + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_mosaic_uses.pdf", paper = "a4r", width = 0, height = 0) #scale estimate is zero for input data
# plot(fire_circle_mn_bam_full_viz, select = 19) + l_dens(type = "cond") + l_rug() + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_soybean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 20) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_forest.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 21) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_temporary_crop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 22) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_perennial_crop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 23) + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_wet_herbaceous.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 24) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_NDVI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 25) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_pop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 26) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_elevation.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 27) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_aspect.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 28) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_TPI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 29) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_circle_mn_bam_full_TRI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_circle_mn_bam_full_viz, select = 30) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()

# January
pdf(
  "Figs/fire_circle_mn_bam_full_Jan.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# February
pdf(
  "Figs/fire_circle_mn_bam_full_Feb.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# March
pdf(
  "Figs/fire_circle_mn_bam_full_Mar.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# April
pdf(
  "Figs/fire_circle_mn_bam_full_Apr.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# May
pdf(
  "Figs/fire_circle_mn_bam_full_May.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# June
pdf(
  "Figs/fire_circle_mn_bam_full_Jun.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# July
pdf(
  "Figs/fire_circle_mn_bam_full_Jul.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# August
pdf(
  "Figs/fire_circle_mn_bam_full_Aug.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# September
pdf(
  "Figs/fire_circle_mn_bam_full_Sep.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# October
pdf(
  "Figs/fire_circle_mn_bam_full_Oct.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# November
pdf(
  "Figs/fire_circle_mn_bam_full_Nov.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# December
pdf(
  "Figs/fire_circle_mn_bam_full_Dec.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# Diagnostics
pdf(
  "Figs/fire_circle_mn_bam_full_check.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
check(fire_circle_mn_bam_full_viz)
dev.off()

## Mean Fire Scar Core Size (lsm_c_core_mn) ----------------------------

head(fire_metrics_df)
print(tail(fire_metrics_df))
skim(fire_metrics_df$lsm_c_core_mn)

# Check distribution

pdf("Figs/histogram_lsm_c_core_mn.pdf", paper = "a4", width = 0, height = 0)
hist(log(fire_metrics_df$lsm_c_core_mn + 1))
dev.off()

fire_core_mn_bam_tweedie <- bam(
  lsm_c_core_mn ~ s(t),
  data = fire_metrics_df,
  family = Tweedie(p = 1.5, link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_core_mn_bam_tweedie)


fire_core_mn_bam_scat <- bam(
  lsm_c_core_mn ~ s(t),
  data = fire_metrics_df,
  family = scat(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_core_mn_bam_scat)


fire_core_mn_bam_gaussian <- bam(
  lsm_c_core_mn ~ s(t),
  data = fire_metrics_df,
  family = gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_core_mn_bam_gaussian)

fire_core_mn_bam_loggaussian <- bam(
  log1p(lsm_c_core_mn) ~ s(t),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_core_mn_bam_loggaussian)

AIC(
  fire_core_mn_bam_tweedie,
  fire_core_mn_bam_scat,
  fire_core_mn_bam_gaussian,
  fire_core_mn_bam_loggaussian
)

bbmle::AICctab(
  fire_core_mn_bam_tweedie,
  fire_core_mn_bam_scat,
  fire_core_mn_bam_gaussian,
  fire_core_mn_bam_loggaussian
)


pdf(
  "Figs/gamcheck_loggaussian_lsm_c_core_mn.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
par(mfrow = c(2, 2))
check(getViz(fire_core_mn_bam_loggaussian))
dev.off()

# Null model
fire_core_mn_bam_null <- bam(
  log1p(lsm_c_core_mn) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)

summary(fire_core_mn_bam_null)

# Full model
sort(fire_core_mn_boruta$finalDecision)
sort(TentativeRoughFix(fire_core_mn_boruta)$finalDecision)


fire_core_mn_bam_full <- bam(
  log1p(lsm_c_core_mn) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(solrad_mean, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_core_mn_bam_full)
saveRDS(fire_core_mn_bam_full, "Output/fire_core_mn_bam_full.rds")

# Compare models
anova(fire_core_mn_bam_null, fire_core_mn_bam_full, test = "Chisq")

AIC(fire_core_mn_bam_null, fire_core_mn_bam_full)
bbmle::AICctab(fire_core_mn_bam_null, fire_core_mn_bam_full)


# Plots

# Full model
fire_core_mn_bam_full <- readRDS("Output/fire_core_mn_bam_full.rds")
fire_core_mn_bam_full_viz <- getViz(fire_core_mn_bam_full)

pdf("Figs/fire_core_mn_bam_full_all.pdf", paper = "a4r", width = 0, height = 0)
print(
  plot(fire_core_mn_bam_full_viz, allTerms = TRUE) +
    l_fitLine() +
    l_ciLine() +
    labs(y = "") +
    # Rotate x-axis text 45 degrees and align it to the right
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  pages = 1
)
dev.off()

pl <- plot(fire_core_mn_bam_full_viz, allTerms = T) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine() +
  theme_get() +
  labs(title = NULL)

pdf(
  "Figs/fire_core_mn_bam_full_all_contours.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
print(pl, pages = 1)
dev.off()

pdf(
  "Figs/fire_core_mn_bam_full_trend.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_core_mn_bam_full_viz, select = 2) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_core_mn_bam_full_TSLF.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_core_mn_bam_full_viz, select = 3) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_core_mn_bam_full_solrad_mean.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_core_mn_bam_full_viz, select = 4) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

# January
pdf("Figs/fire_core_mn_bam_full_Jan.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# February
pdf("Figs/fire_core_mn_bam_full_Feb.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# March
pdf("Figs/fire_core_mn_bam_full_Mar.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# April
pdf("Figs/fire_core_mn_bam_full_Apr.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# May
pdf("Figs/fire_core_mn_bam_full_May.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# June
pdf("Figs/fire_core_mn_bam_full_Jun.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# July
pdf("Figs/fire_core_mn_bam_full_Jul.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# August
pdf("Figs/fire_core_mn_bam_full_Aug.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# September
pdf("Figs/fire_core_mn_bam_full_Sep.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# October
pdf("Figs/fire_core_mn_bam_full_Oct.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# November
pdf("Figs/fire_core_mn_bam_full_Nov.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# December
pdf("Figs/fire_core_mn_bam_full_Dec.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# Diagnostics
pdf(
  "Figs/fire_core_mn_bam_full_check.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
check(fire_core_mn_bam_full_viz)
dev.off()

## Fire Scar Cohesion Variation (lsm_c_dcore_sd) ----------------------------

head(fire_metrics_df)
print(tail(fire_metrics_df))
skim(fire_metrics_df$lsm_c_dcore_sd)

# Check distribution

pdf("Figs/histogram_lsm_c_dcore_sd.pdf", paper = "a4", width = 0, height = 0)
hist(log(fire_metrics_df$lsm_c_dcore_sd + 1))
dev.off()

fire_dcore_sd_bam_tweedie <- bam(
  lsm_c_dcore_sd ~ s(t),
  data = fire_metrics_df,
  family = Tweedie(p = 1.5, link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_dcore_sd_bam_tweedie)


fire_dcore_sd_bam_scat <- bam(
  lsm_c_dcore_sd ~ s(t),
  data = fire_metrics_df,
  family = scat(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_dcore_sd_bam_scat)


fire_dcore_sd_bam_gaussian <- bam(
  lsm_c_dcore_sd ~ s(t),
  data = fire_metrics_df,
  family = gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_dcore_sd_bam_gaussian)

fire_dcore_sd_bam_loggaussian <- bam(
  log1p(lsm_c_dcore_sd) ~ s(t),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_dcore_sd_bam_loggaussian)


AIC(
  fire_dcore_sd_bam_tweedie,
  fire_dcore_sd_bam_scat,
  fire_dcore_sd_bam_gaussian,
  fire_dcore_sd_bam_loggaussian
)

bbmle::AICctab(
  fire_dcore_sd_bam_tweedie,
  fire_dcore_sd_bam_scat,
  fire_dcore_sd_bam_gaussian,
  fire_dcore_sd_bam_loggaussian
)


pdf(
  "Figs/gamcheck_loggaussian_lsm_c_dcore_sd.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
par(mfrow = c(2, 2))
check(getViz(fire_dcore_sd_bam_loggaussian))
dev.off()

# Null model
fire_dcore_sd_bam_null <- bam(
  log1p(lsm_c_dcore_sd) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)

summary(fire_dcore_sd_bam_null)

# Full model
sort(TentativeRoughFix(fire_dcore_sd_boruta)$finalDecision)

fire_dcore_sd_bam_full <- bam(
  log1p(lsm_c_dcore_sd) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(savanna, bs = "cr") +
    s(grassland, bs = "cr") +
    s(pasture, bs = "cr") +
    s(soybean, bs = "cr") +
    s(forest, bs = "cr") +
    s(temporary_crop, bs = "cr") +
    s(wet_herbaceous, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_dcore_sd_bam_full)
saveRDS(fire_dcore_sd_bam_full, "Output/fire_dcore_sd_bam_full.rds")

# Simpler model (lulc only native and anthropic areas)
fire_dcore_sd_bam_simple <- bam(
  log1p(lsm_c_dcore_sd) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(native_area, bs = "cr") +
    s(anthropic_area, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_dcore_sd_bam_simple)

# Compare models
anova(fire_dcore_sd_bam_null, fire_dcore_sd_bam_full, test = "Chisq")

AIC(fire_dcore_sd_bam_null, fire_dcore_sd_bam_full, fire_dcore_sd_bam_simple)
bbmle::AICctab(
  fire_dcore_sd_bam_null,
  fire_dcore_sd_bam_full,
  fire_dcore_sd_bam_simple
)

# Plots

# Full model
fire_dcore_sd_bam_full <- readRDS("Output/fire_dcore_sd_bam_full.rds")

fire_dcore_sd_bam_full_viz <- getViz(fire_dcore_sd_bam_full)

pdf("Figs/fire_dcore_sd_bam_full_all.pdf", paper = "a4r", width = 0, height = 0)
print(
  plot(fire_dcore_sd_bam_full_viz, allTerms = TRUE) +
    l_fitLine() +
    l_ciLine() +
    labs(y = "") +
    # Rotate x-axis text 45 degrees and align it to the right
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  pages = 1
)
dev.off()

pl <- plot(fire_dcore_sd_bam_full_viz, allTerms = T) +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine() +
  theme_get() +
  labs(title = NULL)

pdf(
  "Figs/fire_dcore_sd_bam_full_all_contours.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
print(pl, pages = 1)
dev.off()

pdf(
  "Figs/fire_dcore_sd_bam_full_trend.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_dcore_sd_bam_full_viz, select = 2, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf(
  "Figs/fire_dcore_sd_bam_full_TSLF.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
plot(fire_dcore_sd_bam_full_viz, select = 3, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

# pdf("Figs/fire_dcore_sd_bam_full_dewpoint_temp_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 4, type = "response") + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_temp_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 5) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_pot_evap_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 6) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_precip_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 7) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_solrad_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 8) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_total_evap_max.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 9) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_wind_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 10) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_water_soil_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 11) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_water_soil_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 12) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_savanna.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 13) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
#
# pdf("Figs/fire_dcore_sd_bam_full_mosaic_uses.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 14) + l_fitLine() + l_ciLine() # scale estimate is zero for input data
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_forest.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 15) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_temporary_crop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 16) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_NDVI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 17) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_pop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 18) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_elevation.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 19) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_dcore_sd_bam_full_TRI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_dcore_sd_bam_full_viz, select = 20) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# table(fire_metrics_df[fire_metrics_df$month==2, "t"])

# January
pdf("Figs/fire_dcore_sd_bam_full_Jan.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# February
pdf("Figs/fire_dcore_sd_bam_full_Feb.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# March
pdf("Figs/fire_dcore_sd_bam_full_Mar.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# April
pdf("Figs/fire_dcore_sd_bam_full_Apr.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# May
pdf("Figs/fire_dcore_sd_bam_full_May.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# June
pdf("Figs/fire_dcore_sd_bam_full_Jun.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# July
pdf("Figs/fire_dcore_sd_bam_full_Jul.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# August
pdf("Figs/fire_dcore_sd_bam_full_Aug.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# September
pdf("Figs/fire_dcore_sd_bam_full_Sep.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# October
pdf("Figs/fire_dcore_sd_bam_full_Oct.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# November
pdf("Figs/fire_dcore_sd_bam_full_Nov.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# December
pdf("Figs/fire_dcore_sd_bam_full_Dec.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# Diagnostics
pdf(
  "Figs/fire_dcore_sd_bam_full_check.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
check(fire_dcore_sd_bam_full_viz)
dev.off()

## Fire Scar Isolation Variation (lsm_c_enn_sd) ----------------------------

head(fire_metrics_df)
print(tail(fire_metrics_df))
skim(fire_metrics_df$lsm_c_enn_sd)

# Check distribution

pdf("Figs/histogram_lsm_c_enn_sd.pdf", paper = "a4", width = 0, height = 0)
hist(fire_metrics_df$lsm_c_enn_sd)
dev.off()

fire_enn_sd_bam_tweedie <- bam(
  lsm_c_enn_sd ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = Tweedie(p = 1.5, link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_enn_sd_bam_tweedie)


fire_enn_sd_bam_scat <- bam(
  lsm_c_enn_sd ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = scat(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_enn_sd_bam_scat)


fire_enn_sd_bam_gaussian <- bam(
  lsm_c_enn_sd ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_enn_sd_bam_gaussian)

fire_enn_sd_bam_loggaussian <- bam(
  log1p(lsm_c_enn_sd) ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_enn_sd_bam_loggaussian)


AIC(
  fire_enn_sd_bam_tweedie,
  fire_enn_sd_bam_scat,
  fire_enn_sd_bam_gaussian,
  fire_enn_sd_bam_loggaussian
)

bbmle::AICctab(
  fire_enn_sd_bam_tweedie,
  fire_enn_sd_bam_scat,
  fire_enn_sd_bam_gaussian,
  fire_enn_sd_bam_loggaussian
)


pdf(
  "Figs/gamcheck_loggaussian_lsm_c_enn_sd.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
par(mfrow = c(2, 2))
check(getViz(fire_enn_sd_bam_loggaussian))
dev.off()

# Null model
fire_enn_sd_bam_null <- bam(
  log1p(lsm_c_enn_sd) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)

summary(fire_enn_sd_bam_null)

# Full model
sort(TentativeRoughFix(fire_enn_sd_boruta)$finalDecision)

fire_enn_sd_bam_full <- bam(
  log1p(lsm_c_enn_sd) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(savanna, bs = "cr") +
    s(grassland, bs = "cr") +
    s(pasture, bs = "cr") +
    s(sugar_cane, bs = "cr") +
    s(mosaic_uses, bs = "cr") +
    s(urban, bs = "cr") +
    s(soybean, bs = "cr") +
    s(non_vegetated, bs = "cr") +
    s(forest, bs = "cr") +
    s(temporary_crop, bs = "cr") +
    s(perennial_crop, bs = "cr") +
    s(water, bs = "cr") +
    s(wet_herbaceous, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(aspect, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_enn_sd_bam_full)
saveRDS(fire_enn_sd_bam_full, "Output/fire_enn_sd_bam_full.rds")

# Simpler model (lulc only native and anthropic areas)

fire_enn_sd_bam_simple <- bam(
  log1p(lsm_c_enn_sd) ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(native_area, bs = "cr") +
    s(anthropic_area, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(aspect, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_enn_sd_bam_simple)

# Compare models
anova(fire_enn_sd_bam_null, fire_enn_sd_bam_full, test = "Chisq")

AIC(fire_enn_sd_bam_null, fire_enn_sd_bam_full, fire_enn_sd_bam_simple)
bbmle::AICctab(
  fire_enn_sd_bam_null,
  fire_enn_sd_bam_full,
  fire_enn_sd_bam_simple
)

# Plots

# Full model
fire_enn_sd_bam_full <- readRDS("Output/fire_enn_sd_bam_full.rds")
fire_enn_sd_bam_full_viz <- getViz(fire_enn_sd_bam_full)

pdf("Figs/fire_enn_sd_bam_full_all.pdf", paper = "a4r", width = 0, height = 0)
print(
  plot(fire_enn_sd_bam_full_viz, allTerms = TRUE) +
    l_fitLine() +
    l_ciLine() +
    labs(y = "") +
    # Rotate x-axis text 45 degrees and align it to the right
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  pages = 1
)
dev.off()

pl <- plot(fire_enn_sd_bam_full_viz, allTerms = T) +
  theme_get() +
  labs(title = NULL)

pdf(
  "Figs/fire_enn_sd_bam_full_all_contours.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
print(pl, pages = 1)
dev.off()

pdf("Figs/fire_enn_sd_bam_full_trend.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_enn_sd_bam_full_viz, select = 2, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_enn_sd_bam_full_TSLF.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_enn_sd_bam_full_viz, select = 3, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

# pdf("Figs/fire_enn_sd_bam_full_dewpoint_temp_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 4, type = "response") + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_temp_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 5) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_pot_evap_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 6) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_precip_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 7) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_solrad_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 8) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_total_evap_max.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 9) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_wind_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 10) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_wind_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 11) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_water_soil_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 12) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_water_soil_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 13) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_savanna.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 14) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_forest_plantation.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 15) + l_fitLine() + l_ciLine() # scale estimate is zero for input data
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_grassland.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 16) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_pasture.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 17) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_sugar_cane.pdf", paper = "a4r", width = 0, height = 0) #scale estimate is zero for input data
# plot(fire_enn_sd_bam_full_viz, select = 18) + l_rug() + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_mosaic_uses.pdf", paper = "a4r", width = 0, height = 0) #scale estimate is zero for input data
# plot(fire_enn_sd_bam_full_viz, select = 19) + l_dens(type = "cond") + l_rug() + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_urban.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 20)  + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_soybean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 21) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_non_vegetated.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 22) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_forest.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 23) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_wet_herbaceous.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 24) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_NDVI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 25) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_pop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 26) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_elevation.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 27) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_aspect.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 28) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_TPI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 29) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_enn_sd_bam_full_TRI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_enn_sd_bam_full_viz, select = 30) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()

# January
pdf("Figs/fire_enn_sd_bam_full_Jan.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# February
pdf("Figs/fire_enn_sd_bam_full_Feb.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# March
pdf("Figs/fire_enn_sd_bam_full_Mar.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# April
pdf("Figs/fire_enn_sd_bam_full_Apr.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# May
pdf("Figs/fire_enn_sd_bam_full_May.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# June
pdf("Figs/fire_enn_sd_bam_full_Jun.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# July
pdf("Figs/fire_enn_sd_bam_full_Jul.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# August
pdf("Figs/fire_enn_sd_bam_full_Aug.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# September
pdf("Figs/fire_enn_sd_bam_full_Sep.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# October
pdf("Figs/fire_enn_sd_bam_full_Oct.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# November
pdf("Figs/fire_enn_sd_bam_full_Nov.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# December
pdf("Figs/fire_enn_sd_bam_full_Dec.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# Diagnostics
pdf("Figs/fire_enn_sd_bam_full_check.pdf", paper = "a4r", width = 0, height = 0)
check(fire_enn_sd_bam_full_viz)
dev.off()

## Fire Scar Density (lsm_c_pd) ----------------------------

head(fire_metrics_df)
print(tail(fire_metrics_df))
skim(fire_metrics_df$lsm_c_pd)

# Check distribution

pdf("Figs/histogram_lsm_c_pd.pdf", paper = "a4", width = 0, height = 0)
hist(log(fire_metrics_df$lsm_c_pd + 1))
dev.off()

fire_pd_bam_tweedie <- bam(
  lsm_c_pd ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = Tweedie(p = 1.5, link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_pd_bam_tweedie)


fire_pd_bam_scat <- bam(
  lsm_c_pd ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = scat(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_pd_bam_scat)

fire_pd_bam_gamma <- bam(
  lsm_c_pd ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = Gamma(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_pd_bam_gamma)

fire_pd_bam_gaussian <- bam(
  lsm_c_pd ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_pd_bam_gaussian)

fire_pd_bam_invgaussian <- bam(
  lsm_c_pd ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = inverse.gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_pd_bam_invgaussian)

fire_pd_bam_loggaussian <- bam(
  log1p(lsm_c_pd) ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_pd_bam_loggaussian)

AIC(
  fire_pd_bam_tweedie,
  fire_pd_bam_scat,
  fire_pd_bam_gamma,
  fire_pd_bam_gaussian,
  fire_pd_bam_invgaussian,
  fire_pd_bam_loggaussian
)

bbmle::AICctab(
  fire_pd_bam_tweedie,
  fire_pd_bam_scat,
  fire_pd_bam_gamma,
  fire_pd_bam_gaussian,
  fire_pd_bam_invgaussian,
  fire_pd_bam_loggaussian
)

pdf(
  "Figs/gamcheck_invgaussian_lsm_c_pd.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
par(mfrow = c(2, 2))
check(getViz(fire_pd_bam_invgaussian))
dev.off()

# Null model
fire_pd_bam_null <- bam(
  lsm_c_pd ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15),
  data = fire_metrics_df,
  family = inverse.gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)

summary(fire_pd_bam_null)

# Full model
sort(TentativeRoughFix(fire_pd_boruta)$finalDecision)

fire_pd_bam_full <- bam(
  lsm_c_pd ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(savanna, bs = "cr") +
    s(forest_plantation, bs = "cr") +
    s(grassland, bs = "cr") +
    s(pasture, bs = "cr") +
    s(sugar_cane, bs = "cr") +
    s(mosaic_uses, bs = "cr") +
    s(non_vegetated, bs = "cr") +
    s(forest, bs = "cr") +
    s(temporary_crop, bs = "cr") +
    s(perennial_crop, bs = "cr") +
    s(water, bs = "cr") +
    s(wet_herbaceous, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(aspect, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = inverse.gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_pd_bam_full)
saveRDS(fire_pd_bam_full, "Output/fire_pd_bam_full.rds")

# Simpler model (lulc only native and anthropic areas)

fire_pd_bam_simple <- bam(
  lsm_c_pd ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(native_area, bs = "cr") +
    s(anthropic_area, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(aspect, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = inverse.gaussian(link = "log"),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_pd_bam_simple)

# Compare models
anova(fire_pd_bam_null, fire_pd_bam_full, test = "Chisq")

AIC(fire_pd_bam_null, fire_pd_bam_full, fire_pd_bam_simple)
bbmle::AICtab(fire_pd_bam_null, fire_pd_bam_full, fire_pd_bam_simple)

# Plots

# Full model
fire_pd_bam_full <- readRDS("Output/fire_pd_bam_full.rds")
fire_pd_bam_full_viz <- getViz(fire_pd_bam_full)

pdf("Figs/fire_pd_bam_full_all.pdf", paper = "a4r", width = 0, height = 0)
print(
  plot(fire_pd_bam_full_viz, allTerms = TRUE) +
    l_fitLine() +
    l_ciLine() +
    labs(y = "") +
    # Rotate x-axis text 45 degrees and align it to the right
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  pages = 1
)
dev.off()

pl <- plot(fire_pd_bam_full_viz, allTerms = T) +
  theme_get() +
  labs(title = NULL)

pdf(
  "Figs/fire_pd_bam_full_all_contours.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
print(pl, pages = 1)
dev.off()

pdf("Figs/fire_pd_bam_full_trend.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_pd_bam_full_viz, select = 2, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_pd_bam_full_TSLF.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_pd_bam_full_viz, select = 3, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

# pdf("Figs/fire_pd_bam_full_dewpoint_temp_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 4, type = "response") + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_temp_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 5) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_pot_evap_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 6) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_precip_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 7) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_solrad_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 8) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_total_evap_max.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 9) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_wind_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 10) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_water_soil_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 11) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_water_soil_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 12) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_savanna.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 13) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_grassland.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 14) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_mosaic_uses.pdf", paper = "a4r", width = 0, height = 0) #scale estimate is zero for input data
# plot(fire_pd_bam_full_viz, select = 15) + l_dens(type = "cond") + l_rug() + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_non_vegetated.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 16) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_forest.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 17) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_temporary_crop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 18) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_wet_herbaceous.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 19) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_NDVI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 20) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_pop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 21) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_elevation.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 22) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_aspect.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 23) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_TPI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 24) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_pd_bam_full_TRI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_pd_bam_full_viz, select = 25) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()

# January
pdf("Figs/fire_pd_bam_full_Jan.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# February
pdf("Figs/fire_pd_bam_full_Feb.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# March
pdf("Figs/fire_pd_bam_full_Mar.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# April
pdf("Figs/fire_pd_bam_full_Apr.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# May
pdf("Figs/fire_pd_bam_full_May.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# June
pdf("Figs/fire_pd_bam_full_Jun.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# July
pdf("Figs/fire_pd_bam_full_Jul.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# August
pdf("Figs/fire_pd_bam_full_Aug.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# September
pdf("Figs/fire_pd_bam_full_Sep.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# October
pdf("Figs/fire_pd_bam_full_Oct.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# November
pdf("Figs/fire_pd_bam_full_Nov.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# December
pdf("Figs/fire_pd_bam_full_Dec.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# Diagnostics
pdf("Figs/fire_pd_bam_full_check.pdf", paper = "a4r", width = 0, height = 0)
check(fire_pd_bam_full_viz)
dev.off()

## Fire Scar Fragmentation (lsm_c_split) ----------------------------
head(fire_metrics_df)
print(tail(fire_metrics_df))
skim(fire_metrics_df$lsm_c_split)
# Transform for squared meters for numeric optimization
fire_metrics_df$lsm_c_split <- fire_metrics_df$lsm_c_split / 1000000
fire_metrics_df$lsm_c_split <- log(fire_metrics_df$lsm_c_split)
# Check distribution

pdf("Figs/histogram_lsm_c_split.pdf", paper = "a4", width = 0, height = 0)
hist(fire_metrics_df$lsm_c_split)
dev.off()

# fire_split_bam_tweedie <- bam(lsm_c_split ~ s(t, bs = "cr"),
#                            data = fire_metrics_df,
#                            family = Tweedie(p = 1.5),
#                            discrete = T,
#                            method = "fREML",
#                            nthreads = parallel::detectCores()-2,
#                            samfrac = 0.1)
# summary(fire_split_bam_tweedie)

fire_split_bam_scat <- bam(
  lsm_c_split ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = scat(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_split_bam_scat)

# fire_split_bam_gamma <- bam(lsm_c_split ~ s(t, bs = "cr"),
#                          data = fire_metrics_df,
#                          family = Gamma(),
#                          discrete = T,
#                          method = "fREML",
#                          nthreads = parallel::detectCores()-2,
#                          samfrac = 0.1)
# summary(fire_split_bam_gamma)

fire_split_bam_gaussian <- bam(
  lsm_c_split ~ s(t, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1
)
summary(fire_split_bam_gaussian)

AIC(fire_split_bam_scat, fire_split_bam_gaussian)

bbmle::AICtab(fire_split_bam_scat, fire_split_bam_gaussian)


pdf(
  "Figs/gamcheck_gaussian_lsm_c_split.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
par(mfrow = c(2, 2))
check(getViz(fire_split_bam_gaussian))
dev.off()

# Null model
fire_split_bam_null <- bam(
  lsm_c_split ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)

summary(fire_split_bam_null)

pdf("Figs/gamcheck_null_lsm_c_split.pdf", paper = "a4", width = 0, height = 0)
par(mfrow = c(2, 2))
check(getViz(fire_split_bam_null))
dev.off()

# Full model
sort(fire_split_boruta$finalDecision)

fire_split_bam_full <- bam(
  lsm_c_split ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(savanna, bs = "cr") +
    s(grassland, bs = "cr") +
    s(pasture, bs = "cr") +
    s(sugar_cane, bs = "cr") +
    s(mosaic_uses, bs = "cr") +
    s(urban, bs = "cr") +
    s(soybean, bs = "cr") +
    s(non_vegetated, bs = "cr") +
    s(forest, bs = "cr") +
    s(temporary_crop, bs = "cr") +
    s(perennial_crop, bs = "cr") +
    s(water, bs = "cr") +
    s(wet_herbaceous, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_split_bam_full)
saveRDS(fire_split_bam_full, "Output/fire_split_bam_full.rds")

# Simpler model (lulc only native and anthropic areas)

fire_split_bam_simple <- bam(
  lsm_c_split ~
    # Spatially Varying Seasonality
    te(x, y, month, d = c(2, 1), bs = c("tp", "cc"), k = c(50, 12)) +
    # Global Trend
    s(t, bs = "cr", k = 20) +
    # Fuel Accumulation Curve (The Process)
    # k=15 allows for the "Growth -> Saturation" shape
    s(months_since_fire, bs = "cr", k = 15) +
    s(dewpoint_temp_sd, bs = "cr") +
    s(pot_evap_sd, bs = "cr") +
    s(precip_min, bs = "cr") +
    s(solrad_mean, bs = "cr") +
    s(total_evap_max, bs = "cr") +
    s(wind_min, bs = "cr") +
    s(wind_sd, bs = "cr") +
    s(water_soil_mean, bs = "cr") +
    s(water_soil_sd, bs = "cr") +
    s(native_area, bs = "cr") +
    s(anthropic_area, bs = "cr") +
    s(NDVI, bs = "cr") +
    s(pop, bs = "cr") +
    s(elevation, bs = "cr") +
    s(TPI, bs = "cr") +
    s(TRI, bs = "cr"),
  data = fire_metrics_df,
  family = gaussian(),
  discrete = T,
  method = "fREML",
  nthreads = parallel::detectCores() - 2,
  samfrac = 0.1,
  knots = list(month = c(0.5, 12.5))
)
summary(fire_split_bam_simple)

# Compare models
anova(fire_split_bam_null, fire_split_bam_full, test = "Chisq")

AIC(fire_split_bam_null, fire_split_bam_full, fire_split_bam_simple)
bbmle::AICctab(fire_split_bam_null, fire_split_bam_full, fire_split_bam_simple)

# Plots

# Full model
fire_split_bam_full <- readRDS("Output/fire_split_bam_full.rds")
fire_split_bam_full_viz <- getViz(fire_split_bam_full)

pdf("Figs/fire_split_bam_full_all.pdf", paper = "a4r", width = 0, height = 0)
print(
  plot(fire_split_bam_full_viz, allTerms = TRUE) +
    l_fitLine() +
    l_ciLine() +
    labs(y = "") +
    # Rotate x-axis text 45 degrees and align it to the right
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  pages = 1
)
dev.off()

pl <- plot(fire_split_bam_full_viz, allTerms = T) +
  theme_get() +
  labs(title = NULL)

pdf(
  "Figs/fire_split_bam_full_all_contours.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
print(pl, pages = 1)
dev.off()

pdf("Figs/fire_split_bam_full_trend.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_split_bam_full_viz, select = 2, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

pdf("Figs/fire_split_bam_full_TSLF.pdf", paper = "a4r", width = 0, height = 0)
plot(fire_split_bam_full_viz, select = 3, type = "response") +
  l_dens(type = "cond") +
  l_fitLine() +
  l_ciLine()
dev.off()

#
# pdf("Figs/fire_split_bam_full_dewpoint_temp_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 4, type = "response") + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_temp_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 5) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_pot_evap_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 6) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_precip_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 7) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_solrad_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 8) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_total_evap_max.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 9) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_wind_min.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 10) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_wind_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 11) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_water_soil_mean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 12) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_water_soil_sd.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 13) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_savanna.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 14) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_forest_plantation.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 15)  + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_grassland.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 16) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_pasture.pdf", paper = "a4r", width = 0, height = 0) #scale estimate is zero for input data
# plot(fire_split_bam_full_viz, select = 17) + l_dens(type = "cond") + l_rug() + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_mosaic_uses.pdf", paper = "a4r", width = 0, height = 0) #scale estimate is zero for input data
# plot(fire_split_bam_full_viz, select = 18) + l_dens(type = "cond") + l_rug() + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_soybean.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 19) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_non_vegetated.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 20) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_forest.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 21) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_perennial_crop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 22)  + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_wet_herbaceous.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 23) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_NDVI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 24) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_pop.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 25) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_elevation.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 26) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_TPI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 27) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()
#
# pdf("Figs/fire_split_bam_full_TRI.pdf", paper = "a4r", width = 0, height = 0)
# plot(fire_split_bam_full_viz, select = 28) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# dev.off()

# January
pdf("Figs/fire_split_bam_full_Jan.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# February
pdf("Figs/fire_split_bam_full_Feb.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# March
pdf("Figs/fire_split_bam_full_Mar.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# April
pdf("Figs/fire_split_bam_full_Apr.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# May
pdf("Figs/fire_split_bam_full_May.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# June
pdf("Figs/fire_split_bam_full_Jun.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# July
pdf("Figs/fire_split_bam_full_Jul.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# August
pdf("Figs/fire_split_bam_full_Aug.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# September
pdf("Figs/fire_split_bam_full_Sep.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# October
pdf("Figs/fire_split_bam_full_Oct.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# November
pdf("Figs/fire_split_bam_full_Nov.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# December
pdf("Figs/fire_split_bam_full_Dec.pdf", paper = "a4r", width = 0, height = 0)
plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour()
dev.off()

# Diagnostics
pdf("Figs/fire_split_bam_full_check.pdf", paper = "a4r", width = 0, height = 0)
check(fire_split_bam_full_viz)
dev.off()


## Final Figures -----------------------------------------------------------

### Partial effects -------------------------------------------------------

#### Fire occurrence -------------------------------------------------------
fire_occ_bam_full <- readRDS("Output/fire_occ_bam_full.rds")
fire_occ_bam_full_viz <- getViz(fire_occ_bam_full)

p1 <- plot(sm(fire_occ_bam_full_viz, 2)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "A", x = "Time", y = "")

p2 <- plot(sm(fire_occ_bam_full_viz, 3)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "B", x = "Months since last fire", y = "")

p3 <- plot(sm(fire_occ_bam_full_viz, 11)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "C", x = "Volumetric soil water", y = "")

p4 <- plot(sm(fire_occ_bam_full_viz, 13)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "D",
    x = "Amount of savannas (ha)",
    y = "Partial effect on fire occurrence"
  )

p5 <- plot(sm(fire_occ_bam_full_viz, 15)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "E", x = "Amount of grasslands (ha)", y = "")

p6 <- plot(sm(fire_occ_bam_full_viz, 21)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "F", x = "Amount of forests (ha)", y = "")

p7 <- plot(sm(fire_occ_bam_full_viz, 25)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "G", x = "Amount of wet\n herbaceous vegetation", y = "")

p8 <- plot(sm(fire_occ_bam_full_viz, 7)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "H", x = "Solar radiation", y = "")

p9 <- plot(sm(fire_occ_bam_full_viz, 12)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "I", x = "Volumetric soil water variation", y = "")

# Final figure
final_fig <- (p1$ggObj | p2$ggObj | p3$ggObj) /
  (p4$ggObj | p5$ggObj | p6$ggObj) /
  (p7$ggObj | p8$ggObj | p9$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/fire_occ_partialeffects.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

#### Fire abundance and shape -------------------------------------------------------
fire_pd_bam_full <- readRDS("Output/fire_pd_bam_full.rds")
fire_pd_bam_full_viz <- getViz(fire_pd_bam_full)

fire_circle_mn_bam_full <- readRDS("Output/fire_circle_mn_bam_full.rds")
fire_circle_mn_bam_full_viz <- getViz(fire_circle_mn_bam_full)


p1 <- plot(sm(fire_pd_bam_full_viz, 2)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "A", x = "Time", y = "")

p2 <- plot(sm(fire_circle_mn_bam_full_viz, 2)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "B", x = "Time", y = "")

p3 <- plot(sm(fire_pd_bam_full_viz, 3)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "C", x = "Months since last fire", y = "")

p4 <- plot(sm(fire_circle_mn_bam_full_viz, 3)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "D", x = "Months since last fire", y = "")

p5 <- plot(sm(fire_pd_bam_full_viz, 5)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "E",
    x = "Potential evaporation variation",
    y = "Partial effect on fire density"
  )

p6 <- plot(sm(fire_circle_mn_bam_full_viz, 5)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "D",
    x = "Potential evaporation variation",
    y = "Partial effect on fire shape regularity"
  )

p7 <- plot(sm(fire_pd_bam_full_viz, 26)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "E", x = "Population density", y = "")

p8 <- plot(sm(fire_circle_mn_bam_full_viz, 22)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "F", x = "Amount of temporary crops (ha)", y = "")


# Final figure
final_fig <- (p1$ggObj | p2$ggObj) /
  (p3$ggObj | p4$ggObj) /
  (p5$ggObj | p6$ggObj) /
  (p7$ggObj | p8$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/fire_abund_shape_partialeffects.pdf",
  final_fig,
  width = 8.27,
  height = 11.69,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

#### Fire configuration -------------------------------------------------------
fire_split_bam_full <- readRDS("Output/fire_split_bam_full.rds")
fire_split_bam_full_viz <- getViz(fire_split_bam_full)

fire_enn_sd_bam_full <- readRDS("Output/fire_enn_sd_bam_full.rds")
fire_enn_sd_bam_full_viz <- getViz(fire_enn_sd_bam_full)

p1 <- plot(sm(fire_split_bam_full_viz, 2)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "A", x = "Time", y = "")

p2 <- plot(sm(fire_enn_sd_bam_full_viz, 2)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "B", x = "Time", y = "")

p3 <- plot(sm(fire_split_bam_full_viz, 3)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "C", x = "Months since last fire", y = "")

p4 <- plot(sm(fire_enn_sd_bam_full_viz, 3)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "D", x = "Months since last fire", y = "")

p5 <- plot(sm(fire_split_bam_full_viz, 5)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "E",
    x = "Potential evaporation variation",
    y = "Partial effect on fire fragmentation"
  )

p6 <- plot(sm(fire_enn_sd_bam_full_viz, 5)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "D",
    x = "Potential evaporation variation",
    y = "Partial effect on fire isolation variation"
  )

p7 <- plot(sm(fire_split_bam_full_viz, 7)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "E", x = "Solar radiation", y = "")

p8 <- plot(sm(fire_enn_sd_bam_full_viz, 28)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "F", x = "Elevation (m)", y = "")


# Final figure
final_fig <- (p1$ggObj | p2$ggObj) /
  (p3$ggObj | p4$ggObj) /
  (p5$ggObj | p6$ggObj) /
  (p7$ggObj | p8$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/fire_config_partialeffects.pdf",
  final_fig,
  width = 8.27,
  height = 11.69,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

#### Fire interior -------------------------------------------------------
fire_core_mn_bam_full <- readRDS("Output/fire_core_mn_bam_full.rds")
fire_core_mn_bam_full_viz <- getViz(fire_core_mn_bam_full)

fire_cai_sd_bam_full <- readRDS("Output/fire_cai_sd_bam_full.rds")
fire_cai_sd_bam_full_viz <- getViz(fire_cai_sd_bam_full)

fire_dcore_sd_bam_full <- readRDS("Output/fire_dcore_sd_bam_full.rds")
fire_dcore_sd_bam_full_viz <- getViz(fire_dcore_sd_bam_full)

p1 <- plot(sm(fire_core_mn_bam_full_viz, 2)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "A", x = "Time", y = "")

p2 <- plot(sm(fire_cai_sd_bam_full_viz, 2)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "B", x = "Time", y = "")

p3 <- plot(sm(fire_dcore_sd_bam_full_viz, 2)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "C", x = "Time", y = "")

p4 <- plot(sm(fire_core_mn_bam_full_viz, 3)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "D",
    x = "Months since last fire",
    y = "Partial effect on fire core size"
  )

p5 <- plot(sm(fire_cai_sd_bam_full_viz, 3)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "E",
    x = "Months since last fire",
    y = "Partial effect on fire integrity variation"
  )

p6 <- plot(sm(fire_dcore_sd_bam_full_viz, 3)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "F",
    x = "Months since last fire",
    y = "Partial effect on fire cohesion variation"
  )

p7 <- plot(sm(fire_core_mn_bam_full_viz, 4)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "G", x = "Solar radiation", y = "")

p8 <- plot(sm(fire_cai_sd_bam_full_viz, 5)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "H", x = "Potential evaporation variation", y = "")

p9 <- plot(sm(fire_dcore_sd_bam_full_viz, 11)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "I", x = "Volumetric soil water variation", y = "")

p10 <- plot(sm(fire_cai_sd_bam_full_viz, 7)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "J", x = "Solar radiation", y = "")

p11 <- plot(sm(fire_dcore_sd_bam_full_viz, 20)) +
  l_fitLine(color = "blue", linewidth = 0.5) +
  l_ciPoly(alpha = 0.8) +
  theme_minimal() +
  labs(title = "K", x = "Population density", y = "")


# Final figure
final_fig <- (p1$ggObj | p2$ggObj | p3$ggObj) /
  (p4$ggObj | p5$ggObj | p6$ggObj) /
  (p7$ggObj | p8$ggObj | p9$ggObj) /
  (plot_spacer() | p10$ggObj | p11$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/fire_interior_partialeffects.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

### Partial effect Maps --------------------------------------------------------------------

# Load the Cerrado geometry
Cerrado <- readRDS("Data/Cerrado.rds")
print(Cerrado)

# --- Define a metric projection for Brazil ---
# SIRGAS 2000 / Brazil Polyconic (EPSG:5880) is a good choice for the whole country.
metric_crs <- "EPSG:5880"

# Select only the geometry object from "Cerrado"
Cerrado <- Cerrado %>% dplyr::select(geometry)
Cerrado <- sf::st_transform(Cerrado, metric_crs)

states <- geobr::read_state()
states <- sf::st_transform(states, metric_crs)

x_lim <- st_bbox(Cerrado$geometry)[c(1, 3)]
y_lim <- st_bbox(Cerrado$geometry)[c(2, 4)]

pad <- 0.02
x_pad <- diff(x_lim) * pad
y_pad <- diff(y_lim) * pad

x_lim <- x_lim + c(-x_pad, x_pad)
y_lim <- y_lim + c(-y_pad, y_pad)

#### Fire occurrence (fire_occ) ------------------------

p1 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "A", subtitle = "Jan") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p2 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "B", subtitle = "Feb") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p3 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "C", subtitle = "Mar") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p4 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "D", subtitle = "Apr") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p5 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "E", subtitle = "May") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p6 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "F", subtitle = "Jun") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p7 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "G", subtitle = "Jul") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p8 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "H", subtitle = "Aug") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p9 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "I", subtitle = "Sep") +
  theme(plot.subtitle = element_text(hjust = 0.6))

p10 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "J", subtitle = "Oct") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p11 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "K", subtitle = "Nov") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p12 <- plot(sm(fire_occ_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3.5, 2.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "L", subtitle = "Dec") +
  theme(plot.subtitle = element_text(hjust = 0.5))

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank()
)


p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p3 <- p3 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean
p8 <- p8 + theme_map_clean
p9 <- p9 + theme_map_clean
p10 <- p10 + theme_map_clean
p11 <- p11 + theme_map_clean
p12 <- p12 + theme_map_clean

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")
# p8 <- p8 + theme(legend.position = "none")
p9 <- p9 + theme(legend.position = "none")
p10 <- p10 + theme(legend.position = "none")
p11 <- p11 + theme(legend.position = "none")
p12 <- p12 + theme(legend.position = "none")

# Final figure
final_fig <- (p1$ggObj | p2$ggObj | p3$ggObj | p4$ggObj) /
  (p5$ggObj | p6$ggObj | p7$ggObj | p8$ggObj) /
  (p9$ggObj | p10$ggObj | p11$ggObj | p12$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/bam_fire_occ_months.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

#### Fire Scar Integrity Variation (lsm_c_cai_sd) ------------------------

p1 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "A", subtitle = "Jan") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p2 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "B", subtitle = "Feb") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p3 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "C", subtitle = "Mar") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p4 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "D", subtitle = "Apr") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p5 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "E", subtitle = "May") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p6 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "F", subtitle = "Jun") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p7 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "G", subtitle = "Jul") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p8 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "H", subtitle = "Aug") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p9 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "I", subtitle = "Sep") +
  theme(plot.subtitle = element_text(hjust = 0.6))

p10 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "J", subtitle = "Oct") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p11 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "K", subtitle = "Nov") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p12 <- plot(sm(fire_cai_sd_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.75),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "L", subtitle = "Dec") +
  theme(plot.subtitle = element_text(hjust = 0.5))

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank()
)


p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p3 <- p3 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean
p8 <- p8 + theme_map_clean
p9 <- p9 + theme_map_clean
p10 <- p10 + theme_map_clean
p11 <- p11 + theme_map_clean
p12 <- p12 + theme_map_clean

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")
# p8 <- p8 + theme(legend.position = "none")
p9 <- p9 + theme(legend.position = "none")
p10 <- p10 + theme(legend.position = "none")
p11 <- p11 + theme(legend.position = "none")
p12 <- p12 + theme(legend.position = "none")

# Final figure
final_fig <- (p1$ggObj | p2$ggObj | p3$ggObj | p4$ggObj) /
  (p5$ggObj | p6$ggObj | p7$ggObj | p8$ggObj) /
  (p9$ggObj | p10$ggObj | p11$ggObj | p12$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/bam_fire_cai_sd_months.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

#### Fire Scar Shape Regularity (lsm_c_circle_mn) ------------------------

p1 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "A", subtitle = "Jan") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p2 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "B", subtitle = "Feb") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p3 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "C", subtitle = "Mar") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p4 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "D", subtitle = "Apr") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p5 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "E", subtitle = "May") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p6 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "F", subtitle = "Jun") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p7 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "G", subtitle = "Jul") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p8 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "H", subtitle = "Aug") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p9 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "I", subtitle = "Sep") +
  theme(plot.subtitle = element_text(hjust = 0.6))

p10 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "J", subtitle = "Oct") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p11 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "K", subtitle = "Nov") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p12 <- plot(sm(fire_circle_mn_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.03, 0.025),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "L", subtitle = "Dec") +
  theme(plot.subtitle = element_text(hjust = 0.5))

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank()
)


p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p3 <- p3 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean
p8 <- p8 + theme_map_clean
p9 <- p9 + theme_map_clean
p10 <- p10 + theme_map_clean
p11 <- p11 + theme_map_clean
p12 <- p12 + theme_map_clean

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")
# p8 <- p8 + theme(legend.position = "none")
p9 <- p9 + theme(legend.position = "none")
p10 <- p10 + theme(legend.position = "none")
p11 <- p11 + theme(legend.position = "none")
p12 <- p12 + theme(legend.position = "none")

# Final figure
final_fig <- (p1$ggObj | p2$ggObj | p3$ggObj | p4$ggObj) /
  (p5$ggObj | p6$ggObj | p7$ggObj | p8$ggObj) /
  (p9$ggObj | p10$ggObj | p11$ggObj | p12$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/bam_fire_circle_mn_months.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

#### Mean Fire Scar Core Size (lsm_c_core_mn) ------------------------

p1 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "A", subtitle = "Jan") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p2 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "B", subtitle = "Feb") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p3 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "C", subtitle = "Mar") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p4 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "D", subtitle = "Apr") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p5 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "E", subtitle = "May") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p6 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "F", subtitle = "Jun") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p7 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "G", subtitle = "Jul") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p8 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "H", subtitle = "Aug") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p9 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "I", subtitle = "Sep") +
  theme(plot.subtitle = element_text(hjust = 0.6))

p10 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "J", subtitle = "Oct") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p11 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "K", subtitle = "Nov") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p12 <- plot(sm(fire_core_mn_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-1.5, 1.5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "L", subtitle = "Dec") +
  theme(plot.subtitle = element_text(hjust = 0.5))

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank()
)


p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p3 <- p3 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean
p8 <- p8 + theme_map_clean
p9 <- p9 + theme_map_clean
p10 <- p10 + theme_map_clean
p11 <- p11 + theme_map_clean
p12 <- p12 + theme_map_clean

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")
# p8 <- p8 + theme(legend.position = "none")
p9 <- p9 + theme(legend.position = "none")
p10 <- p10 + theme(legend.position = "none")
p11 <- p11 + theme(legend.position = "none")
p12 <- p12 + theme(legend.position = "none")

# Final figure
final_fig <- (p1$ggObj | p2$ggObj | p3$ggObj | p4$ggObj) /
  (p5$ggObj | p6$ggObj | p7$ggObj | p8$ggObj) /
  (p9$ggObj | p10$ggObj | p11$ggObj | p12$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/bam_fire_core_mn_months.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

#### Internal Fire Scar Cohesion (lsm_c_dcore_sd) ------------------------

p1 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "A", subtitle = "Jan") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p2 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "B", subtitle = "Feb") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p3 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "C", subtitle = "Mar") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p4 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "D", subtitle = "Apr") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p5 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "E", subtitle = "May") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p6 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "F", subtitle = "Jun") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p7 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "G", subtitle = "Jul") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p8 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "H", subtitle = "Aug") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p9 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "I", subtitle = "Sep") +
  theme(plot.subtitle = element_text(hjust = 0.6))

p10 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "J", subtitle = "Oct") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p11 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "K", subtitle = "Nov") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p12 <- plot(sm(fire_dcore_sd_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-0.4, 0.6),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "L", subtitle = "Dec") +
  theme(plot.subtitle = element_text(hjust = 0.5))

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank()
)


p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p3 <- p3 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean
p8 <- p8 + theme_map_clean
p9 <- p9 + theme_map_clean
p10 <- p10 + theme_map_clean
p11 <- p11 + theme_map_clean
p12 <- p12 + theme_map_clean

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")
# p8 <- p8 + theme(legend.position = "none")
p9 <- p9 + theme(legend.position = "none")
p10 <- p10 + theme(legend.position = "none")
p11 <- p11 + theme(legend.position = "none")
p12 <- p12 + theme(legend.position = "none")

# Final figure
final_fig <- (p1$ggObj | p2$ggObj | p3$ggObj | p4$ggObj) /
  (p5$ggObj | p6$ggObj | p7$ggObj | p8$ggObj) /
  (p9$ggObj | p10$ggObj | p11$ggObj | p12$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/bam_fire_dcore_sd_months.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

#### Fire Scar Isolation Variation (lsm_c_enn_sd) ----------------------------

p1 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "A", subtitle = "Jan") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p2 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "B", subtitle = "Feb") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p3 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "C", subtitle = "Mar") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p4 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "D", subtitle = "Apr") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p5 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "E", subtitle = "May") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p6 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "F", subtitle = "Jun") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p7 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "G", subtitle = "Jul") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p8 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "H", subtitle = "Aug") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p9 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "I", subtitle = "Sep") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p10 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "J", subtitle = "Oct") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p11 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "K", subtitle = "Nov") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p12 <- plot(sm(fire_enn_sd_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-3, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "L", subtitle = "Dec") +
  theme(plot.subtitle = element_text(hjust = 0.5))

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank()
)


p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p3 <- p3 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean
p8 <- p8 + theme_map_clean
p9 <- p9 + theme_map_clean
p10 <- p10 + theme_map_clean
p11 <- p11 + theme_map_clean
p12 <- p12 + theme_map_clean

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")
# p8 <- p8 + theme(legend.position = "none")
p9 <- p9 + theme(legend.position = "none")
p10 <- p10 + theme(legend.position = "none")
p11 <- p11 + theme(legend.position = "none")
p12 <- p12 + theme(legend.position = "none")

# Final figure
final_fig <- (p1$ggObj | p2$ggObj | p3$ggObj | p4$ggObj) /
  (p5$ggObj | p6$ggObj | p7$ggObj | p8$ggObj) /
  (p9$ggObj | p10$ggObj | p11$ggObj | p12$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/bam_fire_enn_sd_months.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)


#### Fire Scar Density (lsm_c_pd) -----------------------------------
p1 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "A", subtitle = "Jan") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p2 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "B", subtitle = "Feb") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p3 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "C", subtitle = "Mar") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p4 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "D", subtitle = "Apr") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p5 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "E", subtitle = "May") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p6 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "F", subtitle = "Jun") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p7 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "G", subtitle = "Jul") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p8 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "H", subtitle = "Aug") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p9 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "I", subtitle = "Sep") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p10 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "J", subtitle = "Oct") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p11 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "K", subtitle = "Nov") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p12 <- plot(sm(fire_pd_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.2, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-2, 3),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "L", subtitle = "Dec") +
  theme(plot.subtitle = element_text(hjust = 0.5))

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank()
)


p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p3 <- p3 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean
p8 <- p8 + theme_map_clean
p9 <- p9 + theme_map_clean
p10 <- p10 + theme_map_clean
p11 <- p11 + theme_map_clean
p12 <- p12 + theme_map_clean

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")
# p8 <- p8 + theme(legend.position = "none")
p9 <- p9 + theme(legend.position = "none")
p10 <- p10 + theme(legend.position = "none")
p11 <- p11 + theme(legend.position = "none")
p12 <- p12 + theme(legend.position = "none")

# Final figure
final_fig <- (p1$ggObj | p2$ggObj | p3$ggObj | p4$ggObj) /
  (p5$ggObj | p6$ggObj | p7$ggObj | p8$ggObj) /
  (p9$ggObj | p10$ggObj | p11$ggObj | p12$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/bam_fire_pd_months.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)


#### Fire Scar Fragmentation (lsm_c_split) -----------------------------------
p1 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 1)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-9, 5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "A", subtitle = "Jan") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p2 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 2)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-9, 5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "B", subtitle = "Feb") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p3 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 3)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-9, 5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "C", subtitle = "Mar") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p4 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 4)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-9, 5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "D", subtitle = "Apr") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p5 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 5)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-9, 5),
    na.value = "transparent",
    name = "Partial effect"
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "E", subtitle = "May") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p6 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 6)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-9, 5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "F", subtitle = "Jun") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p7 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 7)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-9, 5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "G", subtitle = "Jul") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p8 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 8)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-9, 5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "H", subtitle = "Aug") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p9 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 9)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-9, 5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "I", subtitle = "Sep") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p10 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 10)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-7, 5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "J", subtitle = "Oct") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p11 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 11)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-7, 5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "K", subtitle = "Nov") +
  theme(plot.subtitle = element_text(hjust = 0.5))

p12 <- plot(sm(fire_split_bam_full_viz, 1), fix = c("month" = 12)) +
  l_fitRaster() +
  l_fitContour(linewidth = 0.5, colour = "white") +
  geom_sf(data = Cerrado, inherit.aes = F, fill = NA, colour = "red") +
  # geom_sf(data = states, inherit.aes = F, fill = NA)+
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_viridis_c(
    limits = c(-7, 5),
    na.value = "transparent",
    name = "Partial effect" # or "grey90", "white", etc.
  ) +
  theme_minimal() +
  labs(x = "", y = "", title = "L", subtitle = "Dec") +
  theme(plot.subtitle = element_text(hjust = 0.5))

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank()
)


p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p3 <- p3 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean
p8 <- p8 + theme_map_clean
p9 <- p9 + theme_map_clean
p10 <- p10 + theme_map_clean
p11 <- p11 + theme_map_clean
p12 <- p12 + theme_map_clean

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")
# p8 <- p8 + theme(legend.position = "none")
p9 <- p9 + theme(legend.position = "none")
p10 <- p10 + theme(legend.position = "none")
p11 <- p11 + theme(legend.position = "none")
p12 <- p12 + theme(legend.position = "none")

# Final figure
final_fig <- (p1$ggObj | p2$ggObj | p3$ggObj | p4$ggObj) /
  (p5$ggObj | p6$ggObj | p7$ggObj | p8$ggObj) /
  (p9$ggObj | p10$ggObj | p11$ggObj | p12$ggObj)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/bam_fire_split_months.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

## Save summary tables -----------------------------------------------------

### Fire occurrence ---------------------------------------------------------

# 1. Create the Summary Object (Calculates P-values once)
sum_obj <- summary(fire_occ_bam_full)

# 2. Extract Parametric Coefficients (Linear terms)
# Convert matrix to data.frame and move row names to a column
p_table <- as.data.frame(sum_obj$p.table)
p_table$Term <- rownames(p_table)
p_table <- p_table[, c("Term", "Estimate", "Std. Error", "z value", "Pr(>|z|)")] # Reorder
names(p_table) <- c("Term", "Estimate", "Std_Error", "z_value", "P_value")

# 3. Extract Smooth Terms (Non-linear s(...) terms)
s_table <- as.data.frame(sum_obj$s.table)
s_table$Term <- rownames(s_table)
s_table <- s_table[, c("Term", "edf", "Ref.df", "Chi.sq", "p-value")] # Reorder
names(s_table) <- c("Term", "EDF", "Ref_DF", "Chisq_stat", "P_value")

# Optional: Print to check
print(p_table)
print(s_table)

# 4. Save to CSV
write.csv(
  p_table,
  "Output/results/bam_fire_occ_parametric_summary.csv",
  row.names = FALSE
)
write.csv(
  s_table,
  "Output/results/bam_fire_occ_smooth_summary.csv",
  row.names = FALSE
)

### Fire Scar Integrity Variation (lsm_c_cai_sd) ----------------------------

# 1. Create the Summary Object (Calculates P-values once)
sum_obj <- summary(fire_cai_sd_bam_full)

# 2. Extract Parametric Coefficients (Linear terms)
# Convert matrix to data.frame and move row names to a column
p_table <- as.data.frame(sum_obj$p.table)
p_table$Term <- rownames(p_table)
p_table <- p_table[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")] # Reorder
names(p_table) <- c("Term", "Estimate", "Std_Error", "t_value", "P_value")

# 3. Extract Smooth Terms (Non-linear s(...) terms)
s_table <- as.data.frame(sum_obj$s.table)
s_table$Term <- rownames(s_table)
s_table <- s_table[, c("Term", "edf", "Ref.df", "F", "p-value")] # Reorder
names(s_table) <- c("Term", "EDF", "Ref_DF", "F_stat", "P_value")

# Optional: Print to check
print(p_table)
print(s_table)

# 4. Save to CSV
write.csv(
  p_table,
  "Output/results/bam_fire_cai_sd_parametric_summary.csv",
  row.names = FALSE
)
write.csv(
  s_table,
  "Output/results/bam_fire_cai_sd_smooth_summary.csv",
  row.names = FALSE
)

### Fire Scar Shape Regularity (lsm_c_circle_mn) ----------------------------

# 1. Create the Summary Object (Calculates P-values once)
sum_obj <- summary(fire_circle_mn_bam_full)

# 2. Extract Parametric Coefficients (Linear terms)
# Convert matrix to data.frame and move row names to a column
p_table <- as.data.frame(sum_obj$p.table)
p_table$Term <- rownames(p_table)
p_table <- p_table[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")] # Reorder
names(p_table) <- c("Term", "Estimate", "Std_Error", "t_value", "P_value")

# 3. Extract Smooth Terms (Non-linear s(...) terms)
s_table <- as.data.frame(sum_obj$s.table)
s_table$Term <- rownames(s_table)
s_table <- s_table[, c("Term", "edf", "Ref.df", "F", "p-value")] # Reorder
names(s_table) <- c("Term", "EDF", "Ref_DF", "F_stat", "P_value")

# Optional: Print to check
print(p_table)
print(s_table)

# 4. Save to CSV
write.csv(
  p_table,
  "Output/results/bam_fire_circle_mn_parametric_summary.csv",
  row.names = FALSE
)
write.csv(
  s_table,
  "Output/results/bam_fire_circle_mn_smooth_summary.csv",
  row.names = FALSE
)

### Mean Fire Scar Core Size (lsm_c_core_mn) ----------------------------

# 1. Create the Summary Object (Calculates P-values once)
sum_obj <- summary(fire_core_mn_bam_full)

# 2. Extract Parametric Coefficients (Linear terms)
# Convert matrix to data.frame and move row names to a column
p_table <- as.data.frame(sum_obj$p.table)
p_table$Term <- rownames(p_table)
p_table <- p_table[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")] # Reorder
names(p_table) <- c("Term", "Estimate", "Std_Error", "t_value", "P_value")

# 3. Extract Smooth Terms (Non-linear s(...) terms)
s_table <- as.data.frame(sum_obj$s.table)
s_table$Term <- rownames(s_table)
s_table <- s_table[, c("Term", "edf", "Ref.df", "F", "p-value")] # Reorder
names(s_table) <- c("Term", "EDF", "Ref_DF", "F_stat", "P_value")

# Optional: Print to check
print(p_table)
print(s_table)

# 4. Save to CSV
write.csv(
  p_table,
  "Output/results/bam_fire_core_mn_parametric_summary.csv",
  row.names = FALSE
)
write.csv(
  s_table,
  "Output/results/bam_fire_core_mn_smooth_summary.csv",
  row.names = FALSE
)

### Internal Fire Scar Heterogeneity (lsm_c_dcore_sd) ----------------------------

# 1. Create the Summary Object (Calculates P-values once)
sum_obj <- summary(fire_dcore_sd_bam_full)

# 2. Extract Parametric Coefficients (Linear terms)
# Convert matrix to data.frame and move row names to a column
p_table <- as.data.frame(sum_obj$p.table)
p_table$Term <- rownames(p_table)
p_table <- p_table[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")] # Reorder
names(p_table) <- c("Term", "Estimate", "Std_Error", "t_value", "P_value")

# 3. Extract Smooth Terms (Non-linear s(...) terms)
s_table <- as.data.frame(sum_obj$s.table)
s_table$Term <- rownames(s_table)
s_table <- s_table[, c("Term", "edf", "Ref.df", "F", "p-value")] # Reorder
names(s_table) <- c("Term", "EDF", "Ref_DF", "F_stat", "P_value")

# Optional: Print to check
print(p_table)
print(s_table)

# 4. Save to CSV
write.csv(
  p_table,
  "Output/results/bam_fire_dcore_sd_parametric_summary.csv",
  row.names = FALSE
)
write.csv(
  s_table,
  "Output/results/bam_fire_dcore_sd_smooth_summary.csv",
  row.names = FALSE
)

### Fire Scar Isolation Variation (lsm_c_enn_sd) ----------------------------

# 1. Create the Summary Object (Calculates P-values once)
sum_obj <- summary(fire_enn_sd_bam_full)

# 2. Extract Parametric Coefficients (Linear terms)
# Convert matrix to data.frame and move row names to a column
p_table <- as.data.frame(sum_obj$p.table)
p_table$Term <- rownames(p_table)
p_table <- p_table[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")] # Reorder
names(p_table) <- c("Term", "Estimate", "Std_Error", "t_value", "P_value")

# 3. Extract Smooth Terms (Non-linear s(...) terms)
s_table <- as.data.frame(sum_obj$s.table)
s_table$Term <- rownames(s_table)
s_table <- s_table[, c("Term", "edf", "Ref.df", "F", "p-value")] # Reorder
names(s_table) <- c("Term", "EDF", "Ref_DF", "F_stat", "P_value")

# Optional: Print to check
print(p_table)
print(s_table)

# 4. Save to CSV
write.csv(
  p_table,
  "Output/results/bam_fire_enn_sd_parametric_summary.csv",
  row.names = FALSE
)
write.csv(
  s_table,
  "Output/results/bam_fire_enn_sd_smooth_summary.csv",
  row.names = FALSE
)

### Fire Scar Density (lsm_c_pd) ----------------------------

# 1. Create the Summary Object (Calculates P-values once)
sum_obj <- summary(fire_pd_bam_full)

# 2. Extract Parametric Coefficients (Linear terms)
# Convert matrix to data.frame and move row names to a column
p_table <- as.data.frame(sum_obj$p.table)
p_table$Term <- rownames(p_table)
p_table <- p_table[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")] # Reorder
names(p_table) <- c("Term", "Estimate", "Std_Error", "t_value", "P_value")

# 3. Extract Smooth Terms (Non-linear s(...) terms)
s_table <- as.data.frame(sum_obj$s.table)
s_table$Term <- rownames(s_table)
s_table <- s_table[, c("Term", "edf", "Ref.df", "F", "p-value")] # Reorder
names(s_table) <- c("Term", "EDF", "Ref_DF", "F_stat", "P_value")

# Optional: Print to check
print(p_table)
print(s_table)

# 4. Save to CSV
write.csv(
  p_table,
  "Output/results/bam_fire_pd_parametric_summary.csv",
  row.names = FALSE
)
write.csv(s_table, "results/bam_fire_pd_smooth_summary.csv", row.names = FALSE)

### Fire Scar Fragmentation (lsm_c_split) ----------------------------

# 1. Create the Summary Object (Calculates P-values once)
sum_obj <- summary(fire_split_bam_full)

# 2. Extract Parametric Coefficients (Linear terms)
# Convert matrix to data.frame and move row names to a column
p_table <- as.data.frame(sum_obj$p.table)
p_table$Term <- rownames(p_table)
p_table <- p_table[, c("Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")] # Reorder
names(p_table) <- c("Term", "Estimate", "Std_Error", "t_value", "P_value")

# 3. Extract Smooth Terms (Non-linear s(...) terms)
s_table <- as.data.frame(sum_obj$s.table)
s_table$Term <- rownames(s_table)
s_table <- s_table[, c("Term", "edf", "Ref.df", "F", "p-value")] # Reorder
names(s_table) <- c("Term", "EDF", "Ref_DF", "F_stat", "P_value")

# Optional: Print to check
print(p_table)
print(s_table)

# 4. Save to CSV
write.csv(
  p_table,
  "Output/results/bam_fire_split_parametric_summary.csv",
  row.names = FALSE
)
write.csv(
  s_table,
  "Output/results/bam_fire_split_smooth_summary.csv",
  row.names = FALSE
)


## Models' evaluation ------------------------------------------------------

### Fire occurrence (binary) ------------------------------------------------------------

set.seed(123)

# 1) Load data and fitted model

fire_occ_df <- readRDS("Data/fire_occ_df_modeling.rds")


# 1. Get the starting year and month (requires collecting a small part)
start_info <- fire_occ_df %>%
  summarise(start_year = min(year), .by = NULL) %>% # Need .by=NULL or similar in older dplyr?
  as.data.table() # Collect the minimum year
start_year <- start_info$start_year
start_month <- fire_occ_df %>%
  filter(year == start_year) %>%
  summarise(start_month = min(month), .by = NULL) %>%
  as.data.table() # Collect the minimum month for the start year
start_month <- start_month$start_month

# 2. Add the 't' column using mutate
fire_occ_df <- fire_occ_df %>%
  mutate(t = (year - start_year) * 12 + (month - start_month) + 1)

# 1. Create the lazy table
fire_lazy <- lazy_dt(fire_occ_df)

# Define your saturation point (e.g., 4 years)
saturation_point <- 48

# 2. Perform the calculation pipeline
fire_dt_with_tslf <- fire_lazy %>%
  group_by(plot_id) %>%
  arrange(t) %>%
  # --- THE FIX ---
  # We define the interval based on the LAG of fire occurrence.
  # This means if a fire happens at t=5, the new group starts at t=6.
  mutate(
    lag_fire = lag(fire_occurrence, default = 0),
    fire_interval_id = cumsum(lag_fire)
  ) %>%
  group_by(plot_id, fire_interval_id) %>%
  mutate(
    # Now, if a fire happens at month 36, this value will be 36.
    # It resets to 1 in the next month.
    months_since_fire = row_number()
  ) %>%
  # --- Saturation Logic (Same as before) ---
  ungroup() %>%
  mutate(
    # Handle pre-history (assume saturated if we haven't seen a previous fire)
    months_since_fire = if_else(
      fire_interval_id == 0,
      as.numeric(saturation_point),
      as.numeric(months_since_fire)
    ),
    # Cap the maximum
    months_since_fire = if_else(
      months_since_fire > saturation_point,
      as.numeric(saturation_point),
      months_since_fire
    )
  ) %>%
  dplyr::select(-lag_fire, -fire_interval_id) %>% # Clean up helper columns
  as_tibble()

# Check the result
glimpse(fire_dt_with_tslf)

# To see the actual result, collect it
fire_occ_df <- fire_dt_with_tslf %>% as_tibble()
head(fire_occ_df)
print(tail(fire_occ_df))

fire_occ_mod <- fire_occ_bam_full

fire_occ_sf <- st_as_sf(
  fire_occ_df,
  coords = c("x", "y"),
  crs = 5880 # any projected CRS is fine since units are meters
)

form_occ <- fire_occ_mod$formula
family_occ <- family(fire_occ_mod)

# 2) Create spatial folds

cv_spautocor <- cv_spatial_autocor(
  x = sample_n(fire_occ_sf, 100000),
  column = "fire_occurrence"
)
cv_spautocor$range_table

# block size in meters (match ecological scale)
block_size <- 200000 # 200 km

fire_occ_df$block_x <- floor(fire_occ_df$x / block_size)
fire_occ_df$block_y <- floor(fire_occ_df$y / block_size)

fire_occ_df$block_id <- interaction(
  fire_occ_df$block_x,
  fire_occ_df$block_y,
  drop = TRUE
)

set.seed(123)

unique_blocks <- unique(fire_occ_df$block_id)
n_folds <- 5

block_to_fold <- data.frame(
  block_id = unique_blocks,
  fold = sample(rep(1:n_folds, length.out = length(unique_blocks)))
)

fire_occ_df <- fire_occ_df |>
  dplyr::left_join(block_to_fold, by = "block_id")

tss_from_pred <- function(obs, pred) {
  roc_obj <- roc(obs, pred, quiet = TRUE)
  best <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("sensitivity", "specificity")
  )
  best["sensitivity"] + best["specificity"] - 1
}

# 3) Cross-validation loop

cv_results <- lapply(sort(unique(fire_occ_df$fold)), function(f) {
  train_df <- fire_occ_df %>% filter(fold != f)
  test_df <- fire_occ_df %>% filter(fold == f)

  mod <- bam(
    formula = form_occ,
    data = train_df,
    family = family_occ,
    method = "fREML",
    discrete = TRUE,
    nthreads = parallel::detectCores() - 2,
    samfrac = 0.1,
    knots = list(month = c(0.5, 12.5))
  )

  pred <- predict(mod, test_df, type = "response")
  obs <- test_df$fire_occurrence

  roc_obj <- roc(obs, pred, quiet = TRUE)

  tibble(
    fold = f,
    AUC = as.numeric(roc_obj$auc),
    Brier = mean((pred - obs)^2),
    TSS = tss_from_pred(obs, pred)
  )
})

cv_occ_summary <- bind_rows(cv_results)
print(cv_occ_summary)

write.csv(
  cv_occ_summary,
  "Output/results_model_evaluation/cv_occ_summary.csv",
  row.names = FALSE
)

# 4) Final summary

cv_occ_summary %>%
  summarise(
    mean_AUC = mean(AUC),
    sd_AUC = sd(AUC),
    mean_Brier = mean(Brier),
    TSS = mean(TSS$sensitivity)
  )

# # A tibble: 1 × 4
# mean_AUC  sd_AUC mean_Brier   TSS
# <dbl>   <dbl>      <dbl> <dbl>
#   1    0.851 0.00773      0.140 0.541

### Continuous Fire metrics ------------------------------------------------------------
cv_continuous_metric <- function(
  df,
  fitted_model,
  response_var,
  fold_col = "fold"
) {
  form <- fitted_model$formula
  fam <- family(fitted_model)

  res <- lapply(sort(unique(df[[fold_col]])), function(f) {
    train_df <- df %>% filter(.data[[fold_col]] != f)
    test_df <- df %>% filter(.data[[fold_col]] == f)

    mod <- bam(
      formula = form,
      data = train_df,
      family = fam,
      method = "fREML",
      discrete = TRUE,
      nthreads = parallel::detectCores() - 2,
      samfrac = 0.1,
      knots = list(month = c(0.5, 12.5))
    )

    pred <- predict(mod, test_df, type = "response")
    obs <- test_df[[response_var]]

    tibble(
      fold = f,
      RMSE = sqrt(mean((pred - obs)^2)),
      pseudoR2 = 1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2)
    )
  })

  bind_rows(res)
}

fire_metrics_df <- readRDS("Data/fire_metrics_df_modeling.rds")

# 1. Get the starting year and month (requires collecting a small part)
start_info <- fire_metrics_df %>%
  summarise(start_year = min(year), .by = NULL) %>% # Need .by=NULL or similar in older dplyr?
  as.data.table() # Collect the minimum year
start_year <- start_info$start_year
start_month <- fire_metrics_df %>%
  filter(year == start_year) %>%
  summarise(start_month = min(month), .by = NULL) %>%
  as.data.table() # Collect the minimum month for the start year
start_month <- start_month$start_month

# 2. Add the 't' column using mutate
fire_metrics_df <- fire_metrics_df %>%
  mutate(t = (year - start_year) * 12 + (month - start_month) + 1)

# Select only keys and the variable we want to transfer to save memory
tslf_lookup <- fire_occ_df %>%
  dplyr::select(plot_id, t, months_since_fire)

fire_metrics_df <- fire_metrics_df %>%
  inner_join(tslf_lookup, by = c("plot_id", "t"))

# Check the (lazy) result structure
print(fire_metrics_df)

# To see the actual result, collect it
fire_metrics_df <- fire_metrics_df %>% as_tibble()
head(fire_metrics_df)
print(tail(fire_metrics_df))


# block size in meters (match ecological scale)
block_size <- 200000 # 200 km

fire_metrics_df$block_x <- floor(fire_metrics_df$x / block_size)
fire_metrics_df$block_y <- floor(fire_metrics_df$y / block_size)

fire_metrics_df$block_id <- interaction(
  fire_metrics_df$block_x,
  fire_metrics_df$block_y,
  drop = TRUE
)

fire_metrics_df$lsm_c_split <- log(fire_metrics_df$lsm_c_split / 1000000)

set.seed(123)

unique_blocks <- unique(fire_metrics_df$block_id)
n_folds <- 5

block_to_fold <- data.frame(
  block_id = unique_blocks,
  fold = sample(rep(1:n_folds, length.out = length(unique_blocks)))
)

fire_metrics_df <- fire_metrics_df |>
  dplyr::left_join(block_to_fold, by = "block_id")

metrics_models <- list(
  cai_sd = fire_cai_sd_bam_full,
  circle_mn = fire_circle_mn_bam_full,
  core_mn = fire_core_mn_bam_full,
  dcore_sd = fire_dcore_sd_bam_full,
  enn_sd = fire_enn_sd_bam_full,
  pd = fire_pd_bam_full,
  split = fire_split_bam_full
)


metrics_results <- lapply(names(metrics_models), function(m) {
  res <- cv_continuous_metric(
    df = fire_metrics_df,
    fitted_model = metrics_models[[m]],
    response_var = paste0("lsm_c_", m)
  )

  res %>%
    summarise(
      metric = m,
      mean_RMSE = mean(RMSE),
      mean_pseudoR2 = mean(pseudoR2)
    )
})

cv_metrics_summary <- bind_rows(metrics_results)
print(cv_metrics_summary)

# > print(cv_metrics_summary)
# # A tibble: 7 x 3
# metric    mean_RMSE mean_pseudoR2
# <chr>         <dbl>         <dbl>
#   1 cai_sd       13.5        -0.593
# 2 circle_mn     0.120      -1.44
# 3 core_mn      42.1        -0.00591
# 4 dcore_sd      1.31        0.0303
# 5 enn_sd     1083.         -0.415
# 6 pd            0.857       0.174
# 7 split         4.15        0.299

write.csv(
  cv_metrics_summary,
  "Output/results_model_evaluation/cv_metrics_summary.csv",
  row.names = FALSE
)


## Uncertainty maps --------------------------------------------------------

### Fire occurrence (binary) ------------------------------------------------------------

cv_pred_list <- vector("list", n_folds)

for (f in 1:n_folds) {
  train_df <- fire_occ_df |> filter(fold != f)
  test_df <- fire_occ_df |> filter(fold == f)

  mod <- bam(
    formula = fire_occ_bam_full$formula,
    data = train_df,
    family = family(fire_occ_bam_full),
    method = "fREML",
    discrete = TRUE,
    nthreads = parallel::detectCores() - 2,
    samfrac = 0.1,
    knots = list(month = c(0.5, 12.5))
  )

  pred <- predict(
    mod,
    test_df,
    type = "response",
    block.size = 300000,
    discrete = TRUE,
    n.threads = parallel::detectCores() - 2
  )

  cv_pred_list[[f]] <- test_df |>
    transmute(
      x,
      y,
      fold = f,
      pred = pred
    )
}

cv_preds <- bind_rows(cv_pred_list)

saveRDS(cv_preds, "Output/cv_fire_occ_predictions.rds")

pred_summary <- cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

pred_summary <- pred_summary |>
  mutate(
    cv_pred = sd_pred / mean_pred
  )

pdf("Figs/fire_occ_cvfold_sd_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = sd_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction SD") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire occurrence predictions",
    subtitle = "Standard deviation across spatial CV folds"
  )
dev.off()

pdf("Figs/fire_occ_cvfold_cv_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = cv_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction CV") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire occurrence predictions",
    subtitle = "Coefficient of variation across spatial CV folds"
  )
dev.off()

### Fire Scar Integrity Variation (lsm_c_cai_sd) ----------------------------

cv_pred_list <- vector("list", n_folds)

for (f in 1:n_folds) {
  train_df <- fire_metrics_df |> filter(fold != f)
  test_df <- fire_metrics_df |> filter(fold == f)

  mod <- bam(
    formula = fire_cai_sd_bam_full$formula,
    data = train_df,
    family = family(fire_cai_sd_bam_full),
    method = "fREML",
    discrete = TRUE,
    nthreads = parallel::detectCores() - 2,
    samfrac = 0.1,
    knots = list(month = c(0.5, 12.5))
  )

  pred <- predict(
    mod,
    test_df,
    type = "response",
    block.size = 300000,
    discrete = TRUE,
    n.threads = parallel::detectCores() - 2
  )

  cv_pred_list[[f]] <- test_df |>
    transmute(
      x,
      y,
      fold = f,
      pred = pred
    )
}

cv_preds <- bind_rows(cv_pred_list)

saveRDS(cv_preds, "Output/cv_fire_cai_sd_predictions.rds")

pred_summary <- cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

pred_summary <- pred_summary |>
  mutate(
    cv_pred = sd_pred / mean_pred
  )

pdf("Figs/fire_cai_sd_cvfold_sd_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = sd_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction SD") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire scar integrity variation predictions",
    subtitle = "Standard deviation across spatial CV folds"
  )
dev.off()

pdf("Figs/fire_cai_sd_cvfold_cv_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = cv_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction CV") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire scar integrity variation predictions",
    subtitle = "Coefficient of variation across spatial CV folds"
  )
dev.off()

### Fire Scar Shape Regularity (lsm_c_circle_mn) ----------------------------

cv_pred_list <- vector("list", n_folds)

for (f in 1:n_folds) {
  train_df <- fire_metrics_df |> filter(fold != f)
  test_df <- fire_metrics_df |> filter(fold == f)

  mod <- bam(
    formula = fire_circle_mn_bam_full$formula,
    data = train_df,
    family = family(fire_circle_mn_bam_full),
    method = "fREML",
    discrete = TRUE,
    nthreads = parallel::detectCores() - 2,
    samfrac = 0.1,
    knots = list(month = c(0.5, 12.5))
  )

  pred <- predict(
    mod,
    test_df,
    type = "response",
    block.size = 300000,
    discrete = TRUE,
    n.threads = parallel::detectCores() - 2
  )

  cv_pred_list[[f]] <- test_df |>
    transmute(
      x,
      y,
      fold = f,
      pred = pred
    )
}

cv_preds <- bind_rows(cv_pred_list)

saveRDS(cv_preds, "Output/cv_fire_circle_mn_predictions.rds")

pred_summary <- cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

pred_summary <- pred_summary |>
  mutate(
    cv_pred = sd_pred / mean_pred
  )

pdf(
  "Figs/fire_circle_mn_cvfold_sd_pred.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
ggplot(pred_summary, aes(x = x, y = y, color = sd_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction SD") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire shape regularity predictions",
    subtitle = "Standard deviation across spatial CV folds"
  )
dev.off()

pdf(
  "Figs/fire_circle_mn_cvfold_cv_pred.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
ggplot(pred_summary, aes(x = x, y = y, color = cv_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction CV") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire shape regularity predictions",
    subtitle = "Coefficient of variation across spatial CV folds"
  )
dev.off()

### Mean Fire Scar Core Size (lsm_c_core_mn) ----------------------------

cv_pred_list <- vector("list", n_folds)

for (f in 1:n_folds) {
  train_df <- fire_metrics_df |> filter(fold != f)
  test_df <- fire_metrics_df |> filter(fold == f)

  mod <- bam(
    formula = fire_core_mn_bam_full$formula,
    data = train_df,
    family = family(fire_core_mn_bam_full),
    method = "fREML",
    discrete = TRUE,
    nthreads = parallel::detectCores() - 2,
    samfrac = 0.1,
    knots = list(month = c(0.5, 12.5))
  )

  pred <- predict(
    mod,
    test_df,
    type = "response",
    block.size = 300000,
    discrete = TRUE,
    n.threads = parallel::detectCores() - 2
  )

  cv_pred_list[[f]] <- test_df |>
    transmute(
      x,
      y,
      fold = f,
      pred = pred
    )
}

cv_preds <- bind_rows(cv_pred_list)

saveRDS(cv_preds, "Output/cv_fire_core_mn_predictions.rds")

pred_summary <- cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

pred_summary <- pred_summary |>
  mutate(
    cv_pred = sd_pred / mean_pred
  )

pdf("Figs/fire_core_mn_cvfold_sd_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = sd_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction SD") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire core size predictions",
    subtitle = "Standard deviation across spatial CV folds"
  )
dev.off()

pdf("Figs/fire_core_mn_cvfold_cv_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = cv_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction CV") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire core size predictions",
    subtitle = "Coefficient of variation across spatial CV folds"
  )
dev.off()

### Internal Fire Scar Heterogeneity (lsm_c_dcore_sd) ----------------------------

cv_pred_list <- vector("list", n_folds)

for (f in 1:n_folds) {
  train_df <- fire_metrics_df |> filter(fold != f)
  test_df <- fire_metrics_df |> filter(fold == f)

  mod <- bam(
    formula = fire_dcore_sd_bam_full$formula,
    data = train_df,
    family = family(fire_dcore_sd_bam_full),
    method = "fREML",
    discrete = TRUE,
    nthreads = parallel::detectCores() - 2,
    samfrac = 0.1,
    knots = list(month = c(0.5, 12.5))
  )

  pred <- predict(
    mod,
    test_df,
    type = "response",
    block.size = 300000,
    discrete = TRUE,
    n.threads = parallel::detectCores() - 2
  )

  cv_pred_list[[f]] <- test_df |>
    transmute(
      x,
      y,
      fold = f,
      pred = pred
    )
}

cv_preds <- bind_rows(cv_pred_list)

saveRDS(cv_preds, "Output/cv_fire_dcore_sd_predictions.rds")

pred_summary <- cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

pred_summary <- pred_summary |>
  mutate(
    cv_pred = sd_pred / mean_pred
  )

pdf(
  "Figs/fire_dcore_sd_cvfold_sd_pred.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
ggplot(pred_summary, aes(x = x, y = y, color = sd_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction SD") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in internal fire heterogeneity predictions",
    subtitle = "Standard deviation across spatial CV folds"
  )
dev.off()

pdf(
  "Figs/fire_dcore_sd_cvfold_cv_pred.pdf",
  paper = "a4",
  width = 0,
  height = 0
)
ggplot(pred_summary, aes(x = x, y = y, color = cv_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction CV") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in internal fire heterogeneity predictions",
    subtitle = "Coefficient of variation across spatial CV folds"
  )
dev.off()

### Fire Scar Isolation Variation (lsm_c_enn_sd) ----------------------------

cv_pred_list <- vector("list", n_folds)

for (f in 1:n_folds) {
  train_df <- fire_metrics_df |> filter(fold != f)
  test_df <- fire_metrics_df |> filter(fold == f)

  mod <- bam(
    formula = fire_enn_sd_bam_full$formula,
    data = train_df,
    family = family(fire_enn_sd_bam_full),
    method = "fREML",
    discrete = TRUE,
    nthreads = parallel::detectCores() - 2,
    samfrac = 0.1,
    knots = list(month = c(0.5, 12.5))
  )

  pred <- predict(
    mod,
    test_df,
    type = "response",
    block.size = 300000,
    discrete = TRUE,
    n.threads = parallel::detectCores() - 2
  )

  cv_pred_list[[f]] <- test_df |>
    transmute(
      x,
      y,
      fold = f,
      pred = pred
    )
}

cv_preds <- bind_rows(cv_pred_list)

saveRDS(cv_preds, "Output/cv_fire_enn_sd_predictions.rds")

pred_summary <- cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

pred_summary <- pred_summary |>
  mutate(
    cv_pred = sd_pred / mean_pred
  )

pdf("Figs/fire_enn_sd_cvfold_sd_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = sd_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction SD") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire isolation variation predictions",
    subtitle = "Standard deviation across spatial CV folds"
  )
dev.off()

pdf("Figs/fire_enn_sd_cvfold_cv_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = cv_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction CV") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire isolation variation predictions",
    subtitle = "Coefficient of variation across spatial CV folds"
  )
dev.off()

### Fire Scar Density (lsm_c_pd) ----------------------------

cv_pred_list <- vector("list", n_folds)

for (f in 1:n_folds) {
  train_df <- fire_metrics_df |> filter(fold != f)
  test_df <- fire_metrics_df |> filter(fold == f)

  mod <- bam(
    formula = fire_pd_bam_full$formula,
    data = train_df,
    family = family(fire_pd_bam_full),
    method = "fREML",
    discrete = TRUE,
    nthreads = parallel::detectCores() - 2,
    samfrac = 0.1,
    knots = list(month = c(0.5, 12.5))
  )

  pred <- predict(
    mod,
    test_df,
    type = "response",
    block.size = 300000,
    discrete = TRUE,
    n.threads = parallel::detectCores() - 2
  )

  cv_pred_list[[f]] <- test_df |>
    transmute(
      x,
      y,
      fold = f,
      pred = pred
    )
}

cv_preds <- bind_rows(cv_pred_list)

saveRDS(cv_preds, "Output/cv_fire_pd_predictions.rds")

pred_summary <- cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

pred_summary <- pred_summary |>
  mutate(
    cv_pred = sd_pred / mean_pred
  )

pdf("Figs/fire_pd_cvfold_sd_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = sd_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction SD") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire scar density predictions",
    subtitle = "Standard deviation across spatial CV folds"
  )
dev.off()

pdf("Figs/fire_pd_cvfold_cv_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = cv_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction CV") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire scar density predictions",
    subtitle = "Coefficient of variation across spatial CV folds"
  )
dev.off()

### Fire Scar Fragmentation (lsm_c_split) ----------------------------

cv_pred_list <- vector("list", n_folds)

for (f in 1:n_folds) {
  train_df <- fire_metrics_df |> filter(fold != f)
  test_df <- fire_metrics_df |> filter(fold == f)

  mod <- bam(
    formula = fire_split_bam_full$formula,
    data = train_df,
    family = family(fire_split_bam_full),
    method = "fREML",
    discrete = TRUE,
    nthreads = parallel::detectCores() - 2,
    samfrac = 0.1,
    knots = list(month = c(0.5, 12.5))
  )

  pred <- predict(
    mod,
    test_df,
    type = "response",
    block.size = 300000,
    discrete = TRUE,
    n.threads = parallel::detectCores() - 2
  )

  cv_pred_list[[f]] <- test_df |>
    transmute(
      x,
      y,
      fold = f,
      pred = pred
    )
}

cv_preds <- bind_rows(cv_pred_list)

saveRDS(cv_preds, "Output/cv_fire_split_predictions.rds")

pred_summary <- cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

pred_summary <- pred_summary |>
  mutate(
    cv_pred = sd_pred / mean_pred
  )

pdf("Figs/fire_split_cvfold_sd_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = sd_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction SD") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire fragmentation predictions",
    subtitle = "Standard deviation across spatial CV folds"
  )
dev.off()

pdf("Figs/fire_split_cvfold_cv_pred.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_summary, aes(x = x, y = y, color = cv_pred)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(name = "Prediction CV") +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Spatial uncertainty in fire fragmentation predictions",
    subtitle = "Coefficient of variation across spatial CV folds"
  )
dev.off()


### Final Figure ------------------------------------------------------------
fire_occ_cv_preds <- readRDS("Output/cv_fire_occ_predictions.rds")

fire_occ_pred_summary <- fire_occ_cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

p1 <- ggplot(fire_occ_pred_summary, aes(x = x, y = y, fill = sd_pred)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Prediction SD",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "A",
    subtitle = "Fire occurrence probability"
  )

fire_cai_sd_cv_preds <- readRDS("Output/cv_fire_cai_sd_predictions.rds")

fire_cai_sd_pred_summary <- fire_cai_sd_cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

p2 <- ggplot(fire_cai_sd_pred_summary, aes(x = x, y = y, fill = sd_pred)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Prediction SD",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "G",
    subtitle = "Fire integrity variation"
  )

fire_circle_mn_cv_preds <- readRDS("Output/cv_fire_circle_mn_predictions.rds")

fire_circle_mn_pred_summary <- fire_circle_mn_cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

p3 <- ggplot(fire_circle_mn_pred_summary, aes(x = x, y = y, fill = sd_pred)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Prediction SD",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "E",
    subtitle = "Fire shape regularity"
  )


fire_core_mn_cv_preds <- readRDS("Output/cv_fire_core_mn_predictions.rds")

fire_core_mn_pred_summary <- fire_core_mn_cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

p4 <- ggplot(fire_core_mn_pred_summary, aes(x = x, y = y, fill = sd_pred)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Prediction SD",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "F",
    subtitle = "Mean fire core size"
  )

fire_dcore_sd_cv_preds <- readRDS("Output/cv_fire_dcore_sd_predictions.rds")

fire_dcore_sd_pred_summary <- fire_dcore_sd_cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

p5 <- ggplot(fire_dcore_sd_pred_summary, aes(x = x, y = y, fill = sd_pred)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Prediction SD",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "H",
    subtitle = "Fire cohesion variation"
  )

fire_enn_sd_cv_preds <- readRDS("Output/cv_fire_enn_sd_predictions.rds")

fire_enn_sd_pred_summary <- fire_enn_sd_cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

p6 <- ggplot(fire_enn_sd_pred_summary, aes(x = x, y = y, fill = sd_pred)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Prediction SD",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "D",
    subtitle = "Fire isolation variation"
  )

fire_pd_cv_preds <- readRDS("Output/cv_fire_pd_predictions.rds")

fire_pd_pred_summary <- fire_pd_cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

p7 <- ggplot(fire_pd_pred_summary, aes(x = x, y = y, fill = sd_pred)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Prediction SD",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "B",
    subtitle = "Fire density"
  )

fire_split_cv_preds <- readRDS("Output/cv_fire_split_predictions.rds")

fire_split_pred_summary <- fire_split_cv_preds |>
  group_by(x, y) |>
  summarise(
    mean_pred = mean(pred),
    sd_pred = sd(pred),
    .groups = "drop"
  )

p8 <- ggplot(fire_split_pred_summary, aes(x = x, y = y, fill = sd_pred)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Prediction SD",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "C",
    subtitle = "Fire fragmentation"
  )

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  plot.subtitle = element_text(hjust = 0.5)
)

p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p3 <- p3 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean
p8 <- p8 + theme_map_clean


# Final figure
final_fig <-
  (p1 | p7) /
  (p8 | p6) /
  (p3 | p4) /
  (p2 | p5)


final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/fire_sp_cv_fold_figure.pdf",
  final_fig,
  width = 8.27,
  height = 11.69,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

# Cluster analysis --------------------------------------------------------

fire_occ_df <- readRDS(
  "Data/final_fire_occ_data.rds"
)
fire_occ_df <- lazy_dt(fire_occ_df)

fire_occ_m <- fire_occ_df %>%
  group_by(plot_id, longitude, latitude, month) %>%
  summarise(fire_occurrence = mean(fire_occurrence, na.rm = T)) %>%
  as_tibble()

skim(fire_occ_m)

wide_fire_occ_m <- fire_occ_m %>%
  pivot_wider(
    id_cols = c(plot_id), # Columns to keep as identifiers
    names_from = month, # Column containing the new column names (jan, feb...)
    values_from = fire_occurrence # Column containing the values
  )

skim(wide_fire_occ_m)

fire_metrics_df <- readRDS(
  "Data/final_clean_data.rds"
)

fire_metrics_df <- dtplyr::lazy_dt(fire_metrics_df, immutable = FALSE)

# Standardize names and variables
fire_metrics_df <- dplyr::rename(fire_metrics_df, x = longitude)
fire_metrics_df <- dplyr::rename(fire_metrics_df, y = latitude)

fire_metrics_m <- fire_metrics_df %>%
  group_by(plot_id, x, y) %>%
  summarise(
    lsm_c_cai_sd = mean(lsm_c_cai_sd, na.rm = T),
    lsm_c_circle_mn = mean(lsm_c_circle_mn, na.rm = T),
    lsm_c_core_mn = mean(lsm_c_core_mn, na.rm = T),
    lsm_c_dcore_sd = mean(lsm_c_dcore_sd, na.rm = T),
    lsm_c_enn_sd = mean(lsm_c_enn_sd, na.rm = T),
    lsm_c_pd = mean(lsm_c_pd, na.rm = T),
    lsm_c_split = mean(lsm_c_split, na.rm = T)
  ) %>%
  as_tibble()

skim(fire_metrics_m)

fire_m <- left_join(wide_fire_occ_m, fire_metrics_m, by = c("plot_id"))
skim(fire_m)

pdf("Figs/fire_occurrence_jan.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `1`, fill = `1`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_feb.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `2`, fill = `2`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_mar.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `3`, fill = `3`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_apr.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `4`, fill = `4`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_may.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `5`, fill = `5`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_jun.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `6`, fill = `6`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_jul.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `7`, fill = `7`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_aug.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `8`, fill = `8`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_sep.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `9`, fill = `9`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_oct.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `10`, fill = `10`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_nov.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `11`, fill = `11`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_occurrence_dec.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = `12`, fill = `12`)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_cai_sd.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = lsm_c_cai_sd, fill = lsm_c_cai_sd)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_circle_mn.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = lsm_c_circle_mn, fill = lsm_c_circle_mn)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_core_mn.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = lsm_c_core_mn, fill = lsm_c_core_mn)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_dcore_sd.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = lsm_c_dcore_sd, fill = lsm_c_dcore_sd)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_enn_sd.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = lsm_c_enn_sd, fill = lsm_c_enn_sd)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_pd.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = lsm_c_pd, fill = lsm_c_pd)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

pdf("Figs/fire_split.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_m, aes(x = x, y = y, z = lsm_c_split, fill = lsm_c_split)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")
dev.off()

fire_m <- na.omit(fire_m)
dim(fire_m)

fire_m <- fire_m[, -c(1, 14, 15)]

fire_m <- as.matrix(scale(fire_m))

# identify BICs for different models
dataBIC <- mclustBIC(fire_m)
print(summary(dataBIC))

pdf("Figs/mclust_BIC.pdf", paper = "a4", width = 0, height = 0)
plot(dataBIC)
dev.off()

mod <- Mclust(
  fire_m, # data for the cluster model
  G = 5 # BIC index for model to be built
)
summary(mod)
mod[["parameters"]][["mean"]] # mean values of clusters

saveRDS(mod, "Output/Mclust_model.rds")
mod <- readRDS("Output/Mclust_model.rds")
# Run on your original Mclust model (mod)
combi <- clustCombi(mod)
summary(combi)

pdf("Figs/mclustcombi_entPlot.pdf", paper = "a4r", width = 0, height = 0)
plot(combi, what = "entropy") # Look for the "elbow" in entropy
dev.off()

pdf("Figs/mclustcombi.pdf", paper = "a4r", width = 0, height = 0)
plot(combi, what = "classification") # Look for the "elbow" in entropy
dev.off()

pdf("Figs/mclustcombi_tree.pdf", paper = "a4r", width = 0, height = 0)
plot(combi, what = "tree") # Look for the "elbow" in entropy
dev.off()

ICL <- mclustICL(fire_m)
plot(ICL)

LRT <- mclustBootstrapLRT(fire_m, modelName = "VVV")
print(LRT)
plot(LRT, G = 5)

ModPred <- predict.Mclust(mod, fire_m) # prediction
fire_m <- left_join(wide_fire_occ_m, fire_metrics_m, by = c("plot_id"))
fire_m <- na.omit(fire_m)

pred_df <- data.frame(fire_m, class = ModPred$classification)

pdf("Figs/mclust_map_class.pdf", paper = "a4", width = 0, height = 0)
ggplot(pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis(option = "turbo")
dev.off()

saveRDS(pred_df, "Output/predictors_fire_regimes_df.rds")

# Pivot longer for plotting
long_pred_df <- pred_df %>%
  pivot_longer(
    cols = c(
      X1,
      X2,
      X3,
      X4,
      X5,
      X6,
      X7,
      X8,
      X9,
      X10,
      X11,
      X12,
      lsm_c_cai_sd,
      lsm_c_circle_mn,
      lsm_c_core_mn,
      lsm_c_dcore_sd,
      lsm_c_enn_sd,
      lsm_c_pd,
      lsm_c_split
    ), # Add all vars
    names_to = "Metric",
    values_to = "Value"
  )

# Plot
pdf("Figs/boxplot_fire_classes.pdf", paper = "a4", width = 0, height = 0)
ggplot(
  long_pred_df,
  aes(x = as.factor(class), y = Value, fill = as.factor(class))
) +
  geom_boxplot(outlier.shape = NA) + # Hide extreme outliers for clarity
  facet_wrap(~Metric, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = turbo(5)) +
  labs(title = "Fire Regime Profiles")
dev.off()


ModPred$z

pdf("Figs/mclust_density_class.pdf", paper = "a4", width = 0, height = 0)
plot(density(ModPred$z[, 1]))
plot(density(ModPred$z[, 2]))
plot(density(ModPred$z[, 3]))
plot(density(ModPred$z[, 4]))
plot(density(ModPred$z[, 5]))
dev.off()

# Convert the matrix to a dataframe and add IDs
z_df <- as.data.frame(ModPred$z)
colnames(z_df) <- paste0("Regime_", 1:5) # Rename columns nicely
z_df$plot_id <- pred_df$plot_id # Assuming column 1 is plot_id

# 2. Assign the "Hard" Class and Uncertainty for sorting
z_df$Assigned_Class <- factor(ModPred$classification)
z_df$Max_Prob <- apply(ModPred$z, 1, max)

# 3. Pivot to Long Format for ggplot
z_long <- z_df %>%
  pivot_longer(
    cols = starts_with("Regime"),
    names_to = "Regime_Prob",
    values_to = "Probability"
  )

# 4. Create the Ordering
# We sort pixels by: 1. Their assigned class, 2. Their confidence (Max_Prob)
# This puts the "purest" pixels at the center of each block and uncertain ones at edges
plot_order <- z_df %>%
  arrange(Assigned_Class, desc(Max_Prob)) %>%
  pull(plot_id)

z_long$plot_id <- factor(z_long$plot_id, levels = plot_order)

# 5. The Plot
pdf("Figs/mclust_uncertainty_class.pdf", paper = "a4r", width = 0, height = 0)

ggplot(z_long, aes(x = plot_id, y = Probability, fill = Regime_Prob)) +
  geom_col(width = 1, position = "stack") + # No gaps between bars
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title = "Fire Regime Classification Structure",
    subtitle = "Each vertical slice is one pixel. Solid blocks = High certainty. Blends = Uncertainty.",
    x = "Pixels (Ordered by Regime)",
    y = "Probability of Membership"
  ) +
  theme_void() + # Remove axes clutter for 18k points
  theme(legend.position = "bottom")
dev.off()

# 1. Calculate Uncertainty
# mclust defines uncertainty as 1 - probability of the most likely group
uncertainty_df <- data.frame(
  x = pred_df[, "x"], # Assuming you have x/y in your matrix or joined object
  y = pred_df[, "y"],
  Uncertainty = 1 - apply(ModPred$z, 1, max)
)

# 2. Plot Map
pdf(
  "Figs/mclust_uncertainty_map_class.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
ggplot(uncertainty_df, aes(x = x, y = y, fill = Uncertainty)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno") + # Dark = Certain, Bright = Uncertain
  coord_fixed() +
  theme_minimal() +
  labs(title = "Map of Classification Uncertainty")
dev.off()

moddr <- MclustDR(mod)
summary(moddr)

pdf("Figs/mclust_fire_boundaries.pdf", paper = "a4", width = 0, height = 0)
plot(moddr, what = "boundaries", ngrid = 200)
dev.off()


fire_m <- fire_m[, -c(1, 14, 15)]

fire_m <- as.matrix(scale(fire_m))
modda <- MclustDA(fire_m, ModPred$classification, modelNames = "VVV")
summary(modda)

# Classes    n     % Model G
# 1 5515 30.63   VVV 3
# 2 3272 18.17   VVV 4
# 3 2844 15.80   VVV 5
# 4 1058  5.88   XXX 1
# 5 5314 29.52   VVV 4
# Classification error = 0.0605
# Brier score          = 0.0459

modcvda <- cvMclustDA(modda, nfold = 5)
unlist(modcvda[3:6])
# ce       se.ce       brier    se.brier
# 0.072376826 0.002218831 0.054618392 0.001134801

moddr.da <- MclustDR(modda)
summary(moddr.da)

dr_coords <- predict(moddr, newdata = mod$data)$dir
loadings_da <- moddr$std.basis

load_df <- data.frame(
  metric = rownames(loadings_da),
  Dir1 = loadings_da[, 1],
  Dir2 = loadings_da[, 2]
)

arrow_scale <- 8 # adjust for aesthetics

load_df <- load_df |>
  mutate(
    xend = Dir1 * arrow_scale,
    yend = Dir2 * arrow_scale
  )
load_df$metric[load_df$metric == "lsm_c_cai_sd"] <- "integ var"
load_df$metric[load_df$metric == "lsm_c_circle_mn"] <- "shape reg"
load_df$metric[load_df$metric == "lsm_c_core_mn"] <- "core size"
load_df$metric[load_df$metric == "lsm_c_dcore_sd"] <- "cohes var"
load_df$metric[load_df$metric == "lsm_c_enn_sd"] <- "isol var"
load_df$metric[load_df$metric == "lsm_c_pd"] <- "dens"
load_df$metric[load_df$metric == "lsm_c_split"] <- "frag"


plot_df <- data.frame(
  Dir1 = dr_coords[, 1],
  Dir2 = dr_coords[, 2],
  class = factor(mod$classification),
  uncertainty = mod$uncertainty
)

round(prop.table(moddr$evalues), 4) * 100

pdf("Figs/mclustDR_fire_biplot.pdf", paper = "a4r", width = 0, height = 0)
ggplot(plot_df, aes(Dir1, Dir2)) +
  geom_point(
    aes(color = class, alpha = 1 - uncertainty),
    size = 0.8
  ) +
  stat_ellipse(
    aes(color = class),
    type = "norm",
    linewidth = 0.9
  ) +
  scale_alpha(range = c(0.1, 0.3), guide = "none") +
  scale_color_viridis_d(option = "turbo", name = "Fire regime") +
  theme_bw() +
  geom_segment(
    data = load_df,
    aes(
      x = 0,
      y = 0,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.15, "cm")),
    color = "red",
    alpha = 0.5,
    linewidth = 0.5
  ) +
  geom_text(
    data = load_df,
    aes(
      x = xend,
      y = yend,
      label = metric
    ),
    inherit.aes = FALSE,
    size = 3,
    color = "red",
    hjust = 0.5,
    vjust = -0.5
  ) +
  labs(
    x = "Discriminant axis 1 (66.84%)",
    y = "Discriminant axis 2 (24.52%)"
  )
dev.off()

pdf("Figs/mclustDR_fire_boundaries.pdf", paper = "a4", width = 0, height = 0)
plot(moddr.da, what = "boundaries", ngrid = 200)
dev.off()

pdf("Figs/mclustDR_fire_densities.pdf", paper = "a4", width = 0, height = 0)
plot(moddr.da, what = "contour")
dev.off()

# PCA
pca_fire <- prcomp(fire_m[, -c(1, 14, 15)], scale. = T)
summary(pca_fire)

pdf("Figs/pca_fire_class.pdf", paper = "a4", width = 0, height = 0)
autoplot(
  pca_fire,
  data = data.frame(classification = as.factor(ModPred$classification)),
  colour = "classification",
  alpha = 0.25,
  frame = TRUE,
  frame.type = "norm",
  loadings = TRUE,
  loadings.label = TRUE
) +
  scale_color_manual(values = turbo(5)) +
  scale_fill_manual(values = turbo(5))
dev.off()


## Final figures ------------------------------------------------------------

### Figure 1
p1 <- ggplot(pred_df, aes(x = x, y = y, fill = as.factor(class))) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis_d(option = "turbo", name = "Fire regime")

p2 <- ggplot(plot_df, aes(Dir1, Dir2)) +
  geom_point(
    aes(color = class, alpha = 1 - uncertainty),
    size = 0.8
  ) +
  stat_ellipse(
    aes(color = class),
    type = "norm",
    linewidth = 0.9
  ) +
  scale_alpha(range = c(0.1, 0.3), guide = "none") +
  scale_color_viridis_d(option = "turbo", name = "Fire regime") +
  theme_bw() +
  geom_segment(
    data = load_df,
    aes(
      x = 0,
      y = 0,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.15, "cm")),
    color = "red",
    alpha = 0.5,
    linewidth = 0.5
  ) +
  geom_text(
    data = load_df,
    aes(
      x = xend,
      y = yend,
      label = metric
    ),
    inherit.aes = FALSE,
    size = 3,
    color = "red",
    hjust = 0.5,
    vjust = -0.5
  ) +
  labs(
    x = "Discriminant axis 1 (66.84%)",
    y = "Discriminant axis 2 (24.52%)"
  )

# Supress redundancies

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank()
)


p1 <- p1 + theme_map_clean


# Tag panels
p1 <- p1 + labs(tag = "A")
p2 <- p2 + labs(tag = "B")


# Final figure
final_fig <- (p1 / p2)


final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/fire_regimes_da.pdf",
  final_fig,
  width = 8.27,
  height = 11.69,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

### Figure 2

fire_m <- left_join(wide_fire_occ_m, fire_metrics_m, by = c("plot_id"))
skim(fire_m)

p1 <- ggplot(fire_m, aes(x = x, y = y, z = `1`, fill = `1`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Jan") +
  theme_bw()

p2 <- ggplot(fire_m, aes(x = x, y = y, z = `2`, fill = `2`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Feb") +
  theme_bw()

p3 <- ggplot(fire_m, aes(x = x, y = y, z = `3`, fill = `3`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Mar") +
  theme_bw()

p4 <- ggplot(fire_m, aes(x = x, y = y, z = `4`, fill = `4`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Apr") +
  theme_bw()

p5 <- ggplot(fire_m, aes(x = x, y = y, z = `5`, fill = `5`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "May") +
  theme_bw()

p6 <- ggplot(fire_m, aes(x = x, y = y, z = `6`, fill = `6`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Jun") +
  theme_bw()

p7 <- ggplot(fire_m, aes(x = x, y = y, z = `7`, fill = `7`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Jul") +
  theme_bw()

p8 <- ggplot(fire_m, aes(x = x, y = y, z = `8`, fill = `8`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Aug") +
  theme_bw()

p9 <- ggplot(fire_m, aes(x = x, y = y, z = `9`, fill = `9`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Sep") +
  theme_bw()

p10 <- ggplot(fire_m, aes(x = x, y = y, z = `10`, fill = `10`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Oct") +
  theme_bw()

p11 <- ggplot(fire_m, aes(x = x, y = y, z = `11`, fill = `11`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Nov") +
  theme_bw()

p12 <- ggplot(fire_m, aes(x = x, y = y, z = `12`, fill = `12`)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    limits = c(0, 1),
    name = "Fire occurrence"
  ) +
  labs(title = "Dec") +
  theme_bw()

p13 <- ggplot(
  fire_m,
  aes(x = x, y = y, z = lsm_c_cai_sd, fill = lsm_c_cai_sd)
) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    name = "Integrity variation",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  theme_bw()

p14 <- ggplot(
  fire_m,
  aes(x = x, y = y, z = lsm_c_circle_mn, fill = lsm_c_circle_mn)
) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    name = "Shape regularity",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  theme_bw()

p15 <- ggplot(
  fire_m,
  aes(x = x, y = y, z = lsm_c_core_mn, fill = lsm_c_core_mn)
) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    name = "Core size\n(ha)",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  theme_bw()

p16 <- ggplot(
  fire_m,
  aes(x = x, y = y, z = lsm_c_dcore_sd, fill = lsm_c_dcore_sd)
) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    name = "Cohesion variation",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  theme_bw()

p17 <- ggplot(
  fire_m,
  aes(x = x, y = y, z = lsm_c_enn_sd, fill = lsm_c_enn_sd)
) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    name = "Isolation variation",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  theme_bw()

p18 <- ggplot(fire_m, aes(x = x, y = y, z = lsm_c_pd, fill = lsm_c_pd)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    name = "Density\n(n° of fires/100 ha)",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  theme_bw()

p19 <- ggplot(fire_m, aes(x = x, y = y, z = lsm_c_split, fill = lsm_c_split)) +
  geom_tile() +
  scale_fill_viridis_c(
    trans = "log1p",
    name = expression(
      atop("Fragmentation", "(×10"^9 * ")")
    ),
    labels = scales::label_number(scale = 1e-9, suffix = ""),
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 8, hjust = 0, lineheight = 0.1),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  theme_bw()

# Supress redundancies

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  plot.title = element_text(hjust = 0.5)
)


p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p3 <- p3 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean
p8 <- p8 + theme_map_clean
p9 <- p9 + theme_map_clean
p10 <- p10 + theme_map_clean
p11 <- p11 + theme_map_clean
p12 <- p12 + theme_map_clean
p13 <- p13 + theme_map_clean
p14 <- p14 + theme_map_clean
p15 <- p15 + theme_map_clean
p16 <- p16 + theme_map_clean
p17 <- p17 + theme_map_clean
p18 <- p18 + theme_map_clean
p19 <- p19 + theme_map_clean

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")
# p8 <- p8 + theme(legend.position = "none")
p9 <- p9 + theme(legend.position = "none")
p10 <- p10 + theme(legend.position = "none")
p11 <- p11 + theme(legend.position = "none")
p12 <- p12 + theme(legend.position = "none")


# Final figure
final_fig <-
  (p1 | p2 | p3 | p4) /
  (p5 | p6 | p7 | p8) /
  (p9 | p10 | p11 | p12)


final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/fire_occ_month_figure.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

# Continuous fire metrics

p18 <- p18 + labs(tag = "A")
p14 <- p14 + labs(tag = "B")
p17 <- p17 + labs(tag = "C")
p19 <- p19 + labs(tag = "D")
p15 <- p15 + labs(tag = "E")
p13 <- p13 + labs(tag = "F")
p16 <- p16 + labs(tag = "G")

final_fig <-
  (p18 | p14 | p17 | p19) /
  (p15 | p13 | p16 | plot_spacer())

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/fire_metrics_figure.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

# Random forest to test fire regime predictability ------------------------

# Boruta
rf_data <- as.data.frame(pred_df[, -c(1, 14, 15)]) %>%
  mutate(class = as.factor(class))

set.seed(42)
boruta_fire <- Boruta(
  class ~ .,
  data = rf_data,
  doTrace = 2
)

# 3. Fix "Tentative" decisions
# Sometimes Boruta is unsure about a variable. This forces a decision.
final_boruta <- TentativeRoughFix(boruta_fire)

print(final_boruta)

# Extract importance history
imp_history <- attStats(final_boruta)
imp_history$Variable <- rownames(imp_history)

# Filter only Confirmed variables
confirmed_vars <- imp_history %>%
  filter(decision == "Confirmed") %>%
  arrange(desc(meanImp))

# Plot
pdf("Figs/boruta_fire_classes.pdf", paper = "a4r", width = 0, height = 0)
ggplot(confirmed_vars, aes(x = reorder(Variable, meanImp), y = meanImp)) +
  geom_col(fill = "forestgreen") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Confirmed Drivers of Fire Regimes (Boruta Analysis)",
    y = "Importance (Z-Score vs Shadow Max)",
    x = ""
  )
dev.off()

# Heatmaps--------------------------------------------------------------------
# 1. Calculate Mean Z-Score for each variable per Regime
# (Z-score helps compare variables with different units like Temperature vs Density)
profile_data <- rf_data %>%
  mutate(across(where(is.numeric), scale)) %>% # Standardize vars
  group_by(class) %>%
  summarise(across(where(is.numeric), mean)) %>%
  pivot_longer(-class, names_to = "Variable", values_to = "Mean_Z_Score")

correct_order <- c(
  # Months in order
  paste0("X", 1:12),
  # Landscape metrics
  "lsm_c_pd",
  "lsm_c_split",
  "lsm_c_circle_mn",
  "lsm_c_enn_sd",
  "lsm_c_core_mn",
  "lsm_c_cai_sd",
  "lsm_c_dcore_sd"
)

profile_data_ordered <- profile_data %>%
  mutate(Variable = factor(Variable, levels = correct_order))

# 2. Plot Heatmap
pdf("Figs/heatmap_fire_classes.pdf", paper = "a4r", width = 0, height = 0)
ggplot(
  profile_data_ordered,
  aes(x = class, y = Variable, fill = Mean_Z_Score)
) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme_minimal() +
  labs(
    title = "Fire Regime Fingerprints",
    subtitle = "Red = Above Average, Blue = Below Average",
    x = "Fire Regime",
    y = ""
  )
dev.off()

mclust_parameters <- as.data.frame(t(mod$parameters$mean)) %>%
  mutate(class = factor(1:n())) %>%
  pivot_longer(-class, names_to = "Variable", values_to = "Parameter")

correct_order <- c(
  # Months in order
  paste0(1:12),
  # Landscape metrics
  "lsm_c_pd",
  "lsm_c_split",
  "lsm_c_circle_mn",
  "lsm_c_enn_sd",
  "lsm_c_core_mn",
  "lsm_c_cai_sd",
  "lsm_c_dcore_sd"
)
mclust_parameters <- mclust_parameters %>%
  mutate(Variable = factor(Variable, levels = correct_order))

pdf(
  "Figs/heatmap_fire_mclust_parameter.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
ggplot(mclust_parameters, aes(x = class, y = Variable, fill = Parameter)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme_minimal() +
  labs(
    title = "Fire Regime Fingerprints",
    subtitle = "Red = Above Average, Blue = Below Average",
    x = "Fire Regime",
    y = ""
  )
dev.off()

# Land use for each fire regime -------------------------------------------

fire_occ_df <- readRDS("Data/fire_occ_df_modeling.rds")
names(fire_occ_df[, c(1:5, 17:34)])
land_df <- fire_occ_df[, c(1:5, 17:34)]
land_df <- lazy_dt(land_df)

land_m <- land_df %>%
  group_by(plot_id, x, y) %>%
  summarise(
    savanna = mean(savanna, na.rm = T),
    grassland = mean(grassland, na.rm = T),
    pasture = mean(pasture, na.rm = T),
    sugar_cane = mean(sugar_cane, na.rm = T),
    mosaic_uses = mean(mosaic_uses, na.rm = T),
    urban = mean(urban, na.rm = T),
    soybean = mean(soybean, na.rm = T),
    non_vegetated = mean(non_vegetated, na.rm = T),
    forest = mean(forest, na.rm = T),
    temporary_crop = mean(temporary_crop, na.rm = T),
    perennial_crop = mean(perennial_crop, na.rm = T),
    water = mean(water, na.rm = T),
    wet_herbaceous = mean(wet_herbaceous, na.rm = T),
    native_area = mean(native_area, na.rm = T),
    anthropic_area = mean(anthropic_area, na.rm = T),
    NDVI = mean(NDVI, na.rm = T),
    pop = mean(pop, na.rm = T)
  ) %>%
  as_tibble()

land_use_fire_regimes <- left_join(pred_df, land_m, by = c("plot_id", "x", "y"))
skim(land_use_fire_regimes)

profile_land_data <- land_use_fire_regimes %>%
  dplyr::select(
    class,
    savanna,
    grassland,
    pasture,
    sugar_cane,
    mosaic_uses,
    urban,
    soybean,
    non_vegetated,
    forest,
    temporary_crop,
    perennial_crop,
    water,
    wet_herbaceous,
    native_area,
    anthropic_area,
    NDVI,
    pop
  ) %>%
  mutate(class = as.factor(class)) %>%
  mutate(across(where(is_double), scale)) %>% # Standardize vars
  group_by(class) %>%
  summarise(across(where(is_double), mean, na.rm = T)) %>%
  pivot_longer(-class, names_to = "Variable", values_to = "Mean_Z_Score")

correct_order <- c(
  "NDVI",
  "pop",
  "native_area",
  "anthropic_area",
  "savanna",
  "grassland",
  "forest",
  "wet_herbaceous",
  "water",
  "pasture",
  "soybean",
  "sugar_cane",
  "mosaic_uses",
  "temporary_crop",
  "perennial_crop",
  "urban",
  "non_vegetated"
)

profile_land_data_ordered <- profile_land_data %>%
  mutate(Variable = factor(Variable, levels = correct_order))

# 2. Plot Heatmap
pdf("Figs/heatmap_fire_classes_land.pdf", paper = "a4r", width = 0, height = 0)
ggplot(
  profile_land_data_ordered,
  aes(x = class, y = Variable, fill = Mean_Z_Score)
) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme_minimal() +
  labs(
    title = "Fire Regime Fingerprints",
    subtitle = "Red = Above Average, Blue = Below Average",
    x = "Fire Regime",
    y = ""
  )
dev.off()

# Decadal cluster analysis ------------------------------------------------

fire_occ_decade_m <- fire_occ_df %>%
  group_by(plot_id, longitude, latitude, decade, month) %>%
  summarise(fire_occurrence = mean(fire_occurrence, na.rm = T)) %>%
  as_tibble()

skim(fire_occ_decade_m)

wide_fire_occ_1stdec_m <- fire_occ_decade_m %>%
  filter(decade == "1985-1994") %>%
  pivot_wider(
    id_cols = c(plot_id), # Columns to keep as identifiers
    names_from = month, # Column containing the new column names (jan, feb...)
    values_from = fire_occurrence # Column containing the values
  )

# Fire metrics
# Define your breaks and labels first (just like in your example)
breaks <- c(1984, 1994, 2004, 2014, 2024)
labels <- c("1985-1994", "1995-2004", "2005-2014", "2015-2024")

# Update the lazy_dt object
fire_metrics_df <- fire_metrics_df %>%
  mutate(
    year_num = as.integer(year),
    decade = cut(year_num, breaks = breaks, labels = labels)
  )

fire_metrics_decade_m <- fire_metrics_df %>%
  group_by(plot_id, x, y, decade) %>%
  summarise(
    lsm_c_cai_sd = mean(lsm_c_cai_sd, na.rm = T),
    lsm_c_circle_mn = mean(lsm_c_circle_mn, na.rm = T),
    lsm_c_core_mn = mean(lsm_c_core_mn, na.rm = T),
    lsm_c_dcore_sd = mean(lsm_c_dcore_sd, na.rm = T),
    lsm_c_enn_sd = mean(lsm_c_enn_sd, na.rm = T),
    lsm_c_pd = mean(lsm_c_pd, na.rm = T),
    lsm_c_split = mean(lsm_c_split, na.rm = T)
  ) %>%
  as_tibble()

skim(fire_metrics_decade_m)


## 1st decade --------------------------------------------------------------

fire_metrics_1stdec_m <- fire_metrics_decade_m %>%
  filter(decade == "1985-1994")

fire_1stdec_m <- left_join(
  wide_fire_occ_1stdec_m,
  fire_metrics_1stdec_m,
  by = c("plot_id")
)
skim(fire_1stdec_m)

fire_1stdec_m <- na.omit(fire_1stdec_m)
dim(fire_1stdec_m)

fire_1stdec_m <- fire_1stdec_m[, -c(1, 14, 15, 16)]

fire_m_raw <- left_join(wide_fire_occ_m, fire_metrics_m, by = c("plot_id"))
skim(fire_m_raw)

fire_m_raw <- na.omit(fire_m_raw)
dim(fire_m_raw)

fire_m_raw <- fire_m_raw[, -c(1, 14, 15)]

# 1. Recover scaling parameters from the GLOBAL data
# (Assuming 'fire_m_raw' is your unscaled matrix from the global analysis)
global_means <- colMeans(fire_m_raw)
global_sds <- apply(fire_m_raw, 2, sd)

# We manually subtract global mean and divide by global SD
fire_1stdec_m <- t((t(fire_1stdec_m) - global_means) / global_sds)

# 'mod' is your fitted Mclust object from the global analysis
fire_1stdec_pred <- predict(mod, newdata = fire_1stdec_m)

fire_1stdec_pred_df <- left_join(
  wide_fire_occ_1stdec_m,
  fire_metrics_1stdec_m,
  by = c("plot_id")
)

fire_1stdec_pred_df <- na.omit(fire_1stdec_pred_df)

fire_1stdec_pred_df <- data.frame(
  fire_1stdec_pred_df,
  class = fire_1stdec_pred$classification
)
skim(fire_1stdec_pred_df)

pdf("Figs/mclust_map_class_1stdec.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_1stdec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  scale_fill_viridis(option = "turbo")
dev.off()


ModPred$z

pdf("Figs/mclust_density_class_1stdec.pdf", paper = "a4", width = 0, height = 0)
plot(density(fire_1stdec_pred$z[, 1]))
plot(density(fire_1stdec_pred$z[, 2]))
plot(density(fire_1stdec_pred$z[, 3]))
plot(density(fire_1stdec_pred$z[, 4]))
plot(density(fire_1stdec_pred$z[, 5]))
dev.off()

# Convert the matrix to a dataframe and add IDs
z_df <- as.data.frame(fire_1stdec_pred$z)
colnames(z_df) <- paste0("Regime_", 1:5) # Rename columns nicely
z_df$plot_id <- fire_1stdec_pred_df$plot_id # Assuming column 1 is plot_id

# 2. Assign the "Hard" Class and Uncertainty for sorting
z_df$Assigned_Class <- factor(fire_1stdec_pred$classification)
z_df$Max_Prob <- apply(fire_1stdec_pred$z, 1, max)

# 3. Pivot to Long Format for ggplot
z_long <- z_df %>%
  pivot_longer(
    cols = starts_with("Regime"),
    names_to = "Regime_Prob",
    values_to = "Probability"
  )

# 4. Create the Ordering
# We sort pixels by: 1. Their assigned class, 2. Their confidence (Max_Prob)
# This puts the "purest" pixels at the center of each block and uncertain ones at edges
plot_order <- z_df %>%
  arrange(Assigned_Class, desc(Max_Prob)) %>%
  pull(plot_id)

z_long$plot_id <- factor(z_long$plot_id, levels = plot_order)

# 5. The Plot
pdf(
  "Figs/mclust_uncertainty_class_1stdec.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)

ggplot(z_long, aes(x = plot_id, y = Probability, fill = Regime_Prob)) +
  geom_col(width = 1, position = "stack") + # No gaps between bars
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title = "Fire Regime Classification Structure",
    subtitle = "Each vertical slice is one pixel. Solid blocks = High certainty. Blends = Uncertainty.",
    x = "Pixels (Ordered by Regime)",
    y = "Probability of Membership"
  ) +
  theme_void() + # Remove axes clutter for 18k points
  theme(legend.position = "bottom")
dev.off()

# 1. Calculate Uncertainty
# mclust defines uncertainty as 1 - probability of the most likely group
uncertainty_df <- data.frame(
  x = fire_1stdec_pred_df[, "x"], # Assuming you have x/y in your matrix or joined object
  y = fire_1stdec_pred_df[, "y"],
  Uncertainty = 1 - apply(fire_1stdec_pred$z, 1, max)
)

# 2. Plot Map
pdf(
  "Figs/mclust_uncertainty_map_class_1stdec.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
ggplot(uncertainty_df, aes(x = x, y = y, fill = Uncertainty)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno") + # Dark = Certain, Bright = Uncertain
  coord_fixed() +
  theme_minimal() +
  labs(title = "Map of Classification Uncertainty")
dev.off()

## 2nd decade --------------------------------------------------------------
wide_fire_occ_2nddec_m <- fire_occ_decade_m %>%
  filter(decade == "1995-2004") %>%
  pivot_wider(
    id_cols = c(plot_id), # Columns to keep as identifiers
    names_from = month, # Column containing the new column names (jan, feb...)
    values_from = fire_occurrence # Column containing the values
  )

fire_metrics_2nddec_m <- fire_metrics_decade_m %>%
  filter(decade == "1995-2004")

fire_2nddec_m <- left_join(
  wide_fire_occ_2nddec_m,
  fire_metrics_2nddec_m,
  by = c("plot_id")
)
skim(fire_2nddec_m)

fire_2nddec_m <- na.omit(fire_2nddec_m)
dim(fire_2nddec_m)

fire_2nddec_m <- fire_2nddec_m[, -c(1, 14, 15, 16)]

fire_m_raw <- left_join(wide_fire_occ_m, fire_metrics_m, by = c("plot_id"))
skim(fire_m_raw)

fire_m_raw <- na.omit(fire_m_raw)
dim(fire_m_raw)

fire_m_raw <- fire_m_raw[, -c(1, 14, 15)]

# 1. Recover scaling parameters from the GLOBAL data
# (Assuming 'fire_m_raw' is your unscaled matrix from the global analysis)
global_means <- colMeans(fire_m_raw)
global_sds <- apply(fire_m_raw, 2, sd)

# We manually subtract global mean and divide by global SD
fire_2nddec_m <- t((t(fire_2nddec_m) - global_means) / global_sds)

# 'mod' is your fitted Mclust object from the global analysis
fire_2nddec_pred <- predict(mod, newdata = fire_2nddec_m)

fire_2nddec_pred_df <- left_join(
  wide_fire_occ_2nddec_m,
  fire_metrics_2nddec_m,
  by = c("plot_id")
)

fire_2nddec_pred_df <- na.omit(fire_2nddec_pred_df)

fire_2nddec_pred_df <- data.frame(
  fire_2nddec_pred_df,
  class = fire_2nddec_pred$classification
)
skim(fire_2nddec_pred_df)

pdf("Figs/mclust_map_class_2nddec.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_2nddec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  scale_fill_viridis(option = "turbo")
dev.off()


# Convert the matrix to a dataframe and add IDs
z_df <- as.data.frame(fire_2nddec_pred$z)
colnames(z_df) <- paste0("Regime_", 1:5) # Rename columns nicely
z_df$plot_id <- fire_2nddec_pred_df$plot_id # Assuming column 1 is plot_id

# 2. Assign the "Hard" Class and Uncertainty for sorting
z_df$Assigned_Class <- factor(fire_2nddec_pred$classification)
z_df$Max_Prob <- apply(fire_2nddec_pred$z, 1, max)

# 3. Pivot to Long Format for ggplot
z_long <- z_df %>%
  pivot_longer(
    cols = starts_with("Regime"),
    names_to = "Regime_Prob",
    values_to = "Probability"
  )

# 4. Create the Ordering
# We sort pixels by: 1. Their assigned class, 2. Their confidence (Max_Prob)
# This puts the "purest" pixels at the center of each block and uncertain ones at edges
plot_order <- z_df %>%
  arrange(Assigned_Class, desc(Max_Prob)) %>%
  pull(plot_id)

z_long$plot_id <- factor(z_long$plot_id, levels = plot_order)

# 5. The Plot
pdf(
  "Figs/mclust_uncertainty_class_2nddec.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)

ggplot(z_long, aes(x = plot_id, y = Probability, fill = Regime_Prob)) +
  geom_col(width = 1, position = "stack") + # No gaps between bars
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title = "Fire Regime Classification Structure",
    subtitle = "Each vertical slice is one pixel. Solid blocks = High certainty. Blends = Uncertainty.",
    x = "Pixels (Ordered by Regime)",
    y = "Probability of Membership"
  ) +
  theme_void() + # Remove axes clutter for 18k points
  theme(legend.position = "bottom")
dev.off()

# 1. Calculate Uncertainty
# mclust defines uncertainty as 1 - probability of the most likely group
uncertainty_df <- data.frame(
  x = fire_2nddec_pred_df[, "x"], # Assuming you have x/y in your matrix or joined object
  y = fire_2nddec_pred_df[, "y"],
  Uncertainty = 1 - apply(fire_2nddec_pred$z, 1, max)
)

# 2. Plot Map
pdf(
  "Figs/mclust_uncertainty_map_class_2nddec.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
ggplot(uncertainty_df, aes(x = x, y = y, fill = Uncertainty)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno") + # Dark = Certain, Bright = Uncertain
  coord_fixed() +
  theme_minimal() +
  labs(title = "Map of Classification Uncertainty")
dev.off()

## 3rd decade --------------------------------------------------------------
wide_fire_occ_3rddec_m <- fire_occ_decade_m %>%
  filter(decade == "2005-2014") %>%
  pivot_wider(
    id_cols = c(plot_id), # Columns to keep as identifiers
    names_from = month, # Column containing the new column names (jan, feb...)
    values_from = fire_occurrence # Column containing the values
  )

fire_metrics_3rddec_m <- fire_metrics_decade_m %>%
  filter(decade == "2005-2014")

fire_3rddec_m <- left_join(
  wide_fire_occ_3rddec_m,
  fire_metrics_3rddec_m,
  by = c("plot_id")
)
skim(fire_3rddec_m)

fire_3rddec_m <- na.omit(fire_3rddec_m)
dim(fire_3rddec_m)

fire_3rddec_m <- fire_3rddec_m[, -c(1, 14, 15, 16)]

fire_m_raw <- left_join(wide_fire_occ_m, fire_metrics_m, by = c("plot_id"))
skim(fire_m_raw)

fire_m_raw <- na.omit(fire_m_raw)
dim(fire_m_raw)

fire_m_raw <- fire_m_raw[, -c(1, 14, 15)]

# 1. Recover scaling parameters from the GLOBAL data
# (Assuming 'fire_m_raw' is your unscaled matrix from the global analysis)
global_means <- colMeans(fire_m_raw)
global_sds <- apply(fire_m_raw, 2, sd)

# We manually subtract global mean and divide by global SD
fire_3rddec_m <- t((t(fire_3rddec_m) - global_means) / global_sds)

# 'mod' is your fitted Mclust object from the global analysis
fire_3rddec_pred <- predict(mod, newdata = fire_3rddec_m)

fire_3rddec_pred_df <- left_join(
  wide_fire_occ_3rddec_m,
  fire_metrics_3rddec_m,
  by = c("plot_id")
)

fire_3rddec_pred_df <- na.omit(fire_3rddec_pred_df)

fire_3rddec_pred_df <- data.frame(
  fire_3rddec_pred_df,
  class = fire_3rddec_pred$classification
)
skim(fire_3rddec_pred_df)

pdf("Figs/mclust_map_class_3rddec.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_3rddec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  scale_fill_viridis(option = "turbo")
dev.off()


# Convert the matrix to a dataframe and add IDs
z_df <- as.data.frame(fire_3rddec_pred$z)
colnames(z_df) <- paste0("Regime_", 1:5) # Rename columns nicely
z_df$plot_id <- fire_3rddec_pred_df$plot_id # Assuming column 1 is plot_id

# 2. Assign the "Hard" Class and Uncertainty for sorting
z_df$Assigned_Class <- factor(fire_3rddec_pred$classification)
z_df$Max_Prob <- apply(fire_3rddec_pred$z, 1, max)

# 3. Pivot to Long Format for ggplot
z_long <- z_df %>%
  pivot_longer(
    cols = starts_with("Regime"),
    names_to = "Regime_Prob",
    values_to = "Probability"
  )

# 4. Create the Ordering
# We sort pixels by: 1. Their assigned class, 2. Their confidence (Max_Prob)
# This puts the "purest" pixels at the center of each block and uncertain ones at edges
plot_order <- z_df %>%
  arrange(Assigned_Class, desc(Max_Prob)) %>%
  pull(plot_id)

z_long$plot_id <- factor(z_long$plot_id, levels = plot_order)

# 5. The Plot
pdf(
  "Figs/mclust_uncertainty_class_3rddec.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)

ggplot(z_long, aes(x = plot_id, y = Probability, fill = Regime_Prob)) +
  geom_col(width = 1, position = "stack") + # No gaps between bars
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title = "Fire Regime Classification Structure",
    subtitle = "Each vertical slice is one pixel. Solid blocks = High certainty. Blends = Uncertainty.",
    x = "Pixels (Ordered by Regime)",
    y = "Probability of Membership"
  ) +
  theme_void() + # Remove axes clutter for 18k points
  theme(legend.position = "bottom")
dev.off()

# 1. Calculate Uncertainty
# mclust defines uncertainty as 1 - probability of the most likely group
uncertainty_df <- data.frame(
  x = fire_3rddec_pred_df[, "x"], # Assuming you have x/y in your matrix or joined object
  y = fire_3rddec_pred_df[, "y"],
  Uncertainty = 1 - apply(fire_3rddec_pred$z, 1, max)
)

# 2. Plot Map
pdf(
  "Figs/mclust_uncertainty_map_class_3rddec.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
ggplot(uncertainty_df, aes(x = x, y = y, fill = Uncertainty)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno") + # Dark = Certain, Bright = Uncertain
  coord_fixed() +
  theme_minimal() +
  labs(title = "Map of Classification Uncertainty")
dev.off()

## 4th decade --------------------------------------------------------------
wide_fire_occ_4thdec_m <- fire_occ_decade_m %>%
  filter(decade == "2015-2024") %>%
  pivot_wider(
    id_cols = c(plot_id), # Columns to keep as identifiers
    names_from = month, # Column containing the new column names (jan, feb...)
    values_from = fire_occurrence # Column containing the values
  )

fire_metrics_4thdec_m <- fire_metrics_decade_m %>%
  filter(decade == "2015-2024")

fire_4thdec_m <- left_join(
  wide_fire_occ_4thdec_m,
  fire_metrics_4thdec_m,
  by = c("plot_id")
)
skim(fire_4thdec_m)

fire_4thdec_m <- na.omit(fire_4thdec_m)
dim(fire_4thdec_m)

fire_4thdec_m <- fire_4thdec_m[, -c(1, 14, 15, 16)]

fire_m_raw <- left_join(wide_fire_occ_m, fire_metrics_m, by = c("plot_id"))
skim(fire_m_raw)

fire_m_raw <- na.omit(fire_m_raw)
dim(fire_m_raw)

fire_m_raw <- fire_m_raw[, -c(1, 14, 15)]

# 1. Recover scaling parameters from the GLOBAL data
# (Assuming 'fire_m_raw' is your unscaled matrix from the global analysis)
global_means <- colMeans(fire_m_raw)
global_sds <- apply(fire_m_raw, 2, sd)

# We manually subtract global mean and divide by global SD
fire_4thdec_m <- t((t(fire_4thdec_m) - global_means) / global_sds)

# 'mod' is your fitted Mclust object from the global analysis
fire_4thdec_pred <- predict(mod, newdata = fire_4thdec_m)

fire_4thdec_pred_df <- left_join(
  wide_fire_occ_4thdec_m,
  fire_metrics_4thdec_m,
  by = c("plot_id")
)

fire_4thdec_pred_df <- na.omit(fire_4thdec_pred_df)

fire_4thdec_pred_df <- data.frame(
  fire_4thdec_pred_df,
  class = fire_4thdec_pred$classification
)
skim(fire_4thdec_pred_df)

pdf("Figs/mclust_map_class_4thdec.pdf", paper = "a4", width = 0, height = 0)
ggplot(fire_4thdec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  scale_fill_viridis(option = "turbo")
dev.off()


# Convert the matrix to a dataframe and add IDs
z_df <- as.data.frame(fire_4thdec_pred$z)
colnames(z_df) <- paste0("Regime_", 1:5) # Rename columns nicely
z_df$plot_id <- fire_4thdec_pred_df$plot_id # Assuming column 1 is plot_id

# 2. Assign the "Hard" Class and Uncertainty for sorting
z_df$Assigned_Class <- factor(fire_4thdec_pred$classification)
z_df$Max_Prob <- apply(fire_4thdec_pred$z, 1, max)

# 3. Pivot to Long Format for ggplot
z_long <- z_df %>%
  pivot_longer(
    cols = starts_with("Regime"),
    names_to = "Regime_Prob",
    values_to = "Probability"
  )

# 4. Create the Ordering
# We sort pixels by: 1. Their assigned class, 2. Their confidence (Max_Prob)
# This puts the "purest" pixels at the center of each block and uncertain ones at edges
plot_order <- z_df %>%
  arrange(Assigned_Class, desc(Max_Prob)) %>%
  pull(plot_id)

z_long$plot_id <- factor(z_long$plot_id, levels = plot_order)

# 5. The Plot
pdf(
  "Figs/mclust_uncertainty_class_4thdec.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)

ggplot(z_long, aes(x = plot_id, y = Probability, fill = Regime_Prob)) +
  geom_col(width = 1, position = "stack") + # No gaps between bars
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title = "Fire Regime Classification Structure",
    subtitle = "Each vertical slice is one pixel. Solid blocks = High certainty. Blends = Uncertainty.",
    x = "Pixels (Ordered by Regime)",
    y = "Probability of Membership"
  ) +
  theme_void() + # Remove axes clutter for 18k points
  theme(legend.position = "bottom")
dev.off()

# 1. Calculate Uncertainty
# mclust defines uncertainty as 1 - probability of the most likely group
uncertainty_df <- data.frame(
  x = fire_4thdec_pred_df[, "x"], # Assuming you have x/y in your matrix or joined object
  y = fire_4thdec_pred_df[, "y"],
  Uncertainty = 1 - apply(fire_4thdec_pred$z, 1, max)
)

# 2. Plot Map
pdf(
  "Figs/mclust_uncertainty_map_class_4thdec.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
ggplot(uncertainty_df, aes(x = x, y = y, fill = Uncertainty)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno") + # Dark = Certain, Bright = Uncertain
  coord_fixed() +
  theme_minimal() +
  labs(title = "Map of Classification Uncertainty")
dev.off()


# Fire regimes transitions ------------------------------------------------
fire_regimes_decade <- rbind(
  fire_1stdec_pred_df,
  fire_2nddec_pred_df,
  fire_3rddec_pred_df,
  fire_4thdec_pred_df
)


fire_regimes_decade$Uncertainty <- c(
  1 - apply(fire_1stdec_pred$z, 1, max),
  1 - apply(fire_2nddec_pred$z, 1, max),
  1 - apply(fire_3rddec_pred$z, 1, max),
  1 - apply(fire_4thdec_pred$z, 1, max)
)

# 1. Pivot Wider to get a "History" for each pixel
regime_trajectories <- fire_regimes_decade %>%
  filter(Uncertainty < 0.1) %>%
  dplyr::select(plot_id, x, y, decade, class) %>%
  # Use names_prefix to make columns standard (e.g., "D_1", "D_2"...)
  # Assuming decades are ordered factors/numeric
  pivot_wider(
    names_from = decade,
    values_from = class,
    names_prefix = "Decade_"
  ) %>%
  # Filter out pixels that might be missing data for some decades
  drop_na()


## Fire regime stability ---------------------------------------------------

# 2. Calculate Stability
stability_df <- regime_trajectories %>%
  rowwise() %>% # Compute row-by-row
  mutate(
    # Count unique regimes across the 4 decades columns
    n_regimes = n_distinct(c_across(starts_with("Decade_"))),

    # Define Stability Class
    Stability_Status = case_when(
      n_regimes == 1 ~ "Stable",
      n_regimes == 2 ~ "Shifted Once/Oscillated",
      n_regimes > 2 ~ "Volatile"
    ),

    # Identify WHICH regime was the stable one (useful for mapping)
    Stable_Regime = ifelse(
      n_regimes == 1,
      as.character(`Decade_1985-1994`),
      "Unstable"
    )
  ) %>%
  ungroup()

# 3. Calculate % of Landscape that is Stable
stability_stats <- stability_df %>%
  dplyr::count(Stability_Status) %>%
  mutate(Percentage = n / sum(n) * 100)

print(stability_stats)

# # A tibble: 3 × 3
# Stability_Status            n Percentage
# <chr>                   <int>      <dbl>
#   1 Shifted Once/Oscillated  4164       45.9
# 2 Stable                   3483       38.4
# 3 Volatile                 1424       15.7

# 4. Aggregate data for the plot (Grouping by flow path makes plotting 100x faster)
sankey_data <- regime_trajectories %>%
  group_by(
    `Decade_1985-1994`,
    `Decade_1995-2004`,
    `Decade_2005-2014`,
    `Decade_2015-2024`
  ) %>%
  dplyr::count() # This weights the flows by number of pixels

# 5. Plot
pdf(
  "Figs/fire_regime_transitions_sankey.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
ggplot(
  sankey_data,
  aes(
    axis1 = `Decade_1985-1994`,
    axis2 = `Decade_1995-2004`,
    axis3 = `Decade_2005-2014`,
    axis4 = `Decade_2015-2024`,
    y = n
  )
) +
  geom_alluvium(aes(fill = factor(`Decade_1985-1994`)), width = 1 / 12) + # Flow colored by start point
  geom_stratum(width = 1 / 12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("1985-1994", "1995-2004", "2005-2014", "2015-2024")
  ) +
  scale_fill_viridis_d(option = "turbo", name = "Start Regime") +
  theme_minimal() +
  labs(title = "Fire Regime Transitions (1985-2024)", y = "Number of Cells")
dev.off()


head(stability_df)
glimpse(stability_df)
unique(stability_df$Stability_Status)

# 7. Plot the Map
pdf("Figs/fire_regime_stability_map.pdf", paper = "a4", width = 0, height = 0)
ggplot(stability_df, aes(x = x, y = y, fill = Stability_Status)) +
  geom_tile() +
  scale_fill_manual(values = viridis(3), name = "Regime stability") +
  coord_fixed() +
  theme_bw() +
  labs(x = "", y = "")
dev.off()

unique(stability_df$Stable_Regime)

stability_df$Stable_Regime <- factor(
  stability_df$Stable_Regime,
  levels = c("1", "2", "3", "4", "5", "Unstable")
)

# 7. Plot the Map
pdf("Figs/fire_regime_stable_map.pdf", paper = "a4", width = 0, height = 0)
ggplot(stability_df, aes(x = x, y = y, fill = Stable_Regime)) +
  geom_tile() +
  scale_fill_manual(
    values = c(turbo(5), viridis(2)[2]),
    name = "Stable regime"
  ) +
  coord_fixed() +
  theme_bw() +
  labs(x = "", y = "")
dev.off()

## Formal tests ------------------------------------------------------------

### Fire regimes differ across decades? -------------------------------------

# Contingency class × decade

tab_class_decade <- table(
  fire_regimes_decade$class,
  fire_regimes_decade$decade
)

#
prop.table(tab_class_decade, margin = 2)

# Observed G statistic
G_obs <- GTest(tab_class_decade)
G_obs$statistic
# 1138.249

# Permutation test
set.seed(123)

n_perm <- 999
G_perm <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  perm_decade <- sample(fire_regimes_decade$decade)
  tab_perm <- table(fire_regimes_decade$class, perm_decade)
  G_perm[i] <- GTest(tab_perm)$statistic
}

mean(G_perm)
# 11.98678
p_perm <- mean(G_perm >= G_obs$statistic)
print(p_perm)
# 0

pdf(
  "Figs/mosaicplot_dec_fire_regimes.pdf",
  paper = "a4r",
  width = 0,
  height = 0
)
mosaicplot(t(tab_class_decade), shade = T, type = "pearson", main = "")
dev.off()

vcd::mosaic(tab_class_decade, shade = T, main = "")

### Transitions are non-random? ---------------------------------------------

transition_matrix <- function(df, d1, d2) {
  tab <- df %>%
    filter(decade %in% c(d1, d2)) %>%
    dplyr::select(plot_id, decade, class) %>%
    pivot_wider(names_from = decade, values_from = class) %>%
    drop_na()

  table(tab[[d1]], tab[[d2]])
}


#### 1st - 2nd decades -------------------------------------------------------

T_1_2 <- transition_matrix(
  fire_regimes_decade,
  "1985-1994",
  "1995-2004"
)

print(T_1_2)

g_test_1_2 <- GTest(T_1_2)
g_test_1_2

set.seed(123)

n_perm <- 999
g_perm <- numeric(n_perm)

# Destination probabilities under null
p_to <- colSums(T_1_2) / sum(T_1_2)

for (i in seq_len(n_perm)) {
  T_perm <- matrix(
    0,
    nrow = nrow(T_1_2),
    ncol = ncol(T_1_2),
    dimnames = dimnames(T_1_2)
  )

  for (r in seq_len(nrow(T_1_2))) {
    n_r <- sum(T_1_2[r, ]) # row total (fixed)
    T_perm[r, ] <- rmultinom(1, size = n_r, prob = p_to)
  }

  g_perm[i] <- suppressWarnings(GTest(T_perm)$statistic)
}

g_obs <- GTest(T_1_2)$statistic
print(g_obs)
# 12099.77

mean(g_perm)
# 15.97476

p_perm_g <- mean(g_perm >= g_obs)
p_perm_g
# 0

# Expected transition matrix under the null
E <- outer(
  rowSums(T_1_2),
  colSums(T_1_2) / sum(T_1_2)
)

# Standardized residuals (this is the key post-hoc tool)
std_resid <- (T_1_2 - E) / sqrt(E)

which(abs(std_resid) > 2, arr.ind = TRUE)

rownames(std_resid)
colnames(std_resid)

posthoc_transitions <- as.data.frame(as.table(std_resid)) |>
  rename(from = Var1, to = Var2, z = Freq) |>
  arrange(desc(abs(z)))

print(posthoc_transitions)

write.csv(
  posthoc_transitions,
  "Output/results/posthoc_dec_1_2_transitions.csv",
  row.names = FALSE
)

#### 2nd - 3rd decades -------------------------------------------------------

T_2_3 <- transition_matrix(
  fire_regimes_decade,
  "1995-2004",
  "2005-2014"
)


g_test_2_3 <- GTest(T_2_3)
g_test_2_3

set.seed(123)

n_perm <- 999
g_perm <- numeric(n_perm)

# Destination probabilities under null
p_to <- colSums(T_2_3) / sum(T_2_3)

for (i in seq_len(n_perm)) {
  T_perm <- matrix(
    0,
    nrow = nrow(T_2_3),
    ncol = ncol(T_2_3),
    dimnames = dimnames(T_2_3)
  )

  for (r in seq_len(nrow(T_2_3))) {
    n_r <- sum(T_2_3[r, ]) # row total (fixed)
    T_perm[r, ] <- rmultinom(1, size = n_r, prob = p_to)
  }

  g_perm[i] <- suppressWarnings(GTest(T_perm)$statistic)
}

g_obs <- GTest(T_2_3)$statistic
print(g_obs)
# 12653.12

mean(g_perm)
# 15.90031

p_perm_g <- mean(g_perm >= g_obs)
p_perm_g
# 0

# Expected transition matrix under the null
E <- outer(
  rowSums(T_2_3),
  colSums(T_2_3) / sum(T_2_3)
)

# Standardized residuals (this is the key post-hoc tool)
std_resid <- (T_2_3 - E) / sqrt(E)

which(abs(std_resid) > 2, arr.ind = TRUE)

rownames(std_resid)
colnames(std_resid)

posthoc_transitions <- as.data.frame(as.table(std_resid)) |>
  rename(from = Var1, to = Var2, z = Freq) |>
  arrange(desc(abs(z)))

print(posthoc_transitions)

write.csv(
  posthoc_transitions,
  "Output/results/posthoc_dec_2_3_transitions.csv",
  row.names = FALSE
)


#### 3rd - 4th decades -------------------------------------------------------

T_3_4 <- transition_matrix(
  fire_regimes_decade,
  "2005-2014",
  "2015-2024"
)


g_test_3_4 <- GTest(T_3_4)
g_test_3_4

set.seed(123)

n_perm <- 999
g_perm <- numeric(n_perm)

# Destination probabilities under null
p_to <- colSums(T_3_4) / sum(T_3_4)

for (i in seq_len(n_perm)) {
  T_perm <- matrix(
    0,
    nrow = nrow(T_3_4),
    ncol = ncol(T_3_4),
    dimnames = dimnames(T_3_4)
  )

  for (r in seq_len(nrow(T_3_4))) {
    n_r <- sum(T_3_4[r, ]) # row total (fixed)
    T_perm[r, ] <- rmultinom(1, size = n_r, prob = p_to)
  }

  g_perm[i] <- suppressWarnings(GTest(T_perm)$statistic)
}

g_obs <- GTest(T_3_4)$statistic
print(g_obs)
# 12877.98

mean(g_perm)
# 15.8865

p_perm_g <- mean(g_perm >= g_obs)
p_perm_g
# 0

# Expected transition matrix under the null
E <- outer(
  rowSums(T_3_4),
  colSums(T_3_4) / sum(T_3_4)
)

# Standardized residuals (this is the key post-hoc tool)
std_resid <- (T_3_4 - E) / sqrt(E)

which(abs(std_resid) > 2, arr.ind = TRUE)

rownames(std_resid)
colnames(std_resid)

posthoc_transitions <- as.data.frame(as.table(std_resid)) |>
  rename(from = Var1, to = Var2, z = Freq) |>
  arrange(desc(abs(z)))

print(posthoc_transitions)

write.csv(
  posthoc_transitions,
  "Output/results/posthoc_dec_3_4_transitions.csv",
  row.names = FALSE
)


### Assoc plots

pdf("Figs/assocplot_by_decade.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
vcd::assoc(T_1_2, shade = TRUE, main = "1985–1994 → 1995–2004")
vcd::assoc(T_2_3, shade = TRUE, main = "1995–2004 → 2005–2014")
vcd::assoc(T_3_4, shade = TRUE, main = "2005–2014 → 2015–2024")
dev.off()


#### Transitions are symmetrical/asymmetrical?

##### 1st - 2nd decades -------------------------------------------------------

transitions <- read.csv("Output/results/posthoc_dec_1_2_transitions.csv")


symmetry_summary <- transitions %>% # full table, no |z| filter
  filter(from != to) %>%
  mutate(
    r1 = pmin(from, to),
    r2 = pmax(from, to),
    dir = ifelse(from == r1, "forward", "backward")
  ) %>%
  dplyr::select(r1, r2, dir, z) %>%
  pivot_wider(names_from = dir, values_from = z) %>%
  mutate(
    delta_Z = abs(forward - backward),
    symmetry = case_when(
      delta_Z < 1 ~ "Approximately symmetric",
      delta_Z < 3 ~ "Moderately asymmetric",
      TRUE ~ "Strongly asymmetric"
    )
  )


print(symmetry_summary)

write.csv(
  symmetry_summary,
  "Output/results/posthoc_dec_1_2_symmetrey.csv",
  row.names = FALSE
)

symmetry_summary %>%
  dplyr::count(symmetry)

#   symmetry                    n
#   <chr>                   <int>
# 1 Approximately symmetric     4
# 2 Moderately asymmetric       5
# 3 Strongly asymmetric         1

##### 2nd - 3rd decades -------------------------------------------------------

transitions <- read.csv("Output/results/posthoc_dec_2_3_transitions.csv")

symmetry_summary <- transitions %>% # full table, no |z| filter
  filter(from != to) %>%
  mutate(
    r1 = pmin(from, to),
    r2 = pmax(from, to),
    dir = ifelse(from == r1, "forward", "backward")
  ) %>%
  dplyr::select(r1, r2, dir, z) %>%
  pivot_wider(names_from = dir, values_from = z) %>%
  mutate(
    delta_Z = abs(forward - backward),
    symmetry = case_when(
      delta_Z < 1 ~ "Approximately symmetric",
      delta_Z < 3 ~ "Moderately asymmetric",
      TRUE ~ "Strongly asymmetric"
    )
  )


print(symmetry_summary)

write.csv(
  symmetry_summary,
  "Output/results/posthoc_dec_2_3_symmetrey.csv",
  row.names = FALSE
)

symmetry_summary %>%
  dplyr::count(symmetry)

# A tibble: 3 × 2
#   symmetry                    n
#   <chr>                   <int>
# 1 Approximately symmetric     2
# 2 Moderately asymmetric       2
# 3 Strongly asymmetric         6

##### 3rd-4th decades -------------------------------------------------------

transitions <- read.csv("Output/results/posthoc_dec_3_4_transitions.csv")

symmetry_summary <- transitions %>% # full table, no |z| filter
  filter(from != to) %>%
  mutate(
    r1 = pmin(from, to),
    r2 = pmax(from, to),
    dir = ifelse(from == r1, "forward", "backward")
  ) %>%
  dplyr::select(r1, r2, dir, z) %>%
  pivot_wider(names_from = dir, values_from = z) %>%
  mutate(
    delta_Z = abs(forward - backward),
    symmetry = case_when(
      delta_Z < 1 ~ "Approximately symmetric",
      delta_Z < 3 ~ "Moderately asymmetric",
      TRUE ~ "Strongly asymmetric"
    )
  )


print(symmetry_summary)

write.csv(
  symmetry_summary,
  "Output/results/posthoc_dec_3_4_symmetrey.csv",
  row.names = FALSE
)

symmetry_summary %>%
  dplyr::count(symmetry)

# A tibble: 2 × 2
#   symmetry                  n
#   <chr>                 <int>
# 1 Moderately asymmetric     4
# 2 Strongly asymmetric       6

### Persistence/Stability and directionality ------------------------------------------

#### 1st - 2nd decades -------------------------------------------------------

# Observed persistence
P_obs <- prop.table(T_1_2, 1) # normaliza por linha
diag(P_obs)

# Persistence under randomness
P_exp <- colSums(T_1_2) / sum(T_1_2)
P_exp

persistence_df <- data.frame(
  class = rownames(T_1_2),
  observed = diag(P_obs),
  expected = colSums(T_1_2) / sum(T_1_2)
) |>
  mutate(excess = observed - expected)

print(persistence_df)

# Permutation-based CI
persist_perm <- matrix(NA, nrow = n_perm, ncol = nrow(T_1_2))

for (i in seq_len(n_perm)) {
  T_perm <- matrix(0, nrow = nrow(T_1_2), ncol = ncol(T_1_2))
  for (r in seq_len(nrow(T_1_2))) {
    T_perm[r, ] <- rmultinom(1, size = sum(T_1_2[r, ]), prob = p_to)
  }
  persist_perm[i, ] <- diag(prop.table(T_perm, 1))
}

persist_ci <- apply(persist_perm, 2, quantile, probs = c(0.025, 0.975))
print(t(persist_ci))
persistence_df <- persistence_df %>%
  mutate(expected.lower = persist_ci[1, ], expected.upper = persist_ci[2, ])

print(persistence_df)

# class  observed   expected    excess expected.lower expected.upper
# 1     1 0.4775811 0.17564124 0.3019399     0.16253687     0.18849558
# 2     2 0.5226572 0.23073439 0.2919228     0.21778173     0.24341637
# 3     3 0.6096572 0.13006059 0.4795966     0.11758201     0.14229635
# 4     4 0.4617886 0.04048468 0.4213039     0.02601626     0.05691057
# 5     5 0.7175981 0.42307910 0.2945190     0.41170482     0.43476516

#### 2nd - 3rd decades -------------------------------------------------------

# Observed persistence
P_obs <- prop.table(T_2_3, 1) # normaliza por linha
diag(P_obs)

# Persistence under randomness
P_exp <- colSums(T_2_3) / sum(T_2_3)
P_exp

persistence_df <- data.frame(
  class = rownames(T_2_3),
  observed = diag(P_obs),
  expected = colSums(T_2_3) / sum(T_2_3)
) |>
  mutate(excess = observed - expected)

print(persistence_df)

# Permutation-based CI
persist_perm <- matrix(NA, nrow = n_perm, ncol = nrow(T_2_3))

for (i in seq_len(n_perm)) {
  T_perm <- matrix(0, nrow = nrow(T_2_3), ncol = ncol(T_2_3))
  for (r in seq_len(nrow(T_2_3))) {
    T_perm[r, ] <- rmultinom(1, size = sum(T_2_3[r, ]), prob = p_to)
  }
  persist_perm[i, ] <- diag(prop.table(T_perm, 1))
}

persist_ci <- apply(persist_perm, 2, quantile, probs = c(0.025, 0.975))
print(t(persist_ci))
persistence_df <- persistence_df %>%
  mutate(expected.lower = persist_ci[1, ], expected.upper = persist_ci[2, ])

print(persistence_df)

# class  observed   expected    excess expected.lower expected.upper
# 1     1 0.5447342 0.20732265 0.3374116     0.20666887      0.2360680
# 2     2 0.5182806 0.24691076 0.2713699     0.21294466      0.2401309
# 3     3 0.7370483 0.15869565 0.5783527     0.12276883      0.1497606
# 4     4 0.5525903 0.06824943 0.4843408     0.06907378      0.1145997
# 5     5 0.5798634 0.31882151 0.2610419     0.31422547      0.3356607

#### 3rd - 4th decades -------------------------------------------------------

# Observed persistence
P_obs <- prop.table(T_3_4, 1) # normaliza por linha
diag(P_obs)

# Persistence under randomness
P_exp <- colSums(T_3_4) / sum(T_3_4)
P_exp

persistence_df <- data.frame(
  class = rownames(T_3_4),
  observed = diag(P_obs),
  expected = colSums(T_3_4) / sum(T_3_4)
) |>
  mutate(excess = observed - expected)

print(persistence_df)

# Permutation-based CI
persist_perm <- matrix(NA, nrow = n_perm, ncol = nrow(T_3_4))

for (i in seq_len(n_perm)) {
  T_perm <- matrix(0, nrow = nrow(T_3_4), ncol = ncol(T_3_4))
  for (r in seq_len(nrow(T_3_4))) {
    T_perm[r, ] <- rmultinom(1, size = sum(T_3_4[r, ]), prob = p_to)
  }
  persist_perm[i, ] <- diag(prop.table(T_perm, 1))
}

persist_ci <- apply(persist_perm, 2, quantile, probs = c(0.025, 0.975))
print(t(persist_ci))
persistence_df <- persistence_df %>%
  mutate(expected.lower = persist_ci[1, ], expected.upper = persist_ci[2, ])

print(persistence_df)

# class  observed   expected    excess expected.lower expected.upper
# 1     1 0.5189946 0.22153487 0.2974597     0.20792631      0.2345044
# 2     2 0.4579615 0.22649548 0.2314660     0.21418506      0.2400188
# 3     3 0.6229272 0.13697111 0.4859561     0.12472963      0.1499640
# 4     4 0.5237173 0.09057485 0.4331425     0.07545983      0.1084221
# 5     5 0.6195008 0.32442369 0.2950771     0.31261447      0.3366942

## Final Figure ------------------------------------------------------------

p1 <- ggplot(stability_df, aes(x = x, y = y, fill = Stability_Status)) +
  geom_tile() +
  scale_fill_manual(values = viridis(3), name = "Regime stability") +
  coord_fixed() +
  theme_bw() +
  labs(x = "", y = "")

p2 <- ggplot(stability_df, aes(x = x, y = y, fill = Stable_Regime)) +
  geom_tile() +
  scale_fill_manual(
    values = c(turbo(5), viridis(2)[2]),
    name = "Stable regime"
  ) +
  coord_fixed() +
  theme_bw() +
  labs(x = "", y = "")

p3 <- ggplot(
  sankey_data,
  aes(
    axis1 = `Decade_1985-1994`,
    axis2 = `Decade_1995-2004`,
    axis3 = `Decade_2005-2014`,
    axis4 = `Decade_2015-2024`,
    y = n
  )
) +
  geom_alluvium(aes(fill = factor(`Decade_1985-1994`)), width = 1 / 12) + # Flow colored by start point
  geom_stratum(width = 1 / 12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("1985-1994", "1995-2004", "2005-2014", "2015-2024")
  ) +
  scale_fill_viridis_d(option = "turbo", name = "Start Regime") +
  theme_minimal() +
  labs(y = "Number of Cells")

p4 <- ggplot(fire_1stdec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis(option = "turbo") +
  labs(x = "", y = "", title = "1985-1994") +
  theme(plot.title = element_text(hjust = 0.5))

p5 <- ggplot(fire_2nddec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis(option = "turbo") +
  labs(x = "", y = "", title = "1995-2004") +
  theme(plot.title = element_text(hjust = 0.5))

p6 <- ggplot(fire_3rddec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis(option = "turbo") +
  labs(x = "", y = "", title = "2005-2014") +
  theme(plot.title = element_text(hjust = 0.5))

p7 <- ggplot(fire_4thdec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis(option = "turbo") +
  labs(x = "", y = "", title = "2015-2024") +
  theme(plot.title = element_text(hjust = 0.5))

# Supress redundancies

theme_map_clean <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank()
)


p1 <- p1 + theme_map_clean
p2 <- p2 + theme_map_clean
p4 <- p4 + theme_map_clean
p5 <- p5 + theme_map_clean
p6 <- p6 + theme_map_clean
p7 <- p7 + theme_map_clean

p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
p7 <- p7 + theme(legend.position = "none")

p3 <- p3 +
  theme(
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10)
  )

# Tag panels
p4 <- p4 + labs(tag = "A")
p5 <- p5 + labs(tag = "B")
p6 <- p6 + labs(tag = "C")
p7 <- p7 + labs(tag = "D")

p1 <- p1 + labs(tag = "E")
p2 <- p2 + labs(tag = "F")

p3 <- p3 + labs(tag = "G")


decade_row <- (p4 | p5 | p6 | p7)

# Final figure
final_fig <-
  decade_row /
  (p1 | p2) /
  p3


final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/fire_regime_dynamics_figure.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

# Land use and cover PCA -------------------------------------------------------

fire_occ_df <- readRDS(
  "Data/final_fire_occ_data.rds"
)
fire_occ_df <- lazy_dt(fire_occ_df)

fire_occ_m <- fire_occ_df %>%
  group_by(plot_id, longitude, latitude, month) %>%
  summarise(fire_occurrence = mean(fire_occurrence, na.rm = T)) %>%
  as_tibble()

skim(fire_occ_m)

wide_fire_occ_m <- fire_occ_m %>%
  pivot_wider(
    id_cols = c(plot_id), # Columns to keep as identifiers
    names_from = month, # Column containing the new column names (jan, feb...)
    values_from = fire_occurrence # Column containing the values
  )

skim(wide_fire_occ_m)


fire_metrics_df <- readRDS(
  "Data/final_clean_data.rds"
)
fire_metrics_df <- lazy_dt(fire_metrics_df)

# Standardize names and variables
fire_metrics_df <- dplyr::rename(fire_metrics_df, x = longitude)
fire_metrics_df <- dplyr::rename(fire_metrics_df, y = latitude)

fire_metrics_m <- fire_metrics_df %>%
  group_by(plot_id, x, y) %>%
  summarise(
    lsm_c_cai_sd = mean(lsm_c_cai_sd, na.rm = T),
    lsm_c_circle_mn = mean(lsm_c_circle_mn, na.rm = T),
    lsm_c_core_mn = mean(lsm_c_core_mn, na.rm = T),
    lsm_c_dcore_sd = mean(lsm_c_dcore_sd, na.rm = T),
    lsm_c_enn_sd = mean(lsm_c_enn_sd, na.rm = T),
    lsm_c_pd = mean(lsm_c_pd, na.rm = T),
    lsm_c_split = mean(lsm_c_split, na.rm = T)
  ) %>%
  as_tibble()

skim(fire_metrics_m)

mod <- readRDS("Output/Mclust_model.rds")
mod$data

ModPred <- predict.Mclust(mod) # prediction
fire_m <- left_join(wide_fire_occ_m, fire_metrics_m, by = c("plot_id"))
fire_m <- na.omit(fire_m)

pred_df <- data.frame(fire_m, class = ModPred$classification)

fire_occ_df <- readRDS("Data/fire_occ_df_modeling.rds")
names(fire_occ_df[, c(1:5, 17:34)])
land_df <- fire_occ_df[, c(1:5, 17:34)]
land_df <- lazy_dt(land_df)

land_m <- land_df %>%
  group_by(plot_id, x, y) %>%
  summarise(
    savanna = mean(savanna, na.rm = T),
    grassland = mean(grassland, na.rm = T),
    pasture = mean(pasture, na.rm = T),
    sugar_cane = mean(sugar_cane, na.rm = T),
    mosaic_uses = mean(mosaic_uses, na.rm = T),
    urban = mean(urban, na.rm = T),
    soybean = mean(soybean, na.rm = T),
    non_vegetated = mean(non_vegetated, na.rm = T),
    forest = mean(forest, na.rm = T),
    temporary_crop = mean(temporary_crop, na.rm = T),
    perennial_crop = mean(perennial_crop, na.rm = T),
    water = mean(water, na.rm = T),
    wet_herbaceous = mean(wet_herbaceous, na.rm = T),
    native_area = mean(native_area, na.rm = T),
    anthropic_area = mean(anthropic_area, na.rm = T),
    NDVI = mean(NDVI, na.rm = T),
    pop = mean(pop, na.rm = T)
  ) %>%
  as_tibble()


land_use_fire_regimes <- left_join(pred_df, land_m, by = c("plot_id", "x", "y"))
skim(land_use_fire_regimes)
land_use_fire_regimes <- na.omit(land_use_fire_regimes)

ggplot(land_use_fire_regimes, aes(x = x, y = y, z = savanna, fill = savanna)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")

# 2. Run PCA on land use data (Exclude ID and Class)
pca_land <- prcomp(land_use_fire_regimes[, c(24:40)], scale. = TRUE)

land_use_fire_regimes$class <- as.factor(land_use_fire_regimes$class)

# --- 1. PREPARE DATA ---

# A. Extract PCA Coordinates (The raw data)
# 'fortify' automatically joins the PC scores with your original data
pca_data <- fortify(pca_land, data = land_use_fire_regimes)

# B. Calculate Centroids (The "Mean" dots you want)
centroids <- pca_data %>%
  group_by(class) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

# C. Extract Loadings (The Arrows)
# We need to scale them so they are visible on the plot
loadings <- as.data.frame(pca_land$rotation)
loadings$Variable <- rownames(loadings)
scale_factor <- 5 # Tweak this if arrows are too short/long
loadings$PC1 <- loadings$PC1 * scale_factor
loadings$PC2 <- loadings$PC2 * scale_factor


# --- 2. BUILD THE PLOT ---
pdf("Figs/PCA_lulc_fire.pdf", paper = "a4r", width = 0, height = 0)

ggplot() +

  # LAYER 1: Ellipses (95% Confidence Regions)
  # We use the raw 'pca_data' here, but 'geom="polygon"' draws shapes, not points.
  stat_ellipse(
    data = pca_data,
    aes(x = PC1, y = PC2, fill = class, color = class),
    geom = "polygon",
    alpha = 0.2,
    level = 0.95,
    show.legend = TRUE
  ) +

  # LAYER 2: Centroids (The Big Dots)
  # We use the 'centroids' dataframe, so only 5 dots appear.
  geom_point(
    data = centroids,
    aes(x = PC1, y = PC2, fill = class),
    size = 6,
    shape = 21,
    color = "black",
    stroke = 1.5
  ) +

  # LAYER 3: Arrows (Loadings)
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "grey40"
  ) +

  # LAYER 4: Arrow Labels
  geom_text(
    data = loadings,
    aes(x = PC1, y = PC2, label = Variable),
    color = "black",
    size = 3.5,
    fontface = "bold",
    vjust = -0.5
  ) +

  # STYLING
  scale_fill_viridis_d(option = "turbo", name = "Fire Regime") +
  scale_color_viridis_d(option = "turbo", name = "Fire Regime") +
  theme_bw() +
  labs(
    title = "Land use and cover of Fire Regimes",
    subtitle = "Centroids (means) and 95% Confidence Ellipses",
    x = paste0(
      "PC1 (",
      round(summary(pca_land)$importance[2, 1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca_land)$importance[2, 2] * 100, 1),
      "%)"
    )
  )
dev.off()

# Climate and physical PCA -------------------------------------------------------

names(fire_occ_df[, c(1:5, 8:16, 35:38)])
env_df <- fire_occ_df[, c(1:5, 8:16, 35:38)]
env_df <- lazy_dt(env_df)

env_m <- env_df %>%
  group_by(plot_id, x, y) %>%
  summarise(
    dewpoint_temp_sd = mean(dewpoint_temp_sd, na.rm = T),
    pot_evap_sd = mean(pot_evap_sd, na.rm = T),
    precip_min = mean(precip_min, na.rm = T),
    solrad_mean = mean(solrad_mean, na.rm = T),
    total_evap_max = mean(total_evap_max, na.rm = T),
    wind_min = mean(wind_min, na.rm = T),
    wind_sd = mean(wind_sd, na.rm = T),
    water_soil_mean = mean(water_soil_mean, na.rm = T),
    water_soil_sd = mean(water_soil_sd, na.rm = T),
    elevation = mean(elevation, na.rm = T),
    aspect = mean(aspect, na.rm = T),
    TPI = mean(TPI, na.rm = T),
    TRI = mean(TRI, na.rm = T)
  ) %>%
  as_tibble()


env_fire_regimes <- left_join(pred_df, env_m, by = c("plot_id", "x", "y"))
skim(env_fire_regimes)
env_fire_regimes <- na.omit(env_fire_regimes)

ggplot(env_fire_regimes, aes(x = x, y = y, z = TPI, fill = TPI)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log1p")

# 2. Run PCA on land use data (Exclude ID and Class)
pca_env <- prcomp(env_fire_regimes[, c(24:36)], scale. = TRUE)

env_fire_regimes$class <- as.factor(env_fire_regimes$class)

# --- 1. PREPARE DATA ---

# A. Extract PCA Coordinates (The raw data)
# 'fortify' automatically joins the PC scores with your original data
pca_data <- fortify(pca_env, data = env_fire_regimes)

# B. Calculate Centroids (The "Mean" dots you want)
centroids <- pca_data %>%
  group_by(class) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

# C. Extract Loadings (The Arrows)
# We need to scale them so they are visible on the plot
loadings <- as.data.frame(pca_env$rotation)
loadings$Variable <- rownames(loadings)
scale_factor <- 5 # Tweak this if arrows are too short/long
loadings$PC1 <- loadings$PC1 * scale_factor
loadings$PC2 <- loadings$PC2 * scale_factor


# --- 2. BUILD THE PLOT ---
pdf("Figs/PCA_env_fire.pdf", paper = "a4r", width = 0, height = 0)

ggplot() +

  # LAYER 1: Ellipses (95% Confidence Regions)
  # We use the raw 'pca_data' here, but 'geom="polygon"' draws shapes, not points.
  stat_ellipse(
    data = pca_data,
    aes(x = PC1, y = PC2, fill = class, color = class),
    geom = "polygon",
    alpha = 0.2,
    level = 0.95,
    show.legend = TRUE
  ) +

  # LAYER 2: Centroids (The Big Dots)
  # We use the 'centroids' dataframe, so only 5 dots appear.
  geom_point(
    data = centroids,
    aes(x = PC1, y = PC2, fill = class),
    size = 6,
    shape = 21,
    color = "black",
    stroke = 1.5
  ) +

  # LAYER 3: Arrows (Loadings)
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "grey40"
  ) +

  # LAYER 4: Arrow Labels
  geom_text(
    data = loadings,
    aes(x = PC1, y = PC2, label = Variable),
    color = "black",
    size = 3.5,
    fontface = "bold",
    vjust = -0.5
  ) +

  # STYLING
  scale_fill_viridis_d(option = "turbo", name = "Fire Regime") +
  scale_color_viridis_d(option = "turbo", name = "Fire Regime") +
  theme_bw() +
  labs(
    title = "Environmental niche of Fire Regimes",
    subtitle = "Centroids (means) and 95% Confidence Ellipses",
    x = paste0(
      "PC1 (",
      round(summary(pca_land)$importance[2, 1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca_land)$importance[2, 2] * 100, 1),
      "%)"
    )
  )
dev.off()


# Spatially varying trends ------------------------------------------------

## Fire occurrence ---------------------------------------------------------

fire_occ_df <- readRDS("Data/fire_occ_df_modeling.rds")

# --- 1. Prepare Data with dtplyr ---

# Convert to a lazy_dt object (High performance, dplyr syntax)
fire_occ_df <- lazy_dt(fire_occ_df)

start_info <- fire_occ_df %>%
  summarise(start_year = min(year), .by = NULL) %>% # Need .by=NULL or similar in older dplyr?
  as.data.table() # Collect the minimum year
start_year <- start_info$start_year
start_month <- fire_occ_df %>%
  filter(year == start_year) %>%
  summarise(start_month = min(month), .by = NULL) %>%
  as.data.table() # Collect the minimum month for the start year
start_month <- start_month$start_month

# 2. Add the 't' column using mutate
fire_occ_df <- fire_occ_df %>%
  mutate(t = (year - start_year) * 12 + (month - start_month) + 1)

# Create Harmonic Terms using standard dplyr syntax
# dtplyr translates this to fast data.table code automatically
fire_occ_df <- fire_occ_df %>%
  mutate(
    sin_term = sin(2 * pi * month / 12),
    cos_term = cos(2 * pi * month / 12)
  ) %>%
  as.data.table()
# Collect the prepared data into a distinct data.table for the loop
# (We need the actual data in memory to split it into chunks)

# --- 2. Define the Efficient Solver Function ---

# A generic solver for ANY family (Binomial, Poisson, Gamma, etc.)
get_harmonic_glm <- function(y, t, sin_t, cos_t, fam_obj) {
  # 1. Validation (Check for NAs or static data)
  # For Binomial: need 0s and 1s. For Gamma: need > 0.
  clean_idx <- which(!is.na(y))
  if (length(clean_idx) < 10) {
    return(list(trend = NA_real_, pval = NA_real_))
  }

  y_clean <- y[clean_idx]

  # Check for zero variance (e.g., all 0s or all 1s) - GLM will crash otherwise
  if (var(y_clean) == 0) {
    return(list(trend = 0, pval = 1))
  }

  # 2. Prepare Design Matrix
  # fastglm requires a Matrix object, not a dataframe
  X <- cbind(1, t[clean_idx], sin_t[clean_idx], cos_t[clean_idx])

  # 3. Fit the Model (The Iterative Step)
  # method=2 uses a highly optimized step-halving solver
  fit <- try(
    fastglm(x = X, y = y_clean, family = fam_obj, method = 2),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(list(trend = NA_real_, pval = NA_real_))
  }

  # 4. Extract Stats
  # fastglm returns coefficients and standard errors directly
  coef_trend <- fit$coefficients[2] # 2nd coef is Time
  se_trend <- fit$se[2]

  # Calculate P-value (Wald Test)
  # For GLMs, we use Z-score (Normal distribution approximation)
  z_stat <- coef_trend / se_trend
  p_val <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)

  return(list(trend = coef_trend, pval = p_val))
}

# --- 3. Run Parallel Loop with dtplyr Syntax ---

# Setup Parallel Backend
n_cores <- parallel::detectCores() - 2
plan(multisession, workers = n_cores)


# Split IDs into chunks (memory safety)
all_plots <- unique(fire_occ_df$plot_id)
chunks <- split(all_plots, ceiling(seq_along(all_plots) / 5000)) # 5000 pixels/chunk

my_family <- binomial(link = "logit")

results_df <- foreach(
  chunk = chunks,
  .combine = bind_rows,
  .packages = c("dtplyr", "dplyr", "data.table", "fastglm")
) %dopar%
  {
    chunk_data <- lazy_dt(fire_occ_df[plot_id %in% chunk])

    chunk_res_dt <- chunk_data %>%
      group_by(plot_id) %>%
      summarise(
        # Pass the family object to the function
        stats = list(get_harmonic_glm(
          fire_occurrence,
          t,
          sin_term,
          cos_term,
          my_family
        ))
      ) %>%
      as_tibble()

    final_chunk <- chunk_res_dt %>%
      mutate(
        Trend = unlist(lapply(stats, `[[`, "trend")),
        P_Value = unlist(lapply(stats, `[[`, "pval"))
      ) %>%
      dplyr::select(-stats)

    return(final_chunk)
  }

# --- 4. Mapping ---

# Merge coordinates back
coords <- fire_occ_df %>%
  lazy_dt() %>%
  dplyr::select(plot_id, x, y) %>%
  distinct() %>%
  as_tibble()

final_map_data <- inner_join(results_df, coords, by = "plot_id")

# Plot
# 1. Determine robust limits (e.g., 1st and 99th percentile)
limits <- quantile(final_map_data$Trend, probs = c(0.01, 0.99), na.rm = TRUE)
max_abs <- max(abs(limits)) # Make it symmetric around 0

x_lim <- st_bbox(Cerrado$geometry)[c(1, 3)]
y_lim <- st_bbox(Cerrado$geometry)[c(2, 4)]

pdf("Figs/fire_occ_sptrend.pdf", paper = "a4r", width = 0, height = 0)
(p1 <- ggplot(final_map_data, aes(x = x, y = y, fill = Trend)) +
  geom_sf(data = Cerrado, inherit.aes = F, fill = gray(0.9, 0.5), color = NA) +
  geom_tile(aes(alpha = -log10(P_Value))) +
  scale_fill_gradient2(
    low = viridis(2)[1],
    mid = gray(0.9, 0.5),
    high = viridis(2)[2],
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    na.value = NA,
    name = "Fire occurrence",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  # Alpha Scale (Control the fade)
  scale_alpha_continuous(
    range = c(0, 1), # 0 = Invisible, 1 = Solid
    limits = c(0, 3), # Cap confidence at P = 0.001 (10^-3)
    oob = scales::squish, # Anything P < 0.001 stays fully opaque
    name = "Confidence\n(-log10 P)"
  ) +
  theme_minimal() +
  labs(title = "A", x = "", y = "") +
  theme(axis.text = element_blank()) +
  guides(alpha = "none"))
dev.off()

## Fire metrics ---------------------------------------------------------
fire_metrics_df <- readRDS("Data/fire_metrics_df_modeling.rds")

fire_metrics_df <- lazy_dt(fire_metrics_df)

start_info <- fire_metrics_df %>%
  summarise(start_year = min(year), .by = NULL) %>% # Need .by=NULL or similar in older dplyr?
  as.data.table() # Collect the minimum year
start_year <- start_info$start_year
start_month <- fire_metrics_df %>%
  filter(year == start_year) %>%
  summarise(start_month = min(month), .by = NULL) %>%
  as.data.table() # Collect the minimum month for the start year
start_month <- start_month$start_month

# 2. Add the 't' column using mutate
fire_metrics_df <- fire_metrics_df %>%
  mutate(t = (year - start_year) * 12 + (month - start_month) + 1)

# Create Harmonic Terms using standard dplyr syntax
# dtplyr translates this to fast data.table code automatically
fire_metrics_df <- fire_metrics_df %>%
  mutate(
    sin_term = sin(2 * pi * month / 12),
    cos_term = cos(2 * pi * month / 12),
    lsm_c_cai_sd = log1p(lsm_c_cai_sd),
    lsm_c_circle_mn = log1p(lsm_c_circle_mn),
    lsm_c_core_mn = log1p(lsm_c_core_mn),
    lsm_c_dcore_sd = log1p(lsm_c_dcore_sd),
    lsm_c_enn_sd = log1p(lsm_c_enn_sd),
    lsm_c_pd = log1p(lsm_c_pd),
    lsm_c_split = log1p(lsm_c_split)
  ) %>%
  as.data.table()
# Collect the prepared data into a distinct data.table for the loop
# (We need the actual data in memory to split it into chunks)

# --- 2. Define the Efficient Solver Function ---

# This function fits the harmonic model using fast Linear Algebra
get_harmonic_stats <- function(y, t, sin_t, cos_t) {
  # 1. Validation: Need enough data points
  if (sum(!is.na(y)) < 10) {
    return(list(trend = NA_real_, pval = NA_real_))
  }

  # 2. Prepare Matrix (Time + Seasonality)
  X <- cbind(1, t, sin_t, cos_t)
  good <- complete.cases(X, y)

  if (sum(good) < 10) {
    return(list(trend = NA_real_, pval = NA_real_))
  }

  # 3. Fit Linear Model (Using .lm.fit for raw speed)
  # This skips all the overhead of the standard lm() function
  fit <- .lm.fit(X[good, , drop = FALSE], y[good])

  # 4. Extract Coefficients
  coef_trend <- fit$coefficients[2] # The 2nd coef is 't' (Trend)

  # 5. Calculate Statistics (Manual calculation is faster than summary())
  res <- fit$residuals
  n <- length(res)
  p <- ncol(X)
  rss <- sum(res^2)

  if (rss <= 0 || (n - p) <= 0) {
    return(list(trend = coef_trend, pval = NA_real_))
  }

  sigma2 <- rss / (n - p)
  XtX_inv <- try(solve(crossprod(X[good, ])), silent = TRUE)

  if (inherits(XtX_inv, "try-error")) {
    return(list(trend = coef_trend, pval = NA_real_))
  }

  se_trend <- sqrt(sigma2 * XtX_inv[2, 2])
  t_stat <- coef_trend / se_trend
  p_val <- 2 * pt(abs(t_stat), df = n - p, lower.tail = FALSE)

  return(list(trend = coef_trend, pval = p_val))
}

### Fire Scar Integrity Variation (lsm_c_cai_sd) ----------------------------
# --- 3. Run Parallel Loop with dtplyr Syntax ---

# Setup Parallel Backend
n_cores <- parallel::detectCores() - 2
plan(multisession, workers = n_cores)

# Split IDs into chunks (memory safety)
all_plots <- unique(fire_metrics_df$plot_id)
chunks <- split(all_plots, ceiling(seq_along(all_plots) / 5000)) # 5000 pixels/chunk

results_df <- foreach(
  chunk = chunks,
  .combine = bind_rows,
  .packages = c("dtplyr", "dplyr", "data.table")
) %dopar%
  {
    # 1. Filter data for this chunk
    # dt_ready is a data.table, so we wrap it in lazy_dt() to use dplyr syntax
    chunk_data <- lazy_dt(fire_metrics_df[plot_id %in% chunk])

    # 2. Apply Solver Pixel-by-Pixel
    chunk_res <- chunk_data %>%
      group_by(plot_id) %>%
      summarise(
        # We return a list column to keep trend and pval together (runs regression once)
        stats = list(get_harmonic_stats(lsm_c_cai_sd, t, sin_term, cos_term))
      ) %>%
      # Unpack the list column into separate columns
      mutate(
        Trend = unlist(lapply(stats, `[[`, "trend")),
        P_Value = unlist(lapply(stats, `[[`, "pval"))
      ) %>%
      dplyr::select(-stats) %>% # Remove the temporary list column
      as_tibble() # Convert to tibble to return from loop

    return(chunk_res)
  }

# --- 4. Mapping ---

# Merge coordinates back
coords <- fire_metrics_df %>%
  lazy_dt() %>%
  dplyr::select(plot_id, x, y) %>%
  distinct() %>%
  as_tibble()

final_map_data <- inner_join(results_df, coords, by = "plot_id")

# Plot

pdf("Figs/fire_cai_sd_sptrend.pdf", paper = "a4r", width = 0, height = 0)
(p2 <- ggplot(final_map_data, aes(x = x, y = y, fill = exp(Trend) - 1)) +
  geom_sf(data = Cerrado, inherit.aes = F, fill = gray(0.9, 0.5), color = NA) +
  geom_tile(aes(alpha = -log10(P_Value))) +
  scale_fill_gradient2(
    low = viridis(2)[1],
    mid = gray(0.9, 0.5),
    high = viridis(2)[2],
    midpoint = 0,
    limits = quantile(
      exp(final_map_data$Trend) - 1,
      c(0.001, 0.999),
      na.rm = TRUE
    ),
    na.value = NA,
    name = "Integrity\nvariation",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  # Alpha Scale (Control the fade)
  scale_alpha_continuous(
    range = c(0, 1), # 0 = Invisible, 1 = Solid
    limits = c(0, 3), # Cap confidence at P = 0.001 (10^-3)
    oob = scales::squish, # Anything P < 0.001 stays fully opaque
    name = "Confidence\n(-log10 P)"
  ) +
  theme_minimal() +
  labs(title = "G", x = "", y = "") +
  theme(axis.text = element_blank()) +
  guides(alpha = "none"))
dev.off()

### Fire Scar Shape Regularity (lsm_c_circle_mn) ----------------------------
# --- 3. Run Parallel Loop with dtplyr Syntax ---

# Setup Parallel Backend
n_cores <- parallel::detectCores() - 2
plan(multisession, workers = n_cores)

# Split IDs into chunks (memory safety)
all_plots <- unique(fire_metrics_df$plot_id)
chunks <- split(all_plots, ceiling(seq_along(all_plots) / 5000)) # 5000 pixels/chunk

results_df <- foreach(
  chunk = chunks,
  .combine = bind_rows,
  .packages = c("dtplyr", "dplyr", "data.table")
) %dopar%
  {
    # 1. Filter data for this chunk
    # dt_ready is a data.table, so we wrap it in lazy_dt() to use dplyr syntax
    chunk_data <- lazy_dt(fire_metrics_df[plot_id %in% chunk])

    # 2. Apply Solver Pixel-by-Pixel
    chunk_res <- chunk_data %>%
      group_by(plot_id) %>%
      summarise(
        # We return a list column to keep trend and pval together (runs regression once)
        stats = list(get_harmonic_stats(lsm_c_circle_mn, t, sin_term, cos_term))
      ) %>%
      # Unpack the list column into separate columns
      mutate(
        Trend = unlist(lapply(stats, `[[`, "trend")),
        P_Value = unlist(lapply(stats, `[[`, "pval"))
      ) %>%
      dplyr::select(-stats) %>% # Remove the temporary list column
      as_tibble() # Convert to tibble to return from loop

    return(chunk_res)
  }

# --- 4. Mapping ---

# Merge coordinates back
coords <- fire_metrics_df %>%
  lazy_dt() %>%
  dplyr::select(plot_id, x, y) %>%
  distinct() %>%
  as_tibble()

final_map_data <- inner_join(results_df, coords, by = "plot_id")

# Plot

pdf("Figs/fire_circle_mn_sptrend.pdf", paper = "a4r", width = 0, height = 0)
(p3 <- ggplot(final_map_data, aes(x = x, y = y, fill = exp(Trend) - 1)) +
  geom_sf(data = Cerrado, inherit.aes = F, fill = gray(0.9, 0.5), color = NA) +
  geom_tile(aes(alpha = -log10(P_Value))) +
  scale_fill_gradient2(
    low = viridis(2)[1],
    mid = gray(0.9, 0.5),
    high = viridis(2)[2],
    midpoint = 0,
    limits = quantile(
      exp(final_map_data$Trend) - 1,
      c(0.001, 0.999),
      na.rm = TRUE
    ),
    na.value = NA,
    name = "Shape\nregularity",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  # Alpha Scale (Control the fade)
  scale_alpha_continuous(
    range = c(0, 1), # 0 = Invisible, 1 = Solid
    limits = c(0, 3), # Cap confidence at P = 0.001 (10^-3)
    oob = scales::squish, # Anything P < 0.001 stays fully opaque
    name = "Confidence\n(-log10 P)"
  ) +
  theme_minimal() +
  labs(title = "C", x = "", y = "") +
  theme(axis.text = element_blank()) +
  guides(alpha = "none"))
dev.off()

### Mean Fire Scar Core Size (lsm_c_core_mn) ----------------------------

# --- 3. Run Parallel Loop with dtplyr Syntax ---

# Setup Parallel Backend
n_cores <- parallel::detectCores() - 2
plan(multisession, workers = n_cores)

# Split IDs into chunks (memory safety)
all_plots <- unique(fire_metrics_df$plot_id)
chunks <- split(all_plots, ceiling(seq_along(all_plots) / 5000)) # 5000 pixels/chunk

results_df <- foreach(
  chunk = chunks,
  .combine = bind_rows,
  .packages = c("dtplyr", "dplyr", "data.table")
) %dopar%
  {
    # 1. Filter data for this chunk
    # dt_ready is a data.table, so we wrap it in lazy_dt() to use dplyr syntax
    chunk_data <- lazy_dt(fire_metrics_df[plot_id %in% chunk])

    # 2. Apply Solver Pixel-by-Pixel
    chunk_res <- chunk_data %>%
      group_by(plot_id) %>%
      summarise(
        # We return a list column to keep trend and pval together (runs regression once)
        stats = list(get_harmonic_stats(lsm_c_core_mn, t, sin_term, cos_term))
      ) %>%
      # Unpack the list column into separate columns
      mutate(
        Trend = unlist(lapply(stats, `[[`, "trend")),
        P_Value = unlist(lapply(stats, `[[`, "pval"))
      ) %>%
      dplyr::select(-stats) %>% # Remove the temporary list column
      as_tibble() # Convert to tibble to return from loop

    return(chunk_res)
  }

# --- 4. Mapping ---

# Merge coordinates back
coords <- fire_metrics_df %>%
  lazy_dt() %>%
  dplyr::select(plot_id, x, y) %>%
  distinct() %>%
  as_tibble()

final_map_data <- inner_join(results_df, coords, by = "plot_id")

# Plot

pdf("Figs/fire_core_mn_sptrend.pdf", paper = "a4r", width = 0, height = 0)
(p4 <- ggplot(final_map_data, aes(x = x, y = y, fill = exp(Trend) - 1)) +
  geom_sf(data = Cerrado, inherit.aes = F, fill = gray(0.9, 0.5), color = NA) +
  geom_tile(aes(alpha = -log10(P_Value))) +
  scale_fill_gradient2(
    low = viridis(2)[1],
    mid = gray(0.9, 0.5),
    high = viridis(2)[2],
    midpoint = 0,
    limits = quantile(
      exp(final_map_data$Trend) - 1,
      c(0.001, 0.999),
      na.rm = TRUE
    ),
    na.value = NA,
    name = "Core\nsize",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  # Alpha Scale (Control the fade)
  scale_alpha_continuous(
    range = c(0, 1), # 0 = Invisible, 1 = Solid
    limits = c(0, 3), # Cap confidence at P = 0.001 (10^-3)
    oob = scales::squish, # Anything P < 0.001 stays fully opaque
    name = "Confidence\n(-log10 P)"
  ) +
  theme_minimal() +
  labs(title = "F", x = "", y = "") +
  theme(axis.text = element_blank()) +
  guides(alpha = "none"))
dev.off()

### Internal Fire Scar Heterogeneity (lsm_c_dcore_sd) ----------------------------

# --- 3. Run Parallel Loop with dtplyr Syntax ---

# Setup Parallel Backend
n_cores <- parallel::detectCores() - 2
plan(multisession, workers = n_cores)

# Split IDs into chunks (memory safety)
all_plots <- unique(fire_metrics_df$plot_id)
chunks <- split(all_plots, ceiling(seq_along(all_plots) / 5000)) # 5000 pixels/chunk

results_df <- foreach(
  chunk = chunks,
  .combine = bind_rows,
  .packages = c("dtplyr", "dplyr", "data.table")
) %dopar%
  {
    # 1. Filter data for this chunk
    # dt_ready is a data.table, so we wrap it in lazy_dt() to use dplyr syntax
    chunk_data <- lazy_dt(fire_metrics_df[plot_id %in% chunk])

    # 2. Apply Solver Pixel-by-Pixel
    chunk_res <- chunk_data %>%
      group_by(plot_id) %>%
      summarise(
        # We return a list column to keep trend and pval together (runs regression once)
        stats = list(get_harmonic_stats(lsm_c_dcore_sd, t, sin_term, cos_term))
      ) %>%
      # Unpack the list column into separate columns
      mutate(
        Trend = unlist(lapply(stats, `[[`, "trend")),
        P_Value = unlist(lapply(stats, `[[`, "pval"))
      ) %>%
      dplyr::select(-stats) %>% # Remove the temporary list column
      as_tibble() # Convert to tibble to return from loop

    return(chunk_res)
  }

# --- 4. Mapping ---

# Merge coordinates back
coords <- fire_metrics_df %>%
  lazy_dt() %>%
  dplyr::select(plot_id, x, y) %>%
  distinct() %>%
  as_tibble()

final_map_data <- inner_join(results_df, coords, by = "plot_id")

# Plot

pdf("Figs/fire_dcore_sd_sptrend.pdf", paper = "a4r", width = 0, height = 0)
(p5 <- ggplot(final_map_data, aes(x = x, y = y, fill = exp(Trend) - 1)) +
  geom_sf(data = Cerrado, inherit.aes = F, fill = gray(0.9, 0.5), color = NA) +
  geom_tile(aes(alpha = -log10(P_Value))) +
  scale_fill_gradient2(
    low = viridis(2)[1],
    mid = gray(0.9, 0.5),
    high = viridis(2)[2],
    midpoint = 0,
    limits = quantile(
      exp(final_map_data$Trend) - 1,
      c(0.001, 0.999),
      na.rm = TRUE
    ),
    na.value = NA,
    name = "Cohesion\nvariation",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  # Alpha Scale (Control the fade)
  scale_alpha_continuous(
    range = c(0, 1), # 0 = Invisible, 1 = Solid
    limits = c(0, 3), # Cap confidence at P = 0.001 (10^-3)
    oob = scales::squish, # Anything P < 0.001 stays fully opaque
    name = "Confidence\n(-log10 P)"
  ) +
  theme_minimal() +
  labs(title = "H", x = "", y = "") +
  theme(axis.text = element_blank()) +
  guides(alpha = "none"))
dev.off()

### Fire Scar Isolation Variation (lsm_c_enn_sd) ----------------------------

# --- 3. Run Parallel Loop with dtplyr Syntax ---

# Setup Parallel Backend
n_cores <- parallel::detectCores() - 2
plan(multisession, workers = n_cores)

# Split IDs into chunks (memory safety)
all_plots <- unique(fire_metrics_df$plot_id)
chunks <- split(all_plots, ceiling(seq_along(all_plots) / 5000)) # 5000 pixels/chunk

results_df <- foreach(
  chunk = chunks,
  .combine = bind_rows,
  .packages = c("dtplyr", "dplyr", "data.table")
) %dopar%
  {
    # 1. Filter data for this chunk
    # dt_ready is a data.table, so we wrap it in lazy_dt() to use dplyr syntax
    chunk_data <- lazy_dt(fire_metrics_df[plot_id %in% chunk])

    # 2. Apply Solver Pixel-by-Pixel
    chunk_res <- chunk_data %>%
      group_by(plot_id) %>%
      summarise(
        # We return a list column to keep trend and pval together (runs regression once)
        stats = list(get_harmonic_stats(lsm_c_enn_sd, t, sin_term, cos_term))
      ) %>%
      # Unpack the list column into separate columns
      mutate(
        Trend = unlist(lapply(stats, `[[`, "trend")),
        P_Value = unlist(lapply(stats, `[[`, "pval"))
      ) %>%
      dplyr::select(-stats) %>% # Remove the temporary list column
      as_tibble() # Convert to tibble to return from loop

    return(chunk_res)
  }

# --- 4. Mapping ---

# Merge coordinates back
coords <- fire_metrics_df %>%
  lazy_dt() %>%
  dplyr::select(plot_id, x, y) %>%
  distinct() %>%
  as_tibble()

final_map_data <- inner_join(results_df, coords, by = "plot_id")

# Plot

pdf("Figs/fire_enn_sd_sptrend.pdf", paper = "a4r", width = 0, height = 0)
(p6 <- ggplot(final_map_data, aes(x = x, y = y, fill = exp(Trend) - 1)) +
  geom_sf(data = Cerrado, inherit.aes = F, fill = gray(0.9, 0.5), color = NA) +
  geom_tile(aes(alpha = -log10(P_Value))) +
  scale_fill_gradient2(
    low = viridis(2)[1],
    mid = gray(0.9, 0.5),
    high = viridis(2)[2],
    midpoint = 0,
    limits = quantile(
      exp(final_map_data$Trend) - 1,
      c(0.001, 0.999),
      na.rm = TRUE
    ),
    na.value = NA,
    name = "Isolation\nvariation",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  # Alpha Scale (Control the fade)
  scale_alpha_continuous(
    range = c(0, 1), # 0 = Invisible, 1 = Solid
    limits = c(0, 3), # Cap confidence at P = 0.001 (10^-3)
    oob = scales::squish, # Anything P < 0.001 stays fully opaque
    name = "Confidence\n(-log10 P)"
  ) +
  theme_minimal() +
  labs(title = "D", x = "", y = "") +
  theme(axis.text = element_blank()) +
  guides(alpha = "none"))
dev.off()

### Fire Scar Density (lsm_c_pd) ----------------------------

# --- 3. Run Parallel Loop with dtplyr Syntax ---

# Setup Parallel Backend
n_cores <- parallel::detectCores() - 2
plan(multisession, workers = n_cores)

# Split IDs into chunks (memory safety)
all_plots <- unique(fire_metrics_df$plot_id)
chunks <- split(all_plots, ceiling(seq_along(all_plots) / 5000)) # 5000 pixels/chunk

results_df <- foreach(
  chunk = chunks,
  .combine = bind_rows,
  .packages = c("dtplyr", "dplyr", "data.table")
) %dopar%
  {
    # 1. Filter data for this chunk
    # dt_ready is a data.table, so we wrap it in lazy_dt() to use dplyr syntax
    chunk_data <- lazy_dt(fire_metrics_df[plot_id %in% chunk])

    # 2. Apply Solver Pixel-by-Pixel
    chunk_res <- chunk_data %>%
      group_by(plot_id) %>%
      summarise(
        # We return a list column to keep trend and pval together (runs regression once)
        stats = list(get_harmonic_stats(lsm_c_pd, t, sin_term, cos_term))
      ) %>%
      # Unpack the list column into separate columns
      mutate(
        Trend = unlist(lapply(stats, `[[`, "trend")),
        P_Value = unlist(lapply(stats, `[[`, "pval"))
      ) %>%
      dplyr::select(-stats) %>% # Remove the temporary list column
      as_tibble() # Convert to tibble to return from loop

    return(chunk_res)
  }

# --- 4. Mapping ---

# Merge coordinates back
coords <- fire_metrics_df %>%
  lazy_dt() %>%
  dplyr::select(plot_id, x, y) %>%
  distinct() %>%
  as_tibble()

final_map_data <- inner_join(results_df, coords, by = "plot_id")

# Plot

pdf("Figs/fire_pd_sptrend.pdf", paper = "a4r", width = 0, height = 0)
(p7 <- ggplot(final_map_data, aes(x = x, y = y, fill = exp(Trend) - 1)) +
  geom_sf(data = Cerrado, inherit.aes = F, fill = gray(0.9, 0.5), color = NA) +
  geom_tile(aes(alpha = -log10(P_Value))) +
  scale_fill_gradient2(
    low = viridis(2)[1],
    mid = gray(0.9, 0.5),
    high = viridis(2)[2],
    midpoint = 0,
    limits = quantile(
      exp(final_map_data$Trend) - 1,
      c(0.001, 0.999),
      na.rm = TRUE
    ),
    na.value = NA,
    name = "Density",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  # Alpha Scale (Control the fade)
  scale_alpha_continuous(
    range = c(0, 1), # 0 = Invisible, 1 = Solid
    limits = c(0, 3), # Cap confidence at P = 0.001 (10^-3)
    oob = scales::squish, # Anything P < 0.001 stays fully opaque
    name = "Confidence\n(-log10 P)"
  ) +
  theme_minimal() +
  labs(title = "B", x = "", y = "") +
  theme(axis.text = element_blank()) +
  guides(alpha = "none"))
dev.off()

### Fire Scar Fragmentation (lsm_c_split) ----------------------------

# --- 3. Run Parallel Loop with dtplyr Syntax ---

# Setup Parallel Backend
n_cores <- parallel::detectCores() - 2
plan(multisession, workers = n_cores)

# Split IDs into chunks (memory safety)
all_plots <- unique(fire_metrics_df$plot_id)
chunks <- split(all_plots, ceiling(seq_along(all_plots) / 5000)) # 5000 pixels/chunk

results_df <- foreach(
  chunk = chunks,
  .combine = bind_rows,
  .packages = c("dtplyr", "dplyr", "data.table")
) %dopar%
  {
    # 1. Filter data for this chunk
    # dt_ready is a data.table, so we wrap it in lazy_dt() to use dplyr syntax
    chunk_data <- lazy_dt(fire_metrics_df[plot_id %in% chunk])

    # 2. Apply Solver Pixel-by-Pixel
    chunk_res <- chunk_data %>%
      group_by(plot_id) %>%
      summarise(
        # We return a list column to keep trend and pval together (runs regression once)
        stats = list(get_harmonic_stats(lsm_c_split, t, sin_term, cos_term))
      ) %>%
      # Unpack the list column into separate columns
      mutate(
        Trend = unlist(lapply(stats, `[[`, "trend")),
        P_Value = unlist(lapply(stats, `[[`, "pval"))
      ) %>%
      dplyr::select(-stats) %>% # Remove the temporary list column
      as_tibble() # Convert to tibble to return from loop

    return(chunk_res)
  }

# --- 4. Mapping ---

# Merge coordinates back
coords <- fire_metrics_df %>%
  lazy_dt() %>%
  dplyr::select(plot_id, x, y) %>%
  distinct() %>%
  as_tibble()

final_map_data <- inner_join(results_df, coords, by = "plot_id")

# Plot
pdf("Figs/fire_split_sptrend.pdf", paper = "a4r", width = 0, height = 0)
(p8 <- ggplot(final_map_data, aes(x = x, y = y, fill = exp(Trend) - 1)) +
  geom_sf(data = Cerrado, inherit.aes = F, fill = gray(0.9, 0.5), color = NA) +
  geom_tile(aes(alpha = -log10(P_Value))) +
  scale_fill_gradient2(
    low = viridis(2)[1],
    mid = gray(0.9, 0.5),
    high = viridis(2)[2],
    midpoint = 0,
    limits = quantile(
      exp(final_map_data$Trend) - 1,
      c(0.001, 0.999),
      na.rm = TRUE
    ),
    na.value = NA,
    name = "Fragmentation",
    breaks = scales::pretty_breaks(n = 4),
    guide = guide_colourbar(
      barwidth = unit(4, "cm"),
      barheight = unit(0.5, "cm"),
      position = "bottom",
      theme = theme(
        legend.spacing.y = unit(8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        legend.title = element_text(size = 9),
        legend.text = element_text(
          size = 8,
          angle = -45,
          hjust = 0.1,
          vjust = 0.5
        )
      )
    )
  ) +
  coord_sf(crs = metric_crs, xlim = x_lim, ylim = y_lim, expand = FALSE) +
  # Alpha Scale (Control the fade)
  scale_alpha_continuous(
    range = c(0, 1), # 0 = Invisible, 1 = Solid
    limits = c(0, 3), # Cap confidence at P = 0.001 (10^-3)
    oob = scales::squish, # Anything P < 0.001 stays fully opaque
    name = "Confidence\n(-log10 P)"
  ) +
  theme_minimal() +
  labs(title = "E", x = "", y = "") +
  theme(axis.text = element_blank()) +
  guides(alpha = "none"))
dev.off()

## Final Figs --------------------------------------------------------

final_fig <-
  (p1 | p7 | p3 | p6) /
  (p8 | p4 | p2 | p5)

final_fig <- final_fig +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(
        t = 10, # top
        r = 10, # right
        b = 10, # bottom
        l = 10, # left
        unit = "mm"
      )
    )
  )

ggsave(
  "Figs/fire_metrics_sptrend_figure.pdf",
  final_fig,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 900,
  device = cairo_pdf
)

# Land use map (2024) --------------------------------------------

# Load the Cerrado geometry
Cerrado <- readRDS("Data/Cerrado.rds")
print(Cerrado)

# --- Define a metric projection for Brazil ---
# SIRGAS 2000 / Brazil Polyconic (EPSG:5880) is a good choice for the whole country.
metric_crs <- "EPSG:5880"

# Select only the geometry object from "Cerrado"
Cerrado <- Cerrado %>% dplyr::select(geometry)
Cerrado <- sf::st_transform(Cerrado, metric_crs)
Cerrado <- st_simplify(Cerrado, dTolerance = 1000)

fire_occ_df <- readRDS("Data/fire_occ_df_modeling.rds")

legend_dt <- read.delim2(
  "Data/Codigos-da-legenda-colecao-10.csv"
)


glimpse(fire_occ_df)


# --- 1. PREPARE THE DATA (Same as before) ---

# A. Create mapping
col_to_id <- tibble(
  col_name = c(
    "forest",
    "savanna",
    "forest_plantation",
    "wet_herbaceous",
    "grassland",
    "pasture",
    "temporary_crop",
    "soybean",
    "sugar_cane",
    "perennial_crop",
    "mosaic_uses",
    "urban",
    "non_vegetated",
    "water"
  ),
  Class_ID = c(3, 4, 9, 11, 12, 15, 19, 39, 20, 36, 21, 24, 22, 33)
)

# B. Process dominant land use
latest_year <- max(fire_occ_df$year)

land_use_map_data <- fire_occ_df %>%
  filter(year == latest_year) %>%
  dplyr::select(plot_id, x, y, all_of(col_to_id$col_name)) %>%
  pivot_longer(
    cols = any_of(col_to_id$col_name),
    names_to = "land_use_col",
    values_to = "area"
  ) %>%
  group_by(plot_id, x, y) %>%
  filter(area == max(area, na.rm = TRUE)) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  left_join(col_to_id, by = c("land_use_col" = "col_name")) %>%
  left_join(legend_dt, by = "Class_ID") %>%
  mutate(Color = as.character(Color))

# --- 2. PREPARE LEGEND & SPATIAL LAYERS ---

# A. Filter Legend (KEY STEP)
# Identify which colors actually exist in the final map data
observed_colors <- unique(land_use_map_data$Color)

# Create a named vector: 'Hex Code' = 'Description'
# We filter the full legend_dt to only rows present in the map
legend_labels <- legend_dt %>%
  filter(as.character(Color) %in% observed_colors) %>%
  dplyr::select(Color, Description) %>%
  deframe() # Converts 2-column df to named vector (Color -> Description)

# Reverse it for scale_fill_identity (needs names as colors)
# We want: breaks = c("#hex1", "#hex2"), labels = c("Forest", "Savanna")
legend_breaks <- names(legend_labels) # The Hex Codes
legend_texts <- unname(legend_labels) # The Descriptions

# B. Spatial Layers
states <- read_state(year = 2020, showProgress = FALSE) %>%
  st_transform(crs = st_crs(Cerrado))

sa <- ne_countries(
  scale = "medium",
  continent = "South America",
  returnclass = "sf"
) %>%
  st_transform(crs = st_crs(Cerrado))

# --- 3. BUILD THE MAP ---

main_map <- ggplot() +
  # 1. Raster Layer
  geom_tile(data = land_use_map_data, aes(x = x, y = y, fill = Color)) +

  # 2. Vector Overlays
  geom_sf(data = states, fill = NA, color = "gray30", size = 0.2) +
  geom_sf(
    data = Cerrado,
    fill = NA,
    color = "black",
    size = 0.15,
    alpha = 0.6
  ) +

  # 3. Dynamic Legend (Filtered)
  scale_fill_identity(
    guide = "legend",
    name = paste("Land Use and Cover in the Cerrado", latest_year),
    breaks = legend_breaks, # Only shows colors present in data
    labels = legend_texts # Uses correct descriptions
  ) +

  # 4. North Arrow & Scale Bar (ggspatial)
  annotation_scale(
    location = "bl", # Bottom Left
    width_hint = 0.2, # Width relative to plot
    style = "ticks",
    line_width = 1,
    height = unit(0.2, "cm"),
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm")
  ) +
  annotation_north_arrow(
    location = "bl", # Bottom Left (above scale)
    which_north = "true",
    pad_x = unit(0.5, "cm"),
    pad_y = unit(1.0, "cm"), # Push it up slightly above scale bar
    style = north_arrow_fancy_orienteering(
      text_col = "black",
      line_col = "black",
      fill = c("white", "black")
    )
  ) +

  # 5. Coordinate System & Grids
  # 'datum = st_crs(Cerrado)' forces the native coordinates (Meters)
  coord_sf(
    datum = st_crs(Cerrado),
    expand = TRUE,
    xlim = st_bbox(Cerrado)[c(1, 3)],
    ylim = st_bbox(Cerrado)[c(2, 4)]
  ) +

  # 6. Formatting Axis Labels (Meters)
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) + # e.g. "6M" or "600km"
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +

  # 7. Theme
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 8, color = "grey30"), # Show grid coords
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.border = element_blank(),
    legend.text = element_text(size = 8),
    panel.grid.major = element_line(color = "gray90", linetype = "dashed") # Styling the grid
  ) +
  labs(x = "Latitude", y = "Longitude")

# --- 4. INSET & COMBINE ---
inset_map <- ggplot() +
  geom_sf(data = sa, fill = "gray95", color = "gray60", size = 0.2) +
  geom_sf(data = states, fill = NA, color = "gray60", size = 0.2) +
  geom_sf(data = Cerrado, fill = "red", color = NA, alpha = 0.5) +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(fill = "white")
  )

# 1. Extract the legend from your main_map
legend <- get_legend(main_map)

# 2. Remove the legend from the main_map so it doesn't print twice
main_map_no_legend <- main_map + theme(legend.position = "none")

# 3. Compose them all together with ggdraw
final_plot <- ggdraw() +
  # Layer 1: The Map (Full screen, no legend)
  draw_plot(main_map_no_legend, x = -0.1) +

  # Layer 2: The Inset Map (Your custom position)
  draw_plot(inset_map, x = 0.75, y = 0.68, width = 0.25, height = 0.3) +

  # Layer 3: The Legend (Manually positioned below the inset)
  # Adjust x and y to slide it exactly where you want
  draw_plot(legend, x = 0.75, y = 0.3, width = 0.2, height = 0.3)

pdf("Figs/lulc_map_2024.pdf", paper = "a4r", width = 0, height = 0)
print(final_plot)
dev.off()
