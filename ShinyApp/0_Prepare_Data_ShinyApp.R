rm(list = ls())

# Prepare data for Shiny App
library(dtplyr)
library(mclust)
library(dplyr)
library(tidyr)
library(skimr)
library(ggplot2)
library(viridis)
library(terra)
library(geobr)
library(sf)

fire_occ_df <- readRDS("Data/final_fire_occ_data.rds")
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


fire_metrics_df <- readRDS("Data/final_clean_data.rds")
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

ggplot(pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  scale_fill_viridis(trans = "log1p", option = "turbo")

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

ggplot(fire_1stdec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  scale_fill_viridis(trans = "log1p", option = "turbo")

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

ggplot(fire_2nddec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  scale_fill_viridis(trans = "log1p", option = "turbo")

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

ggplot(fire_3rddec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  scale_fill_viridis(option = "turbo")

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

ggplot(fire_4thdec_pred_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  scale_fill_viridis(option = "turbo")

save_raster_from_df <- function(df, target_col, filename) {
  # 1. Prepare Data
  # Rounding to 1m precision to prevent floating-point drift
  spatial_data <- df %>%
    select(x, y, all_of(target_col)) %>%
    mutate(x = round(x), y = round(y))

  # 2. Detect Resolution (Distance between centers)
  # Sort unique coords and find the smallest gap
  res_x <- min(diff(sort(unique(spatial_data$x))))
  res_y <- min(diff(sort(unique(spatial_data$y))))

  # 3. Define the Correct Extent (Padding for Centers)
  # If x is the center, the left edge is x - res/2
  min_x <- min(spatial_data$x) - (res_x / 2)
  max_x <- max(spatial_data$x) + (res_x / 2)
  min_y <- min(spatial_data$y) - (res_y / 2)
  max_y <- max(spatial_data$y) + (res_y / 2)

  # 4. Create the Template Raster
  # We specify the extent manually to ensure perfect alignment
  r_template <- rast(
    xmin = min_x,
    xmax = max_x,
    ymin = min_y,
    ymax = max_y,
    resolution = c(res_x, res_y),
    crs = "EPSG:5880"
  )

  # 5. Create Vector Points
  v <- vect(spatial_data, geom = c("x", "y"), crs = "EPSG:5880")

  # 6. Rasterize
  # Using 'near' (Nearest Neighbor) avoids averaging if points slightly drift
  r_final <- rasterize(v, r_template, field = target_col, fun = "first")

  # 7. Save
  writeRaster(r_final, filename, gdal = c("COMPRESS=LZW"), overwrite = TRUE)

  message(paste("Saved:", filename, "| Resolution:", res_x, "m"))
}

# Global Model
save_raster_from_df(
  pred_df,
  "class",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_Global.tif"
)
save_raster_from_df(
  pred_df,
  "X1",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Jan.tif"
)
save_raster_from_df(
  pred_df,
  "X2",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Feb.tif"
)
save_raster_from_df(
  pred_df,
  "X3",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Mar.tif"
)
save_raster_from_df(
  pred_df,
  "X4",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Apr.tif"
)
save_raster_from_df(
  pred_df,
  "X5",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_May.tif"
)
save_raster_from_df(
  pred_df,
  "X6",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Jun.tif"
)
save_raster_from_df(
  pred_df,
  "X7",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Jul.tif"
)
save_raster_from_df(
  pred_df,
  "X8",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Aug.tif"
)
save_raster_from_df(
  pred_df,
  "X9",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Sep.tif"
)
save_raster_from_df(
  pred_df,
  "X10",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Oct.tif"
)
save_raster_from_df(
  pred_df,
  "X11",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Nov.tif"
)
save_raster_from_df(
  pred_df,
  "X12",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Dec.tif"
)
save_raster_from_df(
  pred_df,
  "lsm_c_cai_sd",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Integrity_Var.tif"
)
save_raster_from_df(
  pred_df,
  "lsm_c_circle_mn",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Shape_Reg.tif"
)
save_raster_from_df(
  pred_df,
  "lsm_c_core_mn",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Core_Size.tif"
)

save_raster_from_df(
  pred_df,
  "lsm_c_dcore_sd",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Internal_Heterog.tif"
)
save_raster_from_df(
  pred_df,
  "lsm_c_enn_sd",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Isolation_Var.tif"
)
save_raster_from_df(
  pred_df,
  "lsm_c_pd",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Density.tif"
)
save_raster_from_df(
  pred_df,
  "lsm_c_split",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Fragmentation.tif"
)


# 3. Save the Decades
save_raster_from_df(
  fire_1stdec_pred_df,
  "class",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_1985_1994.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X1",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Jan.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X2",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Feb.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X3",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Mar.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X4",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Apr.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X5",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_May.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X6",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Jun.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X7",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Jul.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X8",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Aug.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X9",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Sep.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X10",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Oct.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X11",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Nov.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "X12",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Freq_Dec.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "lsm_c_cai_sd",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Integrity_Var.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "lsm_c_circle_mn",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Shape_Reg.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "lsm_c_dcore_sd",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Internal_Heterog.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "lsm_c_enn_sd",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Isolation_Var.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "lsm_c_pd",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Density.tif"
)
save_raster_from_df(
  fire_1stdec_pred_df,
  "lsm_c_split",
  "Output/Decadal_Fire_Metrics/Fire1985_1994_Fragmentation.tif"
)

save_raster_from_df(
  fire_2nddec_pred_df,
  "class",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_1995_2004.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X1",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Jan.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X2",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Feb.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X3",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Mar.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X4",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Apr.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X5",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_May.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X6",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Jun.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X7",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Jul.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X8",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Aug.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X9",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Sep.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X10",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Oct.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X11",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Nov.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "X12",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Freq_Dec.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "lsm_c_cai_sd",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Integrity_Var.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "lsm_c_circle_mn",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Shape_Reg.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "lsm_c_dcore_sd",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Internal_Heterog.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "lsm_c_enn_sd",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Isolation_Var.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "lsm_c_pd",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Density.tif"
)
save_raster_from_df(
  fire_2nddec_pred_df,
  "lsm_c_split",
  "Output/Decadal_Fire_Metrics/Fire1995_2004_Fragmentation.tif"
)

save_raster_from_df(
  fire_3rddec_pred_df,
  "class",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_2005_2014.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X1",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Jan.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X2",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Feb.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X3",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Mar.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X4",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Apr.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X5",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_May.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X6",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Jun.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X7",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Jul.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X8",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Aug.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X9",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Sep.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X10",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Oct.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X11",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Nov.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "X12",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Freq_Dec.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "lsm_c_cai_sd",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Integrity_Var.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "lsm_c_circle_mn",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Shape_Reg.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "lsm_c_dcore_sd",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Internal_Heterog.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "lsm_c_enn_sd",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Isolation_Var.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "lsm_c_pd",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Density.tif"
)
save_raster_from_df(
  fire_3rddec_pred_df,
  "lsm_c_split",
  "Output/Decadal_Fire_Metrics/Fire2005_2014_Fragmentation.tif"
)

save_raster_from_df(
  fire_4thdec_pred_df,
  "class",
  "ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_2015_2024.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X1",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Jan.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X2",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Feb.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X3",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Mar.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X4",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Apr.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X5",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_May.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X6",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Jun.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X7",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Jul.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X8",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Aug.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X9",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Sep.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X10",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Oct.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X11",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Nov.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "X12",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Freq_Dec.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "lsm_c_cai_sd",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Integrity_Var.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "lsm_c_circle_mn",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Shape_Reg.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "lsm_c_dcore_sd",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Internal_Heterog.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "lsm_c_enn_sd",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Isolation_Var.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "lsm_c_pd",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Density.tif"
)
save_raster_from_df(
  fire_4thdec_pred_df,
  "lsm_c_split",
  "Output/Decadal_Fire_Metrics/Fire2015_2024_Fragmentation.tif"
)


# 1. Load your Global Raster (To get the Cerrado Extent)
# We need this to crop the vectors so we don't load the whole Amazon
ref_raster <- rast("ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_Global.tif")

# Create a bounding box polygon of your raster
cerrado_bbox <- st_as_sf(as.polygons(ext(ref_raster), crs = "EPSG:5880"))
cerrado_bbox <- st_transform(cerrado_bbox, 4326) # Geobr is in 4674/4326

message("Downloading and processing vectors... this might take a minute.")

# --- 2. Download and Clean Vectors ---

# A. States
states <- read_state(year = 2020, showProgress = FALSE) %>%
  st_transform(4326) %>%
  st_filter(cerrado_bbox) # Keep only states that touch the Cerrado

# B. Municipalities
munis <- read_municipality(year = 2020, showProgress = FALSE) %>%
  st_transform(4326) %>%
  st_filter(cerrado_bbox) %>% # Crucial: Remove non-Cerrado cities
  st_simplify(dTolerance = 0.01) # Simplify geometry for web display speed

# C. Conservation Units
ucs <- read_conservation_units(showProgress = FALSE) %>%
  st_transform(4326) %>%
  st_filter(cerrado_bbox) %>%
  st_simplify(dTolerance = 0.01) %>% # Simplify (approx 1km resolution)
  st_make_valid()

# D. Indigenous Lands
tis <- read_indigenous_land(showProgress = FALSE) %>%
  st_transform(4326) %>%
  st_filter(cerrado_bbox) %>%
  st_simplify(dTolerance = 0.01) %>% # Simplify (approx 1km resolution)
  st_make_valid()

# --- 3. Save as a single lightweight list ---
vector_data <- list(
  states = states,
  munis = munis,
  ucs = ucs,
  tis = tis
)

saveRDS(vector_data, "ShinyApp/CerradoFireRegimes_DSS/vectors_cerrado.rds")
message("Done! 'vectors_cerrado.rds' created. Put this in your Shiny folder.")
