# Load packages -----------------------------------------------------------
devtools::install_github("Nowosad/popgrids")

library(terra)
library(dplyr)
library(usdm)
library(tidyr) # For pivot_wider
library(tidyterra)
library(purrr)
library(data.table)
library(ggplot2)
library(viridis)
library(visdat)
library(data.table)
library(foreach)
library(doFuture) # The backend for foreach to use the future framework
library(exactextractr)
library(sf)
library(popgrids)
library(stars)
library(stringr)
library(elevatr)


# Organizing predictors ---------------------------------------------------

# Land use and land cover -------------------------------------------------
# MapBiomas – Collection 10 of the annual series of Maps of Land Cover and Use of Brazil, accessed on 24 September through the link: https://code.earthengine.google.com/?accept_repo=users/mapbiomas/user-toolkit
data_directory <- "~/Documents/mapbiomas_lulc_10"
all_lulc_files <- list.files(
  data_directory,
  pattern = "\\.tif$",
  full.names = TRUE
)
lulc_timeseries_terra <- rast(all_lulc_files)
print(lulc_timeseries_terra)

legend_dt <- read.delim2(
  "~/Documents/mapbiomas_lulc_10/Codigos-da-legenda-colecao-10.csv"
)

# Load the Cerrado geometry
Cerrado <- readRDS("Cerrado.rds")
print(Cerrado)

# Select only the geometry object from "Cerrado"
Cerrado <- Cerrado %>% dplyr::select(geometry)
Cerrado <- sf::st_transform(Cerrado, crs(lulc_timeseries_terra))

# --- Create the 9km Analysis Grid (in the same metric CRS) ---
cat("Creating 9km analysis grid...\n")
# Load one of the reprojected binary rasters to use as a template
template_raster <- rast("~/Documents/ERA5/precip.grib")[[1]]
print(template_raster)
template_raster <- terra::project(
  template_raster,
  crs(lulc_timeseries_terra),
  threads = T,
  method = "bilinear"
)

template_raster <- crop(template_raster, Cerrado) %>% mask(Cerrado)


grid_9km <- terra::rast(
  extent = ext(template_raster),
  resolution = res(template_raster), # 9000 meters
  crs = crs(template_raster) # Ensures grid and rasters match
)

# This is the key step: convert the grid cells to polygons for sampling
grid_polygons <- terra::as.polygons(grid_9km, dissolve = FALSE)

terra::writeVector(grid_polygons, "grid_polygons.shp")

# Define a path for the yearly output files
output_folder <- "~/Documents/mapbiomas_lulc_10/yearly_lulc_results/"

# Create the directory if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Use n - 1 cores to leave one for our system
registerDoFuture()
num_cores_to_use <- 40
plan(multisession, workers = num_cores_to_use)
cat("Future plan set to 'multisession' with", num_cores_to_use, "workers.\n")

# --- 2. Run the foreach Loop ---
# The result will be a list of the data.tables, just like with future_lapply
results_list <- foreach(
  raster_file = all_lulc_files,
  .packages = c("terra", "data.table", "exactextractr", "sf") # Ensure packages are loaded on each worker
) %dopar%
  {
    # --- The processing logic inside the loop is IDENTICAL to before ---

    # Load and project the grid
    grid_polygons <- vect("grid_polygons.shp") # Use the actual path
    grid_polygons$plot_id <- 1:nrow(grid_polygons)
    grid_polygons_sf <- st_as_sf(grid_polygons)

    # Load the LULC raster
    lulc_raster <- rast(raster_file)

    # --- Use the much simpler and faster built-in "frac" summary ---
    frac_coverage <- exact_extract(
      lulc_raster,
      grid_polygons_sf,
      fun = "frac",
      append_cols = "plot_id"
    )

    # --- Convert the wide output to a long data.table ---
    setDT(frac_coverage)

    long_dt <- melt(
      frac_coverage,
      id.vars = "plot_id",
      variable.name = "code_id",
      value.name = "fraction"
    )

    # Clean up the data
    long_dt[, code_id := as.numeric(gsub("frac_", "", code_id))]
    long_dt <- long_dt[fraction > 0] # Keep only existing classes

    # Add the year
    year <- as.numeric(gsub(".*-(\\d{4})\\.tif$", "\\1", basename(raster_file)))
    long_dt[, year := year]

    # Save the result for the current year
    output_filename <- file.path(
      output_folder,
      paste0("lulc_results_", year, ".csv")
    )
    fwrite(long_dt, file = output_filename)

    # This encourages the worker to release memory it no longer needs.
    gc()

    # Return a success message
    paste("Successfully processed and saved:", basename(raster_file))
  }

# Close the parallel session
plan(sequential)

# Get a list of all the saved yearly result files
saved_files <- list.files(output_folder, pattern = "\\.csv$", full.names = TRUE)

# Read all files and combine them into a single data.table
all_years_long <- rbindlist(lapply(saved_files, fread))

all_years_long <- all_years_long[code_id > 0] # Keep only existing classes

# Reshape from long to wide format
# Each row will be a plot_id-year, and columns will be the LULC classes
lulc_wide <- dcast(
  all_years_long,
  plot_id + year ~ code_id,
  value.var = "fraction",
  fill = 0
) # Fill with 0 for classes not present

# --- Calculate Area (e.g., in Hectares) ---
# MapBiomas has a 30m x 30m resolution. Area of one pixel = 900 sq meters.
# 1 Hectare = 10,000 sq meters.

# Calculate the area of each polygon in square meters
# Load the grid
grid_polygons <- vect("grid_polygons.shp") # Use the actual path
grid_polygons$plot_id <- 1:nrow(grid_polygons)
grid_polygons_sf <- st_as_sf(grid_polygons)

grid_polygons_sf$total_area_sqm <- st_area(grid_polygons_sf)

# Create a data.table with just the plot ID and area in hectares
plot_areas <- data.table(
  plot_id = grid_polygons_sf$plot_id,
  total_area_ha = as.numeric(grid_polygons_sf$total_area_sqm) / 10000
)

# Join the total area to our LULC fraction data
lulc_with_area <- plot_areas[lulc_wide, on = "plot_id"]

# Get the names of the columns that contain the class fractions
code_cols <- as.character(legend_dt$Class_ID)
code_cols_in_data <- intersect(code_cols, names(lulc_with_area))

# Loop through the fraction columns and multiply by the total area
for (col in code_cols_in_data) {
  lulc_with_area[, (col) := get(col) * total_area_ha]
}

summary(lulc_with_area)

# You can now remove the total_area_ha column if you wish
lulc_with_area[, total_area_ha := NULL]

# --- Use the legend to create clean, readable column names ---
# Create a mapping from code to clean name
name_map <- setNames(legend_dt$Description, legend_dt$Class_ID)

# Get the current names (which are codes) and the new names
current_names <- names(lulc_with_area)
new_names <- ifelse(
  current_names %in% names(name_map),
  name_map[current_names],
  current_names
)

# Set the new, clean names
setnames(lulc_with_area, old = current_names, new = new_names)

# our final data now has the correct area in hectares for each class
glimpse(lulc_with_area)
summary(lulc_with_area)

# Remove class without information
lulc_with_area[, `13` := NULL]

# Get the names of only the LULC columns
lulc_cols <- legend_dt$Description
lulc_cols_in_data <- intersect(lulc_cols, names(lulc_with_area))

# Calculate the number of non-zero observations for each class
frequency_counts <- lulc_with_area[,
  lapply(.SD, function(x) sum(x > 0, na.rm = TRUE)),
  .SDcols = lulc_cols_in_data
]

# Reshape the results into a clean table
frequency_table <- melt(
  frequency_counts,
  variable.name = "class_name",
  value.name = "non_zero_count"
)

# Calculate the percentage of total observations
total_obs <- nrow(lulc_with_area)
frequency_table[, percentage_presence := (non_zero_count / total_obs) * 100]

# Order by the least frequent classes
setorder(frequency_table, non_zero_count)

# View the rarest classes
print(frequency_table)

# Identify Persistently Present Classes
num_years <- uniqueN(lulc_with_area$year)
persistence_list <- list()

for (col in lulc_cols_in_data) {
  # For each plot, count how many years the class was present
  years_present_per_plot <- lulc_with_area[,
    .(years_present = sum(get(col) > 0)),
    by = plot_id
  ]

  # Count how many plots had the class present every single year
  persistent_plots_count <- sum(
    years_present_per_plot$years_present == num_years
  )

  persistence_list[[col]] <- data.table(
    class_name = col,
    persistent_plots_count = persistent_plots_count
  )
}

persistence_table <- rbindlist(persistence_list)
setorder(persistence_table, -persistent_plots_count)

print(persistence_table)

# Group by year and check if the sum of the area for each class is > 0
presence_by_year <- lulc_with_area[,
  lapply(.SD, function(x) sum(x, na.rm = TRUE) > 0),
  by = year,
  .SDcols = lulc_cols_in_data
]

# Order the results by year
setorder(presence_by_year, year)

# View the results
print(presence_by_year)
colSums(presence_by_year)

# Simplifying classes

# --- The Fix: Remove whitespace from all column names ---
setnames(
  lulc_with_area,
  old = names(lulc_with_area),
  new = trimws(names(lulc_with_area))
)

# Now, check the cleaned names
names(lulc_with_area)

# --- Define the groups of columns to aggregate ---
forest_cols <- c("Forest Formation", "Floodable Forest", "Mangrove")
non_veg_cols <- c(
  "Photovoltaic Power Plant (beta)",
  "Mining",
  "Beach, Dune and Sand Spot",
  "Other non Vegetated Areas",
  "Hypersaline Tidal Flat",
  "Rocky Outcrop"
)
temp_crop_cols <- c("Other Temporary Crops", "Cotton (beta)", "Rice")
perennial_crop_cols <- c("Other Perennial Crops", "Coffee", "Citrus")
water_cols <- c("River, Lake and Ocean", "Aquaculture")
wet_herbaceous_cols <- c("Wetland", "Herbaceous Sandbank Vegetation")


# --- Perform the aggregations using the robust rowSums(.SD) method ---

# 1. Aggregate Non-Vegetated Areas
lulc_with_area[,
  non_vegetated_areas := rowSums(.SD, na.rm = TRUE),
  .SDcols = non_veg_cols
]

# 2. Aggregate Forest Classes
lulc_with_area[, forest := rowSums(.SD, na.rm = TRUE), .SDcols = forest_cols]

# 3. Aggregate Temporary Crops
lulc_with_area[,
  temporary_crop := rowSums(.SD, na.rm = TRUE),
  .SDcols = temp_crop_cols
]

# 4. Aggregate Perennial Crops
lulc_with_area[,
  perennial_crop := rowSums(.SD, na.rm = TRUE),
  .SDcols = perennial_crop_cols
]

# 5. Aggregate Water Classes
lulc_with_area[, water := rowSums(.SD, na.rm = TRUE), .SDcols = water_cols]

# 6. Aggregate wet vegetation classes
lulc_with_area[,
  wet_herbaceous := rowSums(.SD, na.rm = TRUE),
  .SDcols = wet_herbaceous_cols
]

# --- Now, remove all the original columns at once ---
cols_to_remove <- c(
  non_veg_cols,
  forest_cols,
  temp_crop_cols,
  perennial_crop_cols,
  water_cols,
  wet_herbaceous_cols
)
lulc_with_area[, (cols_to_remove) := NULL]

# Check the new, simplified column names
names(lulc_with_area)

# Get the current names
current_names <- names(lulc_with_area)

# 1. Make all names lowercase
new_names <- tolower(current_names)

# 2. Replace all spaces with underscores
new_names <- gsub(" ", "_", new_names)

# 3. (Optional) Replace multiple underscores with a single one
new_names <- gsub("_{2,}", "_", new_names)

# Apply the new, clean names to the data.table
setnames(lulc_with_area, old = current_names, new = new_names)

# View the new, standardized names
print(names(lulc_with_area))

# --- Step 2: Apply custom shortenings ---
# Get the current names
current_names <- names(lulc_with_area)

# Apply a series of replacements
new_names <- gsub("_formation", "", current_names) # "savanna_formation" -> "savanna"
new_names <- gsub("_of_uses", "_uses", new_names) # "mosaic_of_uses" -> "mosaic_uses"
new_names <- gsub("_areas", "", new_names) # "non_vegetated_areas" -> "non_vegetated"
new_names <- gsub("_area", "", new_names) # "urban_area" -> "urban"

# Apply the new, simplified names to the data.table
setnames(lulc_with_area, old = current_names, new = new_names)

# View the final, simplified names
print(names(lulc_with_area))

# --- Define which classes belong to each meta-category ---
native_classes <- c("forest", "savanna", "wet_herbaceous", "grassland", "water")

anthropic_classes <- c(
  "pasture",
  "temporary_crop",
  "perennial_crop",
  "soybean",
  "sugar_cane",
  "forest_plantation",
  "mosaic_uses",
  "non_vegetated",
  "urban"
)

# Make sure all defined classes are actually in the dataset
native_classes <- intersect(native_classes, names(lulc_with_area))
anthropic_classes <- intersect(anthropic_classes, names(lulc_with_area))

# --- Calculate the total area for each meta-category ---
lulc_with_area[,
  native_area := rowSums(.SD, na.rm = TRUE),
  .SDcols = native_classes
]
lulc_with_area[,
  anthropic_area := rowSums(.SD, na.rm = TRUE),
  .SDcols = anthropic_classes
]

# View the final result with the two new summary columns
glimpse(lulc_with_area[, .(plot_id, year, native_area, anthropic_area)])

summary(lulc_with_area)

saveRDS(
  lulc_with_area,
  "~/Documents/mapbiomas_lulc_10/yearly_lulc_results/lulc_predictors.rds"
)


# NDVI --------------------------------------------------------------------
# Landsat Collection 2 Tier 1 Level 2 Annual NDVI Composite courtesy of the U.S. Geological Survey
# https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_COMPOSITES_C02_T1_L2_ANNUAL_NDVI
data_directory <- "~/Documents/Cerrado_NDVI"
all_ndvi_files <- list.files(
  data_directory,
  pattern = "\\.tif$",
  full.names = TRUE
)
ndvi_timeseries_terra <- rast(all_ndvi_files)
print(ndvi_timeseries_terra)

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

# Load the Cerrado geometry
Cerrado <- readRDS("Cerrado.rds")
print(Cerrado)

# Select only the geometry object from "Cerrado"
Cerrado <- Cerrado %>% dplyr::select(geometry)
Cerrado <- sf::st_transform(Cerrado, metric_crs)
template_raster <- crop(template_raster, Cerrado) %>% mask(Cerrado)


grid_9km <- terra::rast(
  extent = ext(template_raster),
  resolution = res(template_raster), # 9000 meters
  crs = crs(template_raster) # Ensures grid and rasters match
)

print(grid_9km)

# Project and resample to match with grid cells
ndvi_timeseries_terra <- project(ndvi_timeseries_terra, grid_9km, threads = T)
ndvi_timeseries_terra <- resample(
  ndvi_timeseries_terra,
  grid_9km,
  method = "bilinear",
  threads = T
)
print(ndvi_timeseries_terra)

# Rename layers
names(ndvi_timeseries_terra) <- 1985:2024

# Extract data to data.frame
ndvi_timeseries_df <- as.data.frame(
  ndvi_timeseries_terra,
  cells = TRUE,
  xy = TRUE,
  wide = FALSE
)

# Rename columns
names(ndvi_timeseries_df) <- c(
  "plot_id",
  "longitude",
  "latitude",
  "year",
  "NDVI"
)

# Check data
summary(ndvi_timeseries_df)

# Save data
saveRDS(ndvi_timeseries_df, "~/Documents/Cerrado_NDVI/ndvi_predictor.rds")


# Population density ------------------------------------------------------
# Package popgrid
# https://github.com/Nowosad/popgrids
# Load the Cerrado geometry
pop_grid_terra <- rast(pop_grid)
print(pop_grid_terra)

names(pop_grid_terra)
Cerrado <- readRDS("Cerrado.rds")
print(Cerrado)

# Select only the geometry object from "Cerrado"
Cerrado <- Cerrado %>% dplyr::select(geometry)
Cerrado <- sf::st_transform(Cerrado, crs(pop_grid_terra))

pop_grid_cerrado <- crop(pop_grid_terra, Cerrado) %>% mask(Cerrado)
plot(pop_grid_cerrado)

# --- Average 2020 Scenarios ---
# Find all layers that correspond to the year 2020
layers_2020 <- names(pop_grid_cerrado)[str_detect(
  names(pop_grid_cerrado),
  "_2020"
)]
# Calculate the mean of these layers
pop_2020_avg <- mean(pop_grid_cerrado[[layers_2020]])
names(pop_2020_avg) <- "p_2020"

# --- Average 2030 Scenarios ---
# Find all layers that correspond to the year 2030
layers_2030 <- names(pop_grid_cerrado)[str_detect(
  names(pop_grid_cerrado),
  "_2030"
)]
# Calculate the mean of these layers
pop_2030_avg <- mean(pop_grid_cerrado[[layers_2030]])
names(pop_2030_avg) <- "p_2030"

# Create a new stack with the clean, single-year rasters
pop_stack_interp <- c(
  pop_grid_cerrado[["p_1980"]],
  pop_grid_cerrado[["p_1990"]],
  pop_grid_cerrado[["p_2000"]],
  pop_grid_cerrado[["p_2010"]],
  pop_2020_avg,
  pop_2030_avg
)

# Create a list to store the final annual rasters
annual_pop_rasters <- list()

# Get the years from the layer names
key_years <- as.numeric(gsub("p_", "", names(pop_stack_interp)))

# Loop through each time interval to interpolate
for (i in 1:(length(key_years) - 1)) {
  start_year <- key_years[i]
  end_year <- key_years[i + 1]

  start_raster <- pop_stack_interp[[i]]
  end_raster <- pop_stack_interp[[i + 1]]

  # Calculate the total change over the period and the annual change
  interval <- end_year - start_year
  annual_delta <- (end_raster - start_raster) / interval

  # Generate the raster for each year in the interval
  for (year in start_year:(end_year - 1)) {
    step <- year - start_year
    interpolated_raster <- start_raster + (annual_delta * step)
    names(interpolated_raster) <- as.character(year)
    annual_pop_rasters[[as.character(year)]] <- interpolated_raster
  }
}

# Add the very last year (2030)
annual_pop_rasters[["2030"]] <- pop_stack_interp[["p_2030"]]

# Combine the list into a final, multi-layer SpatRaster
final_pop_stack <- rast(annual_pop_rasters)

# Check our final annual stack
print(final_pop_stack)

print(grid_9km)

# Project and resample to match with grid cells
pop_timeseries_terra <- project(final_pop_stack, grid_9km, threads = T)
pop_timeseries_terra <- resample(
  pop_timeseries_terra,
  grid_9km,
  method = "bilinear",
  threads = T
)
print(pop_timeseries_terra)
names(pop_timeseries_terra)
# Rename layers
pop_timeseries_terra <- pop_timeseries_terra[[as.character(1985:2024)]]

# Extract data to data.frame
pop_timeseries_df <- as.data.frame(
  pop_timeseries_terra,
  cells = TRUE,
  xy = TRUE,
  wide = FALSE
)

# Rename columns
names(pop_timeseries_df) <- c("plot_id", "longitude", "latitude", "year", "pop")

# Check data
summary(pop_timeseries_df)

# Save data
saveRDS(pop_timeseries_df, "Modelling/pop_predictor.rds")


# Climate -----------------------------------------------------------------
# Muñoz Sabater, J. (2019): ERA5-Land monthly averaged data from 1950 to present.
# Copernicus Climate Change Service (C3S) Climate Data Store (CDS).
# DOI: 10.24381/cds.68d2bb30 (Accessed on 08-July-2025)
# Read variables
dewpoint_temp <- rast("~/Documents/ERA5/2m_dewpoint_temp.grib")
temp <- rast("~/Documents/ERA5/2m_temp.nc")
pot_evap <- rast("~/Documents/ERA5/pot_evap.grib")
precip <- rast("~/Documents/ERA5/precip.grib")
skin_temp <- rast("~/Documents/ERA5/skin_temp.grib")
solrad <- rast("~/Documents/ERA5/solrad.grib")
total_evap <- rast("~/Documents/ERA5/total_evap.grib")
u_wind <- rast("~/Documents/ERA5/u_wind.grib")
v_wind <- rast("~/Documents/ERA5/v_wind.grib")
water_soil <- rast("~/Documents/ERA5/water_soil_lvl1.grib")

time(temp) <- time(dewpoint_temp)

# Derived variables
# Average wind speed
wind <- (abs(v_wind) + abs(u_wind)) / 2
print(wind)

# Water deficit (total_evap - pot_evap)
water_deficit <- (total_evap - pot_evap)
print(water_deficit)

# Project the rasters
# SIRGAS 2000 / Brazil Polyconic (EPSG:5880) is a good choice for the whole country.
metric_crs <- "EPSG:5880"

dewpoint_temp <- terra::project(
  dewpoint_temp,
  metric_crs,
  threads = T,
  method = "bilinear"
)
temp <- terra::project(temp, metric_crs, threads = T, method = "bilinear")
pot_evap <- terra::project(
  pot_evap,
  metric_crs,
  threads = T,
  method = "bilinear"
)
precip <- terra::project(precip, metric_crs, threads = T, method = "bilinear")
skin_temp <- terra::project(
  skin_temp,
  metric_crs,
  threads = T,
  method = "bilinear"
)
solrad <- terra::project(solrad, metric_crs, threads = T, method = "bilinear")
total_evap <- terra::project(
  total_evap,
  metric_crs,
  threads = T,
  method = "bilinear"
)
wind <- terra::project(wind, metric_crs, threads = T, method = "bilinear")
water_soil <- terra::project(
  water_soil,
  metric_crs,
  threads = T,
  method = "bilinear"
)
water_deficit <- terra::project(
  water_deficit,
  metric_crs,
  threads = T,
  method = "bilinear"
)

# Crop and mask for Cerrado

# Load the Cerrado geometry
Cerrado <- readRDS("Cerrado.rds")
print(Cerrado)

# Select only the geometry object from "Cerrado"
Cerrado <- Cerrado %>% dplyr::select(geometry)
Cerrado <- sf::st_transform(Cerrado, metric_crs)

dewpoint_temp_cerrado <- crop(dewpoint_temp, Cerrado) %>% mask(Cerrado)
temp_cerrado <- crop(temp, Cerrado) %>% mask(Cerrado)
pot_evap_cerrado <- crop(pot_evap, Cerrado) %>% mask(Cerrado)
precip_cerrado <- crop(precip, Cerrado) %>% mask(Cerrado)
skin_temp_cerrado <- crop(skin_temp, Cerrado) %>% mask(Cerrado)
solrad_cerrado <- crop(solrad, Cerrado) %>% mask(Cerrado)
total_evap_cerrado <- crop(total_evap, Cerrado) %>% mask(Cerrado)
wind_cerrado <- crop(wind, Cerrado) %>% mask(Cerrado)
water_soil_cerrado <- crop(water_soil, Cerrado) %>% mask(Cerrado)
water_deficit_cerrado <- crop(water_deficit, Cerrado) %>% mask(Cerrado)

# Standardize wind with other variables
wind_cerrado <- project(wind_cerrado, temp_cerrado, threads = T)

# Create a list of the raster objects to process
# Name the list elements for easy tracking
raster_list <- list(
  dewpoint_temp = dewpoint_temp_cerrado,
  temp = temp_cerrado,
  pot_evap = pot_evap_cerrado,
  precip = precip_cerrado,
  skin_temp = skin_temp_cerrado,
  solrad = solrad_cerrado,
  total_evap = total_evap_cerrado,
  wind = wind_cerrado,
  water_deficit = water_deficit_cerrado,
  water_soil = water_soil_cerrado
)

output_folder <- "~/Documents/ERA5/climate_variable_results/"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

all_variable_summaries <- list()

# Loop through each climate variable
for (var_name in names(raster_list)) {
  cat("Processing and summarizing variable:", var_name, "\n")

  current_raster <- raster_list[[var_name]]

  # Convert raster directly to a LONG data.table
  # This is the key simplification. It's fast and avoids all name/shape issues.
  var_dt <- as.data.frame(
    current_raster,
    cells = TRUE,
    time = TRUE,
    xy = TRUE,
    wide = FALSE
  )
  setDT(var_dt)

  # Add metadata and calculate summaries ---
  # The columns are already named 'cell', 'time', and the variable name.
  setnames(var_dt, old = "values", new = "value") # Rename the value column

  # Add variable name and extract year/month
  var_dt[, variable := var_name]
  var_dt[, year := year(time)]
  var_dt[, month := month(time)]

  # Calculate monthly summaries
  monthly_summary <- var_dt[,
    .(
      mean = mean(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE)
    ),
    by = .(cell, x, y, variable, year, month)
  ]

  all_variable_summaries[[var_name]] <- monthly_summary

  rm(var_dt, monthly_summary)
  gc()
}

saveRDS(
  all_variable_summaries,
  "~/Documents/ERA5/climate_variable_results/list_climate_predictors.rds"
)
lapply(all_variable_summaries, nrow)
lapply(all_variable_summaries, summary)

# Prepare each table for the merge ---
# We will loop through the list and rename the summary columns to be unique
# (e.g., 'mean' becomes 'temp_mean')
prepared_list <- list()
for (var_name in names(all_variable_summaries)) {
  # Get the summary table
  dt <- all_variable_summaries[[var_name]]

  # Rename the summary columns
  setnames(
    dt,
    old = c("mean", "min", "max", "sd"),
    new = paste0(var_name, "_", c("mean", "min", "max", "sd"))
  )

  # We only need the key columns and the newly named summary columns
  key_cols <- c("cell", "x", "y", "year", "month")
  summary_cols <- paste0(var_name, "_", c("mean", "min", "max", "sd"))

  prepared_list[[var_name]] <- dt[, c(key_cols, summary_cols), with = FALSE]
}


# Perform a sequential full outer join
# The Reduce() function will merge the tables one by one.
# all = FALSE ensures that we keep only rows without NAs.

final_climate_dt_complete <- Reduce(
  function(x, y) {
    merge(x, y, by = c("cell", "x", "y", "year", "month"), all = FALSE)
  },
  prepared_list
)

# Rename 'cell' to 'plot_id'
setnames(final_climate_dt_complete, "cell", "plot_id")

# Check the summary of the new, complete dataset
# There should be far fewer (or zero) NAs
summary(final_climate_dt_complete)

# Save results
saveRDS(
  final_climate_dt_complete,
  "~/Documents/ERA5/climate_variable_results/climate_predictors.rds"
)


# Terrain variables -------------------------------------------------------

# Download Elevation Data for the Cerrado ---
# We'll use elevatr to get a digital elevation model (DEM).
# 'z=9' provides a good resolution (~500m) for regional analysis.
# We'll project it directly to our target CRS.
# Load the Cerrado geometry
Cerrado <- readRDS("Cerrado.rds")
print(Cerrado)

# Select only the geometry object from "Cerrado"
Cerrado <- Cerrado %>% dplyr::select(geometry)
Cerrado <- sf::st_transform(Cerrado, metric_crs)
elevation_dem <- get_elev_raster(
  locations = Cerrado,
  z = 9,
  prj = "EPSG:5880",
  clip = "locations"
)

# Calculate Terrain Derivatives
# The terra::terrain() function can calculate everything at once.
terrain_vars <- terrain(
  rast(elevation_dem),
  v = c("slope", "aspect", "TPI", "TRI"),
  unit = "degrees"
)
terrain_vars <- c(rast(elevation_dem), terrain_vars)

names(terrain_vars) <- c("elevation", "slope", "aspect", "TPI", "TRI")

# Align Rasters with our Climate Data Grid
# This is a crucial step to ensure the grids match perfectly.
# We'll resample the terrain rasters to match a template from our climate data.
terrain_aligned <- resample(terrain_vars, temp_cerrado, method = "bilinear")


# Extract the Data
# Now we can use the same as.data.frame() method to extract the terrain data.
terrain_dt <- as.data.frame(terrain_aligned, cells = TRUE, xy = TRUE)
setDT(terrain_dt)
setnames(terrain_dt, "cell", "plot_id")

# View our new terrain predictors, ready to be merged
glimpse(terrain_dt)

# Save results
saveRDS(terrain_dt, "~/Documents/Terrain/terrain_predictors.rds")
