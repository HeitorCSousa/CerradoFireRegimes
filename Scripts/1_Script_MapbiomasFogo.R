# LOAD PACKAGES ---------------------------------------------

library(terra)
library(stars)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)
library(landscapemetrics)
library(usdm)
library(tidyr) # For pivot_wider
library(tidyterra)
library(purrr)
library(data.table)
library(ggplot2)
library(viridis)
library(visdat)
library(ggfortify)
library(corrplot)
# terraOptions(memfrac = 0.9)

# SET UP ENVIRONMENT AND LIST FILES --------------------------------------

data_directory <- "~/Documents/mapbiomas_fogo_4"

# Create a directory to store the merged annual rasters
# This is where we will save our intermediate results
output_directory <- file.path(
  "~/Documents/mapbiomas_fogo_4",
  "MERGED_ANNUAL_RASTERS"
)
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

all_files <- list.files(
  path = data_directory,
  pattern = "\\.tif$",
  full.names = TRUE
)

files_df <- data.frame(filepath = all_files) %>%
  mutate(year = str_extract(basename(filepath), "(?<=-)\\d{4}(?=-)"))


# BATCH PARALLEL PROCESSING (SAVE TO DISK) -------------------------------

unique_years <- unique(files_df$year)

# Define the batch limit
# This is the maximum number of years to process in parallel at one time.
# Start with a conservative number like 8 or 10. You can increase it if you have RAM to spare.
batch_size <- 8

# Create a list of years that still need to be processed
years_to_process <- unique_years[
  !sapply(unique_years, function(year) {
    file.exists(file.path(
      output_directory,
      paste0("merged_fire_", year, ".tif")
    ))
  })
]

if (length(years_to_process) == 0) {
  cat("All years have already been processed!\n")
} else {
  cat(length(years_to_process), "years remaining to process.\n")

  # Create the batches
  batches <- split(
    years_to_process,
    ceiling(seq_along(years_to_process) / batch_size)
  )

  # Set up the parallel cluster
  num_cores <- min(detectCores() - 1, batch_size)
  registerDoParallel(cores = num_cores)
  cat(
    "Processing in batches of up to",
    batch_size,
    "using",
    num_cores,
    "cores.\n"
  )

  # Loop through each batch
  for (i in seq_along(batches)) {
    current_batch <- batches[[i]]
    cat(
      "--- Starting Batch",
      i,
      "of",
      length(batches),
      "--- Years:",
      paste(current_batch, collapse = ", "),
      "---\n"
    )

    # Run foreach on the current small batch
    foreach(
      current_year = current_batch,
      .packages = c("terra", "dplyr", "stringr")
    ) %dopar%
      {
        output_filename <- file.path(
          output_directory,
          paste0("merged_fire_", current_year, ".tif")
        )

        # This check is still useful in case of a crash mid-batch
        if (file.exists(output_filename)) {
          return(NULL)
        }

        cat("Processing year:", current_year, "\n")

        files_for_year <- files_df %>%
          filter(year == current_year) %>%
          pull(filepath)

        # --- More Memory-Efficient Method: VRT ---
        # 1. Create a Virtual Raster (VRT). This is just a pointer file, very low memory.
        virtual_raster <- vrt(files_for_year)

        # 2. Write the VRT to a real file. Terra processes it chunk by chunk.
        writeRaster(virtual_raster, output_filename, overwrite = TRUE)

        return(NULL) # We don't need to return anything
      }
  }

  stopImplicitCluster()
  cat("--- All batches completed! ---\n")
}


# STACK RESULTS FROM DISK ------------------------------------------------

# List all the generated files and stack them
all_output_files <- list.files(
  output_directory,
  pattern = "\\.tif$",
  full.names = TRUE
)
fire_timeseries_terra <- rast(all_output_files)

# Extract years from filenames to ensure correct naming
output_years <- str_extract(basename(all_output_files), "\\d{4}")
names(fire_timeseries_terra) <- output_years

print(fire_timeseries_terra)
plot(fire_timeseries_terra[[1]])

# Load Cerrado geometry
Cerrado <- readRDS("Data/Cerrado.rds")
print(Cerrado)

# Select only Cerrado geometry
Cerrado <- Cerrado %>% dplyr::select(geometry)
Cerrado <- sf::st_transform(Cerrado, crs(fire_timeseries_terra))

fire_timeseries_terra <- crop(fire_timeseries_terra, Cerrado) %>% mask(Cerrado)


# GENERATE BINARY RASTERS ------------------------------------------------

# Assume 'fire_timeseries_terra' is your original SpatRaster

# --- Define a metric projection for Brazil ---
# SIRGAS 2000 / Brazil Polyconic (EPSG:5880) is a good choice for the whole country.
metric_crs <- "EPSG:5880"

# --- Reproject the entire time series to the metric CRS ---
# This may take some time.
cat("Reprojecting time series to metric CRS:", metric_crs, "...\n")
fire_timeseries_metric <- terra::project(
  fire_timeseries_terra,
  metric_crs,
  threads = T,
  method = "near"
)

writeRaster(
  fire_timeseries_metric,
  paste0(
    file.path(
      "~/Documents/mapbiomas_fogo_4",
      "MERGED_ANNUAL_RASTERS",
      "fire_timeseries/"
    ),
    "fire_timeseries_metric.tif"
  ),
  datatype = "INT1U",
  overwrite = T
)

fire_timeseries_metric <- rast(paste0(
  file.path(
    "~/Documents/mapbiomas_fogo_4",
    "MERGED_ANNUAL_RASTERS",
    "fire_timeseries/"
  ),
  "fire_timeseries_metric.tif"
))


# --- Define a small extent to sample from ---
# Get the full extent of the raster
full_ext <- ext(fire_timeseries_metric)

# Calculate the center point
center_x <- (full_ext$xmin + full_ext$xmax) / 2
center_y <- (full_ext$ymin + full_ext$ymax) / 2

# Define a 100km x 100km box around the center
# Note: Using a projected CRS, the units are in meters
sample_extent <- ext(
  center_x - 50000,
  center_x + 50000,
  center_y - 50000,
  center_y + 50000
)

# --- Crop and Extract Values (this is very fast) ---
cat("Cropping a sample block from the raster...\n")
# Crop doesn't read all the data, it just defines the area
raster_block <- crop(fire_timeseries_metric[[1]], sample_extent)

# Now, get all the values from just this small block
cat("Extracting values from the block...\n")
block_values <- values(raster_block)

# --- Check if the values are integers ---
# We do this by seeing if they are equal to themselves when rounded
are_values_integers <- all(block_values == round(block_values), na.rm = TRUE)

if (are_values_integers) {
  cat("The sample block contains only integer values.\n")
} else {
  cat("The sample block contains non-integer (numeric) values.\n")
}

# Load Cerrado geometry
Cerrado <- readRDS("Data/Cerrado.rds")
print(Cerrado)

# Select only Cerrado geometry
Cerrado <- Cerrado %>% dplyr::select(geometry)
Cerrado <- sf::st_transform(Cerrado, crs(fire_timeseries_metric))


## SETUP ------------------------------------------------------------------

# Assume 'fire_timeseries_metric' is your reprojected SpatRaster from the previous steps
# and 'metric_crs' is defined (e.g., "EPSG:5880")

# Directory for the output files
binary_output_dir <- file.path(data_directory, "BINARY_MONTHLY_RASTERS")
if (!dir.exists(binary_output_dir)) {
  dir.create(binary_output_dir, recursive = TRUE)
}


## CREATE AND FILTER TASK LIST --------------------------------------------

# Create a data frame of all possible year-month tasks
all_tasks <- expand.grid(
  year = names(fire_timeseries_metric),
  month = 1:12,
  stringsAsFactors = FALSE
)

# Add a column for the expected output filename for each task
all_tasks <- all_tasks %>%
  mutate(
    binary_filename = file.path(
      binary_output_dir,
      paste0("binary_fire_", year, "_", sprintf("%02d", month), ".tif")
    )
  )

# Filter out tasks that have already been completed
tasks_to_process <- all_tasks %>%
  filter(!file.exists(binary_filename))

if (nrow(tasks_to_process) == 0) {
  cat("All binary rasters have already been generated.\n")
} else {
  cat(nrow(tasks_to_process), "binary rasters remaining to be generated.\n")

  # PARALLEL BATCH PROCESSING ----------------------------------------------

  # Define the batch limit (how many files to process in parallel at once)
  batch_size <- 5

  # Create the list of batches
  batches <- split(
    tasks_to_process,
    ceiling(seq_len(nrow(tasks_to_process)) / batch_size)
  )

  # Set up the parallel cluster
  num_cores <- min(detectCores() - 1, batch_size)
  registerDoParallel(cores = num_cores)
  cat("Processing in", length(batches), "batches using", num_cores, "cores.\n")

  # Loop through each batch
  for (i in seq_along(batches)) {
    current_batch <- batches[[i]]
    cat("--- Starting Batch", i, "of", length(batches), "---\n")

    # The '.packages' argument loads packages on each worker.
    # The 'foreach' loop will process rows of the 'current_batch' data frame.
    foreach(
      j = 1:nrow(current_batch),
      .packages = c("terra")
    ) %dopar%
      {
        # Get task details for the current iteration
        current_year <- current_batch$year[j]
        current_month <- current_batch$month[j]
        binary_filename <- current_batch$binary_filename[j]

        # Select the raster for the current year
        year_raster <- fire_timeseries_metric[[current_year]]

        # Create the binary raster (1 for burn, 0 for no burn)
        binary_raster <- terra::ifel(year_raster == current_month, 1, 0)

        # Only write the file if there is actually burned area
        writeRaster(binary_raster, binary_filename, overwrite = TRUE)

        # We don't need to return anything from the parallel workers
        return(NULL)
      } # End of foreach
  } # End of for loop for batches

  # Stop the parallel cluster to free up resources
  stopImplicitCluster()
  cat("--- All batches completed! ---\n")
}

cat("Part 1 Complete: All binary rasters have been generated and saved.\n")


# CALCULATE METRICS ------------------------------------------------------

options_landscapemetrics(to_disk = TRUE)

# --- Get the list of binary files created in Part 1 ---
binary_output_dir <- file.path(data_directory, "BINARY_MONTHLY_RASTERS")
binary_files <- list.files(
  binary_output_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

# --- Create the 9km Analysis Grid (in the same metric CRS) ---
cat("Creating 9km analysis grid...\n")
# Load one of the reprojected binary rasters to use as a template
template_raster <- rast("~/Documents/ERA5/precip.grib")[[1]]
template_raster <- terra::project(
  template_raster,
  metric_crs,
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

# --- SETUP ---
binary_input_dir <- binary_output_dir
rds_output_dir <- file.path(binary_input_dir, "METRICS_RDS_FILES")
if (!dir.exists(rds_output_dir)) {
  dir.create(rds_output_dir)
}

all_metrics <- list_lsm()

selected_metrics <- all_metrics %>%
  filter(
    # Class-level metrics
    level == "class" |
      # Landscape-level metrics
      (level == "landscape" &
        (type == "diversity metric" | type == "complexity metric"))
  ) %>%
  pull(function_name)

cat("Selected metrics to calculate:\n")
print(selected_metrics)
# Assume 'grid_polygons' is your SpatVector of 9km grid cells

# --- CREATE AND FILTER TASK LIST ---
# Get all binary files
all_files <- list.files(
  binary_input_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

# Create a list of files that still need to be processed
files_to_process <- Filter(
  function(file_path) {
    rds_filename <- file.path(
      rds_output_dir,
      paste0(tools::file_path_sans_ext(basename(file_path)), ".rds")
    )
    !file.exists(rds_filename)
  },
  all_files
)

if (length(files_to_process) == 0) {
  cat("All metrics have already been calculated and saved.\n")
} else {
  cat(length(files_to_process), "files remaining to be processed.\n")

  # --- PARALLEL BATCH PROCESSING ---
  batch_size <- 20
  batches <- split(
    files_to_process,
    ceiling(seq_along(files_to_process) / batch_size)
  )

  num_cores <- min(detectCores() - 2, batch_size)
  registerDoParallel(cores = num_cores)
  cat("Processing in", length(batches), "batches using", num_cores, "cores.\n")

  for (i in seq_along(batches)) {
    current_batch <- batches[[i]]
    cat("--- Starting Batch", i, "of", length(batches), "---\n")

    foreach(
      file_path = current_batch,
      .packages = c("terra", "landscapemetrics")
    ) %dopar%
      {
        filename <- basename(file_path)
        rds_filename <- file.path(
          rds_output_dir,
          paste0(tools::file_path_sans_ext(filename), ".rds")
        )

        # Final check inside the loop in case of restarts
        if (file.exists(rds_filename)) {
          return(NULL)
        }

        binary_raster <- rast(file_path)

        metrics <- sample_lsm(
          landscape = binary_raster,
          y = grid_polygons,
          level = c("class", "landscape"),
          what = selected_metrics,
          transform = FALSE,
          progress = TRUE
        )

        if (nrow(metrics) > 0) {
          current_year <- substr(filename, 13, 16)
          current_month <- as.integer(substr(filename, 18, 19))
          metrics$year <- current_year
          metrics$month <- current_month

          # Save the result for this file to its own RDS file
          saveRDS(metrics, rds_filename)
        }

        return(NULL)
      } # End of foreach
  } # End of for loop for batches

  stopImplicitCluster()
  cat("--- All batches completed! ---\n")
}

# --- Point to the directories ---
rds_output_dir <- file.path(
  "~/Documents/mapbiomas_fogo_4/BINARY_MONTHLY_RASTERS/METRICS_RDS_FILES"
)
# Assume 'grid_polygons' SpatVector from the previous script is available

# --- Read and Combine Existing RDS Files ---
rds_files <- list.files(rds_output_dir, pattern = "\\.rds$", full.names = TRUE)
# Using map_dfr is a robust way to read and bind the files
all_metrics_df <- map_dfr(rds_files, readRDS)

selected_metrics <- all_metrics %>%
  filter(
    # Class-level metrics
    level == "class" |
      # Landscape-level metrics
      (level == "landscape" &
        (type == "diversity metric" | type == "complexity metric"))
  )

selected_metrics$class <- ifelse(selected_metrics$level == "class", 1, NA)

glimpse(all_metrics_df)

setDT(all_metrics_df)
setDT(selected_metrics)

# --- The Filtering Join ---
# This tells data.table to keep only the rows in all_metrics_df
# where the combination of 'level' and 'metric' exists in selected_metrics.
filtered_metrics_dt <- all_metrics_df[
  selected_metrics,
  on = c("level", "metric", "class"),
  nomatch = 0
]

dplyr::distinct(filtered_metrics_dt, level)
dplyr::distinct(filtered_metrics_dt, class)

rm(all_metrics_df)
gc()

# 1. Reshape the data (this step can take some hours)
metrics_wide_sparse <- dcast(
  filtered_metrics_dt,
  plot_id + class + year + month ~ function_name,
  value.var = "value",
  fun.aggregate = mean
)

saveRDS(
  metrics_wide_sparse,
  "~/Documents/mapbiomas_fogo_4/BINARY_MONTHLY_RASTERS/METRICS_RDS_FILES/FILTERED/metrics_wide_sparse.rds"
)

metrics_wide_sparse <- readRDS(
  "~/Documents/mapbiomas_fogo_4/BINARY_MONTHLY_RASTERS/METRICS_RDS_FILES/FILTERED/metrics_wide_sparse.rds"
)

setDT(metrics_wide_sparse)

# Find all rows where the calculation likely failed
error_months <- filtered_metrics_dt[percentage_inside == 0]

# See the list of failed calculations
if (nrow(error_months) > 0) {
  print("Found potential errors in the following plot/date combinations:")
  print(head(error_months[, .(plot_id, year, month)]))
}

# 2. Create a unique lookup table for percentage_inside
#    This table will have one row for each plot_id
percentage_lookup_original <- unique(filtered_metrics_dt[, .(
  plot_id,
  percentage_inside
)])

# Find the plot_ids that have more than one row
percentage_lookup_original[, .N, by = plot_id][N > 1]

unique(percentage_lookup_original$plot_id)
View(error_months[which(plot_id == unique(percentage_lookup_original$plot_id))])

# Remove rows where percentage_inside is 0
# The, group by plot_id and get the mean.
percentage_lookup <- filtered_metrics_dt[
  percentage_inside > 0,
  .(percentage_inside = mean(percentage_inside, na.rm = TRUE)),
  by = plot_id
]

# 3. Join the percentage_inside column back to the wide dataframe
#    This efficiently adds the column using the plot_id as a key
final_metrics <- metrics_wide_sparse[percentage_lookup, on = "plot_id"]

# Check the final result, which now includes the new column
glimpse(final_metrics)

# Landscape metrics are in rows where 'class' is NA
landscape_data <- final_metrics[is.na(class)]
glimpse(landscape_data)

# Class metrics are in rows where 'class' is a number
class_data <- final_metrics[!is.na(class)]
glimpse(class_data)

key_cols <- c("plot_id", "year", "month")

# Identify which columns are ONLY class-level metrics
class_metric_cols <- names(class_data)[str_starts(names(class_data), "lsm_c_")]

class_data_to_merge <- class_data[,
  c(key_cols, class_metric_cols),
  with = FALSE
]

# Identify which columns are ONLY landscape-level metrics
landscape_metric_cols <- names(landscape_data)[str_starts(
  names(landscape_data),
  "lsm_l_"
)]

landscape_data_to_merge <- landscape_data[,
  c(key_cols, landscape_metric_cols),
  with = FALSE
]

merged_data <- landscape_data_to_merge[class_data_to_merge, on = key_cols]

glimpse(merged_data)

grid_polygons$plot_id <- 1:nrow(grid_polygons)

# Get all unique plot_ids from your grid
all_plot_ids <- as.data.table(grid_polygons[, "plot_id"])

# Calculate the centroids of the grid cells
centroids_vec <- centroids(grid_polygons)

# Extract the longitude(x) and latitude(y) coordinates
coords <- as.data.table(crds(centroids_vec))
setnames(coords, c("x", "y"), c("longitude", "latitude"))

# Combine Ids and coordinates into one table
plot_locations_dt <- cbind(all_plot_ids, coords)

# Merge coordinates with the final metrics table
final_data <- plot_locations_dt[merged_data, on = "plot_id"]

# Check
glimpse(final_data)
summary(final_data)

# Clean up memory
rm(filtered_metrics_dt, metrics_wide_sparse, percentage_lookup)
gc()

# Check if there are missing points

# Get the complete list of all plot_ids that should exist
all_ids <- plot_locations_dt$plot_id

# Get all the years and months that were processed
all_years <- unique(final_data$year)
all_months <- unique(final_data$month)

# Create a table with every possible combination of id, year, and month
master_grid <- CJ(plot_id = all_ids, year = all_years, month = all_months)

# Find the missing combinations (Anti-Join)
# Check which rows from our 'master_grid' are NOT present in our final data
actual_data_keys <- final_data[, .(plot_id, year, month)]

missing_combinations <- master_grid[
  !actual_data_keys,
  on = c("plot_id", "year", "month")
]

# Report the results
if (nrow(missing_combinations) == 0) {
  print("Success! The dataset is complete. No missing combinations found.")
} else {
  print(paste(
    "Found",
    nrow(missing_combinations),
    "missing plot/date combinations."
  ))
  print("Here are the first few missing data points:")
  print(head(missing_combinations))
}

unique(missing_combinations$plot_id)

ggplot(
  final_data[final_data$year == 2023 & final_data$month == 9, ],
  aes(x = longitude, y = latitude, color = lsm_c_lpi)
) +
  geom_point(shape = 15, size = 1.5) +
  scale_color_viridis_c(name = "LPI value") +
  labs(
    title = "Spatial Distribution of Largest Patch Index (LPI) - September 2023",
    x = "Longitude (m)",
    y = "Latitude (m)"
  ) +
  theme_minimal() +
  coord_equal()

summary(final_data)
saveRDS(
  final_data,
  "~/Documents/mapbiomas_fogo_4/BINARY_MONTHLY_RASTERS/METRICS_RDS_FILES/FILTERED/final_metrics_data.rds"
)


# Quantify missing data
final_data %>%
  slice(1000000:1001000) %>%
  vis_miss(cluster = TRUE)

setDT(final_data)

# Identify the columns that are sd or cv metrics
variation_cols <- names(final_data)[str_detect(names(final_data), "_sd$|_cv$")]

# Loop through these columns and replace NAs with 0
for (col in variation_cols) {
  final_data[is.na(get(col)), (col) := 0]
}

summary(final_data)

na_percentages <- final_data[,
  lapply(.SD, function(x) sum(is.na(x) | is.nan(x)) / .N * 100),
  .SDcols = selected_metrics$function_name
]


sort(na_percentages, decreasing = TRUE)

na_vector <- unlist(na_percentages)

# Remove unreliable variables and with no variation (pr and prd)
vars_to_remove <- c(
  names(na_vector[na_vector == 100 | na_vector > 5]),
  "lsm_l_pr",
  "lsm_l_prd"
)

print("Removing the following unreliable variables:")
print(vars_to_remove)

# Create a new, clean vector of predictor names
final_metric_vars <- c(
  "plot_id",
  "year",
  "month",
  "latitude",
  "longitude",
  setdiff(selected_metrics$function_name, vars_to_remove)
)

# Clean data
clean_data <- final_data[, ..final_metric_vars]

summary(clean_data)


# VIF ANALYSIS -----------------------------------------------------------

metric_data <- clean_data %>%
  select(where(is.numeric)) %>%
  select(-c(plot_id, month, latitude, longitude)) %>%
  na.omit()

summary(metric_data)

vif_results <- vifstep(metric_data, th = 1.5, size = 10000, method = "spearman")
print(vif_results)

final_clean_data <- exclude(clean_data, vif_results)
summary(final_clean_data)


# Correlation plot
correlation_matrix <- cor(
  final_clean_data[, 6:12],
  use = "pairwise.complete.obs",
  method = "spearman"
)

# Create the correlation plot using ellipses
# A combined plot: ellipses in the upper triangle, numbers in the lower
pdf("Figs/corrplot_firemetrics.pdf", paper = "a4")
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

saveRDS(final_clean_data, "Data/final_clean_data.rds")


# Run the PCA
pca_fire_metrics <- prcomp(
  final_clean_data[, 6:12],
  center = TRUE,
  scale. = TRUE
)
summary(pca_fire_metrics)

# Convert year to numeric
final_clean_data[, year_num := as.integer(year)]

# Define the breaks and labels for the decades
breaks <- c(1984, 1994, 2004, 2014, 2024)
labels <- c("1985-1994", "1995-2004", "2005-2014", "2015-2024")

# Create the decade column
final_clean_data[, decade := cut(year_num, breaks = breaks, labels = labels)]
summary(final_clean_data)

# Create the biplot with ellipses instead of points
pdf("Figs/PCA_firemetrics.pdf", paper = "a4")
autoplot(
  pca_fire_metrics,
  data = final_clean_data[, 6:14], # Provide the data with the grouping variable
  colour = 'decade', # Color the ellipses by decade

  # These arguments create the ellipses
  frame = TRUE,
  frame.type = 'norm',

  # These arguments are for the variable arrows (loadings)
  loadings = TRUE,
  loadings.colour = 'black',
  loadings.label = TRUE,
  loadings.label.colour = 'black',
  loadings.label.size = 3.5
) +
  scale_color_viridis_d(name = "Decade") + # Use a nice color palette
  labs(
    title = "PCA of Fire Landscape Metrics by Decade",
    subtitle = "Ellipses show 95% confidence intervals for each period"
  ) +
  theme_minimal()
dev.off()


cat("Part 2 Complete: Metrics calculated and analyzed.\n")

# Create the fire_occurrence variable ------------------------------------

# It's 1 if total class area > 0, and 0 otherwise.
# We also treat NA as 0, assuming NA means no burned area was found.
# Create a new object with the key columns and the fire_occurrence variable
fire_occurrence_dt <- clean_data[, .(
  plot_id,
  year,
  month,
  fire_occurrence = ifelse(lsm_c_ca > 0 & !is.na(lsm_c_ca), 1, 0)
)]

# Check the new object
print(head(fire_occurrence_dt))
print(table(fire_occurrence_dt$fire_occurrence))

# Join the master grid with the fire occurrence data ---
# This is a left join, keeping all rows from the master_grid.
complete_fire_data <- fire_occurrence_dt[
  master_grid,
  on = c("plot_id", "year", "month")
]

# --- Step 2: Fill the missing values (NAs) with 0 ---
# Where 'fire_occurrence' is NA, it means no fire was recorded, so we set it to 0.
complete_fire_data[is.na(fire_occurrence), fire_occurrence := 0]

# --- Step 3: Verify the result ---
# The table should now show counts for both 0s (no fire) and 1s (fire).
print(table(complete_fire_data$fire_occurrence))

# Merge coordinates with the final metrics table
complete_fire_data <- plot_locations_dt[complete_fire_data, on = "plot_id"]

# Convert year to numeric
complete_fire_data[, year_num := as.integer(year)]

# Define the breaks and labels for the decades
breaks <- c(1984, 1994, 2004, 2014, 2024)
labels <- c("1985-1994", "1995-2004", "2005-2014", "2015-2024")

# Create the decade column
complete_fire_data[, decade := cut(year_num, breaks = breaks, labels = labels)]

summary(complete_fire_data)

saveRDS(complete_fire_data, "Data/final_fire_occ_data.rds")

# AGGREGATE DATA BY DECADE -----------------------------------------------

decadal_means_dt <- final_clean_data[,
  lapply(.SD, mean, na.rm = TRUE),
  by = .(plot_id, longitude, latitude, month, decade),
  .SDcols = names(final_clean_data[, 6:12])
]

summary(decadal_means_dt)

pdf("Figs/decadalmaps_cai_sd.pdf", paper = "a4r", width = 0, height = 0)
ggplot(
  decadal_means_dt,
  aes(x = longitude, y = latitude, fill = lsm_c_cai_sd)
) +
  geom_raster() +
  facet_grid(decade ~ month) +
  scale_fill_viridis_c(
    name = "SD CAI",
    na.value = "transparent",
    transform = "log1p",
    option = "inferno"
  ) +
  coord_equal() +
  labs(
    title = "SD CAI by decade and month",
    x = "Longitude (m)",
    y = "Latitude (m)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8)
  )
dev.off()

pdf("Figs/decadalmaps_circle_mn.pdf", paper = "a4r", width = 0, height = 0)
ggplot(
  decadal_means_dt,
  aes(x = longitude, y = latitude, fill = lsm_c_circle_mn)
) +
  geom_raster() +
  facet_grid(decade ~ month) +
  scale_fill_viridis_c(
    name = "Mean circle",
    na.value = "transparent",
    transform = "log1p",
    option = "inferno"
  ) +
  coord_equal() +
  labs(
    title = "Mean of related circumscribing circle",
    x = "Longitude (m)",
    y = "Latitude (m)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8)
  )
dev.off()

pdf("Figs/decadalmaps_core_mn.pdf", paper = "a4r", width = 0, height = 0)
ggplot(
  decadal_means_dt,
  aes(x = longitude, y = latitude, fill = lsm_c_core_mn)
) +
  geom_raster() +
  facet_grid(decade ~ month) +
  scale_fill_viridis_c(
    name = "Mean core area",
    na.value = "transparent",
    transform = "log1p",
    option = "inferno"
  ) +
  coord_equal() +
  labs(title = "Mean core area", x = "Longitude (m)", y = "Latitude (m)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8)
  )
dev.off()

pdf("Figs/decadalmaps_dcore_sd.pdf", paper = "a4r", width = 0, height = 0)
ggplot(
  decadal_means_dt,
  aes(x = longitude, y = latitude, fill = lsm_c_dcore_sd)
) +
  geom_raster() +
  facet_grid(decade ~ month) +
  scale_fill_viridis_c(
    name = "SD number of disjuct core areas",
    na.value = "transparent",
    transform = "log1p",
    option = "inferno"
  ) +
  coord_equal() +
  labs(
    title = "SD number of disjuct core areas",
    x = "Longitude (m)",
    y = "Latitude (m)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8)
  )
dev.off()

pdf("Figs/decadalmaps_enn_sd.pdf", paper = "a4r", width = 0, height = 0)
ggplot(
  decadal_means_dt,
  aes(x = longitude, y = latitude, fill = lsm_c_enn_sd)
) +
  geom_raster() +
  facet_grid(decade ~ month) +
  scale_fill_viridis_c(
    name = "SD ENN",
    na.value = "transparent",
    transform = "log1p",
    option = "inferno"
  ) +
  coord_equal() +
  labs(title = "SD ENN", x = "Longitude (m)", y = "Latitude (m)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8)
  )
dev.off()

pdf("Figs/decadalmaps_pd.pdf", paper = "a4r", width = 0, height = 0)
ggplot(decadal_means_dt, aes(x = longitude, y = latitude, fill = lsm_c_pd)) +
  geom_raster() +
  facet_grid(decade ~ month) +
  scale_fill_viridis_c(
    name = "Density of fire scars",
    na.value = "transparent",
    transform = "log1p",
    option = "inferno"
  ) +
  coord_equal() +
  labs(
    title = "Density of fire scars",
    x = "Longitude (m)",
    y = "Latitude (m)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8)
  )
dev.off()

pdf("Figs/decadalmaps_split.pdf", paper = "a4r", width = 0, height = 0)
ggplot(decadal_means_dt, aes(x = longitude, y = latitude, fill = lsm_c_split)) +
  geom_raster() +
  facet_grid(decade ~ month) +
  scale_fill_viridis_c(
    name = "Splitting index",
    na.value = "transparent",
    transform = "log1p",
    option = "inferno"
  ) +
  coord_equal() +
  labs(title = "Splitting index", x = "Longitude (m)", y = "Latitude (m)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8)
  )
dev.off()


# Aggregate fire occurrence data by decade
decadal_fire_occ_dt <- complete_fire_data[,
  lapply(.SD, mean, na.rm = TRUE),
  by = .(plot_id, longitude, latitude, month, decade),
  .SDcols = "fire_occurrence"
]

summary(decadal_fire_occ_dt)

pdf("Figs/decadalmaps_fire_occ.pdf", paper = "a4r", width = 0, height = 0)
ggplot(
  decadal_fire_occ_dt,
  aes(x = longitude, y = latitude, fill = fire_occurrence)
) +
  geom_raster() +
  facet_grid(decade ~ month) +
  scale_fill_viridis_c(
    name = "Fire occurrence",
    na.value = "transparent",
    option = "inferno"
  ) +
  coord_equal() +
  labs(
    title = "Fire occurrence by decade and month",
    x = "Longitude (m)",
    y = "Latitude (m)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 8)
  )
dev.off()
