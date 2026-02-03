# Spatio-Temporal Dynamics of Fire Regimes in the Cerrado

**Description:**
This repository contains the complete code and data processing pipeline used to analyze, map, and model fire regimes across the Cerrado biome (1985–2024). The project integrates 40 years of satellite fire data (MapBiomas) with climatic, topographic, and anthropogenic predictors to identify distinct fire regime typologies, assess their stability over time, and project spatial trends.

## Repository Structure

The project is organized into five main directories, reflecting the workflow from raw data processing to the final decision support system.

### 📂 `Scripts/`
Contains the three primary R scripts required to reproduce the analysis. The workflow is sequential:
1.  **`1_Script_MapbiomasFogo.R`**: The core processing engine. It ingests raw MapBiomas Fire rasters, generates binary monthly fire maps, and computes landscape metrics (e.g., density, cohesion, shape) on a 9km hexagonal grid. It also handles data cleaning and VIF analysis.
2.  **`2_Predictors.R`**: The environmental data aggregator. It processes and aligns multi-source predictors (ERA5-Land climate, MapBiomas LULC, NDVI, Population Density, Elevation) to the same spatial grid as the fire metrics.
3.  **`3_Cerrado_SpT_FireRegimes.R`**: The statistical modeling suite. It performs the Model-Based Clustering (to define regimes), runs Big Generalized Additive Models (BAMs) to identify drivers, calculates regime transition probabilities, and analyzes spatio-temporal trends using Random Forests and time-series decomposition.


#### `Scripts/1_Script_MapbiomasFogo.R`

**Description:**
This is the primary data processing pipeline. It transforms raw monthly fire rasters (MapBiomas Fire Collection) into the clean, spatially explicit datasets used for analysis and modeling.

**Workflow:**
1.  **Preprocessing:** Merges monthly GeoTIFFs into annual rasters and reprojects them to a metric CRS (SIRGAS 2000 / Brazil Polyconic, EPSG:5880).
2.  **Binarization:** Converts continuous fire data into binary (burned/unburned) monthly rasters.
3.  **Metric Calculation:** Computes landscape metrics (e.g., `lsm_c_pd`, `lsm_c_split`) on a 9km hexagonal grid using the `landscapemetrics` package. This step utilizes parallel processing and batch handling to manage memory.
4.  **Data Cleaning:** * Filters metrics based on missing data thresholds (>5% missingness removed).
    * Performs Variance Inflation Factor (VIF) analysis to remove collinear variables.
5.  **Exploratory Analysis:** Generates PCA biplots and correlation matrices to visualize relationships between metrics.
6.  **Fire Occurrence:** Generates a binary fire occurrence dataset aligned with the master grid.
7.  **Visualization:** Exports decadal map time-series (PDF format) for visual inspection of metric trends.

**Inputs Required:**
* Raw MapBiomas Fire GeoTIFFs (folder: `~/Documents/mapbiomas_fogo_4` or configured path).
* Cerrado geometry file: `Data/Cerrado.rds`.

**Outputs Generated:**
* **Final Data Objects:**
    * `Data/final_clean_data.rds`: The cleaned dataframe of landscape fire metrics.
    * `Data/final_fire_occ_data.rds`: The binary fire occurrence dataset.
* **Figures:**
    * `Figs/corrplot_firemetrics.pdf`: Correlation matrix of selected metrics.
    * `Figs/PCA_firemetrics.pdf`: Principal Component Analysis showing metric trends by decade.
    * `Figs/decadalmaps_*.pdf`: Spatiotemporal maps for individual metrics (e.g., `cai_sd`, `split`, `pd`).

#### `Scripts/2_Predictors.R`

**Description:**
This script aggregates and standardizes the environmental and anthropogenic predictors used to model fire regimes. It aligns disparate data sources (Land Use, Climate, Topography, Population, and Vegetation Indices) onto the common 9km hexagonal grid defined in the previous step.

**Workflow:**
1.  **Land Use & Land Cover (LULC):**
    * **Source:** MapBiomas Collection 10 (Annual Maps).
    * **Processing:** Calculates the fractional coverage of each LULC class within each 9km grid cell using parallel processing (`foreach` + `exactextractr`).
    * **Simplification:** Aggregates detailed classes into functional groups (e.g., "Forest", "Savanna", "Temporary Crops") and calculates total "Native" vs. "Anthropic" areas per cell/year.
2.  **NDVI (Vegetation Greenness):**
    * **Source:** Landsat Collection 2 Level 2 Annual Composites.
    * **Processing:** Resamples annual NDVI rasters to match the analysis grid resolution and CRS.
3.  **Population Density:**
    * **Source:** `popgrids` package (interpolated global gridded population).
    * **Processing:** Generates annual population density surfaces via interpolation between decadal census years (1985–2030) and extracts mean values per grid cell.
4.  **Climate (ERA5-Land):**
    * **Source:** Copernicus Climate Change Service (ERA5-Land monthly data).
    * **Processing:** Projects variables (Temperature, Precipitation, Evapotranspiration, Wind, Soil Water, etc.) to SIRGAS 2000.
    * **Derivation:** Calculates derived variables like Water Deficit (`total_evap - pot_evap`) and average Wind Speed. Computes monthly statistics (Mean, Min, Max, SD) for every grid cell.
5.  **Terrain:**
    * **Source:** Amazon Web Services Terrain Tiles (via `elevatr`).
    * **Processing:** Derives topographic metrics (Slope, Aspect, TPI, TRI) from the DEM and aligns them to the climate grid.

**Inputs Required:**
* **Raw Data Paths:** The script assumes specific local paths (e.g., `~/Documents/mapbiomas_lulc_10`, `~/Documents/ERA5`). **User Action:** You must adjust the `data_directory` variables at the top of each section to point to your local data storage.
* **Vector Data:** `Data/Cerrado.rds` (Biome boundary).

**Outputs Generated:**
* `Data/lulc_predictors.rds`: Time-series of fractional land cover per plot.
* `Data/ndvi_predictor.rds`: Time-series of mean annual NDVI per plot.
* `Data/pop_predictor.rds`: Time-series of population density per plot.
* `Data/climate_predictors.rds`: Time-series of monthly climatic variables per plot.
* `Data/terrain_predictors.rds`: Static topographic variables per plot.

#### `Scripts/3_Cerrado_SpT_FireRegimes.R`

**Description:**
This script performs the core statistical analysis and modeling of the study. It merges the processed metrics and predictors, runs Generalized Additive Models (GAMs) to identify drivers, classifies fire regimes using model-based clustering, and analyzes temporal trends via time-series decomposition.

**Workflow:**
1.  **Data Merging & Preprocessing:**
    * Integrates all intermediate datasets (`final_clean_data.rds`, `climate_predictors.rds`, `lulc_predictors.rds`, etc.) into a single master dataframe.
    * Standardizes predictor variables (scaling) to ensure comparability in modeling.
    * Handles missing data imputation where appropriate for specific model structures.

2.  **Fire Regime Classification & Validation:**
    * **Clustering:** Uses **Model-Based Clustering (`mclust`)** to identify distinct fire regimes based on the multivariate "fingerprint" of fire metrics (density, size, shape, cohesion, etc.).
    * **Predictability Test:** Trains a **Random Forest classifier** on the identified clusters to assess regime predictability (reporting error rates and Brier scores).
    * **Feature Selection:** Applies the **Boruta algorithm** to identify the most important metrics defining each regime.

3.  **Decadal Transitions & Stability:**
    * **Transition Analysis:** Calculates transition probability matrices to quantify how fire regimes shift between decades (e.g., 1985–1994 to 2015–2024).
    * **Stability Testing:** Performs permutation tests to assess if regime persistence and transitions depart from random expectations.
    * **Visualization:** Generates Sankey diagrams to map the flow of landscape transformation from native to anthropogenic fire regimes.

4.  **Modeling Environmental Drivers (BAMs):**
    * Fits **Big Generalized Additive Models (BAMs)** using the `mgcv` package to predict fire metrics and occurrence.
    * Models non-linear relationships with climatic, topographic, and anthropogenic predictors.
    * Includes spatial smoothers (Gaussian Markov Random Fields) to account for spatial autocorrelation.

5.  **Spatio-Temporal Trend Analysis:**
    * **Time-Series Decomposition:** Applies the **X-13ARIMA-SEATS** method to decompose metric time-series into Seasonal, Trend, and Irregular components.
    * **Spatial Trends:** Quantifies spatially varying trends (e.g., latitudinal gradients) to map where fire metrics are increasing or decreasing across the biome.
    * **Evolutive Seasonality:** Tests for shifts in the timing and intensity of the fire season over the 40-year period.

**Inputs Required:**
* All outputs from `Scripts/1_Script_MapbiomasFogo.R` and `Scripts/2_Predictors.R` (located in the `Data/` folder).

**Outputs Generated:**
* **Models:** Saved `.rds` objects for all fitted BAMs (e.g., `Models/fire_occurrence_bam.rds`) and the Random Forest classifier.
* **Figures:**
    * `Figs/heatmap_fire_mclust_parameter.pdf`: Characterization of the identified fire regimes.
    * `Figs/map_fire_regimes.pdf`: Spatial distribution of regimes across the Cerrado.
    * `Figs/sankey_transition.pdf`: Visualization of regime transitions over time.
    * `Figs/spatial_trends_*.pdf`: Maps showing the strength and direction of trends for each metric.
    * `Figs/partial_effect_*.pdf`: Visual response curves for environmental predictors.
* **Tables:**
    * `Output/results/gam_performance.csv`: Summary of model fit statistics.
    * `Output/results/transition_matrices.csv`: Probability tables of regime shifts.


### 📂 `Data/`
Stores the processed input data and intermediate objects required for modeling.
* **Processed Predictors:** `.rds` files for Climate, LULC, NDVI, Population, and Terrain variables.
* **Final Datasets:** `final_clean_data.rds` (the master dataframe for modeling) and `final_fire_occ_data.rds` (binary occurrence data).
* **Spatial Files:** `Cerrado.rds` (biome boundary) and grid geometries.
* *Note: Raw heavy raster files (e.g., original MapBiomas TIFFs, ERA5 GRIBs) are not included due to size constraints but can be downloaded from their respective sources listed in the scripts.*

### 📂 `Output/`
Contains the statistical and analytical results generated by the scripts.
* **`results/`**: Model performance tables (GAM R²/Deviance), Random Forest error rates, and variable importance scores.
* **`results_model_evaluation/`**: Cross-validation metrics and sensitivity analyses.
* **`Decadal_Fire_Metrics/`**: Aggregated metric summaries used for visualizing long-term trends.
* **`Models/`**: Saved R objects for the fitted BAMs and Random Forest classifiers (allowing users to predict without re-training).

### 📂 `Figs/`
The repository for all static visualizations produced by the analysis.
* **`map_fire_regimes.pdf`**: The final classification map of Cerrado fire regimes.
* **`sankey_transition.pdf`**: Visualizing the flow of regime shifts over four decades.
* **`PCA_firemetrics.pdf`**: Biplots showing the trajectory of fire characteristics through time.
* **`spatial_trends_*.pdf`**: Maps illustrating where specific fire metrics (e.g., frequency, size) are statistically increasing or decreasing.

### 📂 `ShinyApp/`

**Description:**
This directory contains the source code for the **Cerrado Fire Regimes Decision Support System (DSS)**. It includes the data optimization pipeline and the interactive web application used to visualize the study's results.

**Contents:**

* **`0_Prepare_Data_ShinyApp.R`**:
    * **Purpose:** The optimization engine. It takes the heavy, high-resolution outputs from the main analysis (`Scripts/3_Cerrado_SpT_FireRegimes.R`) and transforms them into lightweight, spatially simplified objects suitable for web deployment.
    * **Key Operations:** Aggregates metrics by municipality, simplifies vector geometries (using `rmapshaper` to reduce load times), and pre-calculates summary statistics for the dashboard value boxes.
    * **Output:** Generates the optimized `.rds` files stored in the `CerradoFireRegimes_DSS` folder.

* **`1_ShinyApp_Cerrado_FireRegimes.R`**:
    * **Purpose:** The local development version of the DSS.
    * **Functionality:** Contains the combined `ui` (User Interface) and `server` logic in a single script. It is designed for rapid testing, debugging, and feature iteration on a local machine before deployment.

* **`CerradoFireRegimes_DSS/`**:
    * **Purpose:** The production-ready folder configured for deployment to **shinyapps.io**.
    * **Structure:**
        * `ui.R` & `server.R`: The split application logic required for stable server deployment.
        * `rsconnect/`: Stores web assets (images, CSS styles) for the interface.
    * **Usage:** This is the specific folder targeted when publishing the application to the web.

---

## Reproducibility
To reproduce the analysis:
1.  Clone this repository.
2.  Ensure you have the raw data sources listed in `Scripts/1_Script_MapbiomasFogo.R` and `Scripts/2_Predictors.R`.
3.  Run the scripts in numerical order (1 → 2 → 3).
4.  Launch the Shiny App via `ShinyApp/CerradoFireRegimes_DSS/1_ShinyApp_Cerrado_FireRegimes.R` to visualize the results.

## Citation
If you use this code or data, please cite the associated manuscript and the Zenodo DOI:
* [Insert Manuscript Citation]
* [Insert Zenodo DOI]



