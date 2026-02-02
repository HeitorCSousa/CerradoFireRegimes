library(shiny)
library(leaflet)
library(terra)
library(viridis)
library(sf)
library(geobr)
library(shinyWidgets)
library(dplyr)
library(rsconnect)
library(ggplot2)

# --- 1. CONFIGURATION ---------------------------------------------------------

# Load Vectors (Must be created via prepare_vectors.R first)
vectors <- readRDS("vectors_cerrado.rds")

# File Paths
global_file <- "Fire_Regimes_Global.tif"

decade_files <- c(
  "1985-1994" = "Fire_Regimes_1985_1994.tif",
  "1995-2004" = "Fire_Regimes_1995_2004.tif",
  "2005-2014" = "Fire_Regimes_2005_2014.tif",
  "2015-2024" = "Fire_Regimes_2015_2024.tif"
)

metric_files <- list(
  "Fire frequency in Jan" = "Fire_Freq_Jan.tif",
  "Fire frequency in Feb" = "Fire_Freq_Feb.tif",
  "Fire frequency in Mar" = "Fire_Freq_Mar.tif",
  "Fire frequency in Apr" = "Fire_Freq_Apr.tif",
  "Fire frequency in May" = "Fire_Freq_May.tif",
  "Fire frequency in Jun" = "Fire_Freq_Jun.tif",
  "Fire frequency in Jul" = "Fire_Freq_Jul.tif",
  "Fire frequency in Aug" = "Fire_Freq_Aug.tif",
  "Fire frequency in Sep" = "Fire_Freq_Sep.tif",
  "Fire frequency in Oct" = "Fire_Freq_Oct.tif",
  "Fire frequency in Nov" = "Fire_Freq_Nov.tif",
  "Fire frequency in Dec" = "Fire_Freq_Dec.tif",
  "Fire Integrity Var" = "Fire_Integrity_Var.tif",
  "Fire Core Size" = "Fire_Core_Size.tif",
  "Fire Shape Regularity" = "Fire_Shape_Reg.tif",
  "Fire Internal Heterog" = "Fire_Internal_Heterog.tif",
  "Fire Isolation Var" = "Fire_Isolation_Var.tif",
  "Fire Density" = "Fire_Density.tif",
  "Fire Fragmentation" = "Fire_Fragmentation.tif"
)

# Palettes
regime_pal <- colorFactor(
  palette = viridis(5, option = "turbo"),
  domain = 1:5,
  na.color = "transparent"
)

# Helper to clean names for ID generation
clean_id <- function(x) gsub("[^[:alnum:]]", "_", x)

# --- 2. UI --------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML(
    "
    .leaflet-container { background: #f2f2f2; }
    .leaflet-image-layer { image-rendering: pixelated; }
    .panel-default { margin-bottom: 5px; }
  "
  ))),

  titlePanel("Cerrado Fire Regimes: Decision Support System"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # --- Visualization Controls ---
      h4("1. Visualization"),
      radioGroupButtons(
        "view_mode",
        label = NULL,
        choices = c("Global Model", "Time Series"),
        status = "primary",
        justified = TRUE
      ),

      conditionalPanel(
        condition = "input.view_mode == 'Time Series'",
        sliderTextInput(
          "time_slider",
          "Timeline:",
          choices = names(decade_files),
          selected = names(decade_files)[1],
          animate = animationOptions(interval = 1500),
          grid = TRUE
        )
      ),

      hr(),
      h4("Layer Settings"),
      helpText("Turn on layers in the map to see their opacity sliders here."),
      uiOutput("dynamic_layer_controls"),

      hr(),

      # --- Download Controls ---
      h4("2. Download Data"),

      # NEW: Select which layer to download
      selectInput("dl_layer", "Select Layer to Download:", choices = NULL),

      selectInput(
        "dl_scope",
        "Scope:",
        c(
          "Whole Map",
          "Crop by State",
          "Crop by Municipality",
          "Crop by UC",
          "Crop by TI"
        )
      ),

      conditionalPanel(
        "input.dl_scope == 'Crop by State'",
        selectInput(
          "sel_state",
          "State:",
          choices = sort(unique(vectors$states$name_state))
        )
      ),
      conditionalPanel(
        "input.dl_scope == 'Crop by Municipality'",
        selectInput(
          "sel_state_filter",
          "Filter State:",
          choices = sort(unique(vectors$states$name_state))
        ),
        selectInput("sel_muni", "Municipality:", choices = NULL)
      ),
      conditionalPanel(
        "input.dl_scope == 'Crop by UC'",
        selectizeInput(
          "sel_uc",
          "Unit:",
          choices = sort(unique(vectors$ucs$name_conservation_unit))
        )
      ),
      conditionalPanel(
        "input.dl_scope == 'Crop by TI'",
        selectizeInput(
          "sel_ti",
          "Indigenous Land:",
          choices = sort(unique(vectors$tis$name_indigenous_land))
        )
      ),

      br(),
      downloadButton(
        "download_data",
        "Download GeoTIFF",
        class = "btn-success"
      ),
      hr(),
      plotOutput("legend_plot", height = "150px")
    ),

    mainPanel(
      width = 9,
      leafletOutput("map", height = "85vh")
    )
  )
)
