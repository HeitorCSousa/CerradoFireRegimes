rm(list = ls())

# --- 0. LOAD PACKAGES ---------------------------------------------------------

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
vectors <- readRDS("ShinyApp/CerradoFireRegimes_DSS/vectors_cerrado.rds")

# File Paths
global_file <- "ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_Global.tif"

decade_files <- c(
  "1985-1994" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_1985_1994.tif",
  "1995-2004" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_1995_2004.tif",
  "2005-2014" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_2005_2014.tif",
  "2015-2024" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Regimes_2015_2024.tif"
)

metric_files <- list(
  "Fire frequency in Jan" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Jan.tif",
  "Fire frequency in Feb" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Feb.tif",
  "Fire frequency in Mar" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Mar.tif",
  "Fire frequency in Apr" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Apr.tif",
  "Fire frequency in May" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_May.tif",
  "Fire frequency in Jun" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Jun.tif",
  "Fire frequency in Jul" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Jul.tif",
  "Fire frequency in Aug" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Aug.tif",
  "Fire frequency in Sep" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Sep.tif",
  "Fire frequency in Oct" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Oct.tif",
  "Fire frequency in Nov" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Nov.tif",
  "Fire frequency in Dec" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Freq_Dec.tif",
  "Fire Integrity Var" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Integrity_Var.tif",
  "Fire Core Size" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Core_Size.tif",
  "Fire Shape Regularity" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Shape_Reg.tif",
  "Fire Internal Heterog" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Internal_Heterog.tif",
  "Fire Isolation Var" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Isolation_Var.tif",
  "Fire Density" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Density.tif",
  "Fire Fragmentation" = "ShinyApp/CerradoFireRegimes_DSS/Fire_Fragmentation.tif"
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

# --- 3. SERVER ----------------------------------------------------------------
server <- function(input, output, session) {
  # --- Cache Metric Domains ---
  metric_domains <- reactive({
    doms <- list()
    for (name in names(metric_files)) {
      if (file.exists(metric_files[[name]])) {
        r <- rast(metric_files[[name]])
        rng <- minmax(r)
        doms[[name]] <- c(log1p(rng[1]), log1p(rng[2]))
      }
    }
    return(doms)
  })

  # --- A. Update Inputs ---

  # Update Municipality Dropdown
  observeEvent(input$sel_state_filter, {
    req(vectors$munis)
    filtered <- vectors$munis %>% filter(name_state == input$sel_state_filter)
    updateSelectInput(
      session,
      "sel_muni",
      choices = sort(unique(filtered$name_muni))
    )
  })

  # NEW: Update Download Layer Choices
  # This updates the dropdown to include the current Regime Map + All Metrics
  observe({
    # Determine the current regime name
    regime_name <- if (input$view_mode == "Global Model") {
      "Fire Regimes (Global)"
    } else {
      paste0("Fire Regimes (", input$time_slider, ")")
    }

    # Combine with metrics
    choices <- c(regime_name, names(metric_files))

    updateSelectInput(
      session,
      "dl_layer",
      choices = choices,
      selected = regime_name
    )
  })

  # --- B. Raster Logic ---
  r_display <- reactive({
    file <- if (input$view_mode == "Global Model") {
      global_file
    } else {
      decade_files[[input$time_slider]]
    }
    req(file.exists(file))
    r <- rast(file)
    project(r, "EPSG:3857", method = "near")
  })

  # --- C. Map Initialization ---
  output$map <- renderLeaflet({
    metric_names <- names(metric_files)

    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron, group = "Light") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
      setView(lng = -50, lat = -15, zoom = 5) %>%
      addLayersControl(
        baseGroups = c("Light", "Satellite"),
        overlayGroups = c(
          "Fire Regimes",
          metric_names,
          "States",
          "Municipalities",
          "UCs",
          "TIs"
        ),
        options = layersControlOptions(collapsed = FALSE)
      ) %>%
      hideGroup(c(metric_names, "Municipalities", "UCs", "TIs")) %>%

      addPolygons(
        data = vectors$states,
        fill = FALSE,
        color = "#444",
        weight = 2,
        group = "States"
      ) %>%
      addPolygons(
        data = vectors$munis,
        fill = FALSE,
        color = "#999",
        weight = 1,
        group = "Municipalities"
      ) %>%
      addPolygons(
        data = vectors$ucs,
        fill = FALSE,
        color = "green",
        weight = 2,
        group = "UCs"
      ) %>%
      addPolygons(
        data = vectors$tis,
        fill = FALSE,
        color = "orange",
        weight = 2,
        group = "TIs"
      )
  })

  # --- D. Dynamic UI (Sliders) ---
  output$dynamic_layer_controls <- renderUI({
    req(input$map_groups)
    controls <- list()

    if ("Fire Regimes" %in% input$map_groups) {
      val <- if (!is.null(input$alpha_regime)) input$alpha_regime else 0.8
      controls[[1]] <- div(
        class = "panel panel-default",
        style = "padding:10px;",
        h5("Fire Regimes"),
        sliderInput("alpha_regime", NULL, 0, 1, val, 0.1)
      )
    }

    for (name in names(metric_files)) {
      if (name %in% input$map_groups) {
        id <- paste0("alpha_", clean_id(name))
        val <- if (!is.null(input[[id]])) input[[id]] else 0.7
        controls[[length(controls) + 1]] <- div(
          class = "panel panel-default",
          style = "padding:10px;",
          h5(name),
          sliderInput(id, NULL, 0, 1, val, 0.1)
        )
      }
    }
    do.call(tagList, controls)
  })

  # --- E. Update Layers ---

  # 1. Fire Regimes
  observe({
    req(r_display())
    op <- if (!is.null(input$alpha_regime)) input$alpha_regime else 0.8

    if (is.null(input$map_groups) || "Fire Regimes" %in% input$map_groups) {
      leafletProxy("map") %>%
        clearGroup("Fire Regimes") %>%
        addRasterImage(
          r_display(),
          colors = regime_pal,
          opacity = op,
          group = "Fire Regimes",
          project = FALSE,
          method = "ngb",
          maxBytes = 15 * 1024 * 1024
        )
    } else {
      leafletProxy("map") %>% clearGroup("Fire Regimes")
    }
  })

  # 2. Metrics
  observe({
    map <- leafletProxy("map")

    for (name in names(metric_files)) {
      if (!is.null(input$map_groups) && name %in% input$map_groups) {
        id <- paste0("alpha_", clean_id(name))
        op <- if (!is.null(input[[id]])) input[[id]] else 0.7

        if (file.exists(metric_files[[name]])) {
          r <- project(
            rast(metric_files[[name]]),
            "EPSG:3857",
            method = "bilinear"
          )
          r_log <- log1p(r)
          vals <- values(r_log, mat = FALSE)
          pal <- colorNumeric(
            "viridis",
            range(vals, na.rm = TRUE),
            na.color = "transparent"
          )

          map %>%
            clearGroup(name) %>%
            addRasterImage(
              r_log,
              colors = pal,
              opacity = op,
              group = name,
              project = FALSE,
              maxBytes = 10 * 1024 * 1024
            )
        }
      } else {
        map %>% clearGroup(name)
      }
    }
  })

  # --- F. Legends ---
  observeEvent(input$map_groups, {
    map <- leafletProxy("map") %>% clearControls()

    if ("Fire Regimes" %in% input$map_groups) {
      map %>%
        addLegend(
          pal = regime_pal,
          values = 1:5,
          title = "Regime",
          position = "bottomleft",
          opacity = 1
        )
    }

    domains <- metric_domains()
    for (name in names(metric_files)) {
      if (name %in% input$map_groups && !is.null(domains[[name]])) {
        rng <- domains[[name]]
        pal <- colorNumeric("viridis", domain = rng, na.color = "transparent")
        map %>%
          addLegend(
            pal = pal,
            values = rng,
            title = name,
            position = "bottomleft",
            labFormat = labelFormat(transform = function(x) round(expm1(x), 1))
          )
      }
    }
  })

  # --- G. DOWNLOAD HANDLER (Updated) ---
  output$download_data <- downloadHandler(
    filename = function() {
      # Create clean filename based on selection
      layer_name <- gsub("[^[:alnum:]]", "", input$dl_layer)
      scope_name <- gsub(" ", "", input$dl_scope)
      paste0(layer_name, "_", scope_name, ".tif")
    },

    content = function(file) {
      showNotification("Processing download...", type = "message")

      # 1. IDENTIFY WHICH FILE TO LOAD
      # Check if the user selected a "Regime" map or a "Metric" map
      selected_layer <- input$dl_layer

      file_to_load <- NULL
      data_type <- "FLT4S" # Default for metrics

      # Logic to find the file
      if (grepl("Fire Regimes", selected_layer)) {
        # It's a regime map. Check if global or decade.
        if (grepl("Global", selected_layer)) {
          file_to_load <- global_file
        } else {
          # Extract decade? Or just rely on time_slider if it matches
          # Easier: The user selected "Fire Regimes (1985-1994)".
          # We can pull from decade_files based on input$time_slider
          # BUT: What if they changed the slider after selecting?
          # Safer to match the string or just use current slider if it matches.
          # Let's trust input$time_slider if the string contains it, else default to slider
          file_to_load <- decade_files[[input$time_slider]]
        }
        data_type <- "INT1U" # Optimization for regimes
      } else {
        # It's a metric. Look up in list.
        file_to_load <- metric_files[[selected_layer]]
      }

      # 2. LOAD ORIGINAL RASTER (EPSG:5880)
      req(file_to_load)
      r_orig <- rast(file_to_load)

      # 3. SELECT CUTLINE
      cut_poly <- NULL
      if (input$dl_scope == "Crop by State") {
        cut_poly <- vectors$states %>% filter(name_state == input$sel_state)
      } else if (input$dl_scope == "Crop by Municipality") {
        cut_poly <- vectors$munis %>% filter(name_muni == input$sel_muni)
      } else if (input$dl_scope == "Crop by UC") {
        cut_poly <- vectors$ucs %>%
          filter(name_conservation_unit == input$sel_uc)
      } else if (input$dl_scope == "Crop by TI") {
        cut_poly <- vectors$tis %>% filter(name_indigenous_land == input$sel_ti)
      }

      # 4. CROP & SAVE
      if (is.null(cut_poly)) {
        writeRaster(r_orig, file, datatype = data_type)
      } else {
        cut_poly_proj <- st_transform(cut_poly, crs = crs(r_orig))
        r_crop <- crop(r_orig, cut_poly_proj)
        r_mask <- mask(r_crop, vect(cut_poly_proj))
        writeRaster(r_mask, file, datatype = data_type)
      }
    }
  )

  # Legend Plot
  output$legend_plot <- renderPlot({
    df <- data.frame(Regime = factor(1:5))
    ggplot(df, aes(x = 1, y = Regime, fill = Regime)) +
      geom_tile() +
      scale_fill_viridis_d(option = "turbo", name = "Regime") +
      theme_void() +
      theme(legend.position = "left", legend.text = element_text(size = 12))
  })
}

shinyApp(ui, server)

rsconnect::deployApp(
  '/Users/heito/Documents/GitHub/CerradoFireRegimes/ShinyApp'
)
