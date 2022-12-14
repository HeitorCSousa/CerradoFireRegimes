######################################################
#Script for downloading ERA5 Climate Data for Cerrado#
######################################################

Sys.setenv(TZ = "America/Sao_Paulo")#set time zone

#Load package KrigR 
#devtools::install_github("https://github.com/ErikKusch/KrigR")
library(KrigR)

#Check if other packages are instalee
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

package_vec <- c(
  "tidyr", # for turning rasters into ggplot-dataframes
  "ggplot2", # for plotting
  "viridis", # colour palettes
  "cowplot", # gridding multiple plots
  "rosm", # obtaining satellite maps
  "gimms" # to get some pre-existing data to match in our downscaling
)
sapply(package_vec, install.load.package)

#Before we can proceed, you need to open up an account at the CDS and create 
#an API key by following this link: https://cds.climate.copernicus.eu/api-how-to. 
#This is required so that you may issue download requests through the KrigR package. 

#Register the user ID and API Key as characters 

API_User <- "86619"
API_Key <- "151ac044-f8cf-4ee0-aa00-64da02cf0295"

#A data directory for all of our individual Kriging processes
#A shapefile directory (located within our data directory) for 
#all of the shapefiles that we will use

Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "Data") # folder path for data
Dir.Shapes <- file.path(Dir.Data, "Shapes") # folder path for shapefiles
Dirs <- sapply(c(Dir.Data, Dir.Shapes), function(x) if (!dir.exists(x)) dir.create(x))

#Read and plot Cerrado area
CerradoMask <- readOGR(Dir.Shapes,"Cerrado", verbose = FALSE) # read

windows(10,10)
ggplot() +
  geom_polygon(data = CerradoMask, aes(x = long, y = lat, group = group), colour = "darkred", fill = "black") +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude")

#### STATE MASK (for within-nation borders)
if (!file.exists(file.path(Dir.Shapes, "StateMask.zip"))) { # if not downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces.zip",
                destfile = file.path(Dir.Shapes, "StateMask.zip")
  ) # download cultural vector
  unzip(file.path(Dir.Shapes, "StateMask.zip"), exdir = Dir.Shapes) # unzip data
}

StateMask <- readOGR(Dir.Shapes, "ne_10m_admin_1_states_provinces", verbose = FALSE) # read

windows(20,10)
ggplot() +
  geom_polygon(data = StateMask, aes(x = long, y = lat, group = group), colour = "darkred", fill = "black") +
  coord_sf(xlim = c(-85, -25), ylim = c(-60, 20), expand = FALSE) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude")

#### ECOREGIONS (for ecoregional borders)
if (!file.exists(file.path(Dir.Shapes, "WWF_ecoregions"))) { # if not downloaded yet
  download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
                destfile = file.path(Dir.Shapes, "wwf_ecoregions.zip")
  ) # download regions
  unzip(file.path(Dir.Shapes, "wwf_ecoregions.zip"), exdir = file.path(Dir.Shapes, "WWF_ecoregions")) # unzip data
}

EcoregionMask <- readOGR(file.path(Dir.Shapes, "WWF_ecoregions", "official", "wwf_terr_ecos.shp"), verbose = FALSE) # read

windows(20,10)
ggplot() +
  geom_polygon(data = EcoregionMask, aes(x = long, y = lat, group = group), colour = "darkred", fill = "black") +
  coord_sf(xlim = c(-85, -25), ylim = c(-60, 20), expand = FALSE) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude")


rasterOptions(tmpdir=file.path("D:/Documentos/Raster")) 

####################
#Plotting functions#
####################

#Raw Data#
##########
#This function plots the raw data that we will krig and exports a single plot 
#of all layers of the input raster

Plot_Raw <- function(Raw, Shp = NULL, Dates, Legend = "Air Temperature [K]") {
  Raw_df <- as.data.frame(Raw, xy = TRUE) # turn raster into dataframe
  colnames(Raw_df)[c(-1, -2)] <- Dates # set colnames
  Raw_df <- gather(data = Raw_df, key = Values, value = "value", colnames(Raw_df)[c(-1, -2)]) #  make ggplot-ready
  Raw_plot <- ggplot() + # create a plot
    geom_raster(data = Raw_df, aes(x = x, y = y, fill = value)) + # plot the raw data
    facet_wrap(~Values) + # split raster layers up
    theme_bw() +
    labs(x = "Longitude", y = "Latitude") + # make plot more readable
    scale_fill_gradientn(name = Legend, colours = inferno(100)) # add colour and legend
  if (!is.null(Shp)) { # if a shape has been designated
    Raw_plot <- Raw_plot + geom_polygon(data = Shp, aes(x = long, y = lat, group = group), colour = "black", fill = "NA") # add shape
  }
  return(Raw_plot)
} # export the plot

#Covariates#
############
#This function plots the covariate data we will be using by showing each 
#variable on both resolution levels side-by-side

Plot_Covs <- function(Covs, Shp = NULL) {
  Plots_ls <- as.list(rep(NA, nlayers(Covs[[1]]) * 2)) # create as many plots as there are covariates variables * 2
  Counter <- 1 # counter for list position
  for (Variable in 1:nlayers(Covs[[1]])) { # loop over all covariate variables
    Covs_Iter <- list(Covs[[1]][[Variable]], Covs[[2]][[Variable]]) # extract the data for this variable
    for (Plot in 1:2) { # loop over both resolutions for the current variable
      Cov_df <- as.data.frame(Covs_Iter[[Plot]], xy = TRUE) # turn raster into data frame
      colnames(Cov_df)[3] <- "Values" # assign a column name to the data column
      Plots_ls[[Counter]] <- ggplot() + # create plot
        geom_raster(data = Cov_df, aes(x = x, y = y, fill = Values)) + # plot the covariate data
        theme_bw() +
        labs(x = "Longitude", y = "Latitude") + # make plot more readable
        scale_fill_gradientn(name = names(Covs_Iter[[Plot]]), colours = cividis(100)) + # add colour and legend
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) # reduce margins (for fusing of plots)
      if (!is.null(Shp)) { # if a shape has been designated
        Plots_ls[[Counter]] <- Plots_ls[[Counter]] + geom_polygon(data = Shp, aes(x = long, y = lat, group = group), colour = "black", fill = "NA") # add shape
      }
      Counter <- Counter + 1 # raise list counter
    } # end of resolution loop
  } # end of variable loop
  ggPlot <- plot_grid(plotlist = Plots_ls, ncol = 2, labels = "AUTO") # fuse the plots into one big plot
  return(ggPlot)
} # export the plot


#Kriged products#
#################
#This function plots the Kriging outputs by splitting them up according to 
#whether they are Kriging predictions or the uncertainties associated with them

Plot_Krigs <- function(Krigs, Shp = NULL, Dates, Legend = "Air Temperature [K]") {
  Type_vec <- c("Prediction", "Standard Error") # these are the output types of krigR
  Colours_ls <- list(inferno(100), rev(viridis(100))) # we want separate colours for the types
  Plots_ls <- as.list(NA, NA) # this list will be filled with the output plots
  for (Plot in 1:2) { # loop over both output types
    Krig_df <- as.data.frame(Krigs[[Plot]], xy = TRUE) # turn raster into dataframe
    colnames(Krig_df)[c(-1, -2)] <- paste(Type_vec[Plot], Dates) # set colnames
    Krig_df <- gather(data = Krig_df, key = Values, value = "value", colnames(Krig_df)[c(-1, -2)]) # make ggplot-ready
    Plots_ls[[Plot]] <- ggplot() + # create plot
      geom_raster(data = Krig_df, aes(x = x, y = y, fill = value)) + # plot the kriged data
      facet_wrap(~Values) + # split raster layers up
      theme_bw() +
      labs(x = "Longitude", y = "Latitude") + # make plot more readable
      scale_fill_gradientn(name = Legend, colours = Colours_ls[[Plot]]) + # add colour and legend
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) # reduce margins (for fusing of plots)
    if (!is.null(Shp)) { # if a shape has been designated
      Plots_ls[[Plot]] <- Plots_ls[[Plot]] + geom_polygon(data = Shp, aes(x = long, y = lat, group = group), colour = "black", fill = "NA") # add shape
    }
  } # end of type-loop
  ggPlot <- plot_grid(plotlist = Plots_ls, ncol = 1, labels = "AUTO") # fuse the plots into one big plot
  return(ggPlot)
} # export the plot

##############
#Climate Data#
##############
Extent <- extent(CerradoMask)
windows(10,10)
bmaps.plot(bbox = CerradoMask, type = "AerialWithLabels", quiet = TRUE, progress = "none")

Dir.CerradoExt <- file.path(Dir.Data, "Cerrado_Extent")
#dir.create(Dir.CerradoExt)

Variable_List("era5-land")
#For more info, check:
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means?tab=overview
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview

Cerrado_Raw_temp <- download_ERA(
  Variable = "2m_temperature",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_temp

Cerrado_Raw_skin_temp <- download_ERA(
  Variable = "skin_temperature",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_skin_temp

Cerrado_Raw_temp_lv1 <- download_ERA(
  Variable = "soil_temperature_level_1",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_temp_lv1

Cerrado_Raw_temp_lv2 <- download_ERA(
  Variable = "soil_temperature_level_1",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_temp_lv2

Cerrado_Raw_sol <- download_ERA(
  Variable = "surface_solar_radiation_downwards",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_sol

Cerrado_Raw_water_lv1 <- download_ERA(
  Variable = "volumetric_soil_water_layer_1",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_water_lv1

Cerrado_Raw_water_lv2 <- download_ERA(
  Variable = "volumetric_soil_water_layer_1",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_water_lv2

Cerrado_Raw_water_veg <- download_ERA(
  Variable = "skin_reservoir_content",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_water_veg

Cerrado_Raw_precip <- download_ERA(
  Variable = "total_precipitation",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
  )

Cerrado_Raw_precip

Cerrado_Raw_leaf_high <- download_ERA(
  Variable = "leaf_area_index_high_vegetation",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_leaf_high

Cerrado_Raw_leaf_low <- download_ERA(
  Variable = "leaf_area_index_low_vegetation",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_leaf_low

Cerrado_Raw_total_evap <- download_ERA(
  Variable = "total_evaporation",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_total_evap

Cerrado_Raw_pot_evap <- download_ERA(
  Variable = "potential_evaporation",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_pot_evap

Cerrado_Raw_dewpoint_temp <- download_ERA(
  Variable = "2m_dewpoint_temperature",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_dewpoint_temp

Cerrado_Raw_surf_press <- download_ERA(
  Variable = "surface_pressure",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_surf_press


Cerrado_Raw_wind_u <- download_ERA(
  Variable = "10m_u_component_of_wind",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_wind_u

Cerrado_Raw_wind_v <- download_ERA(
  Variable = "10m_v_component_of_wind",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  Extent = CerradoMask,
  Dir = Dir.CerradoExt,
  API_User = API_User,
  API_Key = API_Key,
  SingularDL = TRUE
)

Cerrado_Raw_wind_v

#Derived variables
#Average wind speed
Cerrado_Raw_wind<-(abs(Cerrado_Raw_wind_v)+abs(Cerrado_Raw_wind_u))/2

#Water dephicit (total_evap - pot_evap)
Cerrado_Raw_total_evap<-stack("Cerrado_total_evap.grd")
Cerrado_Raw_pot_evap<-stack("Cerrado_pot_evap.grd")
Cerrado_Raw_water_deficit<-(Cerrado_Raw_total_evap-Cerrado_Raw_pot_evap)
Cerrado_Raw_water_deficit

#Relative humidity
#Following formulae in https://confluence.ecmwf.int/pages/viewpage.action?pageId=171411214

Rdry <- 287.0597 ; Rvap <- 461.5250
a1 <- 611.21 ; a3 <- 17.502 ; a4 <- 32.19 ; T0 <- 273.16

Esat <- a1*exp(a3*(Cerrado_Raw_dewpoint_temp-T0)/(Cerrado_Raw_dewpoint_temp-a4))

Cerrado_Raw_humidity <- ((Rdry/Rvap)*Esat)/(Cerrado_Raw_surf_press-((1-Rdry/Rvap)*Esat))
Cerrado_Raw_humidity

Cerrado_Raw_RH <- (0.263 * Cerrado_Raw_surf_press *Cerrado_Raw_humidity)/
  exp((17.67*(Cerrado_Raw_temp-T0))/(Cerrado_Raw_temp-29.65))
Cerrado_Raw_RH

writeRaster(Cerrado_Raw_temp,"Cerrado_temp_2m.grd")
writeRaster(Cerrado_Raw_skin_temp,"Cerrado_skin_temp.grd")
writeRaster(Cerrado_Raw_temp_lv1,"Cerrado_temp_lv1.grd")
writeRaster(Cerrado_Raw_temp_lv2,"Cerrado_temp_lv2.grd")
writeRaster(Cerrado_Raw_sol,"Cerrado_sol.grd")
writeRaster(Cerrado_Raw_total_evap,"Cerrado_total_evap.grd")
writeRaster(Cerrado_Raw_pot_evap,"Cerrado_pot_evap.grd")
writeRaster(Cerrado_Raw_water_deficit,"Cerrado_water_deficit.grd")
writeRaster(Cerrado_Raw_leaf_high,"Cerrado_leaf_high.grd")
writeRaster(Cerrado_Raw_leaf_low,"Cerrado_leaf_low.grd")
writeRaster(Cerrado_Raw_water_lv1,"Cerrado_water_lv1.grd")
writeRaster(Cerrado_Raw_water_lv2,"Cerrado_water_lv2.grd")
writeRaster(Cerrado_Raw_precip,"Cerrado_precip.grd")
writeRaster(Cerrado_Raw_wind,"Cerrado_wind.grd")
writeRaster(Cerrado_Raw_humidity,"Cerrado_humidity.grd")
writeRaster(Cerrado_Raw_RH,"Cerrado_RH.grd")
writeRaster(Cerrado_Raw_dewpoint_temp,"Cerrado_dewpoint.grd")
writeRaster(Cerrado_Raw_surf_press,"Cerrado_surf_press.grd")


Raw_df_test <- as.data.frame(Cerrado_Raw_wind, xy = TRUE) # turn raster into dataframe
(colnames(Raw_df_test)[c(-1, -2)]) # set colnames

windows(10,10)
plot(Cerrado_Raw_temp[[1:12]])

#############################
#Covariates for krigin - DEM#
#############################

Cerrado_Raw_temp<-stack("Cerrado_temp_2m.grd")
z1.agreg<-stack("E:/Heitor/BD_SIG/MODIS_Fire/MCD64A1/Bbin.all.Cerrado.MCD64A1.2.5km.grd")

BA.all.Cerrado.LTDR<-stack("E:/Heitor/BD_SIG/AVHRR_LTDR_Fire/Pixel/BA.all.LTDR.Cerrado.tif")
BA.all.Cerrado.LTDR

Covs_ls <- download_DEM(
  Train_ras = Cerrado_Raw_temp,
  Target_res = BA.all.Cerrado.LTDR[[1]],
  Dir = Dir.CerradoExt,
  Shape = CerradoMask,
  Keep_Temporary = TRUE
)

windows(20,10)
Plot_Covs(Covs_ls,CerradoMask)

KrigStart <- Sys.time()
Cerrado_temp_Krig <- krigR(
  Data = Cerrado_Raw_temp, # data we want to krig as a raster object
  Covariates_coarse = Covs_ls[[1]], # training covariate as a raster object
  Covariates_fine = Covs_ls[[2]], # target covariate as a raster object
  Keep_Temporary = FALSE, # we don't want to retain the individually kriged layers on our hard-drive
  Cores = 4, # we want to krig on just one core
  FileName = "Cerrado_Ext.nc", # the file name for our full kriging output
  Dir = Dir.CerradoExt # which directory to save our final input in
)
KrigStop <- Sys.time()
