# CerradoFireRegimes
Repository for the data, code, and products of research aiming to investigate the fire regimes in Cerrado in a large spatio-temporal scale

## Code for analyzes
The code "Script_SpatioTemporal_Fire_Regimes_Cerrado.R" has the main analyzes present in the manuscript. Many of the commented lines are slow to run, but we preferred to maintain them for transparency. However, you can jump these parts and the code will run just fine. We saved many of the outputs as RDS files to save the time for readers.

The code "Script_ERA5.R" downloads the weather (climate) data from the ERA5-land database using package KrigR. You do not need to run this code to perform the main analizes, but we maintained for transparency.

## Data
We provide the data to perform all the analyzes:
- The AVHRR LTDR Fire and MODIS fire products in pixel and grid. We provide the product already gathered and cropped for the Cerrado distribution.
- The predictors to perform the BAMs

## Output
We provide some outputs including the BAMs and Mclust in RDS files. The BAMs predict the fire occurrence probability and extent using weather, vegetation, physical, and anthropogenic variables. The Mclust classified the fire regimes using the monthly aggregated data of fire occurrence and extent from the AVHRR LTDR and MODIS fire products. We also provide the output from these classifications in raster files (for GIS users).

## Animation
The animation directory contains the monthly data (1982-2018) from the AVHRR LTDR Fire products (fire occurrence and extent) cropped for the Cerrado distribution animated in HTML and GIF files. The subdirectories and other files are necessary to run the HTML animations properly. 
