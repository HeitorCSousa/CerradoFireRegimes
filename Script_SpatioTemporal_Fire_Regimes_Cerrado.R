###############
#Load packages#
###############
library(gdalUtils)
library(raster)
library(rgdal)

Sys.setenv(TZ = "America/Sao_Paulo")
rasterOptions(tmpdir=file.path("/Volumes/Extreme SSD/Heitor/BD_SIG/temp"),
              chunksize = 1e+06,maxmemory = 1e+08,
              progress = "text", tmptime = 12,
              todisk = T)



# If you have lots of files then you can make a loop to do all this for you
setwd("Data/AVHRR_LTDR_Fire/Grid")

files<-list.files(pattern = "ESACCI-L4_FIRE-BA-AVHRR-LTDR-", recursive = TRUE)
files

filename <- substr(files,6,11)
filename <- paste0("BA", filename)
filename

BA.all.LTDR<-stack(files)
BA.all.LTDR
windows(10,10)
plot(BA.all.LTDR[[1]])

#########################
#Le shapefile do Cerrado#
#########################
Cerrado<-readOGR("Data/AVHRR_LTDR_Fire/Grid/Cerrado.shp")#Le shapefile
plot(BA.all.LTDR[[9]])
plot(Cerrado,add=T)

#Cortar rasters para Cerrado
BA.all.LTDR.Cerrado<-crop(BA.all.LTDR, extent(Cerrado))
BA.all.LTDR.Cerrado <- mask(BA.all.LTDR.Cerrado, Cerrado)

plot(BA.all.LTDR.Cerrado[[1]])
plot(Cerrado,axes=T,add=T)

writeRaster(BA.all.LTDR.Cerrado,"BA.all.LTDR.Cerrado.tif")
BA.all.LTDR.Cerrado<-stack("BA.all.LTDR.Cerrado.tif")

#Le raster (stack)##########################
#BA.all.LTDR.Cerrado<-stack("BA.all.LTDR.Cerrado.tif")
############################################

# BA.all.LTDR.Cerrado
# crs(BA.all.LTDR.Cerrado)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# Cerrado
# plot(BA.all.LTDR.Cerrado)
# plot(BA.all.LTDR.Cerrado[[1]])
# plot(Cerrado,axes=T,add=T)


#########################################################
#Calcula soma de ?rea queimada por m?s para todo o bioma#
#########################################################
library(forecast)
soma_BA_LTDR<-cellStats(BA.all.LTDR.Cerrado,sum)
soma_BA_LTDR[144]
length(soma_BA_LTDR)
#soma_BA_LTDR<-c(soma_BA_LTDR[1:144],rep(NA,12),soma_BA_LTDR[145:432])


#Analise da serie temporal
ts.soma.BA.LTDR<-ts(soma_BA_LTDR,start=c(1982, 01), end=c(2017, 12), frequency=12)
windows(10,10)
acf(ts.soma.BA.LTDR)
acf(ts.soma.BA.LTDR,lag.max = 20)
pacf(ts.soma.BA.LTDR,lag.max = 20)
ts.soma.BA.LTDR
#ts.soma.BA.LTDR<-ts.soma.BA.LTDR*0.01 #Para transformar em ha multiplicar pelo fator 0.01 (Giglio et al., 2018 - Collection 6 MODIS Burned Area Product User?s Guide Version 1.2)

plot(ts.soma.BA.LTDR,ylab="Área (ha)",xlab="Tempo (anos)")

#Decomposi??o da s?rie temporal
plot(decompose(ts.soma.BA.LTDR))
Queima.decom <- decompose(ts.soma.BA.LTDR,type="multiplicative")
plot(Queima.decom)
Trend <- Queima.decom$trend
Seasonal <- Queima.decom$seasonal
ts.plot(cbind(Trend, Trend * Seasonal), lty = 1:2,ylab="Area (m²)")

decomp<-stl(ts.soma.BA.LTDR,s.window = 13, t.window=13)
summary(decomp)
decomp<-stl(ts.soma.BA.LTDR, s.window = 4,t.window=13)
summary(decomp)
plot(decomp)
autoplot(decomp)
monthplot(ts.soma.BA.LTDR, ylab="Area (m²)", type="h")
monthplot(decomp, choice=c("seasonal"), type="h")
monthplot(decomp, choice="trend", type="h")
monthplot(decomp, choice="remainder", type="h")

fit1 <- StructTS(ts.soma.BA.LTDR, type = "trend")
plot(ts.soma.BA.LTDR)
plot(cbind(fitted(fit1), resids=resid(fit1)))
tsdiag(fit1)

fit2 <- StructTS(ts.soma.BA.LTDR, type = "BSM")
print(fit2)
plot(ts.soma.BA.LTDR)
autoplot(fit2)
monthplot(fit2,type="h")
monthplot(fit2,type="h",choice="slope")
monthplot(fit2,type="h",choice="level")

autoplot(cbind(fitted(fit2), resids=resid(fit2)))
tsdiag(fit2)

fit3 <- ets(ts.soma.BA.LTDR,ic="bic")
summary(fit3)
accuracy(ts.soma.BA.LTDR)
autoplot(fit3)

fit3 %>% forecast(h=240) %>%
  autoplot() 

# Porcentagem da variacao explicada por cada componente da serie temporal
p.sazonal<-(var(decomp$time.series[, 1]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.sazonal

p.tendencia<-(var(decomp$time.series[, 2]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.tendencia

p.residual<-(var(decomp$time.series[, 3]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.residual

p.sazonal+p.tendencia+p.residual
t<-1:length(ts.soma.BA.LTDR)
m1<-lm(soma_BA_LTDR~t)
summary(m1)

#Raster time series - bom para visualizar entre anos x meses
#install.packages("rts",dep=T)
library(rts)

rts.BA.Cerrado.LTDR<-rts(BA.all.LTDR.Cerrado,as.yearmon(time(ts.soma.BA.LTDR)))
rts.BA.Cerrado.LTDR

rts.BA.Cerrado.LTDR.monthly<-apply.yearly(rts.BA.Cerrado.LTDR,mean)
plot(rts.BA.Cerrado.LTDR[[1]])

windows(10,10)
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421))]])#jan
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+1)]])#fev
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+2)]])#mar
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+3)]])#abr
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+4)]])#mai
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+5)]])#jun
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+6)]])#jul
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+7)]])#ago
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+8)]])#set
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+9)]])#out
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+10)]])#nov
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                             289,301,313,325,337,349,361,373,385,397,409,421)+11)]])#dez

########################
#Transformar em bin?rio#
########################
Bbin.all.Cerrado.LTDR<-BA.all.LTDR.Cerrado
Bbin.all.Cerrado.LTDR[Bbin.all.Cerrado.LTDR>0]<-1
Bbin.all.Cerrado.LTDR[Bbin.all.Cerrado.LTDR<0]<-0
writeRaster(Bbin.all.Cerrado.LTDR,"Bbin.all.Cerrado.LTDR.tif",overwrite=T)
plot(Bbin.all.Cerrado.LTDR[[1:12]])
Freq.all.Cerrado.LTDR<-sum(Bbin.all.Cerrado.LTDR) 
plot(Freq.all.Cerrado.LTDR)
plot(Freq.all.Cerrado.LTDR/432)
hist(Freq.all.Cerrado.LTDR,breaks = seq(0,250,5))
summary(Freq.all.Cerrado.LTDR)
abline(v=72,col="red")

setwd("/Volumes/Extreme SSD/Heitor//Doutorado/Analises/Cap1_SpatioTemporal_Fire_Regimes_Cerrado")

#Load rasters
###########################################################
#Our dependent variables (burned area and burned/unburned)#
###########################################################
#ALTDR Fire_cci BA v1.0 products (also called FireCCILT10 for short)
#Binary raster of burned areas (1=burned; 0=unburned)
#Temporal resolution: 1982-2018 (except 1994)
#Spatial resolution: 0.05 degrees (approximately 5 km at the Equator)
Bbin.all.Cerrado.LTDR<-stack("Data/AVHRR_LTDR_Fire/Grid/Bbin.all.Cerrado.LTDR.tif")
Bbin.all.Cerrado.LTDR

#ALTDR Fire_cci BA v1.0 products (also called FireCCILT10 for short)
#Raster of burned areas (0-N=burned area in m²; -1=not observed;-2=unburnable,i.e.water,bare land,urban,etc.)
#Temporal resolution: 1982-2018 (except 1994)
#Spatial resolution: 0.05 degrees (approximately 5 km at the Equator)
BA.all.Cerrado.LTDR<-stack("Data/AVHRR_LTDR_Fire/Grid/BA.all.LTDR.Cerrado.tif")
BA.all.Cerrado.LTDR

#######################
#Independent variables#
#######################
#Climate variables
#From ERA5 Climate Data
#Temporal resolution: 1982-2020
#Spatial resolution: 0.1 degrees (9km)

Cerrado.temp.2m<-stack("Data/Predictors/ERA5/Cerrado_temp_2m.grd")
Cerrado.skin.temp<-stack("Data/Predictors/ERA5/Cerrado_skin_temp.grd")
Cerrado.temp.lv1<-stack("Data/Predictors/ERA5/Cerrado_temp_lv1.grd")
Cerrado.temp.lv2<-stack("Data/Predictors/ERA5/Cerrado_temp_lv2.grd")
Cerrado.sol<-stack("Data/Predictors/ERA5/Cerrado_sol.grd")
Cerrado.water.lv1<-stack("Data/Predictors/ERA5/Cerrado_water_lv1.grd")
Cerrado.water.lv2<-stack("Data/Predictors/ERA5/Cerrado_water_lv2.grd")
Cerrado.precip<-stack("Data/Predictors/ERA5/Cerrado_precip.grd")
Cerrado.wind<-stack("Data/Predictors/ERA5/Cerrado_wind.grd")
Cerrado.RH<-stack("Data/Predictors/ERA5/Cerrado_RH.grd")
Cerrado.total.evap<-stack("Data/Predictors/ERA5/Cerrado_total_evap.grd")
Cerrado.pot.evap<-stack("Data/Predictors/ERA5/Cerrado_pot_evap.grd")
Cerrado.water.deficit<-stack("Data/Predictors/ERA5/Cerrado_water_deficit.grd")

#Much of the next lines are commented because files are large to provide online. We maintained them here for transparency. We can provide such files with prior contact. However, the next analyses are not dependent on these files

# DEM.Cerrado<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/CSR_UFMG/bioma_cerrado_altitude/bioma_cerrado_altitude.tif")
# DEM.Cerrado
# decliv.Cerrado<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/CSR_UFMG/bioma_cerrado_declividade/bioma_cerrado_declividade.tif")
# decliv.Cerrado

######################
#Land use - Mapbiomas#
######################
#Temporal resolution: 1985-2019
#Spatial resolution: 30m (Landsat)
#See legend in /Volumes/Extreme SSD/Heitor/BD_SIG/MapBiomas/Códigos_das_classes_da_legenda_e_paleta_de_cores.pdf
# files<-list.files(path="/Volumes/Extreme SSD/Heitor/BD_SIG/MapBiomas",
#                   pattern = "COLECAO_5_DOWNLOADS_COLECOES_ANUAL_CERRADO_CERRADO-",
#                   recursive = TRUE)
# files
# files<-paste0("/Volumes/Extreme SSD/Heitor/BD_SIG/MapBiomas/",files)
# 
# filename <- substr(files,79,82)
# filename
# 
# cerrado.lu<-stack(files)
# 
# #writeRaster(cerrado.lu,"cerrado.lu.mapbiomas1985_2019.tif")
# 
# rasterOptions(tmpdir=file.path("/Volumes/Extreme SSD/Heitor/BD_SIG/temp"),
#               chunksize = 1e+06,maxmemory = 1e+08,
#               progress = "text", tmptime = 12,
#               todisk = T)
# cerrado.lu<-stack("cerrado.lu.mapbiomas1985_2019.tif")
# 
# beginCluster(4)
# 
# #remotes::install_github("hagc/rasterB")
# require(rasterB)
# cerrado.lu.5km<-list()
# for(i in 2:35){
#   cerrado.lu.5km[i]<-aggregateB(cerrado.lu[[i]],186,fun="modal")
#   print(cerrado.lu.5km)
#   writeRaster(cerrado.lu.5km[[i]],paste0("cerrado.lu.5km",i,".grd"))
# }
# 
# plot(cerrado.lu.5km[[35]])
# cerrado.lu.5km.full<-stack(cerrado.lu.5km)

#writeRaster(cerrado.lu.5km.full,"cerrado.lu.5km.full.grd")
#endCluster()

# cerrado.lu.5km.full<-stack("cerrado.lu.5km.full.grd")
# cerrado.lu.5km.full
# 
# cerrado.lu.25km.full<-aggregate(cerrado.lu.5km.full,5,fun="modal")
# 
# cerrado.lu.25km.full<-resample(cerrado.lu.25km.full,Bbin.all.Cerrado.LTDR,method="ngb",
#                                filename="cerrado.lu.25km.full.grd",overwrite=T)
cerrado.lu.25km.full<-stack("cerrado.lu.25km.full.grd")
windows(10,10)
plot(cerrado.lu.25km.full[[c(1,5,10,15,20,25,30,35)]])

cerrado.lu.df<-as.data.frame(cerrado.lu.25km.full,na.rm=T,xy=T)

library(lubridate)
library(tidyverse)
dates.lu<-c(seq(ym("1985-01"),ym("2019-12"),by="years"))

names(cerrado.lu.df)<-c("x","y",as.character(dates.lu))

cerrado.lu.dflong<-gather(cerrado.lu.df,Date,land_use,-x,-y,factor_key = T)
cerrado.lu.dflong$Year.x<-year(cerrado.lu.dflong$Date)

############
#Vegetation#
############
#Vegetation (IBGE) - 1:250.000

# BrasilVeg <- readOGR(Dir.Shapes,"Brasil_vege_area", verbose = FALSE) # read
# CerradoVeg<-crop(BrasilVeg, extent(CerradoMask))
# shapefile(CerradoVeg,"CerradoVeg.shp")
# CerradoVeg <- readOGR("CerradoVeg.shp", verbose = FALSE) # read
# windows(10,10)
# plot(CerradoVeg)
# 
# df.Cerrado<-as.data.frame(Bbin.all.Cerrado.LTDR[[1]],xy = TRUE)
# completos <- complete.cases(df.Cerrado[,c("Bbin.all.Cerrado.LTDR.1")])
# df.Cerrado <- droplevels(df.Cerrado[completos,])
# xy.Cerrado<-df.Cerrado[1:2]
# coordinates(xy.Cerrado) <- ~ x + y
# proj4string(xy.Cerrado) <- crs(CerradoVeg)
# test <- data.frame(xx=over(xy.Cerrado,CerradoVeg))
# head(test)
# nm_uveg_Cerrado<-as.data.frame(cbind(test$xx.nm_uveg,test$xx.leg_uveg,xy.Cerrado@coords))
# levels(as.factor(nm_uveg_Cerrado$V1))
# head(nm_uveg_Cerrado)
# write.table(nm_uveg_Cerrado,"nm_uveg_Cerrado.txt",sep="\t")
nm_uveg_Cerrado<-read.table("nm_uveg_Cerrado.txt")
head(nm_uveg_Cerrado)
nm_uveg_Cerrado$V2<-as.factor(nm_uveg_Cerrado$V2)
str(nm_uveg_Cerrado)
summary(nm_uveg_Cerrado$V2)
nm_uveg_Cerrado<-cbind(nm_uveg_Cerrado[,3],nm_uveg_Cerrado[,4],nm_uveg_Cerrado[,2])
unique(nm_uveg_Cerrado[,1])
unique(nm_uveg_Cerrado[,2])

head(nm_uveg_Cerrado)

VegCerrado.r<-rasterFromXYZ(nm_uveg_Cerrado,crs = crs(Bbin.all.Cerrado.LTDR))
VegCerrado.r
windows(10,10)
plot(VegCerrado.r)

#####################
#Aboveground Biomass#
#####################
#Maurizio Santoro & Oliver Cartus
#Version 3.0 - 31 May 2018
#Santoro, M., Cartus, O., Mermoz, S., Bouvet, A., Le Toan, T., Carvalhais, N., Rozendaal, D., Herold, M., Avitabile, V., Quegan, S., Carreiras, J., Rauste, Y., Balzter, H., Schmullius, C., Seifert, F.M., 2018, GlobBiomass global above-ground biomass and growing stock volume datasets, available on-line at http://globbiomass.org/products/global-mapping
#Reference system: Lat-long, WGS-84
#Pixel spacing: The GSV and AGB estimates are provided with a pixel spacing of 0.0008888° (roughly
#corresponding to 100 m at the Equator).
#unit: tons/ha i.e., Mg/ha
#Definition: the mass, expressed as oven-dry weight of the woody parts (stem, bark, 
#branches and twigs) of all living trees excluding stump and roots
# biomass.Cerrado<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/N00W060_agb/N00W060_agb.tif")
# windows(10,10)
# plot(biomass.Cerrado)
# plot(CerradoMask,add=T)
# 
# biomass.Cerrado<-crop(biomass.Cerrado, extent(CerradoMask))
# biomass.Cerrado <- mask(biomass.Cerrado, CerradoMask)
# writeRaster(biomass.Cerrado,"Biomass_Cerrado.grd")
# biomass.Cerrado<-raster("Biomass_Cerrado.grd")
# windows(10,10)
# plot(biomass.Cerrado)

###############
#Canopy height#
###############
#Simard et al. (2011)
#1 km spatial resolution
#2005 data
# canopy.height<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/Simard_Pinto_3DGlobalVeg_JGR.tif/Simard_Pinto_3DGlobalVeg_L3C.tif")
# canopy.height.Cerrado<-crop(canopy.height, extent(CerradoMask))
# canopy.height.Cerrado <- mask(canopy.height.Cerrado, CerradoMask)
# writeRaster(canopy.height.Cerrado,"canopy.height.Cerrado.grd")
# canopy.height.Cerrado<-raster("canopy.height.Cerrado.grd")
# 
# windows(10,10)
# plot(canopy.height.Cerrado)

####################
#Tree cover percent#
####################
#Hansen et al. (2013)
#spatial resolution of 1 arc-second per pixel, 
#or approximately 30 meters per pixel at the equator
# tree.cover.1<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/Hansen_GFC_2020/Hansen_GFC-2020-v1.8_treecover2000_00N_050W.tif")
# tree.cover.2<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/Hansen_GFC_2020/Hansen_GFC-2020-v1.8_treecover2000_00N_060W.tif")
# tree.cover.3<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/Hansen_GFC_2020/Hansen_GFC-2020-v1.8_treecover2000_10S_050W.tif")
# tree.cover.4<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/Hansen_GFC_2020/Hansen_GFC-2020-v1.8_treecover2000_10S_060W.tif")
# tree.cover.5<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/Hansen_GFC_2020/Hansen_GFC-2020-v1.8_treecover2000_10S_070W.tif")
# tree.cover.6<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/Hansen_GFC_2020/Hansen_GFC-2020-v1.8_treecover2000_20S_050W.tif")
# tree.cover.7<-raster("/Volumes/Extreme SSD/Heitor/BD_SIG/Hansen_GFC_2020/Hansen_GFC-2020-v1.8_treecover2000_20S_060W.tif")
# 
# tree.cover.full<-mosaic(tree.cover.1,
#                         tree.cover.2,
#                         tree.cover.3,
#                         tree.cover.4,
#                         tree.cover.5,
#                         tree.cover.6,
#                         tree.cover.7,fun=mean)
# beginCluster(4)
# tree.cover.Cerrado<-crop(tree.cover.full, extent(CerradoMask))
# tree.cover.Cerrado <- mask(tree.cover.Cerrado, CerradoMask)
# writeRaster(tree.cover.Cerrado,"tree.cover.Cerrado.grd")
# tree.cover.Cerrado<-raster("tree.cover.Cerrado.grd")
# windows(10,10)
# plot(tree.cover.Cerrado)
#endCluster()

#Roads#
#######
# CerradoRoads <- readOGR(Dir.Shapes,"Roads_Cerrado", verbose = FALSE) # read
# 
# # CerradoRoads.dist<-raster("Cerrado_Roads_Dist.tif")
# # CerradoRoads.dist<-crop(CerradoRoads.dist, extent(CerradoMask))
# # CerradoRoads.dist <- mask(CerradoRoads.dist, CerradoMask)
# # writeRaster(CerradoRoads.dist,"CerradoRoads.dist.grd")
# CerradoRoads.dist<-raster("CerradoRoads.dist.grd")
# 
# windows(10,10)
# plot(CerradoRoads.dist)
#plot(CerradoRoads,add=T)


#Aggregate to match resolutions - target = 25km#
################################################
# DEM.Cerrado
# decliv.Cerrado
# beginCluster(4)
# DEM.Cerrado.25km<-aggregate(DEM.Cerrado,250)
# decliv.Cerrado.25km<-aggregate(decliv.Cerrado,250)
# DEM.Cerrado.25km<-projectRaster(DEM.Cerrado.25km,Bbin.all.Cerrado.LTDR)
# decliv.Cerrado.25km<-projectRaster(decliv.Cerrado.25km,Bbin.all.Cerrado.LTDR)
# DEM.Cerrado.25km
# decliv.Cerrado.25km
# writeRaster(DEM.Cerrado.25km,"DEM.Cerrado.25km.grd")
# writeRaster(decliv.Cerrado.25km,"decliv.Cerrado.25km.grd")
# 
# 
# biomass.Cerrado
# res(Bbin.all.Cerrado.LTDR)/res(biomass.Cerrado)
# biomass.Cerrado.25km<-aggregate(biomass.Cerrado,281)
# biomass.Cerrado.25km<-resample(biomass.Cerrado.25km,Bbin.all.Cerrado.LTDR,
#                               filename="biomass.Cerrado.25km.grd")
# biomass.Cerrado.25km
# 
# canopy.height.Cerrado
# res(Bbin.all.Cerrado.LTDR)/res(canopy.height.Cerrado)
# canopy.height.Cerrado.25km<-aggregate(canopy.height.Cerrado,30)
# canopy.height.Cerrado.25km<-resample(canopy.height.Cerrado.25km,Bbin.all.Cerrado.LTDR,
#                             filename="canopy.height.Cerrado.25km.grd")
# canopy.height.Cerrado.25km
# 
# tree.cover.Cerrado
# res(Bbin.all.Cerrado.LTDR)/res(tree.cover.Cerrado)
# tree.cover.Cerrado.25km<-aggregate(tree.cover.Cerrado,1000)
# tree.cover.Cerrado.25km<-resample(tree.cover.Cerrado.25km,Bbin.all.Cerrado.LTDR,
#                          filename = "tree.cover.Cerrado.25km.grd")
# tree.cover.Cerrado.25km
# 
# CerradoRoads.dist
# res(Bbin.all.Cerrado.LTDR)/res(CerradoRoads.dist)
# CerradoRoads.dist.25km<-aggregate(CerradoRoads.dist,282)
# CerradoRoads.dist.25km<-resample(CerradoRoads.dist.25km,Bbin.all.Cerrado.LTDR,
#                                 filename = "CerradoRoads.dist.25km.grd")
# CerradoRoads.dist.25km
# 
# VegCerrado.r
# VegCerrado.r<-raster::resample(VegCerrado.r,Bbin.all.Cerrado.LTDR,method="ngb",
#                                  filename = "VegCerrado.r.25km.grd",overwrite=T)


DEM.Cerrado.25km<-raster("DEM.Cerrado.25km.grd")
decliv.Cerrado.25km<-raster("decliv.Cerrado.25km.grd")
biomass.Cerrado.25km<-raster("biomass.Cerrado.25km.grd")
canopy.height.Cerrado.25km<-raster("canopy.height.Cerrado.25km.grd")
tree.cover.Cerrado.25km<-raster("tree.cover.Cerrado.25km.grd")
CerradoRoads.dist.25km<-raster("CerradoRoads.dist.25km.grd")
VegCerrado.r<-raster("VegCerrado.r.25km.grd")

##############################################
#Check collinearity problems with covariables#
##############################################
library(usdm)
cov.stack<-stack(DEM.Cerrado.25km,decliv.Cerrado.25km,biomass.Cerrado.25km,canopy.height.Cerrado.25km,
      tree.cover.Cerrado.25km,CerradoRoads.dist.25km,VegCerrado.r)
names(cov.stack)<-c("DEM.Cerrado","decliv.Cerrado","biomass.Cerrado",
                    "canopy.height.Cerrado","tree.cover.Cerrado",
                    "CerradoRoads.dist","VegCerrado.r")
cov.stack
vifcor(cov.stack,th=0.8)

#subsample 5% of pixels and calculate pairwise correlations
library(corrplot)
cor.m<- cor(sampleRandom(cov.stack, size= ncell(cov.stack[[1]]) * 0.05 ), 
            method = "spearman")
cor.m
#plot correlation matrix
windows(10,10)
df <- corrplot(cor.m, method = "number")

#Conclusion: Exclude tree cover from analysis
fixed.covs.df<-as.data.frame(cov.stack,na.rm=T,xy=T)
saveRDS(fixed.covs.df,"fixed.covs.df.rds")


#Aggregate to match resolutions - target = 25km
beginCluster(4)
Cerrado.temp.2m<-resample(Cerrado.temp.2m,Bbin.all.Cerrado.LTDR,
                                  filename = "Cerrado.temp.2m.25km.grd",overwrite=T)
Cerrado.temp.2m

Cerrado.skin.temp<-resample(Cerrado.skin.temp,Bbin.all.Cerrado.LTDR,
                          filename = "Cerrado.skin.temp.25km.grd",overwrite=T)
Cerrado.skin.temp

Cerrado.temp.lv1<-resample(Cerrado.temp.lv1,Bbin.all.Cerrado.LTDR,
                            filename = "Cerrado.temp.lv1.25km.grd",overwrite=T)
Cerrado.temp.lv1

Cerrado.temp.lv2<-resample(Cerrado.temp.lv2,Bbin.all.Cerrado.LTDR,
                           filename = "Cerrado.temp.lv2.25km.grd",overwrite=T)
Cerrado.temp.lv2

Cerrado.sol<-resample(Cerrado.sol,Bbin.all.Cerrado.LTDR,
                           filename = "Cerrado.sol.25km.grd",overwrite=T)
Cerrado.sol

Cerrado.water.lv1<-resample(Cerrado.water.lv1,Bbin.all.Cerrado.LTDR,
                            filename = "Cerrado.water.lv1.25km.grd",overwrite=T)
Cerrado.water.lv1

Cerrado.water.lv2<-resample(Cerrado.water.lv2,Bbin.all.Cerrado.LTDR,
                            filename = "Cerrado.water.lv2.25km.grd",overwrite=T)
Cerrado.water.lv2

Cerrado.precip<-resample(Cerrado.precip,Bbin.all.Cerrado.LTDR,
                            filename = "Cerrado.precip.25km.grd",overwrite=T)
Cerrado.precip

Cerrado.wind<-resample(Cerrado.wind,Bbin.all.Cerrado.LTDR,
                         filename = "Cerrado.wind.25km.grd",overwrite=T)
Cerrado.wind

Cerrado.RH<-resample(Cerrado.RH,Bbin.all.Cerrado.LTDR,
                       filename = "Cerrado.RH.25km.grd",overwrite=T)
Cerrado.RH

Cerrado.total.evap<-resample(Cerrado.total.evap,Bbin.all.Cerrado.LTDR,
                     filename = "Cerrado.total.evap.25km.grd",overwrite=T)
Cerrado.total.evap

Cerrado.pot.evap<-resample(Cerrado.pot.evap,Bbin.all.Cerrado.LTDR,
                             filename = "Cerrado.pot.evap.25km.grd",overwrite=T)
Cerrado.pot.evap

Cerrado.water.deficit<-resample(Cerrado.water.deficit,Bbin.all.Cerrado.LTDR,
                             filename = "Cerrado.water.deficit.25km.grd",overwrite=T)
Cerrado.water.deficit

#Transform to data.frame (wide format)
Cerrado.temp.2m.df<-as.data.frame(Cerrado.temp.2m,na.rm=T,xy=T)
Cerrado.skin.temp.df<-as.data.frame(Cerrado.skin.temp,na.rm=T,xy=T)
Cerrado.temp.lv1.df<-as.data.frame(Cerrado.temp.lv1,na.rm=T,xy=T)
Cerrado.temp.lv2.df<-as.data.frame(Cerrado.temp.lv2,na.rm=T,xy=T)
Cerrado.sol.df<-as.data.frame(Cerrado.sol,na.rm=T,xy=T)
Cerrado.water.lv1.df<-as.data.frame(Cerrado.water.lv1,na.rm=T,xy=T)
Cerrado.water.lv2.df<-as.data.frame(Cerrado.water.lv2,na.rm=T,xy=T)
Cerrado.precip.df<-as.data.frame(Cerrado.precip,na.rm=T,xy=T)
Cerrado.wind.df<-as.data.frame(Cerrado.wind,na.rm=T,xy=T)
Cerrado.RH.df<-as.data.frame(Cerrado.RH,na.rm=T,xy=T)
Cerrado.total.evap.df<-as.data.frame(Cerrado.total.evap,na.rm=T,xy=T)
Cerrado.pot.evap.df<-as.data.frame(Cerrado.pot.evap,na.rm=T,xy=T)
Cerrado.water.deficit.df<-as.data.frame(Cerrado.water.deficit,na.rm=T,xy=T)


endCluster()

#Transform to long format
library("dplyr")
library("tidyr")
library(lubridate)
dates.study<-c(seq(ym("1982-01"),ym("1993-12"),by="months"),
               seq(ym("1995-01"),ym("2018-12"),by="months"))
months.study<-month(dates.study,label = T, abbr = T)
years.study<-year(dates.study)

dates.climate<-c(seq(ym("1982-01"),ym("2020-12"),by="months"))

names(Cerrado.temp.2m.df)<-c(as.character(dates.climate),"x","y")

Cerrado.temp.2m.dflong<-gather(Cerrado.temp.2m.df,Date,temp_2m,-x,-y,factor_key = T)
Cerrado.temp.2m.dflong$Month<-month(Cerrado.temp.2m.dflong$Date,label=T,abbr=T)
Cerrado.temp.2m.dflong$Year<-year(Cerrado.temp.2m.dflong$Date)

Cerrado.skin.temp.dflong<-gather(Cerrado.skin.temp.df,Index,skin_temp,-x,-y,factor_key = T)
Cerrado.temp.lv1.dflong<-gather(Cerrado.temp.lv1.df,Index,temp_lv1,-x,-y,factor_key = T)
Cerrado.temp.lv2.dflong<-gather(Cerrado.temp.lv2.df,Index,temp_lv2,-x,-y,factor_key = T)
Cerrado.sol.dflong<-gather(Cerrado.sol.df,Index,sol,-x,-y,factor_key = T)
# Cerrado.evap.soil.dflong<-gather(Cerrado.evap.soil.df,Index,evap_soil,-x,-y,factor_key = T)
# Cerrado.evap.veg.dflong<-gather(Cerrado.evap.veg.df,Index,evap_veg,-x,-y,factor_key = T)
Cerrado.water.lv1.dflong<-gather(Cerrado.water.lv1.df,Index,water_lv1,-x,-y,factor_key = T)
Cerrado.water.lv2.dflong<-gather(Cerrado.water.lv2.df,Index,water_lv2,-x,-y,factor_key = T)
Cerrado.precip.dflong<-gather(Cerrado.precip.df,Index,precip,-x,-y,factor_key = T)
Cerrado.wind.dflong<-gather(Cerrado.wind.df,Index,wind,-x,-y,factor_key = T)
Cerrado.RH.dflong<-gather(Cerrado.RH.df,Index,RH,-x,-y,factor_key = T)
Cerrado.total.evap.dflong<-gather(Cerrado.total.evap.df,Index,total.evap,-x,-y,factor_key = T)
Cerrado.pot.evap.dflong<-gather(Cerrado.pot.evap.df,Index,pot.evap,-x,-y,factor_key = T)
Cerrado.water.deficit.dflong<-gather(Cerrado.water.deficit.df,Index,water.deficit,-x,-y,factor_key = T)


cerrado.climate.df<-cbind(Cerrado.temp.2m.dflong,
                          Cerrado.skin.temp.dflong[,4],
                          Cerrado.temp.lv1.dflong[,4],
                          Cerrado.temp.lv2.dflong[,4],
                          Cerrado.sol.dflong[,4],
                          Cerrado.water.lv1.dflong[,4],
                          Cerrado.water.lv2.dflong[,4],
                          Cerrado.precip.dflong[,4],
                          Cerrado.wind.dflong[,4],
                          Cerrado.RH.dflong[,4],
                          Cerrado.total.evap.dflong[,4],
                          Cerrado.pot.evap.dflong[,4],
                          Cerrado.water.deficit.dflong[,4])
head(cerrado.climate.df)
names(cerrado.climate.df)<-c("x","y","Date","temp_2m","Month","Year",
                             "skin_temp","temp_lv1","temp_lv2","sol",
                             "water_lv1","water_lv2",
                             "precip","wind",
                             "RH","total.evap","pot.evap","water.deficit")

#Check collinearity problems with climatic variables#
#####################################################

vifstep(cerrado.climate.df[sample(1:length(cerrado.climate.df[,4]),
                                 length(cerrado.climate.df[,4])*0.05),c(4,7:18)])

vifcor(cerrado.climate.df[sample(1:length(cerrado.climate.df[,4]),
                                 length(cerrado.climate.df[,4])*0.05),c(4,7:18)])


cor.climate<- cor(cerrado.climate.df[sample(1:length(cerrado.climate.df[,4]),
                                      length(cerrado.climate.df[,4])*0.05),
                               c(4,7:18)], 
            method = "spearman")
cor.climate
#plot correlation matrix
windows(10,10)
quartz(8,8)
corrplot(cor.climate, method = "number")
write.table(cerrado.climate.df,"cerrado.climate.df.txt",sep="/t")
saveRDS(cerrado.climate.df,"cerrado.climate.df.rds")

#Conclusion: exclude temp_lv2 water_lv2 skin_temp temp_lv1 water.deficit 

#################
#Hovmoller plots#
#################
library("animation")
library("dplyr")
library("ggplot2")
library("gstat")
library("maps")
library(devtools)
#install_github("andrewzm/STRbook")
library("STRbook")
library("grid")
library("gridExtra")

# beginCluster(4)
# 
# Bbin.all.Cerrado.LTDR.df<-as.data.frame(Bbin.all.Cerrado.LTDR,na.rm=T,xy=T)
# 
# BA.all.Cerrado.LTDR.df<-as.data.frame(BA.all.Cerrado.LTDR,na.rm=T,xy=T)
# 
# saveRDS(Bbin.all.Cerrado.LTDR.df,"Bbin.all.Cerrado.LTDR.df.rds")
# 
# saveRDS(BA.all.Cerrado.LTDR.df,"BA.all.Cerrado.LTDR.df.rds")
# 
# Bbin.all.Cerrado.LTDR.df<-readRDS("Bbin.all.Cerrado.LTDR.df.rds")
# 
# BA.all.Cerrado.LTDR.df<-readRDS("BA.all.Cerrado.LTDR.df.rds")
# 
# names(Bbin.all.Cerrado.LTDR.df)<-c(as.character(dates.study),"x","y")
# 
# Bbin.all.Cerrado.LTDR.dflong<-gather(Bbin.all.Cerrado.LTDR.df,Date,Freq_fire,-x,-y,factor_key = T)
# Bbin.all.Cerrado.LTDR.dflong$Month<-month(Bbin.all.Cerrado.LTDR.dflong$Date,label=T,abbr=T)
# Bbin.all.Cerrado.LTDR.dflong$Year<-year(Bbin.all.Cerrado.LTDR.dflong$Date)
# Bbin.all.Cerrado.LTDR.dflong$t<-as.integer(Bbin.all.Cerrado.LTDR.dflong$Date)
# 
# BA.all.Cerrado.LTDR.dflong<-gather(BA.all.Cerrado.LTDR.df,Date,Extent_fire,-x,-y,factor_key = T)
# 
# cerrado.fire.df<-cbind(Bbin.all.Cerrado.LTDR.dflong,
#                         BA.all.Cerrado.LTDR.dflong[,4])
# names(cerrado.fire.df)<-c("x","y","Date","Freq_fire","Month","Year","t","Exten_Fire")
# cerrado.fire.df$Exten_Fire[cerrado.fire.df$Exten_Fire<0]<-NA
# 
# saveRDS(cerrado.fire.df,"cerrado.fire.df.rds")
# endCluster()
cerrado.fire.df<-readRDS("cerrado.fire.df.rds")

extent(BA.all.Cerrado.LTDR)
extent(Bbin.all.Cerrado.LTDR)

spat_Freqfire <- group_by(cerrado.fire.df, y, x) %>%    # group by lon-lat
  summarise(Freq_fire = mean(Freq_fire))     # mean for each lon-lat

spat_Extenfire <- group_by(cerrado.fire.df, y, x) %>%    # group by lon-lat
  summarise(Exten_Fire = mean(Exten_Fire))     # mean for each lon-lat

## ------------------------------------------------------------------------
lat_Freqfire_means <- ggplot(spat_Freqfire) +
  geom_point(aes(y, Freq_fire)) +
  xlab("Latitude (deg)") +
  ylab("Fire probability") + theme_bw()

lon_Freqfire_means <- ggplot(spat_Freqfire) +
  geom_point(aes(x, Freq_fire)) +
  xlab("Longitude (deg)") +
  ylab("Fire probability") + theme_bw()

lat_Extenfire_means <- ggplot(spat_Extenfire) +
  geom_point(aes(y, Exten_Fire)) +
  xlab("Latitude (deg)") +
  ylab("Fire extent") + theme_bw()

lon_Extenfire_means <- ggplot(spat_Extenfire) +
  geom_point(aes(x, Exten_Fire)) +
  xlab("Longitude (deg)") +
  ylab("Fire extent") + theme_bw()

windows(10,10)
quartz(8,8)
lat_Freqfire_means
windows(10,10)
lon_Freqfire_means

windows(10,10)
lat_Extenfire_means
windows(10,10)
lon_Extenfire_means

## ------------------------------------------------------------------------
Freq_fire_av <- group_by(cerrado.fire.df, Date) %>%
  summarise(meanFreqFire = mean(Freq_fire))

windows(10,10)
plot(Freq_fire_av$meanFreqFire,type="l")

Exten_fire_av <- group_by(cerrado.fire.df, Date) %>%
  summarise(meanExtenFire = mean(Exten_Fire,na.rm=T))

windows(10,10)
plot(Exten_fire_av$meanExtenFire,type="l")

###################
#Frequency of fire#
###################
(lim_lat <- extent(Bbin.all.Cerrado.LTDR)[3:4])        # latitude range
(lim_long <- extent(Bbin.all.Cerrado.LTDR)[1:2])           # longitude range
(lim_t <- c(1,nlayers(Bbin.all.Cerrado.LTDR)))           # time range
lat_axis <- seq(lim_lat[1],       # latitude axis
                lim_lat[2],
                length=25)
long_axis <- seq(lim_long[1],       # latitude axis
                 lim_long[2],
                 length=25)
t_axis <- seq(lim_t[1],           # time axis
              lim_t[2],
              length=nlayers(Bbin.all.Cerrado.LTDR))

lat_t_grid <- expand.grid(lat = lat_axis,
                          t = t_axis)

long_t_grid <- expand.grid(long = long_axis,
                           t = t_axis)

## -----------------------------------------------------------

fire_grid <- cerrado.fire.df


## ------------------------------------------------------------------------
Freqfire_lat_Hov <- group_by(fire_grid, y, Date) %>%
  summarise(Freq_fire = mean(Freq_fire))

Freqfire_long_Hov <- group_by(fire_grid, x, Date) %>%
  summarise(Freq_fire = mean(Freq_fire))

Freqfire_lat_Hov$t<-as.integer(Freqfire_lat_Hov$Date)
Freqfire_long_Hov$t<-as.integer(Freqfire_long_Hov$Date)

## ------------------------------------------------------------------------
Hovmoller_lat <- ggplot(Freqfire_lat_Hov) +            # take data
  geom_tile(aes(x = y, y = t, fill = Freq_fire)) + # plot
  fill_scale(name = "Fire probability") +     # add color scale
  scale_y_reverse() +             # rev y scale
  ylab("Time (months)") +     # add y label
  xlab("Latitude (degrees)") +    # add x label
  theme_bw()                      # change theme

Hovmoller_long <- ggplot(Freqfire_long_Hov) +            # take data
  geom_tile(aes(x = x, y = t, fill = Freq_fire)) + # plot
  fill_scale(name = "Fire probability") +     # add color scale
  scale_y_reverse() +             # rev y scale
  ylab("Time (months)") +     # add y label
  xlab("Longitude (degrees)") +    # add x label
  theme_bw()                      # change theme

windows(12,10)
Hovmoller_lat

windows(12,10)
Hovmoller_long

## ------------------------------------------------------------------------

Extenfire_lat_Hov <- group_by(fire_grid, y, Date) %>%
  summarise(Exten_Fire = mean(Exten_Fire))

Extenfire_long_Hov <- group_by(fire_grid, x, Date) %>%
  summarise(Exten_Fire = mean(Exten_Fire))

Extenfire_lat_Hov$t<-as.integer(Extenfire_lat_Hov$Date)
Extenfire_long_Hov$t<-as.integer(Extenfire_long_Hov$Date)

## ------------------------------------------------------------------------
Extenfire_Hovmoller_lat <- ggplot(Extenfire_lat_Hov) +            # take data
  geom_tile(aes(x = y, y = t, fill = Exten_Fire)) + # plot
  fill_scale(name = "Fire extent") +     # add color scale
  scale_y_reverse() +             # rev y scale
  ylab("Time (months)") +     # add y label
  xlab("Latitude (degrees)") +    # add x label
  theme_bw()                      # change theme

Extenfire_Hovmoller_long <- ggplot(Extenfire_long_Hov) +            # take data
  geom_tile(aes(x = x, y = t, fill = Exten_Fire)) + # plot
  fill_scale(name = "Fire extent") +     # add color scale
  scale_y_reverse() +             # rev y scale
  ylab("Time (months)") +     # add y label
  xlab("Longitude (degrees)") +    # add x label
  theme_bw()                      # change theme

windows(12,10)
Extenfire_Hovmoller_lat

windows(12,10)
Extenfire_Hovmoller_long

###########
#rasterVis#
###########

library(rasterVis)
library(dichromat)
library(zoo)
tt <- as.yearmon(dates.study)
Bbin.all.Cerrado.LTDR<-setZ(Bbin.all.Cerrado.LTDR,tt)
names(Bbin.all.Cerrado.LTDR) <- as.character(tt)

BA.all.Cerrado.LTDR<-setZ(BA.all.Cerrado.LTDR,tt)
names(BA.all.Cerrado.LTDR)<-as.character(tt)

Freq.fire.LTDR.Cerrado<-stack("Data/AVHRR_LTDR_Fire/Pixel/Freq.fire.LTDR.Cerrado.grd")
Freq.fire.LTDR.Cerrado<-setZ(Freq.fire.LTDR.Cerrado,
                             seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month'))
names(Freq.fire.LTDR.Cerrado)<-as.character(month.abb)

FireExtent.LTDR.Cerrado<-stack("Data/AVHRR_LTDR_Fire/Pixel/FireExtent.Cerrado.LTDR.monthly.grd")
FireExtent.LTDR.Cerrado<-setZ(FireExtent.LTDR.Cerrado,
                             seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month'))
names(FireExtent.LTDR.Cerrado)<-as.character(month.abb)

levelplot(Bbin.all.Cerrado.LTDR[[1:12]],par.settings=BuRdTheme)

levelplot(Bbin.all.Cerrado.LTDR[[9]],par.settings=BuRdTheme)

levelplot(BA.all.Cerrado.LTDR[[1:12]],par.settings=BuRdTheme)

levelplot(BA.all.Cerrado.LTDR[[9]],par.settings=BuRdTheme)

myTheme <- BuRdTheme()
myTheme$panel.background$col = 'gray' 

#windows(10,10)
quartz(8,8)
levelplot(Freq.fire.LTDR.Cerrado,par.settings=BuRdTheme) + latticeExtra::layer(sp.polygons(Cerrado))

#MODIS
setwd("Data/MODIS_Fire/MCD64A1")
Fire.MCD64A1.month<-stack("Freq.fire.soma.MCD64A1.Cerrado.grd")
setwd("/Users/heito/Documents/CerradoFireRegimes")
Fire.MCD64A1.month<-setZ(Fire.MCD64A1.month,
                              seq(as.Date('2008-01-15'), as.Date('2008-12-15'), 'month'))
names(Fire.MCD64A1.month)<-as.character(month.abb)


quartz(8,8)
levelplot(Fire.MCD64A1.month/20,par.settings=BuRdTheme) + latticeExtra::layer(sp.polygons(Cerrado))

# windows(10,10)
quartz(8,8)

LogFireExtent.LTDR.Cerrado<- log10(FireExtent.LTDR.Cerrado)
LogFireExtent.LTDR.Cerrado[LogFireExtent.LTDR.Cerrado<0] <- NA

levelplot(LogFireExtent.LTDR.Cerrado,par.settings=infernoTheme) + latticeExtra::layer(sp.polygons(Cerrado))

#Land use
cerrado.lu.25km.full[cerrado.lu.25km.full==5]<-3
cerrado.lu.25km.1985<-ratify(cerrado.lu.25km.full[[1]])
cerrado.lu.25km.1986<-ratify(cerrado.lu.25km.full[[2]])
cerrado.lu.25km.1987<-ratify(cerrado.lu.25km.full[[3]])
cerrado.lu.25km.1988<-ratify(cerrado.lu.25km.full[[4]])
cerrado.lu.25km.1989<-ratify(cerrado.lu.25km.full[[5]])
cerrado.lu.25km.1990<-ratify(cerrado.lu.25km.full[[6]])
cerrado.lu.25km.1991<-ratify(cerrado.lu.25km.full[[7]])
cerrado.lu.25km.1992<-ratify(cerrado.lu.25km.full[[8]])
cerrado.lu.25km.1993<-ratify(cerrado.lu.25km.full[[9]])
cerrado.lu.25km.1994<-ratify(cerrado.lu.25km.full[[10]])
cerrado.lu.25km.1995<-ratify(cerrado.lu.25km.full[[11]])
cerrado.lu.25km.1996<-ratify(cerrado.lu.25km.full[[12]])
cerrado.lu.25km.1997<-ratify(cerrado.lu.25km.full[[13]])
cerrado.lu.25km.1998<-ratify(cerrado.lu.25km.full[[14]])
cerrado.lu.25km.1999<-ratify(cerrado.lu.25km.full[[15]])
cerrado.lu.25km.2000<-ratify(cerrado.lu.25km.full[[16]])
cerrado.lu.25km.2001<-ratify(cerrado.lu.25km.full[[17]])
cerrado.lu.25km.2002<-ratify(cerrado.lu.25km.full[[18]])
cerrado.lu.25km.2003<-ratify(cerrado.lu.25km.full[[19]])
cerrado.lu.25km.2004<-ratify(cerrado.lu.25km.full[[20]])
cerrado.lu.25km.2005<-ratify(cerrado.lu.25km.full[[21]])
cerrado.lu.25km.2006<-ratify(cerrado.lu.25km.full[[22]])
cerrado.lu.25km.2007<-ratify(cerrado.lu.25km.full[[23]])
cerrado.lu.25km.2008<-ratify(cerrado.lu.25km.full[[24]])
cerrado.lu.25km.2009<-ratify(cerrado.lu.25km.full[[25]])
cerrado.lu.25km.2010<-ratify(cerrado.lu.25km.full[[26]])
cerrado.lu.25km.2011<-ratify(cerrado.lu.25km.full[[27]])
cerrado.lu.25km.2012<-ratify(cerrado.lu.25km.full[[28]])
cerrado.lu.25km.2013<-ratify(cerrado.lu.25km.full[[29]])
cerrado.lu.25km.2014<-ratify(cerrado.lu.25km.full[[30]])
cerrado.lu.25km.2015<-ratify(cerrado.lu.25km.full[[31]])
cerrado.lu.25km.2016<-ratify(cerrado.lu.25km.full[[32]])
cerrado.lu.25km.2017<-ratify(cerrado.lu.25km.full[[33]])
cerrado.lu.25km.2018<-ratify(cerrado.lu.25km.full[[34]])
cerrado.lu.25km.2019<-ratify(cerrado.lu.25km.full[[35]])

levels(cerrado.lu.25km.1993)[[1]]==levels(cerrado.lu.25km.1999)[[1]]
levels(cerrado.lu.25km.2000)[[1]]==levels(cerrado.lu.25km.2004)[[1]]
levels(cerrado.lu.25km.2005)[[1]]==levels(cerrado.lu.25km.2016)[[1]]
levels(cerrado.lu.25km.2005)[[1]]==levels(cerrado.lu.25km.2019)[[1]]


rat <- levels(cerrado.lu.25km.2019)[[1]]
rat$landcover <- c("Unclassified",
                   "Forest",#3
                   "Savanna",#4
                   "Forest plantation",#9
                   "Grassland",#12
                   "Pasture",#15
                   "Sugar Cane",#20
                   "Dunes",#23
                   "Urban area",#24
                   "Water",#33
                   "Soy bean",#39
                   "Other crops")#41

levels(cerrado.lu.25km.2019) <- rat
levels(cerrado.lu.25km.2018) <- rat
levels(cerrado.lu.25km.2017) <- rat[-c(10),]
levels(cerrado.lu.25km.2016) <- rat[-c(10),]
levels(cerrado.lu.25km.2015) <- rat
levels(cerrado.lu.25km.2014) <- rat
levels(cerrado.lu.25km.2013) <- rat
levels(cerrado.lu.25km.2012) <- rat
levels(cerrado.lu.25km.2011) <- rat
levels(cerrado.lu.25km.2010) <- rat
levels(cerrado.lu.25km.2009) <- rat
levels(cerrado.lu.25km.2008) <- rat
levels(cerrado.lu.25km.2007) <- rat
levels(cerrado.lu.25km.2006) <- rat
levels(cerrado.lu.25km.2005) <- rat
levels(cerrado.lu.25km.2004) <- rat[-c(9),]
levels(cerrado.lu.25km.2003) <- rat[-c(9),]
levels(cerrado.lu.25km.2002) <- rat[-c(9),]
levels(cerrado.lu.25km.2001) <- rat[-c(9),]
levels(cerrado.lu.25km.2000) <- rat[-c(9),]
levels(cerrado.lu.25km.1999) <- rat[-c(9,11),]
levels(cerrado.lu.25km.1998) <- rat[-c(9,11),]
levels(cerrado.lu.25km.1997) <- rat[-c(9,11),]
levels(cerrado.lu.25km.1996) <- rat[-c(9,11),]
levels(cerrado.lu.25km.1995) <- rat[-c(9,11),]
levels(cerrado.lu.25km.1994) <- rat[-c(9,11),]
levels(cerrado.lu.25km.1993) <- rat[-c(9,11),]
levels(cerrado.lu.25km.1992) <- rat[-c(7,9,11),]
levels(cerrado.lu.25km.1991) <- rat[-c(7,9,11),]
levels(cerrado.lu.25km.1990) <- rat[-c(7,9,11),]
levels(cerrado.lu.25km.1989) <- rat[-c(7,9,11),]
levels(cerrado.lu.25km.1988) <- rat[-c(4,7,9,11),]
levels(cerrado.lu.25km.1987) <- rat[-c(7,9,11),]
levels(cerrado.lu.25km.1986) <- rat[-c(7,9,11),]
levels(cerrado.lu.25km.1985) <- rat[-c(7,9,11),]

#Years 2005 to 2015, 2018, and 2019
quartz(8,8)
levelplot(cerrado.lu.25km.2014, main="2014",
          col.regions=c("NA",#Unclassified(0)
                                              "#006400",#Forest(3)
                                              "#32CD32",#Savanna(4)
                                              "#935132",#Forest plantation(9)
                                              "#B8AF4F",#Grassland(12)
                                              "#FFD966",#Pasture(15)
                                              "#C27BA0",#Sugar Cane(20)
                                              "#DD7E6B",#Dunes(23)
                                              "#aa0000",#Urban area(24)
                                              "#0000FF",#Water(33)
                                              "#c59ff4",#Soy bean(39)
                                              "#e787f8"))#Other crops(41)

#Years 2016 and 2017
quartz(8,8)
levelplot(cerrado.lu.25km.2017, main="2017",
          col.regions=c("NA",#Unclassified(0)
                                              "#006400",#Forest(3)
                                              "#32CD32",#Savanna(4)
                                              "#935132",#Forest plantation(9)
                                              "#B8AF4F",#Grassland(12)
                                              "#FFD966",#Pasture(15)
                                              "#C27BA0",#Sugar Cane(20)
                                              "#DD7E6B",#Dunes(23)
                                              "#aa0000",#Urban area(24)
                                              #"#0000FF",#Water(33)
                                              "#c59ff4",#Soy bean(39)
                                              "#e787f8"))#Other crops(41)

#Years 2000 to 2004
quartz(8,8)
levelplot(cerrado.lu.25km.2004, main="2004",
          col.regions=c("NA",#Unclassified(0)
                                              "#006400",#Forest(3)
                                              "#32CD32",#Savanna(4)
                                              "#935132",#Forest plantation(9)
                                              "#B8AF4F",#Grassland(12)
                                              "#FFD966",#Pasture(15)
                                              "#C27BA0",#Sugar Cane(20)
                                              "#DD7E6B",#Dunes(23)
                                              #"#aa0000",#Urban area(24)
                                              "#0000FF",#Water(33)
                                              "#c59ff4",#Soy bean(39)
                                              "#e787f8"))#Other crops(41)

#Years 1993 to 1999
quartz(8,8)
levelplot(cerrado.lu.25km.1999, main="1999",
          col.regions=c("NA",#Unclassified(0)
                                              "#006400",#Forest(3)
                                              "#32CD32",#Savanna(4)
                                              "#935132",#Forest plantation(9)
                                              "#B8AF4F",#Grassland(12)
                                              "#FFD966",#Pasture(15)
                                              "#C27BA0",#Sugar Cane(20)
                                              "#DD7E6B",#Dunes(23)
                                              #"#aa0000",#Urban area(24)
                                              "#0000FF",#Water(33)
                                              #"#c59ff4",#Soy bean(39)
                                              "#e787f8"))#Other crops(41)

#Years 1985 to 1992
quartz(8,8)
levelplot(cerrado.lu.25km.1992, main="1992",
          col.regions=c("NA",#Unclassified(0)
                                              "#006400",#Forest(3)
                                              "#32CD32",#Savanna(4)
                                              "#935132",#Forest plantation(9)
                                              "#B8AF4F",#Grassland(12)
                                              "#FFD966",#Pasture(15)
                                              #"#C27BA0",#Sugar Cane(20)
                                              "#DD7E6B",#Dunes(23)
                                              #"#aa0000",#Urban area(24)
                                              "#0000FF",#Water(33)
                                              #"#c59ff4",#Soy bean(39)
                                              "#e787f8"))#Other crops(41)

#From 5 to 5 years
quartz(8,8)
levelplot(cerrado.lu.25km.1985, main = "1985",
          col.regions=c("NA",#Unclassified(0)
                        "#006400",#Forest(3)
                        "#32CD32",#Savanna(4)
                        "#935132",#Forest plantation(9)
                        "#B8AF4F",#Grassland(12)
                        "#FFD966",#Pasture(15)
                        #"#C27BA0",#Sugar Cane(20)
                        "#DD7E6B",#Dunes(23)
                        #"#aa0000",#Urban area(24)
                        "#0000FF",#Water(33)
                        #"#c59ff4",#Soy bean(39)
                        "#e787f8"))#Other crops(41)

levelplot(cerrado.lu.25km.1990, main = "1990",
          col.regions=c("NA",#Unclassified(0)
                        "#006400",#Forest(3)
                        "#32CD32",#Savanna(4)
                        "#935132",#Forest plantation(9)
                        "#B8AF4F",#Grassland(12)
                        "#FFD966",#Pasture(15)
                        #"#C27BA0",#Sugar Cane(20)
                        "#DD7E6B",#Dunes(23)
                        #"#aa0000",#Urban area(24)
                        "#0000FF",#Water(33)
                        #"#c59ff4",#Soy bean(39)
                        "#e787f8"))#Other crops(41)

levelplot(cerrado.lu.25km.1995, main = "1995",
          col.regions=c("NA",#Unclassified(0)
                        "#006400",#Forest(3)
                        "#32CD32",#Savanna(4)
                        "#935132",#Forest plantation(9)
                        "#B8AF4F",#Grassland(12)
                        "#FFD966",#Pasture(15)
                        "#C27BA0",#Sugar Cane(20)
                        "#DD7E6B",#Dunes(23)
                        #"#aa0000",#Urban area(24)
                        "#0000FF",#Water(33)
                        #"#c59ff4",#Soy bean(39)
                        "#e787f8"))#Other crops(41)

levelplot(cerrado.lu.25km.2000, main = "2000",
          col.regions=c("NA",#Unclassified(0)
                        "#006400",#Forest(3)
                        "#32CD32",#Savanna(4)
                        "#935132",#Forest plantation(9)
                        "#B8AF4F",#Grassland(12)
                        "#FFD966",#Pasture(15)
                        "#C27BA0",#Sugar Cane(20)
                        "#DD7E6B",#Dunes(23)
                        #"#aa0000",#Urban area(24)
                        "#0000FF",#Water(33)
                        "#c59ff4",#Soy bean(39)
                        "#e787f8"))#Other crops(41)

levelplot(cerrado.lu.25km.2005, main = "2005",
          col.regions=c("NA",#Unclassified(0)
                        "#006400",#Forest(3)
                        "#32CD32",#Savanna(4)
                        "#935132",#Forest plantation(9)
                        "#B8AF4F",#Grassland(12)
                        "#FFD966",#Pasture(15)
                        "#C27BA0",#Sugar Cane(20)
                        "#DD7E6B",#Dunes(23)
                        "#aa0000",#Urban area(24)
                        "#0000FF",#Water(33)
                        "#c59ff4",#Soy bean(39)
                        "#e787f8"))#Other crops(41)

levelplot(cerrado.lu.25km.2010, main = "2010",
          col.regions=c("NA",#Unclassified(0)
                        "#006400",#Forest(3)
                        "#32CD32",#Savanna(4)
                        "#935132",#Forest plantation(9)
                        "#B8AF4F",#Grassland(12)
                        "#FFD966",#Pasture(15)
                        "#C27BA0",#Sugar Cane(20)
                        "#DD7E6B",#Dunes(23)
                        "#aa0000",#Urban area(24)
                        "#0000FF",#Water(33)
                        "#c59ff4",#Soy bean(39)
                        "#e787f8"))#Other crops(41)

levelplot(cerrado.lu.25km.2015, main = "2015",
          col.regions=c("NA",#Unclassified(0)
                        "#006400",#Forest(3)
                        "#32CD32",#Savanna(4)
                        "#935132",#Forest plantation(9)
                        "#B8AF4F",#Grassland(12)
                        "#FFD966",#Pasture(15)
                        "#C27BA0",#Sugar Cane(20)
                        "#DD7E6B",#Dunes(23)
                        "#aa0000",#Urban area(24)
                        "#0000FF",#Water(33)
                        "#c59ff4",#Soy bean(39)
                        "#e787f8"))#Other crops(41)

levelplot(cerrado.lu.25km.2019, main = "2019",
          col.regions=c("NA",#Unclassified(0)
                        "#006400",#Forest(3)
                        "#32CD32",#Savanna(4)
                        "#935132",#Forest plantation(9)
                        "#B8AF4F",#Grassland(12)
                        "#FFD966",#Pasture(15)
                        "#C27BA0",#Sugar Cane(20)
                        "#DD7E6B",#Dunes(23)
                        "#aa0000",#Urban area(24)
                        "#0000FF",#Water(33)
                        "#c59ff4",#Soy bean(39)
                        "#e787f8"))#Other crops(41)



###########
#Modelling#
###########
# 
# cerrado.climate.df<-readRDS("Data/Predictors/cerrado.climate.df.rds")
# fixed.covs.df<-readRDS("fixed.covs.df.rds")
# cerrado.fire.df<-readRDS("cerrado.fire.df.rds")
# 
# cerrado.climate.df<-cerrado.climate.df[,-c(7,8,9,12,18)]#exclude temp_lv2 water_lv2 skin_temp temp_lv1 water.deficit
# cerrado.fire.covs.df<-left_join(cerrado.fire.df,cerrado.climate.df,by=c("x","y","Date"))
# cerrado.fire.covs.df<-left_join(cerrado.fire.covs.df,cerrado.lu.dflong,by=c("x","y","Year.x"))
# cerrado.fire.covs.df<-left_join(cerrado.fire.covs.df,fixed.covs.df[,-7],by=c("x","y"))
# cerrado.fire.covs.df$ID<-interaction(cerrado.fire.covs.df$x,cerrado.fire.covs.df$y)
# cerrado.fire.covs.df<-cerrado.fire.covs.df[,c(27,1:3,5:7,4,8,9,12:18,20:26)]
# 
# saveRDS(cerrado.fire.covs.df,"cerrado.fire.covs.df.rds")
library(mgcv)
cerrado.fire.covs.df<-readRDS("Data/Predictors/cerrado.fire.covs.df.rds")
head(cerrado.fire.covs.df)
tail(cerrado.fire.covs.df)

###########################
#Air temperature - temp_2m#
###########################
cerrado.fire.covs.df$temp_2m<-cerrado.fire.covs.df$temp_2m-273.15

Freq_fire_temp2m <- group_by(cerrado.fire.covs.df, temp_2m) %>%
  summarise(meanFreqFire = mean(Freq_fire))


windows(10,10)
quartz(8,8)
ggplot(Freq_fire_temp2m, aes(temp_2m,meanFreqFire)) +
  geom_point(size=.8, colour="orange", alpha=0.5) +
  stat_smooth(method=bam, formula=y~s(x,bs="cr"), level=0.95, colour="gray30", size=1.25,
              method.args = list(family = "binomial",discrete=T)) +
  labs(x="Air temperature (°C)", y="Fire probability")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(colour="black"))


Exten_fire_temp2m <- group_by(cerrado.fire.covs.df, temp_2m) %>%
  summarise(meanExtenFire = mean(Exten_Fire,na.rm=T))

windows(10,10)
quartz(8,8)
ggplot(Exten_fire_temp2m, aes(temp_2m,meanExtenFire)) +
  geom_point(size=.8, colour="orange", alpha=0.5) +
  stat_smooth(method=bam, formula=y~s(x,bs="cr"), level=0.95, colour="gray30", size=1.25,
              method.args = list(family = "poisson",discrete=T)) +
  labs(x="Air temperature (°C)", y="Fire Extent (m²)")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(colour="black"))

########################
#Precipitation - precip#
########################

Freq_fire_precip <- group_by(cerrado.fire.covs.df, precip) %>%
  summarise(meanFreqFire = mean(Freq_fire))


windows(10,10)
quartz(8,8)

ggplot(Freq_fire_precip, aes(precip,meanFreqFire)) +
  geom_point(size=.8, colour="orange", alpha=0.5) +
  stat_smooth(method=bam, formula=y~s(x,bs="cr"), level=0.95, colour="gray30", size=1.25,
              method.args = list(family = "binomial",discrete=T)) +
  labs(x="Precipitation (m)", y="Fire probability")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(colour="black"))


Exten_fire_precip <- group_by(cerrado.fire.covs.df, precip) %>%
  summarise(meanExtenFire = mean(Exten_Fire,na.rm=T))

windows(10,10)
quartz(8,8)

ggplot(Exten_fire_precip, aes(precip,meanExtenFire)) +
  geom_point(size=.8, colour="orange", alpha=0.5) +
  stat_smooth(method=bam, formula=y~s(x,bs="cr"), level=0.95, colour="gray30", size=1.25,
              method.args = list(family = "poisson",discrete=T)) +
  labs(x="Precipitation (m)", y="Fire Extent (m²)")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(colour="black"))


##################
#Insolation - sol#
##################

Freq_fire_sol <- group_by(cerrado.fire.covs.df, sol) %>%
  summarise(meanFreqFire = mean(Freq_fire))


windows(10,10)
quartz(8,8)

ggplot(Freq_fire_sol, aes(sol,meanFreqFire)) +
  geom_point(size=.8, colour="orange", alpha=0.5) +
  stat_smooth(method=bam, formula=y~s(x,bs="cr"), level=0.95, colour="gray30", size=1.25,
              method.args = list(family = "binomial",discrete=T)) +
  labs(x="Solar radiation (J/m²)", y="Fire probability")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(colour="black"))


Exten_fire_sol <- group_by(cerrado.fire.covs.df, sol) %>%
  summarise(meanExtenFire = mean(Exten_Fire,na.rm=T))

windows(10,10)
quartz(8,8)

ggplot(Exten_fire_sol, aes(sol,meanExtenFire)) +
  geom_point(size=.8, colour="orange", alpha=0.5) +
  stat_smooth(method=bam, formula=y~s(x,bs="cr"), level=0.95, colour="gray30", size=1.25,
              method.args = list(family = "poisson",discrete=T)) +
  labs(x="Solar radiation (J/m²)", y="Fire Extent (m²)")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(colour="black"))

#Space-time object
# xy.Cerrado<-cerrado.fire.covs.df[,2:3]
# coordinates(xy.Cerrado) <- ~ x + y
# proj4string(xy.Cerrado) <- crs(Bbin.all.Cerrado.LTDR)
# xy.Cerrado
# 
# st.fire.df<-cerrado.fire.covs.df[,c(1,4,8)]
# st.fire.df$Date.x<-as.POSIXct(st.fire.df$Date.x)
# 
# FireProb.st<-stConstruct(st.fire.df,
#                          "ID","Date.x",
#                          xy.Cerrado)



##########
#Land use#
##########
table(cerrado.fire.covs.df$land_use,cerrado.fire.covs.df$Year.x)
cerrado.fire.covs.df$land_use[cerrado.fire.covs.df$land_use==5]<-3#Transform Mangrove to forest

table(cerrado.fire.covs.df$land_use,cerrado.fire.covs.df$Year.x)

Probfire_landuse <-group_by(cerrado.fire.covs.df, land_use, Year.x) %>%
  summarise(Prob_Fire_mean = mean(Freq_fire),Prob_Fire_sd = sd(Freq_fire))

Extenfire_landuse <- group_by(cerrado.fire.covs.df, land_use, Year.x) %>%
  summarise(Exten_Fire_mean = mean(Exten_Fire),Exten_Fire_sd = sd(Exten_Fire))

windows(10,10)
quartz(8,8)
par(mar=c(8,4,1,1))
boxplot(Exten_Fire_mean~land_use,data=Extenfire_landuse,
        ylab="Fire extent (m²)",xlab="",las=2,bty="n",
        names = c("Unclassified",
                  "Forest",
                  "Savanna",
                  "Forest plantation",
                  "Grassland",
                  "Pasture",
                  "Sugar Cane",
                  "Urban area",
                  "Water",
                  "Soy bean",
                  "Other crops")
)

windows(10,10)
quartz(8,8)
par(mar=c(8,4,1,1))
boxplot(Prob_Fire_mean~land_use,data=Probfire_landuse,
        ylab="Fire occurrence",xlab="",las=2,bty="n",
        names = c("Unclassified",
                  "Forest",
                  "Savanna",
                  "Forest plantation",
                  "Grassland",
                  "Pasture",
                  "Sugar Cane",
                  "Urban area",
                  "Water",
                  "Soy bean",
                  "Other crops")
)


############
#Vegetation#
############
Probfire_veg <-group_by(cerrado.fire.covs.df, VegCerrado.r, Year.x) %>%
  summarise(Prob_Fire_mean = mean(Freq_fire),Prob_Fire_sd = sd(Freq_fire))

Extenfire_veg <- group_by(cerrado.fire.covs.df, VegCerrado.r, Year.x) %>%
  summarise(Exten_Fire_mean = mean(Exten_Fire),Exten_Fire_sd = sd(Exten_Fire))

table(cerrado.fire.covs.df$VegCerrado.r)

windows(10,10)
quartz(8,8)
par(mar=c(8,4,1,1))
boxplot(Exten_Fire_mean~VegCerrado.r,data=Extenfire_veg,
        ylab="Fire extent (m²)",xlab="",las=2,bty="n"#,
        # names = c("Unclassified",
        #           "Forest",
        #           "Savanna",
        #           "Forest plantation",
        #           "Grassland",
        #           "Pasture",
        #           "Sugar Cane",
        #           "Urban area",
        #           "Water",
        #           "Soy bean",
        #           "Other crops")
)

windows(10,10)
par(mar=c(8,4,1,1))
boxplot(Prob_Fire_mean~VegCerrado.r,data=Probfire_veg,
        ylab="Fire occurrence",xlab="",las=2,bty="n"#,
        # names = c("Unclassified",
        #           "Forest",
        #           "Savanna",
        #           "Forest plantation",
        #           "Grassland",
        #           "Pasture",
        #           "Sugar Cane",
        #           "Urban area",
        #           "Water",
        #           "Soy bean",
        #           "Other crops")
)


###############
#Random Forest#
###############
library(randomForest)
library(Boruta)
library(RRF)
#library(bigrf)
options(na.action=na.omit)

completos <- complete.cases(cerrado.fire.covs.df)
new.cerrado.fire.covs.df <- droplevels(cerrado.fire.covs.df[completos,])

##################
#Fire probability#
##################

Borutafireprob <- Boruta(as.factor(Freq_fire) ~., data = new.cerrado.fire.covs.df[,-c(1:7,9)], doTrace = 2)
Borutafireprob
saveRDS(Borutafireprob, "Borutafireprob.rds")
plot(Borutafireprob)
quartz(8,12)
plot(Borutafireprob, xlab = "", xaxt = "n",bty="n")
lz<-lapply(1:ncol(Borutafireprob$ImpHistory),function(i)
  Borutafireprob$ImpHistory[is.finite(Borutafireprob$ImpHistory[,i]),i])
names(lz) <- colnames(Borutafireprob$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(Borutafireprob$ImpHistory), cex.axis = 0.7)

quartz(8,8)
plotImpHistory(Borutafireprob,bty="n",ylim=c(0,150))

attStats(Borutafireprob)

# Ordinary RF
rf.fireprob <- randomForest(as.factor(Freq_fire) ~ .,
                             proximity = F,
                             importance = T,
                             keep.forest = F,
                             ntree = 1000,
                            replace = F,
                            sampsize = 0.05*nrow(cerrado.fire.covs.df),
                             data = new.cerrado.fire.covs.df[,-c(1:7,9)],
                             flagReg = 0,
                            na.action=na.exclude)
rf.fireprob
rf.fireprob$importance

# Guided and regularized RF
rrf.fireprob <- RRF(new.cerrado.fire.covs.df[,-c(1:9)],
                    as.factor(new.cerrado.fire.covs.df[, 8]),
                     flagReg = 0,
                    replace = F,
                    sampsize = 0.05*nrow(new.cerrado.fire.covs.df))

(impRF.fireprob <- rrf.fireprob$importance)
impRF.fireprob <- impRF.fireprob[, "MeanDecreaseGini"]
imp.fireprob <- impRF.fireprob/(max(impRF.fireprob)) # normalize the importance score
gamma <- 0.5
coefReg.fireprob <- (1 - gamma) + gamma * imp.fireprob # weighted average

grrf.fireprob <- RRF(new.cerrado.fire.covs.df[,-c(1:9)],
                     as.factor(new.cerrado.fire.covs.df[, 8]),
                      coefReg = coefReg.fireprob,
                      flagReg = 1,
                      importance = T,
                      ntree = 1000,
                     replace = F,
                     sampsize = 0.05*nrow(new.cerrado.fire.covs.df))
RRF::importance(grrf.fireprob)
grrf.fireprob$confusion # matriz de confusao
str(grrf.fireprob)
tail(grrf.fireprob$err.rate ) # taxa de erro da predicao (OOB)
grrf.fireprob$feaSet # variaveis selecionadas
names(new.cerrado.fire.covs.df[,-c(1:9)])[grrf.fireprob$feaSet] # nomes das variaveis selecionadas

#windows(t="Classification Error - SC", h = 8, w = 12)
quartz(8,12)
plot(grrf.fireprob, bty = "n", main = NULL)

#windows(t="GRRF Variable Importance - LF", h = 8, w = 12)
quartz(8,12)
varImpPlot(grrf.fireprob, main = NULL)


# 10-fold cross-validation
set.seed(29062018)
step <-  1 - (1 / 10)
rrfcv.fireprob <- rrfcv(new.cerrado.fire.covs.df[,-c(1:9)],
                        as.factor(new.cerrado.fire.covs.df[, 8]),
                         cv.fold = 10,
                         step = step)
#windows(8, 8, t = "rrfcv.LF")
quartz("rrfcv.LF", 8, 8)
with(rrfcv.fireprob, plot(n.var, error.cv,
                           type = "o",
                           cex = 1.5,
                           xaxp = c(1, 12, 11),
                           xlab = "Number of predictors",
                           ylab = "CV Mean Error Rate",
                           las = 1,
                           bty = "n"))


# 100 replicates of 10-fold cross-validation
trials <- 100
rrfcv.100.fireprob <- foreach(icount(trials), .packages = "RRF") %dopar% {
  rrfcv(cerrado.fire.covs.df[,-c(1:9)],
        as.factor(cerrado.fire.covs.df[, 8]),
        cv.fold = 10,
        step	 = step)
}
(error.cv.fireprob <- sapply(rrfcv.100.fireprob, "[[", "error.cv"))

#windows (w = 8, h = 8, t = "Matplot SF")
quartz(w = 12, h = 8, "Matplot LF")
matplot(rrfcv.100.fireprob[[1]]$n.var,
        cbind(error.cv.fireprob, rowMeans(error.cv.fireprob)),
        type = "l",
        lwd = c(rep(0.5, ncol(error.cv.fireprob)), 5),
        col = c(rep(rgb(69, 117, 180, 100, max = 255), ncol(error.cv.fireprob)), rgb(215, 48, 39, max = 255)),
        lty = 1,
        bty = "n",
        las = 1,
        xaxp = c(1, 3, 2),
        xlim=c(1, 3),
        ylim = c(0.9, 1),
        xlab = "Number of Predictors",
        ylab = "CV Mean Error Rate"
)  

rowMeans(error.cv.fireprob)

#############
#Fire extent#
#############
Extent_Fire_Cerrado.df <- new.cerrado.fire.covs.df
Extent_Fire_Cerrado.df$Exten_Fire <- log10(Extent_Fire_Cerrado.df$Exten_Fire)
Extent_Fire_Cerrado.df <- Extent_Fire_Cerrado.df[is.finite(Extent_Fire_Cerrado.df$Exten_Fire),]
summary(Extent_Fire_Cerrado.df$Exten_Fire)
boxplot(Extent_Fire_Cerrado.df$Exten_Fire)
Extent_Fire_Cerrado.df <- Extent_Fire_Cerrado.df[Extent_Fire_Cerrado.df$Exten_Fire>2,]


Borutafireext <- Boruta(Exten_Fire ~., data = Extent_Fire_Cerrado.df[,-c(1:8)], doTrace = 2)
Borutafireext
saveRDS(Borutafireext, "Borutafireext.rds")
plot(Borutafireext)
quartz(8,12)
plot(Borutafireext, xlab = "", xaxt = "n",bty="n")
lz<-lapply(1:ncol(Borutafireext$ImpHistory),function(i)
  Borutafireext$ImpHistory[is.finite(Borutafireext$ImpHistory[,i]),i])
names(lz) <- colnames(Borutafireext$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(Borutafireext$ImpHistory), cex.axis = 0.7)

quartz(8,8)
plotImpHistory(Borutafireext,bty="n",ylim=c(0,150))

attStats(Borutafireext)

# Ordinary RF
rf.fireExten <- randomForest(Exten_Fire ~ .,
                            proximity = F,
                            importance = T,
                            keep.forest = F,
                            ntree = 1000,
                            replace = F,
                            sampsize = 0.05*nrow(new.cerrado.fire.covs.df),
                            data = new.cerrado.fire.covs.df[,-c(1:8)],
                            flagReg = 0,
                            na.action=na.exclude)
rf.fireExten
rf.fireExten$importance

# Guided and regularized RF
rrf.fireExten <- RRF(new.cerrado.fire.covs.df[,-c(1:9)],
                    new.cerrado.fire.covs.df[, 9],
                    flagReg = 0,
                    replace = F,
                    sampsize = 0.05*nrow(new.cerrado.fire.covs.df))

(impRF.fireExten <- rrf.fireExten$importance)
impRF.fireExten <- impRF.fireExten[, "MeanDecreaseGini"]
imp.fireExten <- impRF.fireExten/(max(impRF.fireExten)) # normalize the importance score
gamma <- 0.5
coefReg.fireExten <- (1 - gamma) + gamma * imp.fireExten # weighted average

grrf.fireExten <- RRF(new.cerrado.fire.covs.df[,-c(1:9)],
                     new.cerrado.fire.covs.df[, 9],
                     coefReg = coefReg.fireExten,
                     flagReg = 1,
                     importance = T,
                     ntree = 1000,
                     replace = F,
                     sampsize = 0.05*nrow(new.cerrado.fire.covs.df))
RRF::importance(grrf.fireExten)
grrf.fireExten$confusion # matriz de confusao
str(grrf.fireExten)
tail(grrf.fireExten$err.rate ) # taxa de erro da predicao (OOB)
grrf.fireExten$feaSet # variaveis selecionadas
names(new.cerrado.fire.covs.df[,-c(1:9)])[grrf.fireExten$feaSet] # nomes das variaveis selecionadas

quartz(t="Classification Error - SC", h = 8, w = 12)
quartz(8,12)
plot(grrf.fireExten, bty = "n", main = NULL)

quartz( h = 8, w = 12)
varImpPlot(grrf.fireExten, main = NULL)


# 10-fold cross-validation
set.seed(29062018)
step <-  1 - (1 / 10)
rrfcv.fireExten <- rrfcv(new.cerrado.fire.covs.df[,-c(1:8)],
                        as.factor(new.cerrado.fire.covs.df[, 9]),
                        cv.fold = 10,
                        step = step)
quartz(8, 8, t = "rrfcv.LF")
#quartz("rrfcv.LF", 8, 8)
with(rrfcv.fireExten, plot(n.var, error.cv,
                          type = "o",
                          cex = 1.5,
                          xaxp = c(1, 12, 11),
                          xlab = "Number of predictors",
                          ylab = "CV Mean Error Rate",
                          las = 1,
                          bty = "n"))


# 100 replicates of 10-fold cross-validation
trials <- 100
rrfcv.100.fireExten <- foreach(icount(trials), .packages = "RRF") %dopar% {
  rrfcv(cerrado.fire.covs.df[,-c(1:9)],
        as.factor(cerrado.fire.covs.df[, 8]),
        cv.fold = 10,
        step	 = step)
}
(error.cv.fireExten <- sapply(rrfcv.100.fireExten, "[[", "error.cv"))

quartz(w = 8, h = 8, t = "Matplot SF")
#quartz(w = 12, h = 8, "Matplot LF")
matplot(rrfcv.100.fireExten[[1]]$n.var,
        cbind(error.cv.fireExten, rowMeans(error.cv.fireExten)),
        type = "l",
        lwd = c(rep(0.5, ncol(error.cv.fireExten)), 5),
        col = c(rep(rgb(69, 117, 180, 100, max = 255), ncol(error.cv.fireExten)), rgb(215, 48, 39, max = 255)),
        lty = 1,
        bty = "n",
        las = 1,
        xaxp = c(1, 3, 2),
        xlim=c(1, 3),
        ylim = c(0.9, 1),
        xlab = "Number of Predictors",
        ylab = "CV Mean Error Rate"
)  

rowMeans(error.cv.fireExten)

#####
#BAM#
#####

library(mgcv)

##################
#Fire probability#
##################
new.cerrado.fire.covs.df$land_use[new.cerrado.fire.covs.df$land_use==5]<-3#Mangrove to forest

bam0.FireProb <- bam(Freq_fire~te(x,y,t,bs="cr"),
            data = new.cerrado.fire.covs.df,family="binomial",discrete = T)
summary(bam0.FireProb)
quartz(8,8)
plot(bam0.FireProb)

bam.full.FireProb <- bam(Freq_fire~te(x,y,t,bs="cr")+
                           s(RH,bs="cr")+
                           s(total.evap,bs="cr")+
                           s(pot.evap,bs="cr")+
                  s(precip,bs="cr")+
                  s(temp_2m,bs="cr")+
                  s(sol,bs="cr")+
                  s(wind,bs="cr")+
                  s(water_lv1,bs="cr")+
                  s(decliv.Cerrado,bs="cr"),
                data = new.cerrado.fire.covs.df,family="binomial",discrete=T)

summary(bam.full.FireProb)

bam.full.FireProb2 <- bam(Freq_fire~te(x,y,t,bs="cr")+
                            s(RH,bs="cr")+
                            s(total.evap,bs="cr")+
                            s(pot.evap,bs="cr")+
                          s(precip,bs="cr")+
                            s(temp_2m,bs="cr")+
                            s(sol,bs="cr")+
                            s(wind,bs="cr")+
                            s(water_lv1,bs="cr"),
                          data = new.cerrado.fire.covs.df,family="binomial",discrete=T)
summary(bam.full.FireProb2)

bam.full.FireProb3 <- bam(Freq_fire~te(x,y,t,bs="cr",d=c(2,1))+
                            s(RH,bs="cr")+
                            s(total.evap,bs="cr")+
                            s(pot.evap,bs="cr")+
                          s(precip,bs="cr")+
                            s(temp_2m,bs="cr")+
                            s(sol,bs="cr")+
                            s(wind,bs="cr")+
                            s(water_lv1,bs="cr")+
                            s(decliv.Cerrado,bs="cr"),
                          data = new.cerrado.fire.covs.df,family="binomial",discrete=T,na.action=na.fail)
summary(bam.full.FireProb3)

bam.full.FireProb4 <- bam(Freq_fire~te(x,y,t,bs="cr",d=c(2,1))+
                            s(RH,bs="cr")+
                            s(total.evap,bs="cr")+
                            s(pot.evap,bs="cr")+
                            s(precip,bs="cr")+
                            s(temp_2m,bs="cr")+
                            s(sol,bs="cr")+
                            s(wind,bs="cr")+
                            s(water_lv1,bs="cr")+
                            s(decliv.Cerrado,bs="cr")+
                            as.factor(land_use)+
                            as.factor(VegCerrado.r),
                          data = new.cerrado.fire.covs.df,family="binomial",discrete=T,na.action=na.fail)
summary(bam.full.FireProb4)


bam.full.FireProb5 <- bam(Freq_fire~te(x,y,t,bs="cr",d=c(2,1))+
                            s(RH,bs="cr")+
                            s(total.evap,bs="cr")+
                            s(pot.evap,bs="cr")+
                            s(precip,bs="cr")+
                            s(temp_2m,bs="cr")+
                            s(sol,bs="cr")+
                            s(wind,bs="cr")+
                            s(water_lv1,bs="cr")+
                            s(decliv.Cerrado,bs="cr")+
                            s(DEM.Cerrado)+
                            s(biomass.Cerrado)+
                            s(canopy.height.Cerrado)+
                            s(CerradoRoads.dist)+
                            as.factor(land_use)+
                            as.factor(VegCerrado.r),
                          data = new.cerrado.fire.covs.df,family="binomial",discrete=T,na.action=na.fail)

summary(bam.full.FireProb5)


bam.trend.FireProb <- bam(Freq_fire~te(x,y,bs="cr")+t,
                          data = new.cerrado.fire.covs.df,family="binomial",discrete=T)

summary(bam.trend.FireProb)
plot(bam.trend.FireProb)
plot(bam.trend.FireProb,all.terms=T)


AIC(bam.full.FireProb,bam.full.FireProb2,bam.full.FireProb3,bam.full.FireProb4)

plot(bam.full.FireProb3,all.terms=T,pages=1)
plot(bam.full.FireProb5,all.terms=T,pages=1)


#windows(10,10)
quartz(8,8)
plot(bam.trend.FireProb,all.terms=T,pages=1)




# 
# saveRDS(bam0.FireProb,"bam0.FireProb.rds")
# saveRDS(bam.full.FireProb,"bam.full.FireProb.rds")
# saveRDS(bam.full.FireProb2,"bam.full.FireProb2.rds")
# saveRDS(bam.full.FireProb3,"bam.full.FireProb3.rds")
# 
# bam0.FireProb<-readRDS("bam0.FireProb.rds")
# bam.full.FireProb<- readRDS("bam.full.FireProb.rds")
# bam.full.FireProb2<- readRDS("bam.full.FireProb2.rds")
# bam.full.FireProb3<- readRDS("bam.full.FireProb3.rds")

anova(bam0.FireProb,bam.full.FireProb,test="Chisq")
anova(bam.full.FireProb,bam.full.FireProb2,test="Chisq")
anova(bam.full.FireProb2,bam.full.FireProb3,test="Chisq")
anova(bam.full.FireProb3,bam.full.FireProb4,test="Chisq")


#############
#Final plots#
#############
summary(lm_RH <- lm(RH~Year.x,data=new.cerrado.fire.covs.df))
summary(lm_temp <- lm(temp_2m~Year.x,data=new.cerrado.fire.covs.df))
summary(lm_potevap<- lm(pot.evap~Year.x,data=new.cerrado.fire.covs.df))
summary(lm_totalevap<- lm(total.evap~Year.x,data=new.cerrado.fire.covs.df))
summary(lm_precip<- lm(precip*10000~Year.x,data=new.cerrado.fire.covs.df))
summary(lm_wind<- lm(wind~Year.x,data=new.cerrado.fire.covs.df))
summary(lm_sol<- lm(sol~Year.x,data=new.cerrado.fire.covs.df))
summary(lm_water<- lm(water_lv1~Year.x,data=new.cerrado.fire.covs.df))

summary(lm_fire<- lm(Freq_fire~Year.x,data=new.cerrado.fire.covs.df))
summary(lm_fire2<- lm(Exten_Fire~Year.x,data=Extent_Fire_Cerrado.df))


quartz(8,12)
ggplot(new.cerrado.fire.covs.df, 
       aes(x = Year.x, y = RH, group = Year.x)) + 
  ggdist::stat_halfeye(
    adjust = 5,
    width = .8,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .1,
    outlier.shape = NA
  ) +
  geom_abline(slope=lm_RH$coefficients[2],intercept = lm_RH$coefficients[1], color="red")+
  coord_cartesian( clip = "off")+
  labs(x="", y="Relative humidity (%)")

quartz(8,12)
ggplot(new.cerrado.fire.covs.df, 
       aes(x = Year.x, y = temp_2m, group = Year.x)) + 
  ggdist::stat_halfeye(
    adjust = 5,
    width = .8,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .1,
    outlier.shape = NA
  ) +
  geom_abline(slope=lm_temp$coefficients[2],intercept = lm_temp$coefficients[1], color="red")+
  coord_cartesian( clip = "off")+
  labs(x="", y="Air temperature (°C)")

quartz(8,12)
ggplot(new.cerrado.fire.covs.df, 
       aes(x = Year.x, y = pot.evap, group = Year.x)) + 
  ggdist::stat_halfeye(
    adjust = 5,
    width = .8,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .1,
    outlier.shape = NA
  ) +
  geom_abline(slope=lm_potevap$coefficients[2],intercept = lm_potevap$coefficients[1], color="red")+
  coord_cartesian( clip = "off",ylim = c(-.02,0))+
  labs(x="", y="Potential evaporation (m)")

quartz(8,12)
ggplot(new.cerrado.fire.covs.df, 
       aes(x = Year.x, y = total.evap, group = Year.x)) + 
  ggdist::stat_halfeye(
    adjust = 5,
    width = .8,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .1,
    outlier.shape = NA
  ) +
  geom_abline(slope=lm_totalevap$coefficients[2],intercept = lm_totalevap$coefficients[1], color="red")+
  coord_cartesian( clip = "off")+
  labs(x="", y="Evaporation (m)")

quartz(8,12)
ggplot(new.cerrado.fire.covs.df, 
       aes(x = Year.x, y = precip*10000, group = Year.x)) + 
  ggdist::stat_halfeye(
    adjust = 5,
    width = .8,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .1,
    outlier.shape = NA
  ) +
  geom_abline(slope=lm_precip$coefficients[2],intercept = lm_precip$coefficients[1], color="red")+
  coord_cartesian( clip = "off",ylim=c(0,300))+
  labs(x="", y="Precipitation (mm)")

quartz(8,12)
ggplot(new.cerrado.fire.covs.df, 
       aes(x = Year.x, y = wind, group = Year.x)) + 
  ggdist::stat_halfeye(
    adjust = 5,
    width = .8,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .1,
    outlier.shape = NA
  ) +
  geom_abline(slope=lm_wind$coefficients[2],intercept = lm_wind$coefficients[1], color="red")+
  coord_cartesian( clip = "off",ylim=c(0,3))+
  labs(x="", y="Wind speed (m/s)")

quartz(8,12)
ggplot(new.cerrado.fire.covs.df, 
       aes(x = Year.x, y = sol, group = Year.x)) + 
  ggdist::stat_halfeye(
    adjust = 5,
    width = .8,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .1,
    outlier.shape = NA
  ) +
  geom_abline(slope=lm_sol$coefficients[2],intercept = lm_sol$coefficients[1], color="red")+
  coord_cartesian( clip = "off")+
  labs(x="", y="Insolation (J/m2)")

quartz(8,12)
ggplot(new.cerrado.fire.covs.df, 
       aes(x = Year.x, y = water_lv1, group = Year.x)) + 
  ggdist::stat_halfeye(
    adjust = 5,
    width = .8,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .1,
    outlier.shape = NA
  ) +
  geom_abline(slope=lm_water$coefficients[2],intercept = lm_water$coefficients[1], color="red")+
  coord_cartesian( clip = "off")+
  labs(x="", y="Volumetric soil water (m3/m3)")

quartz(8,12)
ggplot(new.cerrado.fire.covs.df, 
       aes(x = Year.x, y = Freq_fire, group = Year.x)) + 
  ggdist::stat_halfeye(
    adjust = 5,
    width = .8,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .1,
    outlier.shape = NA
  ) +
  geom_abline(slope=lm_fire$coefficients[2],intercept = lm_fire$coefficients[1], color="red")+
  coord_cartesian( clip = "off",ylim=c(0,1))+
  labs(x="", y="Fire probability")

quartz(8,12)
ggplot(Extent_Fire_Cerrado.df, 
       aes(x = Year.x, y = Exten_Fire, group = Year.x)) + 
  ggdist::stat_halfeye(
    adjust = 5,
    width = .8,
    .width = 0,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .1,
    outlier.shape = NA
  ) +
  geom_abline(slope=lm_fire2$coefficients[2],intercept = lm_fire2$coefficients[1], color="red")+
  coord_cartesian( clip = "off",ylim=c(1,9))+
  labs(x="", y="log10(Fire extent +1)")


summary(bam.full.FireProb3)

library(mgcViz)

bam.trend.FireProb.viz <- getViz(bam.trend.FireProb)

print(plot(bam.trend.FireProb.viz, allTerms = T), pages = 1) # Calls print.plotGam()

quartz(12,10)
pl <- plot(bam.trend.FireProb.viz, allTerms = T) + l_dens2D("cond") + l_fitContour() +  l_points() + l_fitLine(linetype = 3) +
  l_ciLine(colour = 2) + l_ciBar() + l_fitPoints(size = 1, col = 2) + 
  theme_get() + labs(title = NULL)
print(pl,pages=1)

quartz(8,8)
plot(sm(bam.trend.FireProb.viz, 1)) + l_fitRaster() + l_fitContour() + 
  labs(title = NULL,x = "Longitude (m)",y="Latitude (m)",z="Abundance")

#Full model
bam.full.FireProb.viz <- getViz(bam.full.FireProb3)

print(plot(bam.full.FireProb.viz, allTerms = T), pages = 1) # Calls print.plotGam()

quartz(12,10)
pl <- plot(bam.full.FireProb.viz, allTerms = T) + l_dens2D("joint") + l_fitContour() +  l_points() + l_fitLine(linetype = 3) +
  l_ciLine(colour = 2) + l_ciBar() + l_fitPoints(size = 1, col = 2) + 
  theme_get() + labs(title = NULL)
print(pl,pages=1)

quartz(8,8)
plot(bam.full.FireProb.viz, select = 2) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireProb.viz,x="RH"))

plot(bam.full.FireProb.viz, select = 3) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireProb.viz,x="total.evap"))

plot(bam.full.FireProb.viz, select = 4) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireProb.viz,x="pot.evap"))

plot(bam.full.FireProb.viz, select = 5) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireProb.viz,x="precip"))

plot(bam.full.FireProb.viz, select = 6) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireProb.viz,x="temp_2m"))

plot(bam.full.FireProb.viz, select = 7) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireProb.viz,x="sol"))

plot(bam.full.FireProb.viz, select = 8) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireProb.viz,x="wind"))

plot(bam.full.FireProb.viz, select = 9) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireProb.viz,x="water_lv1"))

plot(bam.full.FireProb.viz, select = 10) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireProb.viz,x="decliv.Cerrado"))


summary(new.cerrado.fire.covs.df[new.cerrado.fire.covs.df$Month.x=="fev",7])
#January
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(37, 421, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#February
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(38, 422, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#March
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(39, 423, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#April
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(40, 424, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#May
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(41, 425, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#June
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(42, 426, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#July
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(43, 427, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#August
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(44, 428, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#September
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(45, 429, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#October
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(46, 430, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#November
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(47, 431, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#December
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireProb.viz, 1), 
                fix = list("t" = round(seq(48, 432, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#############
#Fire Extent#
#############
hist(log10(new.cerrado.fire.covs.df$Exten_Fire+1))
summary(log10(new.cerrado.fire.covs.df$Exten_Fire+1))
bam0.FireExtent <- bam(log(Exten_Fire+1)~te(x,y,bs="cr")+s(t,bs="cr"),
          data = new.cerrado.fire.covs.df,discrete=T)
summary(bam0.FireExtent)

Extent_Fire_Cerrado.df <- new.cerrado.fire.covs.df
Extent_Fire_Cerrado.df$Exten_Fire <- log10(Extent_Fire_Cerrado.df$Exten_Fire)
Extent_Fire_Cerrado.df <- Extent_Fire_Cerrado.df[is.finite(Extent_Fire_Cerrado.df$Exten_Fire),]
summary(Extent_Fire_Cerrado.df$Exten_Fire)
boxplot(Extent_Fire_Cerrado.df$Exten_Fire)
Extent_Fire_Cerrado.df <- Extent_Fire_Cerrado.df[Extent_Fire_Cerrado.df$Exten_Fire>2,]


bam.full.FireExtent <- bam(Exten_Fire~te(x,y,bs="cr")+s(t,bs="cr")+
                s(sol,bs="cr")+
                s(wind,bs="cr")+
                s(DEM.Cerrado,bs="cr")+
                s(decliv.Cerrado,bs="cr")+
                s(canopy.height.Cerrado,bs="cr")+
                s(biomass.Cerrado,bs="cr")+
                s(CerradoRoads.dist,bs="cr")+
                  as.factor(land_use)+
                  as.factor(VegCerrado.r),
              data = Extent_Fire_Cerrado.df,discrete=T)


bam.full.FireExtent2 <- bam(Exten_Fire~te(x,y,bs="cr")+s(t,bs="cr")+
                              s(sol,bs="cr")+
                              s(wind,bs="cr")+
                              s(DEM.Cerrado,bs="cr")+
                              s(decliv.Cerrado,bs="cr")+
                              s(canopy.height.Cerrado,bs="cr")+
                              s(biomass.Cerrado,bs="cr")+
                              s(CerradoRoads.dist,bs="cr")+
                              as.factor(land_use),
                           data = Extent_Fire_Cerrado.df,discrete=T)


bam.full.FireExtent3 <- bam(Exten_Fire~te(y,x,t,bs="cr",d=c(2,1))+
                              s(sol,bs="cr")+
                              s(wind,bs="cr")+
                              s(DEM.Cerrado,bs="cr")+
                              s(decliv.Cerrado,bs="cr")+
                              s(canopy.height.Cerrado,bs="cr")+
                              s(biomass.Cerrado,bs="cr")+
                              s(CerradoRoads.dist,bs="cr")+
                              as.factor(land_use)+
                              as.factor(VegCerrado.r),
                           data = Extent_Fire_Cerrado.df,discrete=T)
Extent_Fire_Cerrado.df$land_use <- as.factor(Extent_Fire_Cerrado.df$land_use)
Extent_Fire_Cerrado.df$VegCerrado.r <- as.factor(Extent_Fire_Cerrado.df$VegCerrado.r)

bam.full.FireExtent4 <- bam(Exten_Fire~te(x,y,t,bs="cr",d=c(2,1))+
                              s(RH,bs="cr")+
                              s(total.evap,bs="cr")+
                              s(pot.evap,bs="cr")+
                              s(precip,bs="cr")+
                              s(temp_2m,bs="cr")+
                              s(sol,bs="cr")+
                              s(wind,bs="cr")+
                              s(water_lv1,bs="cr")+
                              s(DEM.Cerrado,bs="cr")+
                              s(decliv.Cerrado,bs="cr")+
                              s(canopy.height.Cerrado,bs="cr")+
                              s(biomass.Cerrado,bs="cr")+
                              s(CerradoRoads.dist,bs="cr")+
                              land_use+
                              VegCerrado.r,
                            data = Extent_Fire_Cerrado.df,discrete=T)

bam.full.FireExtent5 <- bam(Exten_Fire~te(x,y,t,bs="cr",d=c(2,1))+
                              s(sol,bs="cr")+
                              s(wind,bs="cr")+
                              s(DEM.Cerrado,bs="cr")+
                              s(decliv.Cerrado,bs="cr")+
                              s(canopy.height.Cerrado,bs="cr")+
                              s(biomass.Cerrado,bs="cr")+
                              s(CerradoRoads.dist,bs="cr")+
                              land_use+
                              VegCerrado.r,
                            data = Extent_Fire_Cerrado.df,discrete=T)
summary(bam.full.FireExtent4)
summary(bam.full.FireExtent5)

bam.null.FireExtent5 <- bam(Exten_Fire~te(x,y,t,bs="cr",d=c(2,1)),
                            data = Extent_Fire_Cerrado.df,discrete=T)
summary(bam.null.FireExtent5)

bam.trend.FireExten <- bam(Exten_Fire~te(x,y,bs="cr")+t,
                          data = Extent_Fire_Cerrado.df,discrete=T)

summary(bam.trend.FireExten)
#windows(10,10)
quartz(8,8)
plot(bam.trend.FireExten,all.terms=T,pages=1)

# saveRDS(bam0.FireExtent,"bam0.FireExtent.rds")
# saveRDS(bam.full.FireExtent,"bam.full.FireExtent.rds")
# saveRDS(bam.full.FireExtent2,"bam.full.FireExtent2.rds")
# saveRDS(bam.full.FireExtent3,"bam.full.FireExtent3.rds")
# 
# bam0.FireExtent<-readRDS("bam0.FireExtent.rds")
# bam.full.FireExtent<- readRDS("bam.full.FireExtent.rds")
# bam.full.FireExtent2<- readRDS("bam.full.FireExtent2.rds")
# bam.full.FireExtent3<- readRDS("bam.full.FireExtent3.rds")

AIC(bam.full.FireExtent,bam.full.FireExtent2,bam.full.FireExtent3,bam.full.FireExtent4,bam.full.FireExtent5)

summary(bam.full.FireExtent4)
summary(bam.full.FireExtent5)

plot(bam.full.FireExtent4, pages=1)
quartz(8,12)
plot(bam.full.FireExtent5, pages=1)


anova(bam0.FireExtent,bam.full.FireExtent,test="Chisq")
anova(bam.full.FireExtent,bam.full.FireExtent2,test="Chisq")
anova(bam.full.FireExtent2,bam.full.FireExtent3,test="Chisq")
anova(bam.full.FireExtent3,bam.full.FireExtent4,test="Chisq")
anova(bam.full.FireExtent4,bam.full.FireExtent5,test="Chisq")


#############
#Final plots#
#############
bam.trend.FireExtent.viz <- getViz(bam.trend.FireExten)

quartz(8,12)
print(plot(bam.trend.FireExtent.viz, allTerms = T), pages = 1) # Calls print.plotGam()

quartz(12,10)
pl <- plot(bam.trend.FireExtent.viz, allTerms = T) + l_dens2D("cond") + l_fitContour() +  l_points() + l_fitLine(linetype = 3) +
  l_ciLine(colour = 2) + l_ciBar() + l_fitPoints(size = 1, col = 2) + 
  theme_get() + labs(title = NULL)
print(pl,pages=1)

quartz(8,8)
plot(sm(bam.trend.FireExtent.viz, 1)) + l_fitRaster() + l_fitContour() + 
  labs(title = NULL,x = "Longitude (m)",y="Latitude (m)",z="Abundance")

#Full model
bam.full.FireExtent.viz <- getViz(bam.full.FireExtent5)

quartz(8,8)
print(plot(bam.full.FireExtent.viz, allTerms = T), pages = 1) # Calls print.plotGam()

quartz(12,10)
pl <- plot(bam.full.FireExtent.viz, allTerms = T) + l_dens2D("cond") + l_fitContour() + l_fitLine(linetype = 1) +
  l_ciLine(linetype=2) + l_ciBar() + l_fitPoints(size = 1, col = 2)  + labs(title = NULL)
print(pl,pages=1)

quartz(8,8)
plot(bam.full.FireExtent.viz, select = 2) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireExtent.viz,x="sol"))

plot(bam.full.FireExtent.viz, select = 3) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireExtent.viz,x="wind"))

plot(bam.full.FireExtent.viz, select = 4) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireExtent.viz,x="DEM.Cerrado"))

plot(bam.full.FireExtent.viz, select = 5) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireExtent.viz,x="decliv.Cerrado"))

plot(bam.full.FireExtent.viz, select = 6) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireExtent.viz,x="canopy.height.Cerrado"))

plot(bam.full.FireExtent.viz, select = 7) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireExtent.viz,x="biomass.Cerrado"))

plot(bam.full.FireExtent.viz, select = 8) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
plot(ALE(bam.full.FireExtent.viz,x="CerradoRoads.dist"))

# plot(bam.full.FireExtent.viz, select = 9) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# plot(bam.full.FireExtent.viz, select = 10) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# plot(bam.full.FireExtent.viz, select = 11) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# plot(bam.full.FireExtent.viz, select = 12) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# plot(bam.full.FireExtent.viz, select = 13) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
# plot(bam.full.FireExtent.viz, select = 14) + l_dens(type = "cond") + l_fitLine() + l_ciLine()
quartz(8,12)
plot(bam.full.FireExtent.viz,allTerms = T, select=9)

quartz(8,12)
  ggplot(Extent_Fire_Cerrado.df, 
         aes(x = as.factor(land_use), y = Exten_Fire)) + 
  ggdist::stat_halfeye(adjust = 5,
                       width = .5,
                       .width = 0,
                       justification = -.2,
                       point_colour = NA)

quartz(8,12)
plot(bam.full.FireExtent.viz, select=10)

quartz(8,12)
ggplot(Extent_Fire_Cerrado.df, 
       aes(x = as.factor(VegCerrado.r), y = Exten_Fire)) + 
  ggdist::stat_halfeye(adjust = 5,
                       width = 1,
                       .width = 0,
                       justification = -.1,
                       point_colour = NA)

#January
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(37, 421, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#February
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(38, 422, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#March
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(39, 423, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#April
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(40, 424, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#May
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(41, 425, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#June
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(42, 426, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#July
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(43, 427, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#August
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(44, 428, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#September
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(45, 429, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#October
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(46, 430, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#November
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(47, 431, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()

#December
quartz(12,12)
pl <- plotSlice(x = sm(bam.full.FireExtent.viz, 1), 
                fix = list("t" = round(seq(48, 432, length.out = 33),0)))
pl + l_fitRaster() + l_fitContour()  + l_rug()


# saveRDS(bam.full.FireProb3,"bam.full.FireProb3.rds")
# saveRDS(bam.full.FireExtent5,"bam.full.FireExtent5.rds")


#############################################
## Plot predictions and overlay observations#
#############################################
library(spTimer)
bam.full.FireProb3<- readRDS("Output/bam.full.FireProb3.rds")
bam.full.FireExtent5<- readRDS("Output/bam.full.FireExtent5.rds")

Fire.month.mean <- group_by(new.cerrado.fire.covs.df, x, y, Month.x) %>%
  summarise_all(mean)
Fire.month.mean$land_use<-new.cerrado.fire.covs.df$land_use[new.cerrado.fire.covs.df$Year.x==2018]
#Fire.month.mean$land_use[Fire.month.mean$land_use==5]<-3#Mangrove to Forest
Fire.month.mean<-droplevels.data.frame(Fire.month.mean)
summary(Fire.month.mean)

FireExtent.month.mean <- group_by(Extent_Fire_Cerrado.df, x, y, Month.x) %>%
  summarise_all(mean)

FireExtent.month.mean$land_use <- group_by(Extent_Fire_Cerrado.df, x, y, Month.x) %>%
  summarise(land_use = modal(land_use))

FireExtent.month.mean$VegCerrado.r <- group_by(Extent_Fire_Cerrado.df, x, y, Month.x) %>%
  summarise(VegCerrado.r = modal(VegCerrado.r))


FireExtent.month.mean$land_use <- FireExtent.month.mean$land_use$land_use
FireExtent.month.mean$VegCerrado.r <- FireExtent.month.mean$VegCerrado.r$VegCerrado.r

#Fire.month.mean$land_use[Fire.month.mean$land_use==5]<-3#Mangrove to Forest

FireExtent.month.mean<-droplevels.data.frame(FireExtent.month.mean)

pred.bam.mean.FireProb<-predict(bam.full.FireProb3,
                                Fire.month.mean,
                                type="response",se.fit = T)

pred.bam.mean.FireExtent<-predict(bam.full.FireExtent5,
                                  FireExtent.month.mean,
                                  type="response",se.fit=T)

spT.validation(Fire.month.mean$Freq_fire, pred.bam.mean.FireProb$fit)
spT.validation(FireExtent.month.mean$Exten_Fire, pred.bam.mean.FireExtent$fit)

pred.bam.mean.FireProb.df<-data.frame(ProbFire = pred.bam.mean.FireProb$fit,
                                      ProbFire.se = pred.bam.mean.FireProb$se.fit,
                                      x = Fire.month.mean$x,
                                      y = Fire.month.mean$y,
                                      t = Fire.month.mean$Month.x)

pred.bam.mean.FireExtent.df<-data.frame(ExtentFire = pred.bam.mean.FireExtent$fit,
                                        ExtentFire.se = pred.bam.mean.FireExtent$se.fit,
                                        x = FireExtent.month.mean$x,
                                        y = FireExtent.month.mean$y,
                                        t = FireExtent.month.mean$Month.x)

# saveRDS(pred.bam.mean.FireProb.df,"pred.bam.mean.FireProb.df.rds")
# saveRDS(pred.bam.mean.FireExtent.df,"pred.bam.mean.FireExtent.df.rds")

library(ggplot2)
library(STRbook)
library(tidyr)
library(dplyr)
pred.bam.mean.FireProb.df<-readRDS("Output/pred.bam.mean.FireProb.df.rds")
pred.bam.mean.FireExtent.df<-readRDS("Output/pred.bam.mean.FireExtent.df.rds")

#Fire occurrence
g1.pred.bam.FireProb <- ggplot() +
  geom_raster(data = pred.bam.mean.FireProb.df,
              aes(x, y, fill = pmin(pmax(ProbFire, 0), 1))) +
  facet_wrap(~t, nrow = 3, ncol = 4) +
  fill_scale(limits = c(0, 1),
             name = expression(Y[t])) +
  col_scale(name = "Fire occurrence", limits=c(0, 1)) +
  theme_bw()

windows(10,10)
quartz(8,8)
g1.pred.bam.FireProb

summary(pred.bam.mean.FireProb.df$ProbFire.se)
## Plot prediction standard errors
g2.pred.bam.FireProb <- ggplot() +
  geom_raster(data = pred.bam.mean.FireProb.df,
              aes(x, y, fill = pmin(ProbFire.se, 0.05))) +
  facet_wrap(~t, nrow = 3, ncol = 4) +
  fill_scale(palette = "BrBG",
             limits = c(0, 0.05),
             name = expression(s.e.)) +
  theme_bw()

windows(10,10)
quartz(8,8)
g2.pred.bam.FireProb

summary(pred.bam.mean.FireExtent.df$ExtentFire)

#Fire Extent
g1.pred.bam.FireExtent <- ggplot() +
  geom_raster(data = pred.bam.mean.FireExtent.df,
              aes(x, y, fill = pmin(pmax(ExtentFire, 4), 8))) +
  facet_wrap(~t, nrow = 3, ncol = 4) +
  fill_scale(limits = c(4, 8),
             name = expression(log(Y[t]))) +
  col_scale(name = "Fire extent (m2)") +
  theme_bw()

#windows(10,10)
quartz(8,8)
g1.pred.bam.FireExtent



#windows(10,10)
quartz(8,8)
boxplot(pred.bam.mean.FireExtent.df$ExtentFire.se)
## Plot prediction standard errors
g2.pred.bam.FireExtent <- ggplot() +
  geom_raster(data = pred.bam.mean.FireExtent.df,
              aes(x, y, fill = pmin(ExtentFire.se, 0.15))) +
  facet_wrap(~t, nrow = 3, ncol = 4) +
  fill_scale(palette = "BrBG",
             limits = c(0,.15),
             name = expression(s.e.)) +
  theme_bw()

#windows(10,10)
quartz(8,8)
g2.pred.bam.FireExtent


##################
#Cluster analysis#
##################
library(raster)
library(tidyverse)
library(corrplot)
library(usdm)
library(ggfortify)
setwd("Data/AVHRR_LTDR_Fire/Pixel")

Freq.fire.LTDR.Cerrado <- stack("Freq.fire.LTDR.Cerrado.grd")
FireExtent.LTDR.Cerrado <- stack("FireExtent.Cerrado.LTDR.monthly.grd")

Freq.sd.fire.LTDR.Cerrado <- stack("Freq.sd.fire.LTDR.Cerrado.grd")
FireExtent.sd.LTDR.Cerrado <- stack("FireExtent.sd.Cerrado.LTDR.monthly.grd")

setwd("/Users/heito/Documents/CerradoFireRegimes")

Freq.fire.LTDR.Cerrado.df<-as.data.frame(Freq.fire.LTDR.Cerrado,na.rm=T,xy=T)
FireExtent.LTDR.Cerrado.df<-as.data.frame(FireExtent.LTDR.Cerrado,na.rm=T,xy=T)
Freq.sd.fire.LTDR.Cerrado.df<-as.data.frame(Freq.sd.fire.LTDR.Cerrado,na.rm=T,xy=T)
FireExtent.sd.LTDR.Cerrado.df<-as.data.frame(FireExtent.sd.LTDR.Cerrado,na.rm=T,xy=T)

Freq.fire.LTDR.Cerrado.df$ID <- interaction(Freq.fire.LTDR.Cerrado.df$x,
                                            Freq.fire.LTDR.Cerrado.df$y)

FireExtent.LTDR.Cerrado.df$ID <- interaction(FireExtent.LTDR.Cerrado.df$x,
                                             FireExtent.LTDR.Cerrado.df$y)

Freq.sd.fire.LTDR.Cerrado.df$ID <- interaction(Freq.sd.fire.LTDR.Cerrado.df$x,
                                               Freq.sd.fire.LTDR.Cerrado.df$y)

FireExtent.sd.LTDR.Cerrado.df$ID <- interaction(FireExtent.sd.LTDR.Cerrado.df$x,
                                                FireExtent.sd.LTDR.Cerrado.df$y)

Fire.LTDR.Cerrado.df <- full_join(Freq.fire.LTDR.Cerrado.df,
                                  FireExtent.LTDR.Cerrado.df,
                                  by = c("ID"))
  
Fire.LTDR.Cerrado.df <- full_join(Fire.LTDR.Cerrado.df,
                                  Freq.sd.fire.LTDR.Cerrado.df,
                                  by = c("ID"))

Fire.LTDR.Cerrado.df <- full_join(Fire.LTDR.Cerrado.df,
                                  FireExtent.sd.LTDR.Cerrado.df,
                                  by = c("ID"))


Fire.LTDR.Cerrado.df <- na.omit(Fire.LTDR.Cerrado.df)

Fire.LTDR.Cerrado.df <- Fire.LTDR.Cerrado.df[,-c(16,17,30,31,44,45)]

names(Fire.LTDR.Cerrado.df) <- c("x","y",
                                 "mean.freq.jan",
                                 "mean.freq.feb",
                                 "mean.freq.mar",
                                 "mean.freq.apr",
                                 "mean.freq.may",
                                 "mean.freq.jun",
                                 "mean.freq.jul",
                                 "mean.freq.aug",
                                 "mean.freq.sep",
                                 "mean.freq.oct",
                                 "mean.freq.nov",
                                 "mean.freq.dec",
                                 "ID",
                                 "mean.ba.jan",
                                 "mean.ba.feb",
                                 "mean.ba.mar",
                                 "mean.ba.apr",
                                 "mean.ba.may",
                                 "mean.ba.jun",
                                 "mean.ba.jul",
                                 "mean.ba.aug",
                                 "mean.ba.sep",
                                 "mean.ba.oct",
                                 "mean.ba.nov",
                                 "mean.ba.dec",
                                 "sd.freq.jan",
                                 "sd.freq.feb",
                                 "sd.freq.mar",
                                 "sd.freq.apr",
                                 "sd.freq.may",
                                 "sd.freq.jun",
                                 "sd.freq.jul",
                                 "sd.freq.aug",
                                 "sd.freq.sep",
                                 "sd.freq.oct",
                                 "sd.freq.nov",
                                 "sd.freq.dec",
                                 "sd.ba.jan",
                                 "sd.ba.feb",
                                 "sd.ba.mar",
                                 "sd.ba.apr",
                                 "sd.ba.may",
                                 "sd.ba.jun",
                                 "sd.ba.jul",
                                 "sd.ba.aug",
                                 "sd.ba.sep",
                                 "sd.ba.oct",
                                 "sd.ba.nov",
                                 "sd.ba.dec"
                                 )

cor.fire <- cor(Fire.LTDR.Cerrado.df[,-c(1,2,15)])
quartz(8,12)
corrplot(cor.fire, method = "number", type="upper", number.cex = 0.3)

vifstep(Fire.LTDR.Cerrado.df[,-c(1,2,15)])
vifcor(Fire.LTDR.Cerrado.df[,-c(1,2,15)])

Fire.LTDR.Cerrado.df <- Fire.LTDR.Cerrado.df[,c("x","y",
  "mean.freq.jan",
  "mean.freq.feb",
  "mean.freq.mar",
  "mean.freq.apr",
  "mean.freq.may",
  "mean.freq.jun",
  "mean.freq.jul",
  "mean.freq.aug",
  "mean.freq.sep",
  "mean.freq.oct",
  "mean.freq.nov",
  "mean.freq.dec",
  "ID",
  "mean.ba.jan",
  "mean.ba.feb",
  "mean.ba.mar",
  "mean.ba.apr",
  "mean.ba.may",
  "mean.ba.jun",
  "mean.ba.jul",
  "mean.ba.aug",
  "mean.ba.sep",
  "mean.ba.oct",
  "mean.ba.nov",
  "mean.ba.dec",
  "sd.ba.jan",
  "sd.ba.feb",
  "sd.ba.mar",
  "sd.ba.apr",
  "sd.ba.may",
  "sd.ba.jun",
  "sd.ba.jul",
  "sd.ba.aug",
  "sd.ba.sep",
  "sd.ba.oct",
  "sd.ba.nov",
  "sd.ba.dec"
)]

cor.fire <- cor(Fire.LTDR.Cerrado.df[,-c(1,2,15)])
quartz(8,12)
corrplot(cor.fire, method = "number", type="upper", number.cex = 0.3)

Freq.fire.LTDR.Cerrado
FireExtent.LTDR.Cerrado

FireProbExten.LTDR.Cerrado<-stack(Freq.fire.LTDR.Cerrado,
                                  FireExtent.LTDR.Cerrado,
                                  FireExtent.sd.LTDR.Cerrado)
names(FireProbExten.LTDR.Cerrado)<-c("FireProb.1","FireProb.2","FireProb.3",
                                     "FireProb.4","FireProb.5","FireProb.6",
                                     "FireProb.7","FireProb.8","FireProb.9",
                                     "FireProb.10","FireProb.11","FireProb.12",
                                     "FireExten.1","FireExten.2","FireExten.3",
                                     "FireExten.4","FireExten.5","FireExten.6",
                                     "FireExten.7","FireExten.8","FireExten.9",
                                     "FireExten.10","FireExten.11","FireExten.12",
                                     "sdExten.1","sdExten.2","sdExten.3",
                                     "sdExten.4","sdExten.5","sdExten.6",
                                     "sdExten.7","sdExten.8","sdExten.9",
                                     "sdExten.10","sdExten.11","sdExten.12")
FireProbExten.LTDR.Cerrado


#Using monthly data and including coefficient of variation of fire extent#
##########################################################################
ValsFireMonth_mat <- values(FireProbExten.LTDR.Cerrado)
rownames(ValsFireMonth_mat) <- 1:dim(ValsFireMonth_mat)[1]

scale1d <- function(x,df,na.rm=TRUE){
  (x-mean(df,na.rm=TRUE))/sd(df,na.rm=TRUE)
}

ValsFireMonth_mat <- na.omit(ValsFireMonth_mat)
dim(ValsFireMonth_mat)

ValsFireMonth_mat[,c(1:12)] <- scale1d(ValsFireMonth_mat[,c(1:12)],
                                       ValsFireMonth_mat[,c(1:12)])

ValsFireMonth_mat[,c(13:24)] <- scale1d(ValsFireMonth_mat[,c(13:24)],
                                       ValsFireMonth_mat[,c(13:24)])

ValsFireMonth_mat[,c(25:36)] <- scale1d(ValsFireMonth_mat[,c(25:36)],
                                        ValsFireMonth_mat[,c(25:36)])

library(mclust)

set.seed(123)
dataBIC <- mclustBIC(ValsFireMonth_mat)# identify BICs for different models
print(summary(dataBIC)) # show summary of top-ranking models
#windows(10,10)
quartz(height = 8, width = 8)
plot(dataBIC)

set.seed(123)
mod <- Mclust(ValsFireMonth_mat, G=2)

mod$G # number of groups/clusters in model
summary(mod)
mod[["parameters"]][["mean"]] # mean values of clusters

windows(10,10)
quartz(8,8)
plot(mod, what = "uncertainty")

ModPred <- predict.Mclust(mod, ValsFireMonth_mat) # prediction
Pred_ras <- Freq.fire.LTDR.Cerrado[[1]] # establishing a rediction raster
values(Pred_ras) <- NA # set everything to NA
# set values of prediction raster to corresponding classification according to rowname
values(Pred_ras)[as.numeric(rownames(ValsFireMonth_mat))] <- as.vector(ModPred$classification)
Pred_ras

quartz(8,8)
plot(Pred_ras, # what to plot
     col = c("#fdae61","#abd9e9","#d7191c"), # colours for groups ,"#2c7bb6",
     #colNA = "gray", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2, # horizontal size of legend
     ylab = "Latitude", xlab = "Longitude"
)

library(terra)
Pred_ras_interp <- terra::focal(rast(Pred_ras),9,fun="modal",na.policy="only",na.rm=T)

#windows(10,10)
quartz(8,8)
plot(Pred_ras_interp, # what to plot
     type = "class",
     col = c("#fdae61","#abd9e9") # colours for groups,
     #colNA = "gray", # which colour to assign to NA values
     #legend.shrink=1, # vertical size of legend
     #legend.width=2, # horizontal size of legend
)
plot(Cerrado,add=T)
freq(Pred_ras_interp)

table(values(Pred_ras))

#Uncertainty#
#############
ModPred$z
plot(density(ModPred$z[,1]))
plot(density(ModPred$z[,2]))

#PCA#
#####
pca.fire <- prcomp(Fire.LTDR.Cerrado.df[,-c(1,2,15)], scale. = T)
summary(pca.fire)

quartz(8,8)
autoplot(pca.fire, 
         data = data.frame(classification = as.factor(ModPred$classification)), 
         colour = 'classification', alpha=0.5,
         loadings = TRUE,
         loadings.label = TRUE)+
  scale_color_manual(values = c("#fdae61","#abd9e9"))

moddr <- MclustDR(mod)
summary(moddr)
#plot(mod2dr, what = "pairs")
quartz(8,8)
#windows(10,10)
plot(moddr, what = "boundaries", ngrid = 200)

moddr <- MclustDR(mod,lambda=.5)
summary(moddr)
#plot(mod2dr, what = "pairs")
#windows(10,10)
quartz(8,8)
plot(moddr, what = "boundaries", ngrid = 200)

modda <- MclustDA(ValsFireMonth_mat,ModPred$classification)
summary(modda)

moddr.da <- MclustDR(modda)
summary(moddr.da)
plot(moddr.da, what = "boundaries", ngrid = 200)


plot(modda,"classification")
plot(modda,"error")


#plot(mod2dr, what = "pairs")
#windows(10,10)
quartz(8,8)
plot(moddr, what = "boundaries", ngrid = 200)

# writeRaster(Pred_ras,"FireRegimesCerrado_months.grd",overwrite=T)
# writeRaster(Pred_ras_interp,"FireRegimesCerrado_months_interp.grd",overwrite=T)



Freq.fire.LTDR.Cerrado
FireExtent.LTDR.Cerrado

FireProbExten.LTDR.Cerrado<-stack(Freq.fire.LTDR.Cerrado,
                                  FireExtent.LTDR.Cerrado)
names(FireProbExten.LTDR.Cerrado)<-c("FireProb.1","FireProb.2","FireProb.3",
                                     "FireProb.4","FireProb.5","FireProb.6",
                                     "FireProb.7","FireProb.8","FireProb.9",
                                     "FireProb.10","FireProb.11","FireProb.12",
                                     "FireExten.1","FireExten.2","FireExten.3",
                                     "FireExten.4","FireExten.5","FireExten.6",
                                     "FireExten.7","FireExten.8","FireExten.9",
                                     "FireExten.10","FireExten.11","FireExten.12")
FireProbExten.LTDR.Cerrado

#Using monthly data just with mean fire occurrence (frequency) and extent#
##########################################################################
ValsFireMonth_mat <- values(FireProbExten.LTDR.Cerrado)
rownames(ValsFireMonth_mat) <- 1:dim(ValsFireMonth_mat)[1]

scale1d <- function(x,df,na.rm=TRUE){
  (x-mean(df,na.rm=TRUE))/sd(df,na.rm=TRUE)
}

ValsFireMonth_mat <- na.omit(ValsFireMonth_mat)
dim(ValsFireMonth_mat)

ValsFireMonth_mat[,c(1:12)] <- scale1d(ValsFireMonth_mat[,c(1:12)],
                                       ValsFireMonth_mat[,c(1:12)])

ValsFireMonth_mat[,c(13:24)] <- scale1d(ValsFireMonth_mat[,c(13:24)],
                                       ValsFireMonth_mat[,c(13:24)])

library(mclust)

set.seed(123)
dataBIC <- mclustBIC(ValsFireMonth_mat)# identify BICs for different models
print(summary(dataBIC)) # show summary of top-ranking models
#windows(10,10)
quartz(8,8)
plot(dataBIC)

set.seed(123)
mod <- Mclust(ValsFireMonth_mat,G=2)

mod$G # number of groups/clusters in model
summary(mod)
mod[["parameters"]][["mean"]] # mean values of clusters

# windows(10,10)
# quartz(8,8)
plot(mod)

ModPred <- predict.Mclust(mod, ValsFireMonth_mat) # prediction
Pred_ras <- Freq.fire.LTDR.Cerrado[[1]] # establishing a rediction raster
values(Pred_ras) <- NA # set everything to NA
# set values of prediction raster to corresponding classification according to rowname
values(Pred_ras)[as.numeric(rownames(ValsFireMonth_mat))] <- as.vector(ModPred$classification)
Pred_ras

quartz(8,8)
plot(Pred_ras, # what to plot
     col = c("#2c7bb6","#fdae61"), # "#2c7bb6",colours for groups ,
     #colNA = "gray", # which colour to assign to NA values
     legend.shrink=1, # vertical size of legend
     legend.width=2, # horizontal size of legend
     ylab = "Latitude", xlab = "Longitude"
)

library(terra)
Pred_ras_interp <- terra::focal(rast(Pred_ras),9,fun="modal",na.policy="only",na.rm=T)

#windows(10,10)
quartz(8,8)
plot(Pred_ras_interp, # what to plot
     type = "class",
     col = c("#2c7bb6","#fdae61") # colours for groups"#2c7bb6",
     #colNA = "gray", # which colour to assign to NA values
     #legend.shrink=1, # vertical size of legend
     #legend.width=2, # horizontal size of legend
)
plot(Cerrado,add=T)
freq(Pred_ras_interp)

table(values(Pred_ras))


#Uncertainty#
#############
ModPred$z
quartz(8,8)
plot(density(ModPred$z[ModPred$classification==1,1] -1, from = 0, to = 0.5), 
     col=alpha("#2c7bb6",0.5), bty = "n", ylim = c(0,3.5), bty="n", lwd = 2,
     main = "", xlab = "Uncertainty probability", ylab = "Density")
lines(density(ModPred$z[ModPred$classification==2,2]-1, from = 0, to = 0.5), 
      col=alpha("#fdae61",0.5), lwd = 2)

#PCA#
#####
pca.fire <- prcomp(ValsFireMonth_mat, scale. = T)
summary(pca.fire)

quartz(8,8)
autoplot(pca.fire, 
         data = data.frame(classification = as.factor(ModPred$classification)), 
         colour = 'classification', 
         loadings = TRUE,
         loadings.label = TRUE, alpha=.5,
         frame = TRUE, frame.type = 'norm',)+
  scale_color_manual(values = c("#2c7bb6","#fdae61"))+
  scale_fill_manual(values = c("#2c7bb6","#fdae61"))

#Discriminant analysis
moddr <- MclustDR(mod)
summary(moddr)
#plot(mod2dr, what = "pairs")
quartz(8,8)
#windows(10,10)
plot(moddr, what = "boundaries")
plot(moddr, what = "density")
plot(moddr, what = "evalues")


moddr <- MclustDR(mod,lambda=.5)
summary(moddr)
#plot(mod2dr, what = "pairs")
#windows(10,10)
quartz(8,8)
plot(moddr, what = "boundaries", ngrid = 200)

modda <- MclustDA(ValsFireMonth_mat,ModPred$classification)
summary(modda)
#plot(modda)

moddr.da <- MclustDR(modda)
summary(moddr.da)
quartz(height = 8, width = 8)
plot(moddr.da, "scatterplot", dimens = c(1,2), bty = "n",
     colors = c(alpha("#abd9e9",0.5),alpha("#fdae61",0.5)),
     symbols = c(16,16))

plot(moddr.da, "scatterplot", dimens = c(1,5), 
     colors = c(alpha("#abd9e9",0.5),alpha("#fdae61",0.5)))

writeRaster(Pred_ras,"FireRegimesCerrado_months.grd",overwrite=T)
writeRaster(Pred_ras_interp,"FireRegimesCerrado_months_interp.grd",overwrite=T)


#Using MODIS data#
##################
setwd("Data/MODIS_Fire/MCD64A1")
Fire.MCD64A1.month<-stack("Freq.fire.soma.MCD64A1.Cerrado.grd")
setwd("/Users/heito/Documents/CerradoFireRegimes")


ValsFireMonth_mat3 <- values(Fire.MCD64A1.month)
rownames(ValsFireMonth_mat3) <- 1:dim(ValsFireMonth_mat3)[1]

ValsFireMonth_mat3 <- na.omit(ValsFireMonth_mat3)

cor.fire <- cor(ValsFireMonth_mat3)
quartz(8,12)
corrplot(cor.fire, method = "number", type="upper", number.cex = 1)

dim(ValsFireMonth_mat3)

ValsFireMonth_mat3[,c(1:12)] <- scale1d(ValsFireMonth_mat3[,c(1:12)],
                                        ValsFireMonth_mat3[,c(1:12)])

#Next lines take a long time to run
set.seed(123)
dataBIC3 <- mclustBIC(ValsFireMonth_mat3,G=1:5)# identify BICs for different models
print(summary(dataBIC3)) # show summary of top-ranking models
windows(10,10)
quartz(8,8)
plot(dataBIC3)

dataBIC4 <- mclustBIC(ValsFireMonth_mat3,x=dataBIC3,G=6:9,modelNames="EEV")# identify BICs for different models
dataBIC4<-mclustBICupdate(dataBIC3,dataBIC4)
print(summary(dataBIC4)) # show summary of top-ranking models
windows(10,10)
quartz(8,8)
plot(dataBIC4)

# dataBIC5 <- mclustBIC(ValsFireMonth_mat3,G=10:12,modelNames="EEV")# identify BICs for different models
# dataBIC5<-mclustBICupdate(dataBIC4,dataBIC5)
# print(summary(dataBIC5)) # show summary of top-ranking models
# windows(10,10)
# plot(dataBIC5)

mod3 <- Mclust(ValsFireMonth_mat3,# data for the cluster model
               modelNames="EEV",
               G = 5 # BIC index for model to be built
)

mod4 <- Mclust(ValsFireMonth_mat3,# data for the cluster model
               modelNames="EEV",
               G = 3 # BIC index for model to be built
)

mod5 <- Mclust(ValsFireMonth_mat3,# data for the cluster model
               modelNames="EEV",
               G = 2 # BIC index for model to be built
)

ModPred3 <- predict.Mclust(mod3, ValsFireMonth_mat3) # prediction

plot(density(ModPred3$z[,1]))
plot(density(ModPred3$z[,2]))
plot(density(ModPred3$z[,3]))
plot(density(ModPred3$z[,4]))
plot(density(ModPred3$z[,5]))

ModPred4 <- predict.Mclust(mod4, ValsFireMonth_mat3) # prediction

plot(density(ModPred4$z[,1]))
plot(density(ModPred4$z[,2]))
plot(density(ModPred4$z[,3]))

ModPred5 <- predict.Mclust(mod5, ValsFireMonth_mat3) # prediction

quartz(8,8)
plot(density(ModPred5$z[ModPred5$classification==1,1]-1, from = 0, to = 0.1), 
     col=alpha("#fdae61",.5),lwd=2,ylim=c(0,550), bty="n", xlim = c(0, 0.005),
     main="", xlab = "Uncertainty probability", ylab = "Density")
lines(density(ModPred5$z[ModPred5$classification==2,2]-1, from = 0, to = 0.1), 
      col=alpha("#2c7bb6",0.5), lwd = 2)

#Discriminant analysis
modda.modis <- MclustDA(ValsFireMonth_mat3,ModPred3$classification)
summary(modda.modis)
#plot(modda)

modda.modis2 <- MclustDA(ValsFireMonth_mat3,ModPred4$classification)
summary(modda.modis2)

modda.modis3 <- MclustDA(ValsFireMonth_mat3,ModPred5$classification)
summary(modda.modis3)

moddr.da.modis <- MclustDR(modda.modis3)
summary(moddr.da.modis)
quartz(height = 8, width = 8)
plot(moddr.da.modis, "scatterplot", dimens = c(1,2), bty = "n",
     colors = c(alpha("#fdae61",0.5),alpha("#2c7bb6",0.5)),
     symbols = c(16,16))


saveRDS(mod3,"mod3.mclust.rds")
saveRDS(mod5, "mod5.mclust.rds")
saveRDS(modda.modis3,"modda.modis3.rds")
saveRDS(moddr.da.modis,"moddr.da.modis.rds")

summary(mod3)

mod5 <- readRDS("mod5.mclust.rds")
mod5[["parameters"]][["mean"]] # mean values of clusters
summary(mod5)

saveRDS(ModPred3,"MclustModelPred3.rds")
saveRDS(ModPred5,"MclustModelPred5.rds")

Pred_ras5 <- Fire.MCD64A1.month[[1]] # establishing a rediction raster
values(Pred_ras5) <- NA # set everything to NA
# set values of prediction raster to corresponding classification according to rowname
values(Pred_ras5)[as.numeric(rownames(ValsFireMonth_mat3))] <- as.vector(ModPred5$classification)

Pred_ras5 <- rast(Pred_ras5)

#plot(Fire.MCD64A1.month)

# Pred_ras3 <- rast("FireRegimesCerrado_MODIS_5clusters.grd")

library(viridis)
#windows(10,10)
quartz(8,8)
plot(Pred_ras5, # what to plot
     type = "class",
     col = c("#fdae61","#2c7bb6") # colours for groups,"#ffffbf",,"#d7191c","#abd9e9"
     #colNA = "gray", # which colour to assign to NA values
     #legend.shrink=1, # vertical size of legend
     #legend.width=2, # horizontal size of legend
)

plot(Cerrado,add=T)
freq(Pred_ras5)
Pred_ras5<-as.raster(Pred_ras5)

table(values(Pred_ras5))

#PCA
pca.fire <- prcomp(ValsFireMonth_mat3, scale. = T)
summary(pca.fire)

quartz(8,8)
autoplot(pca.fire, 
         data = data.frame(classification = as.factor(ModPred5$classification)), 
         colour = 'classification', 
         loadings = TRUE,
         frame = TRUE, frame.type = 'norm',
         loadings.label = TRUE, alpha=.5) +
  scale_color_manual(values = c("#fdae61","#2c7bb6"))

# writeRaster(Pred_ras5,"FireRegimesCerrado_MODIS_2clusters.grd",overwrite=T)

