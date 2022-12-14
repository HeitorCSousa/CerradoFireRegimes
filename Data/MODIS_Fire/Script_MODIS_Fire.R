
# install.packages("raster",dep=T)
# install.packages("rgdal",dep=T)
# install.packages("gdalUtils",dep=T)

library(gdalUtils)

library(raster)
library(rgdal)

rasterOptions(tmpdir=file.path("D:/Documentos/Raster")) 
rasterOptions()

# Provides detailed data on hdf4 files but takes ages

gdalinfo("MCD64CMQ.A2000306.006.2018149224217.hdf")

# Tells me what subdatasets are within my hdf4 MODIS files and makes them into a list

sds <- get_subdatasets("MCD64CMQ.A2000306.006.2018149224217.hdf")
sds

# I'm only interested in the first subdataset and I can use gdal_translate to convert it to a .tif
sds.1<-readGDAL("MCD64CMQ.A2000306.006.2018149224217.hdf")
sds.1
gdal_translate(sds[0], dst_dataset = "BA_11_2000.tif")

# Load and plot the new .tif

rast <- raster("BA2000306.tif")
plot(rast)
zoom(rast)

# If you have lots of files then you can make a loop to do all this for you

files <- dir(pattern = ".hdf")
files

filename <- substr(files,11,17)
filename <- paste0("BA", filename, ".tif")
filename

i <- 1

for (i in 1:216){
  sds <- get_subdatasets(files[i])
  gdal_translate(sds[1], dst_dataset = filename[i])
}

BA.all<-stack(filename)
BA.all

#########################
#Reads Cerrado shapefile#
#########################
Cerrado<-readOGR("/Volumes/Extreme SSD/Heitor/BD_SIG/MODIS_FireCerrado")#Le shapefile
plot(Cerrado,axes=T)

#Crop rasters to Cerrado extent

res(BA.all)<-c(0.25,0.25)
BA.all <- setExtent(BA.all, extent(-180, 180, -90, 90))
crs(BA.all)<-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
Cerrado
BA.all

plot(BA.all[[1]])
plot(Cerrado,add=T)
zoom(BA.all[[1]])

BA.all.Cerrado<-crop(BA.all, extent(Cerrado))
BA.all.Cerrado <- mask(BA.all.Cerrado, Cerrado)

plot(BA.all.Cerrado[[1]])
plot(Cerrado,axes=T,add=T)

writeRaster(BA.all.Cerrado,"BA.all.Cerrado.tif")

#Reads raster (stack)##########################
BA.all.Cerrado<-stack("BA.all.Cerrado.tif")
############################################

BA.all.Cerrado
crs(BA.all.Cerrado)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
Cerrado
plot(BA.all.Cerrado)
plot(BA.all.Cerrado[[1]])
plot(Cerrado,axes=T,add=T)

###########################################################
#Calculates sum of burned area by month for Cerrado extent#
###########################################################
soma_BA<-cellStats(BA.all.Cerrado,sum)

#Time series analyses
ts.soma.BA<-ts(soma_BA,start=c(2000, 11), end=c(2018, 10), frequency=12)
ts.soma.BA<-ts.soma.BA*0.01 #Para transformar em ha multiplicar pelo fator 0.01 (Giglio et al., 2018 - Collection 6 MODIS Burned Area Product User?s Guide Version 1.2)

plot(ts.soma.BA,ylab="?rea (ha)",xlab="Tempo (anos)")

#Time series decomposition
plot(decompose(ts.soma.BA))
Queima.decom <- decompose(ts.soma.BA,type="multiplicative")
plot(Queima.decom)
Trend <- Queima.decom$trend
Seasonal <- Queima.decom$seasonal
ts.plot(cbind(Trend, Trend * Seasonal), lty = 1:2,ylab="?rea (ha)")

decomp<-stl(ts.soma.BA, "per")
summary(decomp)
plot(decomp)
monthplot(ts.soma.BA, ylab="?rea (ha)", type="h")
monthplot(decomp, choice="seasonal")
monthplot(decomp, choice="trend")
monthplot(decomp, choice="remainder")

# Percentage of explained variation by each component
p.sazonal<-(var(decomp$time.series[, 1]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.sazonal

p.tendencia<-(var(decomp$time.series[, 2]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.tendencia

p.residual<-(var(decomp$time.series[, 3]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.residual

p.sazonal+p.tendencia+p.residual

###############################################################
#Calculates sum of burned area by month for Conservation Units#
###############################################################
install.packages("psych",dep=T)
library(psych)
UCs.Cerrado<-readOGR(".","UCs_Cerrado")#Le shapefile

soma_BA_UC<-extract(BA.all.Cerrado,UCs.Cerrado,sum)
head(soma_BA_UC)
soma.BA.UC.table<-data.frame(soma_BA_UC,UCs.Cerrado$GRUPO4)
head(soma.BA.UC.table)
tail(soma.BA.UC.table)
str(soma.BA.UC.table)
somatotal.BA.UC.table<-aggregate(soma.BA.UC.table[,c(1:216)], list(soma.BA.UC.table$UCs.Cerrado.GRUPO4), sum, na.rm=T)
str(somatotal.BA.UC.table)

somatotal.BA.UC.table<-as.data.frame(t(somatotal.BA.UC.table[,-1]))
names(somatotal.BA.UC.table)<-c("PI","US")
str(somatotal.BA.UC.table)
head(somatotal.BA.UC.table)
tail(somatotal.BA.UC.table)

sapply(somatotal.BA.UC.table,mean)
boxplot(somatotal.BA.UC.table$PI, somatotal.BA.UC.table$US,ylim=c(0,3e+07))
summary(somatotal.BA.UC.table)

hist(somatotal.BA.UC.table$PI)
hist(log(somatotal.BA.UC.table$PI))
hist(somatotal.BA.UC.table$US)
hist(log(somatotal.BA.UC.table$US))
shapiro.test(log(somatotal.BA.UC.table$US))
shapiro.test(log(somatotal.BA.UC.table$PI))

t.test(somatotal.BA.UC.table$PI, somatotal.BA.UC.table$US,paired=T)
t.test(log(somatotal.BA.UC.table$PI), log(somatotal.BA.UC.table$US),paired=T)
wilcox.test(log(somatotal.BA.UC.table$PI), log(somatotal.BA.UC.table$US),paired=T)

diff.US.PI<-somatotal.BA.UC.table$US-somatotal.BA.UC.table$PI
summary(diff.US.PI)

write.table(somatotal.BA.UC.table,"BA_UCs.txt")

#Time series analyses in Conservation Units of integral protection
ts.soma.BA.PI<-ts(somatotal.BA.UC.table$PI,start=c(2000, 11), end=c(2018, 10), frequency=12)
ts.soma.BA.PI<-ts.soma.BA.PI*0.01 #Para transformar em ha multiplicar pelo fator 0.01 (Giglio et al., 2018 - Collection 6 MODIS Burned Area Product User?s Guide Version 1.2)

plot(ts.soma.BA.PI,ylab="?rea (ha)",xlab="Tempo (anos)")

#Time series decomposition
plot(decompose(ts.soma.BA.PI))
Queima.decom.PI <- decompose(ts.soma.BA.PI,type="multiplicative")
plot(Queima.decom.PI)
Trend.PI <- Queima.decom.PI$trend
Seasonal.PI<- Queima.decom.PI$seasonal
ts.plot(cbind(Trend.PI, Trend.PI * Seasonal.PI), lty = 1:2,ylab="?rea (ha)")

decomp.PI<-stl(ts.soma.BA.PI, "per")
summary(decomp.PI)
plot(decomp.PI)
monthplot(ts.soma.BA.PI, ylab="?rea (ha)", type="h")
monthplot(decomp.PI, choice="seasonal", type="h")
monthplot(decomp.PI, choice="trend", type="h")
monthplot(decomp.PI, choice="remainder", type="h")

# Percentage of explained variation by each component
p.sazonal.PI<-(var(decomp.PI$time.series[, 1]))/(var(decomp.PI$time.series[, 1])+var(decomp.PI$time.series[, 2])+var(decomp.PI$time.series[, 3]))
p.sazonal.PI

p.tendencia.PI<-(var(decomp.PI$time.series[, 2]))/(var(decomp.PI$time.series[, 1])+var(decomp.PI$time.series[, 2])+var(decomp.PI$time.series[, 3]))
p.tendencia.PI

p.residual.PI<-(var(decomp.PI$time.series[, 3]))/(var(decomp.PI$time.series[, 1])+var(decomp.PI$time.series[, 2])+var(decomp.PI$time.series[, 3]))
p.residual.PI

p.sazonal.PI+p.tendencia.PI+p.residual.PI


#Time series analyses for conservation units of sustainable use
ts.soma.BA.US<-ts(somatotal.BA.UC.table$US,start=c(2000, 11), end=c(2018, 10), frequency=12)
ts.soma.BA.US<-ts.soma.BA.US*0.01 #Para transformar em ha multiplicar pelo fator 0.01 (Giglio et al., 2018 - Collection 6 MODIS Burned Area Product User?s Guide Version 1.2)

plot(ts.soma.BA.US,ylab="?rea (ha)",xlab="Tempo (anos)")

#Time series decomposition
plot(decompose(ts.soma.BA.US))
Queima.decom.US <- decompose(ts.soma.BA.US,type="multiplicative")
plot(Queima.decom.US)
Trend.US <- Queima.decom.US$trend
Seasonal.US<- Queima.decom.US$seasonal
ts.plot(cbind(Trend.US, Trend.US * Seasonal.US), lty = 1:2,ylab="?rea (ha)")

decomp.US<-stl(ts.soma.BA.US, "per")
summary(decomp.US)
plot(decomp.US)
monthplot(ts.soma.BA.US, ylab="?rea (ha)", type="h")
monthplot(decomp.US, choice="seasonal", type="h")
monthplot(decomp.US, choice="trend", type="h")
monthplot(decomp.US, choice="remainder", type="h")

# Percentage of explained variation by each component
p.sazonal.US<-(var(decomp.US$time.series[, 1]))/(var(decomp.US$time.series[, 1])+var(decomp.US$time.series[, 2])+var(decomp.US$time.series[, 3]))
p.sazonal.US

p.tendencia.US<-(var(decomp.US$time.series[, 2]))/(var(decomp.US$time.series[, 1])+var(decomp.US$time.series[, 2])+var(decomp.US$time.series[, 3]))
p.tendencia.US

p.residual.US<-(var(decomp.US$time.series[, 3]))/(var(decomp.US$time.series[, 1])+var(decomp.US$time.series[, 2])+var(decomp.US$time.series[, 3]))
p.residual.US

p.sazonal.US+p.tendencia.US+p.residual.US




#Raster time series
install.packages("rts",dep=T)
library(rts)


rts.BA.Cerrado<-rts(BA.all.Cerrado,as.yearmon(time(ts.soma.BA)))
rts.BA.Cerrado

rts.BA.Cerrado.monthly<-apply.yearly(rts.BA.Cerrado,mean)
plot(rts.BA.Cerrado[[1]])

plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+2)]])#jan
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+3)]])#fev
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+4)]])#mar
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+5)]])#abr
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+6)]])#mai
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+7)]])#jun
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+8)]])#jul
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+9)]])#ago
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+10)]])#set
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+11)]])#out
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205))]])#nov
plot(rts.BA.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+1)]])#dez

#####################
#Transform to binary#
#####################
Bbin.all.Cerrado<-BA.all.Cerrado
Bbin.all.Cerrado[Bbin.all.Cerrado>0]<-1
writeRaster(Bbin.all.Cerrado,"Bbin.all.Cerrado.tif")
plot(Bbin.all.Cerrado[[1:12]])
Freq.all.Cerrado<-sum(Bbin.all.Cerrado) 
plot(Freq.all.Cerrado)
plot(Freq.all.Cerrado/216)
hist(Freq.all.Cerrado)
summary(Freq.all.Cerrado)
abline(v=49,col="red")

#January
Freq.jan.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+2)]]) 
plot(Freq.jan.Cerrado)
plot(Freq.jan.Cerrado/18)
hist(Freq.jan.Cerrado)
summary(Freq.jan.Cerrado)
abline(v=1,col="red")

#February
Freq.fev.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+3)]]) 
plot(Freq.fev.Cerrado)
plot(Freq.fev.Cerrado/18)
hist(Freq.fev.Cerrado)
summary(Freq.fev.Cerrado)
abline(v=1,col="red")

#March
Freq.mar.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+4)]]) 
plot(Freq.mar.Cerrado)
plot(Freq.mar.Cerrado/18)
hist(Freq.mar.Cerrado)
summary(Freq.mar.Cerrado)
abline(v=1,col="red")

#April
Freq.abr.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+5)]]) 
plot(Freq.abr.Cerrado)
plot(Freq.abr.Cerrado/18)
hist(Freq.abr.Cerrado)
summary(Freq.abr.Cerrado)
abline(v=1,col="red")

#May
Freq.mai.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+6)]]) 
plot(Freq.mai.Cerrado)
plot(Freq.mai.Cerrado/18)
hist(Freq.mai.Cerrado)
summary(Freq.mai.Cerrado)
abline(v=1,col="red")

#June
Freq.jun.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+7)]]) 
plot(Freq.jun.Cerrado)
plot(Freq.jun.Cerrado/18)
hist(Freq.jun.Cerrado)
summary(Freq.jun.Cerrado)
abline(v=2,col="red")

#July
Freq.jul.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+8)]]) 
plot(Freq.jul.Cerrado)
plot(Freq.jul.Cerrado/18)
hist(Freq.jul.Cerrado)
summary(Freq.jul.Cerrado)
abline(v=4,col="red")

#August
Freq.ago.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+9)]]) 
plot(Freq.ago.Cerrado)
plot(Freq.ago.Cerrado/18)
hist(Freq.ago.Cerrado)
summary(Freq.ago.Cerrado)
abline(v=8,col="red")

#September
Freq.set.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+10)]]) 
plot(Freq.set.Cerrado)
plot(Freq.set.Cerrado/18)
hist(Freq.set.Cerrado)
summary(Freq.set.Cerrado)
abline(v=11,col="red")

#October
Freq.out.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+11)]]) 
plot(Freq.out.Cerrado)
plot(Freq.out.Cerrado/18)
hist(Freq.out.Cerrado)
summary(Freq.out.Cerrado)
abline(v=9,col="red")

#November
Freq.nov.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205))]]) 
plot(Freq.nov.Cerrado)
plot(Freq.nov.Cerrado/18)
hist(Freq.nov.Cerrado)
summary(Freq.nov.Cerrado)
abline(v=3,col="red")

#December
Freq.dez.Cerrado<-sum(Bbin.all.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205)+1)]]) 
plot(Freq.dez.Cerrado)
plot(Freq.dez.Cerrado/18)
hist(Freq.dez.Cerrado)
summary(Freq.dez.Cerrado)
abline(v=1,col="red")

stack
par(mfrow=c(3,4))
plot(stack((Freq.jan.Cerrado/18),(Freq.fev.Cerrado/18),(Freq.mar.Cerrado/18),(Freq.abr.Cerrado/18),(Freq.mai.Cerrado/18),(Freq.jun.Cerrado/18),
(Freq.jul.Cerrado/18),(Freq.ago.Cerrado/18),(Freq.set.Cerrado/18),(Freq.out.Cerrado/18),(Freq.nov.Cerrado/18),(Freq.dez.Cerrado/18)))


#Using data from global fire atlas: https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1642

fire.ign<-stack(c("E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2003.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2004.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2005.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2006.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2007.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2008.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2009.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2010.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2011.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2012.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2013.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2014.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2015.tif",
  "E:/Heitor/BD_SIG/CMS_Global_Fire_Atlas_1642/data/Global_fire_atlas_day_of_burn_yearly_2016.tif"))

fire.ign
Cerrado
plot(fire.ign)

plot(fire.ign[[1]])
plot(Cerrado,add=T,col="blue",lwd=2)

Cerrado.proj<-spTransform(Cerrado,"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

plot(fire.ign[[1]])
plot(Cerrado.proj,add=T,col="blue",lwd=2)

#Crop rasters to Cerrado extent


fire.ign.Cerrado<-crop(fire.ign, extent(Cerrado.proj))
fire.ign.Cerrado <- mask(fire.ign.Cerrado, Cerrado.proj)

plot(fire.ign.Cerrado)
plot(fire.ign.Cerrado[[1]])
plot(Cerrado.proj,add=T,col="NA",lwd=1)

fire.ign.Cerrado.jan<-reclassify(fire.ign.Cerrado,c(1,31,1,31,366,0))
fire.ign.Cerrado.fev<-reclassify(fire.ign.Cerrado,c(1,31,0,31,60,1,60,366,0))
fire.ign.Cerrado.mar<-reclassify(fire.ign.Cerrado,c(1,60,0,60,91,1,91,366,0))
fire.ign.Cerrado.abr<-reclassify(fire.ign.Cerrado,c(1,91,0,91,121,1,121,366,0))
fire.ign.Cerrado.mai<-reclassify(fire.ign.Cerrado,c(1,121,0,121,152,1,152,366,0))
fire.ign.Cerrado.jun<-reclassify(fire.ign.Cerrado,c(1,152,0,152,182,1,182,366,0))
fire.ign.Cerrado.jul<-reclassify(fire.ign.Cerrado,c(1,182,0,182,213,1,213,366,0))
fire.ign.Cerrado.ago<-reclassify(fire.ign.Cerrado,c(1,213,0,213,243,1,243,366,0))
fire.ign.Cerrado.set<-reclassify(fire.ign.Cerrado,c(1,243,0,243,274,1,274,366,0))
fire.ign.Cerrado.out<-reclassify(fire.ign.Cerrado,c(1,274,0,274,304,1,304,366,0))
fire.ign.Cerrado.nov<-reclassify(fire.ign.Cerrado,c(1,304,0,304,335,1,335,366,0))
fire.ign.Cerrado.dez<-reclassify(fire.ign.Cerrado,c(1,335,0,335,366,1))

fire.ign.Cerrado.jan
freq(fire.ign.Cerrado.jan)
plot(fire.ign.Cerrado.jan)
freq.fire.Cerrado.jan<-sum(fire.ign.Cerrado.jan,na.rm=TRUE)
plot(freq.fire.Cerrado.jan)
freq.fire.Cerrado.fev<-sum(fire.ign.Cerrado.fev,na.rm=TRUE)
freq.fire.Cerrado.mar<-sum(fire.ign.Cerrado.mar,na.rm=TRUE)
freq.fire.Cerrado.abr<-sum(fire.ign.Cerrado.abr,na.rm=TRUE)
freq.fire.Cerrado.mai<-sum(fire.ign.Cerrado.mai,na.rm=TRUE)
freq.fire.Cerrado.jun<-sum(fire.ign.Cerrado.jun,na.rm=TRUE)
freq.fire.Cerrado.jul<-sum(fire.ign.Cerrado.jul,na.rm=TRUE)
freq.fire.Cerrado.ago<-sum(fire.ign.Cerrado.ago,na.rm=TRUE)
freq.fire.Cerrado.set<-sum(fire.ign.Cerrado.set,na.rm=TRUE)
freq.fire.Cerrado.out<-sum(fire.ign.Cerrado.out,na.rm=TRUE)
freq.fire.Cerrado.nov<-sum(fire.ign.Cerrado.nov,na.rm=TRUE)
freq.fire.Cerrado.dez<-sum(fire.ign.Cerrado.dez,na.rm=TRUE)

freq.fire.Cerrado<-stack(freq.fire.Cerrado.jan,freq.fire.Cerrado.fev,freq.fire.Cerrado.mar,freq.fire.Cerrado.abr,freq.fire.Cerrado.mai,freq.fire.Cerrado.jun,
freq.fire.Cerrado.jul,freq.fire.Cerrado.ago,freq.fire.Cerrado.set,freq.fire.Cerrado.out,freq.fire.Cerrado.nov,freq.fire.Cerrado.dez)

plot(freq.fire.Cerrado)
freq.fire.Cerrado<-crop(freq.fire.Cerrado, extent(Cerrado.proj))
freq.fire.Cerrado <- mask(freq.fire.Cerrado, Cerrado.proj)
plot(freq.fire.Cerrado[[1]])

writeRaster(freq.fire.Cerrado,"freq.fire.Cerrado.mes.grd")

freq.fire.chuva<-freq.fire.Cerrado[[1]]+freq.fire.Cerrado[[2]]+freq.fire.Cerrado[[11]]+freq.fire.Cerrado[[12]]
freq.fire.precoce<-freq.fire.Cerrado[[3]]+freq.fire.Cerrado[[4]]+freq.fire.Cerrado[[5]]
freq.fire.modal<-freq.fire.Cerrado[[6]]+freq.fire.Cerrado[[7]]
freq.fire.tardia<-freq.fire.Cerrado[[8]]+freq.fire.Cerrado[[9]]+freq.fire.Cerrado[[10]]

writeRaster(freq.fire.chuva,"freq.fire.chuva.grd")
writeRaster(freq.fire.precoce,"freq.fire.precoce.grd")
writeRaster(freq.fire.modal,"freq.fire.modal.grd")
writeRaster(freq.fire.tardia,"freq.fire.tardia.grd")

#################
#AVHRR LTDR Fire#
#################

# If you have lots of files then you can make a loop to do all this for you
setwd("/Volumes/Extreme SSD/Heitor/BD_SIG/AVHRR_LTDR_Fire/Pixel")
files<-list.files(pattern = "-BA-AVHRR-LTDR-fv1.1-BA.tif", recursive = TRUE)
files

filename <- substr(files,87,94)
filename <- paste0("BA", filename)
filename

BA.all.LTDR<-stack(files)
BA.all.LTDR
plot(BA.all.LTDR[[1]])

#########################
#Reads Cerrado shapefile#
#########################
Cerrado<-readOGR(".","Cerrado")#Le shapefile
plot(BA.all.LTDR[[1]])
plot(Cerrado,add=T)

#Crop rasters to Cerrado extent
BA.all.LTDR.Cerrado<-crop(BA.all.LTDR, extent(Cerrado))
BA.all.LTDR.Cerrado <- mask(BA.all.LTDR.Cerrado, Cerrado)

plot(BA.all.LTDR.Cerrado[[1]])
plot(Cerrado,axes=T,add=T)

writeRaster(BA.all.LTDR.Cerrado,"BA.all.LTDR.Cerrado.tif")

#Reads raster (stack)##########################
BA.all.LTDR.Cerrado<-stack("BA.all.LTDR.Cerrado.tif")
############################################

# BA.all.LTDR.Cerrado
# crs(BA.all.LTDR.Cerrado)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# Cerrado
# plot(BA.all.LTDR.Cerrado)
# plot(BA.all.LTDR.Cerrado[[1]])
# plot(Cerrado,axes=T,add=T)


###########################################################
#Calculates sum of burned area by month for Cerrado extent#
###########################################################
soma_BA_LTDR<-cellStats(BA.all.LTDR.Cerrado,sum)
soma_BA_LTDR[144]
length(soma_BA_LTDR)
soma_BA_LTDR<-c(soma_BA_LTDR[1:144],rep(NA,12),soma_BA_LTDR[145:432])


#Time series analysis
ts.soma.BA.LTDR<-ts(soma_BA_LTDR,start=c(1982, 01), end=c(2018, 12), frequency=12)
ts.soma.BA.LTDR
#ts.soma.BA.LTDR<-ts.soma.BA*0.01 #Para transformar em ha multiplicar pelo fator 0.01 (Giglio et al., 2018 - Collection 6 MODIS Burned Area Product User?s Guide Version 1.2)

plot(ts.soma.BA.LTDR,ylab="Área (ha)",xlab="Tempo (anos)")

#Time series decomposition
plot(decompose(ts.soma.BA.LTDR))
Queima.decom <- decompose(ts.soma.BA.LTDR,type="multiplicative")
plot(Queima.decom)
Trend <- Queima.decom$trend
Seasonal <- Queima.decom$seasonal
ts.plot(cbind(Trend, Trend * Seasonal), lty = 1:2,ylab="?rea (ha)")

decomp<-stl(ts.soma.BA.LTDR, "per")
summary(decomp)
plot(decomp)
monthplot(ts.soma.BA.LTDR, ylab="?rea (ha)", type="h")
monthplot(decomp, choice="seasonal")
monthplot(decomp, choice="trend")
monthplot(decomp, choice="remainder")

# Porcentagem da variacao explicada por cada componente da serie temporal
p.sazonal<-(var(decomp$time.series[, 1]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.sazonal

p.tendencia<-(var(decomp$time.series[, 2]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.tendencia

p.residual<-(var(decomp$time.series[, 3]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.residual

p.sazonal+p.tendencia+p.residual

#Raster time series - bom para visualizar entre anos x meses
#install.packages("rts",dep=T)
library(rts)

rts.BA.Cerrado.LTDR<-rts(BA.all.LTDR.Cerrado,as.yearmon(time(ts.soma.BA.LTDR)))
rts.BA.Cerrado.LTDR

rts.BA.Cerrado.LTDR.monthly<-apply.yearly(rts.BA.Cerrado.LTDR,mean)
rts.BA.Cerrado.LTDR.monthly
writeRaster(rts.BA.Cerrado.LTDR.monthly,"FireExtent.Cerrado.LTDR.monthly.grd")
plot(rts.BA.Cerrado.LTDR[[1]])

windows(10,10)
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421))]])#jan

FireExtent.Cerrado.jan<-mean(BA.all.LTDR.Cerrado[[c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)]],na.rm=T)

FireExtent.sd.Cerrado.jan<-cv(BA.all.LTDR.Cerrado[[c(1,13,25,37,49,61,73,85,97,
                                                     109,121,133,145,157,169,181,
                                                     193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)]],
                                                  na.rm=T)

plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+1)]])#fev

FireExtent.Cerrado.fev<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+1)]],na.rm=T)

FireExtent.sd.Cerrado.fev<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,
                                                     109,121,133,145,157,169,181,
                                                     193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,
                                                     373,385,397,409,421)+1)]],na.rm=T)

plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+2)]])#mar

FireExtent.Cerrado.mar<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,
                                                     109,121,133,145,157,169,181,
                                                     193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+2)]],
                             na.rm=T)

FireExtent.sd.Cerrado.mar<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,
                                                     109,121,133,145,157,169,181,
                                                     193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+2)]],
                             na.rm=T)
  
plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+3)]])#abr

FireExtent.Cerrado.abr<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,
                                                     109,121,133,145,157,169,181,
                                                     193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+3)]],
                             na.rm=T)

FireExtent.sd.Cerrado.abr<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,
                                                     109,121,133,145,157,169,181,
                                                     193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,
                                                     373,385,397,409,421)+3)]],
                             na.rm=T)

plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+4)]])#mai

FireExtent.Cerrado.mai<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,
                                                     109,121,133,145,157,169,181,
                                                     193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+4)]],
                             na.rm=T)

FireExtent.sd.Cerrado.mai<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,
                                                     109,121,133,145,157,169,181,
                                                     193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,
                                                   373,385,397,409,421)+4)]],
                             na.rm=T)

plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+5)]])#jun

FireExtent.Cerrado.jun<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+5)]],na.rm=T)

FireExtent.sd.Cerrado.jun<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,
                                                        109,121,133,145,157,169,181,
                                                        193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,
                                                     385,397,409,421)+5)]],na.rm=T)

plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+6)]])#jul

FireExtent.Cerrado.jul<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+6)]],na.rm=T)

FireExtent.sd.Cerrado.jul<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+6)]],na.rm=T)

plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+7)]])#ago

FireExtent.Cerrado.ago<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+7)]],na.rm=T)

FireExtent.sd.Cerrado.ago<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+7)]],na.rm=T)

plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+8)]])#set

FireExtent.Cerrado.set<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+8)]],na.rm=T)

FireExtent.sd.Cerrado.set<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+8)]],na.rm=T)

plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+9)]])#out

FireExtent.Cerrado.out<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+9)]],na.rm=T)

FireExtent.sd.Cerrado.out<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+9)]],na.rm=T)

plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+10)]])#nov

FireExtent.Cerrado.nov<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+10)]],na.rm=T)

FireExtent.sd.Cerrado.nov<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+10)]],na.rm=T)

plot(rts.BA.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                        289,301,313,325,337,349,361,373,385,397,409,421)+11)]])#dez

FireExtent.Cerrado.dez<-mean(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                289,301,313,325,337,349,361,373,385,397,409,421)+11)]],na.rm=T)

FireExtent.sd.Cerrado.dez<-cv(BA.all.LTDR.Cerrado[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+11)]],na.rm=T)

FireExtent.Cerrado.LTDR.monthly<-stack(FireExtent.Cerrado.jan,
                                       FireExtent.Cerrado.fev,
                                       FireExtent.Cerrado.mar,
                                       FireExtent.Cerrado.abr,
                                       FireExtent.Cerrado.mai,
                                       FireExtent.Cerrado.jun,
                                       FireExtent.Cerrado.jul,
                                       FireExtent.Cerrado.ago,
                                       FireExtent.Cerrado.set,
                                       FireExtent.Cerrado.out,
                                       FireExtent.Cerrado.nov,
                                       FireExtent.Cerrado.dez)


FireExtent.Cerrado.LTDR.monthly[FireExtent.Cerrado.LTDR.monthly<0]<-NA
FireExtent.Cerrado.LTDR.monthly

FireExtent.sd.Cerrado.LTDR.monthly<-stack(FireExtent.sd.Cerrado.jan,
                                       FireExtent.sd.Cerrado.fev,
                                       FireExtent.sd.Cerrado.mar,
                                       FireExtent.sd.Cerrado.abr,
                                       FireExtent.sd.Cerrado.mai,
                                       FireExtent.sd.Cerrado.jun,
                                       FireExtent.sd.Cerrado.jul,
                                       FireExtent.sd.Cerrado.ago,
                                       FireExtent.sd.Cerrado.set,
                                       FireExtent.sd.Cerrado.out,
                                       FireExtent.sd.Cerrado.nov,
                                       FireExtent.sd.Cerrado.dez)

plot(FireExtent.sd.Cerrado.LTDR.monthly)
plot(density(FireExtent.sd.Cerrado.LTDR.monthly))

writeRaster(FireExtent.Cerrado.LTDR.monthly,"FireExtent.Cerrado.LTDR.monthly.grd")
writeRaster(FireExtent.sd.Cerrado.LTDR.monthly,"FireExtent.sd.Cerrado.LTDR.monthly.grd")


########################
#Transformar em bin?rio#
########################
Bbin.all.Cerrado.LTDR<-BA.all.LTDR.Cerrado
Bbin.all.Cerrado.LTDR[Bbin.all.Cerrado.LTDR>0]<-1
Bbin.all.Cerrado.LTDR[Bbin.all.Cerrado.LTDR<0]<-0
writeRaster(Bbin.all.Cerrado.LTDR,"Bbin.all.Cerrado.LTDR.tif",overwrite=T)

Bbin.all.Cerrado.LTDR <- stack("Bbin.all.Cerrado.LTDR.tif")

plot(Bbin.all.Cerrado.LTDR[[1:12]])
Freq.all.Cerrado.LTDR<-sum(Bbin.all.Cerrado.LTDR) 
plot(Freq.all.Cerrado.LTDR)
plot(Freq.all.Cerrado.LTDR/432)
hist(Freq.all.Cerrado.LTDR)
summary(Freq.all.Cerrado.LTDR)
abline(v=17,col="red")

#Para janeiro
Freq.LTDR.jan.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)]]) 
plot(Freq.LTDR.jan.Cerrado)
plot(Freq.LTDR.jan.Cerrado/36)
hist(Freq.LTDR.jan.Cerrado)
summary(Freq.LTDR.jan.Cerrado)
abline(v=1,col="red")

Freq.sd.LTDR.jan.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                    289,301,313,325,337,349,361,373,385,397,409,421)]], aszero=T) 

#Para fevereiro
Freq.LTDR.fev.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+1)]]) 
plot(Freq.LTDR.fev.Cerrado)
plot(Freq.LTDR.fev.Cerrado/36)
hist(Freq.LTDR.fev.Cerrado)
summary(Freq.LTDR.fev.Cerrado)
abline(v=1,col="red")

Freq.sd.LTDR.fev.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+1)]], aszero=T) 

#Para mar?o
Freq.LTDR.mar.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+2)]]) 
plot(Freq.LTDR.mar.Cerrado)
plot(Freq.LTDR.mar.Cerrado/36)
hist(Freq.LTDR.mar.Cerrado)
summary(Freq.LTDR.mar.Cerrado)
abline(v=1,col="red")

Freq.sd.LTDR.mar.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+2)]],aszero=T) 

#Para abril
Freq.LTDR.abr.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+3)]]) 
plot(Freq.LTDR.abr.Cerrado)
plot(Freq.LTDR.abr.Cerrado/36)
hist(Freq.LTDR.abr.Cerrado)
summary(Freq.LTDR.abr.Cerrado)
abline(v=1,col="red")

Freq.sd.LTDR.abr.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+3)]],aszero=T) 

#Para maio
Freq.LTDR.mai.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+4)]]) 
plot(Freq.LTDR.mai.Cerrado)
plot(Freq.LTDR.mai.Cerrado/36)
hist(Freq.LTDR.mai.Cerrado)
summary(Freq.LTDR.mai.Cerrado)
abline(v=1,col="red")

Freq.sd.LTDR.mai.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+4)]],aszero = T) 

#Para junho
Freq.LTDR.jun.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+5)]]) 
plot(Freq.LTDR.jun.Cerrado)
plot(Freq.LTDR.jun.Cerrado/36)
hist(Freq.LTDR.jun.Cerrado)
summary(Freq.LTDR.jun.Cerrado)
abline(v=2,col="red")

Freq.sd.LTDR.jun.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+5)]], aszero = T) 

#Para julho
Freq.LTDR.jul.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+6)]]) 
plot(Freq.LTDR.jul.Cerrado)
plot(Freq.LTDR.jul.Cerrado/36)
hist(Freq.LTDR.jul.Cerrado)
summary(Freq.LTDR.jul.Cerrado)
abline(v=4,col="red")

Freq.sd.LTDR.jul.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+6)]], aszero = T) 

#Para agosto
Freq.LTDR.ago.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+7)]]) 
plot(Freq.LTDR.ago.Cerrado)
plot(Freq.LTDR.ago.Cerrado/36)
hist(Freq.LTDR.ago.Cerrado)
summary(Freq.LTDR.ago.Cerrado)
abline(v=8,col="red")

Freq.sd.LTDR.ago.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+7)]], aszero = T) 

#Para setembro
Freq.LTDR.set.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+8)]]) 
plot(Freq.LTDR.set.Cerrado)
plot(Freq.LTDR.set.Cerrado/36)
hist(Freq.LTDR.set.Cerrado)
summary(Freq.LTDR.set.Cerrado)
abline(v=11,col="red")

Freq.sd.LTDR.set.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+8)]], aszero = T) 

#Para outubro
Freq.LTDR.out.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+9)]]) 
plot(Freq.LTDR.out.Cerrado)
plot(Freq.LTDR.out.Cerrado/36)
hist(Freq.LTDR.out.Cerrado)
summary(Freq.LTDR.out.Cerrado)
abline(v=9,col="red")

Freq.sd.LTDR.out.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+9)]], aszero = T) 

#Para novembro
Freq.LTDR.nov.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+10)]]) 
plot(Freq.LTDR.nov.Cerrado)
plot(Freq.LTDR.nov.Cerrado/36)
hist(Freq.LTDR.nov.Cerrado)
summary(Freq.LTDR.nov.Cerrado)
abline(v=3,col="red")

Freq.sd.LTDR.nov.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+10)]], aszero = T)

#Para dezembro
Freq.LTDR.dez.Cerrado<-sum(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                           289,301,313,325,337,349,361,373,385,397,409,421)+11)]]) 
plot(Freq.LTDR.dez.Cerrado)
plot(Freq.LTDR.dez.Cerrado/36)
hist(Freq.LTDR.dez.Cerrado)
summary(Freq.LTDR.dez.Cerrado)
abline(v=1,col="red")

Freq.sd.LTDR.dez.Cerrado<-cv(Bbin.all.Cerrado.LTDR[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229,241,253,265,277,
                                                     289,301,313,325,337,349,361,373,385,397,409,421)+11)]], aszero = T) 

par(mfrow=c(3,4))
plot(stack((Freq.LTDR.jan.Cerrado/36),(Freq.LTDR.fev.Cerrado/36),(Freq.LTDR.mar.Cerrado/36),
           (Freq.LTDR.abr.Cerrado/36),(Freq.LTDR.mai.Cerrado/36),(Freq.LTDR.jun.Cerrado/36),
           (Freq.LTDR.jul.Cerrado/36),(Freq.LTDR.ago.Cerrado/36),(Freq.LTDR.set.Cerrado/36),
           (Freq.LTDR.out.Cerrado/36),(Freq.LTDR.nov.Cerrado/36),(Freq.LTDR.dez.Cerrado/36)))

Freq.fire.LTDR.Cerrado<-stack((Freq.LTDR.jan.Cerrado/36),(Freq.LTDR.fev.Cerrado/36),(Freq.LTDR.mar.Cerrado/36),
                              (Freq.LTDR.abr.Cerrado/36),(Freq.LTDR.mai.Cerrado/36),(Freq.LTDR.jun.Cerrado/36),
                              (Freq.LTDR.jul.Cerrado/36),(Freq.LTDR.ago.Cerrado/36),(Freq.LTDR.set.Cerrado/36),
                              (Freq.LTDR.out.Cerrado/36),(Freq.LTDR.nov.Cerrado/36),(Freq.LTDR.dez.Cerrado/36))

Freq.sd.fire.LTDR.Cerrado<-stack(Freq.sd.LTDR.jan.Cerrado,
                                 Freq.sd.LTDR.fev.Cerrado,
                                 Freq.sd.LTDR.mar.Cerrado,
                                 Freq.sd.LTDR.abr.Cerrado,
                                 Freq.sd.LTDR.mai.Cerrado,
                                 Freq.sd.LTDR.jun.Cerrado,
                                 Freq.sd.LTDR.jul.Cerrado,
                                 Freq.sd.LTDR.ago.Cerrado,
                                 Freq.sd.LTDR.set.Cerrado,
                                 Freq.sd.LTDR.out.Cerrado,
                                 Freq.sd.LTDR.nov.Cerrado,
                                 Freq.sd.LTDR.dez.Cerrado)

plot(Freq.sd.fire.LTDR.Cerrado)
plot(density(Freq.sd.fire.LTDR.Cerrado))

writeRaster(Freq.fire.LTDR.Cerrado,"Freq.fire.LTDR.Cerrado.grd")

writeRaster(Freq.sd.fire.LTDR.Cerrado,"Freq.sd.fire.LTDR.Cerrado.grd")

##############
#MCD641A data#
##############
#Resolution 500 m
#Pixel values:
#>=1 = burned (actually, the value is the julian day when the pixel burned. These values will be set to 1)
#0 = unburned
#-1 = unmapped due to unsufficient data (These values will be set to NA)
#-2 = water (These values will be set to NA)

rasterOptions(tmpdir=file.path("/Users/heito/Documents/raster"),
              chunksize = 1e+06,maxmemory = 1e+09,
            progress = "text",
            todisk = T) 
rasterOptions()

# If you have lots of files then you can make a loop to do all this for you
setwd("/Volumes/Extreme SSD/Heitor/BD_SIG/MODIS_Fire/MCD64A1")

files<-list.files(pattern = "_Burn_Date_doy", recursive = TRUE)
files

filename <- substr(files,26,32)
filename


BA.all.MCD64A1<-stack(files)
names(BA.all.MCD64A1)<-filename
BA.all.MCD64A1
windows(10,10)
plot(BA.all.MCD64A1[[240]])

########################
#Transformar em binario#
########################
Bbin.all.Cerrado.MCD64A1<-BA.all.MCD64A1

beginCluster(7)

z1 <- clusterR(Bbin.all.Cerrado.MCD64A1,reclassify, args=list(c(-Inf,-0.5,NA,0.5,Inf,1)))

save(z1,file="Bbin.all.Cerrado.MCD64A1.RData")

writeRaster(z1,"Bbin.all.Cerrado.MCD64A1.grd",overwrite=T)

z1
names(z1)<-filename

windows(10,10)
plot(z1[[1:12]])
Freq.all.Cerrado.MCD64A1<-sum(z1,na.rm=T) 
plot(Freq.all.Cerrado.MCD64A1)
plot(Freq.all.Cerrado.MCD64A1/nlayers(z1))

#Cortar rasters para Cerrado
setwd("E:/Heitor/BD_SIG/MODIS_Fire")
Cerrado<-readOGR(".","Cerrado")#Le shapefile
setwd("E:/Heitor/BD_SIG/MODIS_Fire/MCD64A1")
Freq.all.Cerrado.MCD64A1.crop<-crop(Freq.all.Cerrado.MCD64A1, extent(Cerrado))
Freq.all.Cerrado.MCD64A1.crop <- mask(Freq.all.Cerrado.MCD64A1.crop, Cerrado)

windows(10,10)
plot(Freq.all.Cerrado.MCD64A1.crop/nlayers(z1))

#save(Freq.all.Cerrado.MCD64A1.crop,file="Freq.all.Cerrado.MCD64A1.RData")

#writeRaster(Freq.all.Cerrado.MCD64A1.crop,"Freq.all.Cerrado.MCD64A1.grd",overwrite=T)

setwd("/Volumes/Extreme SSD/Heitor/BD_SIG/MODIS_Fire/MCD64A1")
z1<-stack("Bbin.all.Cerrado.MCD64A1.grd")
dim(z1)

write("TMPDIR = '<E:/Heitor/BD_SIG/MODIS_Fire/MCD64A1>'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

for(i in 1:nlayers(z1)){
  v<-rasterToPoints(z1[[i]])
  write.table(v,paste0("df",i,".txt"),sep = "\t")
}


hist(Freq.all.Cerrado.MCD64A1)
summary(Freq.all.Cerrado.MCD64A1)
abline(v=17,col="red")

#########################################################
#Calcula soma de ?rea queimada por m?s para todo o bioma#
#########################################################
#soma_BA_MCD64A1<-cellStats(z1,sum)
soma_BA_MCD64A1<-read.table("soma_BA_MCD64A1.txt",h=T)
soma_BA_MCD64A1<-soma_BA_MCD64A1[,2]
str(soma_BA_MCD64A1[144])

#write.table(soma_BA_MCD64A1,"soma_BA_MCD64A1.txt",sep="\t")
soma_BA_MCD64A1 <- read.table("soma_BA_MCD64A1.txt",h=T)
length(soma_BA_MCD64A1)

#Analise da serie temporal
ts.soma.BA.MCD64A1<-ts(soma_BA_MCD64A1$BA,start=c(2001, 01), end=c(2020, 12), frequency=12)
ts.soma.BA.MCD64A1
ts.soma.BA.MCD64A1<-(ts.soma.BA.MCD64A1*500^2)/10000 #Para transformar em ha

#windows(10,10)
months <- rep(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),39)
quartz(8,12)
plot(soma_BA_LTDR/10000,ylab="Area (ha)",xlab="",type="l",bty="n",
     axes=F,ylim=c(0,14000000),xlim = c(1,468))
lines(229:468,(soma_BA_MCD64A1$BA*500^2)/10000,col="red")
axis(1,1:468,xaxp = c(1,468,467),labels = months, las=2, cex.axis = .5)
axis(2)


#Decomposi??o da s?rie temporal
plot(decompose(ts.soma.BA.MCD64A1))
Queima.decom <- decompose(ts.soma.BA.MCD64A1,type="multiplicative")
plot(Queima.decom)
Trend <- Queima.decom$trend
Seasonal <- Queima.decom$seasonal
ts.plot(cbind(Trend, Trend * Seasonal), lty = 1:2,ylab="Área (ha)")

decomp<-stl(ts.soma.BA.MCD64A1, "per")
summary(decomp)
plot(decomp)
monthplot(ts.soma.BA.MCD64A1, ylab="Área (ha)", type="h")
monthplot(decomp, choice="seasonal")
monthplot(decomp, choice="trend")
monthplot(decomp, choice="remainder")

# Porcentagem da variacao explicada por cada componente da serie temporal
p.sazonal<-(var(decomp$time.series[, 1]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.sazonal

p.tendencia<-(var(decomp$time.series[, 2]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.tendencia

p.residual<-(var(decomp$time.series[, 3]))/(var(decomp$time.series[, 1])+var(decomp$time.series[, 2])+var(decomp$time.series[, 3]))
p.residual

p.sazonal+p.tendencia+p.residual

saveRDS(soma_BA_MCD64A1,file="soma_BA_MCD64A1.RData")

#Raster time series - bom para visualizar entre anos x meses
#install.packages("rts",dep=T)
library(rts)
Sys.setenv(TZ = "America/Sao_Paulo")

rts.BA.Cerrado.MCD64A1<-rts(z1,as.yearmon(time(ts.soma.BA.MCD64A1)))
rts.BA.Cerrado.MCD64A1

#rts.BA.Cerrado.MCD64A1.monthly<-apply.yearly(rts.BA.Cerrado.MCD64A1,mean)
plot(rts.BA.Cerrado.MCD64A1[[1]])

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229))]])#jan

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+1)]])#fev

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+2)]])#mar

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+3)]])#abr

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+4)]])#mai

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+5)]])#jun

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+6)]])#jul

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+7)]])#ago

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+8)]])#set

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+9)]])#out

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+10)]])#nov

windows(10,10)
plot(rts.BA.Cerrado.MCD64A1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+11)]])#dez


#Para janeiro
Freq.MCD64A1.jan.Cerrado<-sum(z1[[c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)]],na.rm=T) 
Freq.MCD64A1.jan.Cerrado<-crop(Freq.MCD64A1.jan.Cerrado, extent(Cerrado))
Freq.MCD64A1.jan.Cerrado <- mask(Freq.MCD64A1.jan.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.jan.Cerrado,main="Jan")
# windows(10,10)
# plot(Freq.MCD64A1.jan.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.jan.Cerrado)
summary(Freq.MCD64A1.jan.Cerrado)
#abline(v=1,col="red")

Freq.sd.MCD64A1.jan.Cerrado<-cv(z1[[c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)]],na.rm=T, aszero = T) 

#Para fevereiro
Freq.MCD64A1.fev.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+1)]],na.rm=T) 
Freq.MCD64A1.fev.Cerrado<-crop(Freq.MCD64A1.fev.Cerrado, extent(Cerrado))
Freq.MCD64A1.fev.Cerrado <- mask(Freq.MCD64A1.fev.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.fev.Cerrado,main="Fev")
# windows(10,10)
# plot(Freq.MCD64A1.fev.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.fev.Cerrado)
summary(Freq.MCD64A1.fev.Cerrado)
#abline(v=1,col="red")

Freq.sd.MCD64A1.fev.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+1)]],na.rm=T, aszero = T) 

#Para marco
Freq.MCD64A1.mar.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+2)]],na.rm=T) 
Freq.MCD64A1.mar.Cerrado<-crop(Freq.MCD64A1.mar.Cerrado, extent(Cerrado))
Freq.MCD64A1.mar.Cerrado <- mask(Freq.MCD64A1.mar.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.mar.Cerrado,main="Mar")
# windows(10,10)
# plot(Freq.MCD64A1.mar.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.mar.Cerrado)
summary(Freq.MCD64A1.mar.Cerrado)
#abline(v=1,col="red")

Freq.sd.MCD64A1.mar.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+2)]],na.rm=T, aszero = T) 

#Para abril
Freq.MCD64A1.abr.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+3)]],na.rm=T) 
Freq.MCD64A1.abr.Cerrado<-crop(Freq.MCD64A1.abr.Cerrado, extent(Cerrado))
Freq.MCD64A1.abr.Cerrado <- mask(Freq.MCD64A1.abr.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.abr.Cerrado,main="Abr")
# windows(10,10)
# plot(Freq.MCD64A1.abr.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.abr.Cerrado)
summary(Freq.MCD64A1.abr.Cerrado)
#abline(v=1,col="red")

Freq.sd.MCD64A1.abr.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+3)]],na.rm=T, aszero = T) 

#Para maio
Freq.MCD64A1.mai.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+4)]],na.rm=T) 
Freq.MCD64A1.mai.Cerrado<-crop(Freq.MCD64A1.mai.Cerrado, extent(Cerrado))
Freq.MCD64A1.mai.Cerrado <- mask(Freq.MCD64A1.mai.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.mai.Cerrado,main="Mai")
# windows(10,10)
# plot(Freq.MCD64A1.mai.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.mai.Cerrado)
summary(Freq.MCD64A1.mai.Cerrado)
#abline(v=1,col="red")

Freq.sd.MCD64A1.mai.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+4)]],na.rm=T, aszero = T) 

#Para junho
Freq.MCD64A1.jun.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+5)]],na.rm=T) 
Freq.MCD64A1.jun.Cerrado<-crop(Freq.MCD64A1.jun.Cerrado, extent(Cerrado))
Freq.MCD64A1.jun.Cerrado <- mask(Freq.MCD64A1.jun.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.jun.Cerrado,main="Jun")
# windows(10,10)
# plot(Freq.MCD64A1.jun.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.jun.Cerrado)
summary(Freq.MCD64A1.jun.Cerrado)
#abline(v=2,col="red")

Freq.sd.MCD64A1.jun.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+5)]],na.rm=T, aszero = T) 

#Para julho
Freq.MCD64A1.jul.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+6)]],na.rm=T) 
Freq.MCD64A1.jul.Cerrado<-crop(Freq.MCD64A1.jul.Cerrado, extent(Cerrado))
Freq.MCD64A1.jul.Cerrado <- mask(Freq.MCD64A1.jul.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.jul.Cerrado,main="Jul")
# windows(10,10)
# plot(Freq.MCD64A1.jul.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.jul.Cerrado)
summary(Freq.MCD64A1.jul.Cerrado)
#abline(v=4,col="red")

Freq.sd.MCD64A1.jul.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+6)]],na.rm=T, aszero = T)

#Para agosto
Freq.MCD64A1.ago.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+7)]],na.rm=T) 
Freq.MCD64A1.ago.Cerrado<-crop(Freq.MCD64A1.ago.Cerrado, extent(Cerrado))
Freq.MCD64A1.ago.Cerrado <- mask(Freq.MCD64A1.ago.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.ago.Cerrado,main="Ago")
# windows(10,10)
# plot(Freq.MCD64A1.ago.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.ago.Cerrado)
summary(Freq.MCD64A1.ago.Cerrado)
#abline(v=8,col="red")

Freq.sd.MCD64A1.ago.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+7)]],na.rm=T, aszero = T) 

#Para setembro
Freq.MCD64A1.set.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+8)]],na.rm=T) 
Freq.MCD64A1.set.Cerrado<-crop(Freq.MCD64A1.set.Cerrado, extent(Cerrado))
Freq.MCD64A1.set.Cerrado <- mask(Freq.MCD64A1.set.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.set.Cerrado,main="Set")
# windows(10,10)
# plot(Freq.MCD64A1.set.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.set.Cerrado)
summary(Freq.MCD64A1.set.Cerrado)
#abline(v=11,col="red")

Freq.sd.MCD64A1.set.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+8)]],na.rm=T, aszero = T) 

#Para outubro
Freq.MCD64A1.out.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+9)]],na.rm=T) 
Freq.MCD64A1.out.Cerrado<-crop(Freq.MCD64A1.out.Cerrado, extent(Cerrado))
Freq.MCD64A1.out.Cerrado <- mask(Freq.MCD64A1.out.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.out.Cerrado,main="Out")
# windows(10,10)
# plot(Freq.MCD64A1.out.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.out.Cerrado)
summary(Freq.MCD64A1.out.Cerrado)
#abline(v=9,col="red")

Freq.sd.MCD64A1.out.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+9)]],na.rm=T, aszero = T) 

#Para novembro
Freq.MCD64A1.nov.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+10)]],na.rm=T) 
Freq.MCD64A1.nov.Cerrado<-crop(Freq.MCD64A1.nov.Cerrado, extent(Cerrado))
Freq.MCD64A1.nov.Cerrado <- mask(Freq.MCD64A1.nov.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.nov.Cerrado,main="Nov")
# windows(10,10)
# plot(Freq.MCD64A1.nov.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.nov.Cerrado)
summary(Freq.MCD64A1.nov.Cerrado)
#abline(v=3,col="red")

Freq.sd.MCD64A1.nov.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+10)]],na.rm=T, aszero = T) 

#Para dezembro
Freq.MCD64A1.dez.Cerrado<-sum(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+11)]],na.rm=T) 
Freq.MCD64A1.dez.Cerrado<-crop(Freq.MCD64A1.dez.Cerrado, extent(Cerrado))
Freq.MCD64A1.dez.Cerrado <- mask(Freq.MCD64A1.dez.Cerrado, Cerrado)
windows(10,10)
plot(Freq.MCD64A1.dez.Cerrado,main="Dez")
# windows(10,10)
# plot(Freq.MCD64A1.dez.Cerrado/20)
# windows(10,10)
# hist(Freq.MCD64A1.dez.Cerrado)
summary(Freq.MCD64A1.dez.Cerrado)
#abline(v=1,col="red")

Freq.sd.MCD64A1.dez.Cerrado<-cv(z1[[(c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,169,181,193,205,217,229)+11)]],na.rm=T, aszero = T) 


windows(10,10)
par(mfrow=c(3,4))
plot(stack((Freq.MCD64A1.jan.Cerrado/20),(Freq.MCD64A1.fev.Cerrado/20),(Freq.MCD64A1.mar.Cerrado/20),
           (Freq.MCD64A1.abr.Cerrado/20),(Freq.MCD64A1.mai.Cerrado/20),(Freq.MCD64A1.jun.Cerrado/20),
           (Freq.MCD64A1.jul.Cerrado/20),(Freq.MCD64A1.ago.Cerrado/20),(Freq.MCD64A1.set.Cerrado/20),
           (Freq.MCD64A1.out.Cerrado/20),(Freq.MCD64A1.nov.Cerrado/20),(Freq.MCD64A1.dez.Cerrado/20)))

windows(10,10)
par(mfrow=c(3,4))
plot(stack((Freq.MCD64A1.jan.Cerrado),(Freq.MCD64A1.fev.Cerrado),(Freq.MCD64A1.mar.Cerrado),
           (Freq.MCD64A1.abr.Cerrado),(Freq.MCD64A1.mai.Cerrado),(Freq.MCD64A1.jun.Cerrado),
           (Freq.MCD64A1.jul.Cerrado),(Freq.MCD64A1.ago.Cerrado),(Freq.MCD64A1.set.Cerrado),
           (Freq.MCD64A1.out.Cerrado),(Freq.MCD64A1.nov.Cerrado),(Freq.MCD64A1.dez.Cerrado)))

Freq.fire.soma.MCD64A1.Cerrado<-stack((Freq.MCD64A1.jan.Cerrado),(Freq.MCD64A1.fev.Cerrado),(Freq.MCD64A1.mar.Cerrado),
                              (Freq.MCD64A1.abr.Cerrado),(Freq.MCD64A1.mai.Cerrado),(Freq.MCD64A1.jun.Cerrado),
                              (Freq.MCD64A1.jul.Cerrado),(Freq.MCD64A1.ago.Cerrado),(Freq.MCD64A1.set.Cerrado),
                              (Freq.MCD64A1.out.Cerrado),(Freq.MCD64A1.nov.Cerrado),(Freq.MCD64A1.dez.Cerrado))

writeRaster(Freq.fire.soma.MCD64A1.Cerrado,"Freq.fire.soma.MCD64A1.Cerrado.grd")

Freq.sd.fire.MCD64A1.Cerrado<-stack(Freq.sd.MCD64A1.jan.Cerrado,
                                      Freq.sd.MCD64A1.fev.Cerrado,
                                      Freq.sd.MCD64A1.mar.Cerrado,
                                      Freq.sd.MCD64A1.abr.Cerrado,
                                      Freq.sd.MCD64A1.mai.Cerrado,
                                      Freq.sd.MCD64A1.jun.Cerrado,
                                      Freq.sd.MCD64A1.jul.Cerrado,
                                      Freq.sd.MCD64A1.ago.Cerrado,
                                      Freq.sd.MCD64A1.set.Cerrado,
                                      Freq.sd.MCD64A1.out.Cerrado,
                                      Freq.sd.MCD64A1.nov.Cerrado,
                                      Freq.sd.MCD64A1.dez.Cerrado)

plot(Freq.sd.fire.MCD64A1.Cerrado)
plot(density(Freq.sd.fire.MCD64A1.Cerrado))

writeRaster(Freq.sd.fire.MCD64A1.Cerrado, "Freq.sd.fire.MCD64A1.Cerrado.grd")

Freq.fire.MCD64A1.chuva<-Freq.fire.soma.MCD64A1.Cerrado[[1]]+Freq.fire.soma.MCD64A1.Cerrado[[2]]+Freq.fire.soma.MCD64A1.Cerrado[[11]]+Freq.fire.soma.MCD64A1.Cerrado[[12]]
Freq.fire.MCD64A1.precoce<-Freq.fire.soma.MCD64A1.Cerrado[[3]]+Freq.fire.soma.MCD64A1.Cerrado[[4]]+Freq.fire.soma.MCD64A1.Cerrado[[5]]
Freq.fire.MCD64A1.modal<-Freq.fire.soma.MCD64A1.Cerrado[[6]]+Freq.fire.soma.MCD64A1.Cerrado[[7]]
Freq.fire.MCD64A1.tardia<-Freq.fire.soma.MCD64A1.Cerrado[[8]]+Freq.fire.soma.MCD64A1.Cerrado[[9]]+Freq.fire.soma.MCD64A1.Cerrado[[10]]

writeRaster(Freq.fire.MCD64A1.chuva,"Freq.fire.MCD64A1.chuva.grd")
writeRaster(Freq.fire.MCD64A1.precoce,"Freq.fire.MCD64A1.precoce.grd")
writeRaster(Freq.fire.MCD64A1.modal,"Freq.fire.MCD64A1.modal.grd")
writeRaster(Freq.fire.MCD64A1.tardia,"Freq.fire.MCD64A1.tardia.grd")

#################
#Mudar resolução#
#################
setwd("E:/Heitor/BD_SIG/MODIS_Fire/MCD64A1")
beginCluster(4)
rasterOptions(tmpdir=file.path("D:/Documentos/Raster"),
              progress = "text",
              todisk = T) 
z1.agreg<-raster::aggregate(z1,fact=2,fun=sum,na.rm=T)#Muda para 1km de res e soma os pixels queimados
z1.agreg
writeRaster(z1.agreg,"Bbin.all.Cerrado.MCD64A1.1km.grd",overwrite=T)

z1.agreg.2.5km<-raster::aggregate(z1,fact=5,fun=sum,na.rm=T)#Muda para 2.5km de res e soma os pixels queimados
z1.agreg.2.5km
writeRaster(z1.agreg.2.5km,"Bbin.all.Cerrado.MCD64A1.2.5km.grd",overwrite=T)

z1.agreg.9km<-raster::aggregate(z1,fact=18,fun=sum,na.rm=T)#Muda para 9km de res e soma os pixels queimados
z1.agreg.9km
plot(z1.agreg.9km[[1]])
writeRaster(z1.agreg.9km,"Bbin.all.Cerrado.MCD64A1.9km.grd",overwrite=T)

endCluster()
