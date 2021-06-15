library(sp)
library(devtools)
library(tmap)
library(raster)
library(rgeos)
library(rgdal)
library(mapdata)
library(climate)
library(maptools)
library(dismo)
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
library(tmap)
library("spatstat")
devtools::install_github("ecohealthalliance/fasterize")
library(fasterize)
## open marburg occurance data file##
getwd()
setwd("C:/Users/sj20e051/Documents/biogeography/ecohealth/")
dat <- read.csv(file = "C:/Users/sj20e051/Documents/biogeography/ecohealth/cov_co.csv", header = TRUE)
obs.data <- dat[, c("x", "y")]
plot(obs.data)
## computational grid from the distribution of marburg reservior polygons##
data(wrld_simpl)
res<- readOGR("C:/Users/sj20e051/Documents/biogeography/ecohealth/mammal_cov/TERRESTRIAL_MAMMALS.shp")
e <- extent(res)
rgeos::set_RGEOS_CheckValidity(2L)
res <- crop(wrld_simpl, e)

CP <- as(extent(res), "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(wrld_simpl))

## Clip the map
out <- gIntersection(wrld_simpl, CP, byid=TRUE)


plot(res)

# ratify raster
r <- ratify(ras)

# Create levels
rat <- levels(r)[[1]]
rat$names <- nam_df$nam
rat$IDs <- nam_df$ID
levels(r) <- rat


bbox(res_nip)
e <- extent(res)
plot(ras)
plot(wrld_simpl,  xlim = c(-10, 130),
     ylim = c(-30, 50),
     axes = TRUE)
res_cov <- aggregate(res_cov,dissolve=T)
rgeos::set_RGEOS_CheckValidity(2L)
## Crop to the desired extent, then plot
#-10,130,-30, 47
extend(e, 10)
res_cov<- crop(wrld_simpl, extent(res_cov))
plot(res_cov)
dev.off()
points(obs.data$x, obs.data$y)
proj4string(res_nip)
resolution <- 0.2
r <- raster(res_mvd, resolution = resolution)
(nrow <- nrow(r))
(ncol <- ncol(r))
nrow*ncol
r[] <- 0
tab <- table(cellFromXY(r, p))
tab
r[as.numeric(names(tab))] <- tab
grid <- rasterToPolygons(r)

coordinates(obs.data)=c("x","y")
projection(obs.data)="+proj=longlat +datum=WGS84 +ellps=WGS84"
obs.data@data[,c("x","y")]=coordinates(obs.data)  
## add spatial buffers to the occurance data###
x1 = circles(obs.data, d=10000, lonlat=T)
x2 = circles(obs.data, d=50000, lonlat=T)
x3 = circles(obs.data, d=100000, lonlat=T)

x5 <- polygons(x1)
x10 <- polygons(x2)
x15 <- polygons(x3)

##generation of presence points at spatial buffer of 100km##
p = spsample(x5, 500, type='random', iter=1000)
plot(p)

                              
##generation of pseudo absence points at spatial buffer of distribution of the reservoir##

# Randomly sample points (same number as our observed points)

background <- spsample(res,n= 1000,"random", iter= 5000) 
plot(background)


res(tmin_m)
res(lc_raster)
res(alt)
res(lf_e)
res(pop_e)
lf_e <- resample(lf_raster, lc_raster, method = 'bilinear')
pop_e <- resample(pop_raster, lc_raster, method = 'bilinear')
##climatic covariates##
tmin <- getData(name = "worldclim",var = "tmin", res = 2.5, path = "C:/Users/sj20e051/Documents/biogeography/ecohealth/tmin/")
tmin_m <- mean(tmin)
tmax <- getData(name = "worldclim",var = "tmax", res = 2.5, path = "C:/Users/sj20e051/Documents/biogeography/ecohealth/tmax/")
tmax_m <- mean(tmax)
ppt <- getData(name = "worldclim",var = "prec", res = 2.5, path = "C:/Users/sj20e051/Documents/biogeography/ecohealth/ppt/")
ppt_m <- mean(ppt)
alt <- getData(name = "worldclim",var = "alt", res = 2.5, path = "C:/Users/sj20e051/Documents/biogeography/ecohealth/")
plot(ppt_m)

## For multiple files, could use a for loop
## Input directory
library(ncdf4)
setwd('C:/Users/sj20e051/Documents/biogeography/ecohealth/diseasex/')
files <- list.files(pattern='*.nc', full.names=TRUE)
s <- stack(files)
s_tmin<- writeRaster(s, filename="tmin_s", bylayer=TRUE, format="GTiff")


rastlist<- list.files(path = "C:/Users/sj20e051/Documents/biogeography/ecohealth/tmin/", pattern='*.tif$', all.files=TRUE, full.names=FALSE)
tmin <- raster::stack(paste0("C:/Users/sj20e051/Documents/biogeography/ecohealth/tmin/", rastlist))

setwd('C:/Users/sj20e051/Documents/biogeography/ecohealth/tmax/')
files <- list.files(pattern='*.nc', full.names=TRUE)
s <- stack(files)
s_tmax<- writeRaster(s, filename="tmax_s", bylayer=TRUE, format="GTiff")


rastlist<- list.files(path = "C:/Users/sj20e051/Documents/biogeography/ecohealth/tmax/", pattern='*.tif$', all.files=TRUE, full.names=FALSE)
tmax <- raster::stack(paste0("C:/Users/sj20e051/Documents/biogeography/ecohealth/tmax/", rastlist))

setwd('C:/Users/sj20e051/Documents/biogeography/ecohealth/ppt/')
files <- list.files(pattern='*.nc', full.names=TRUE)
s <- stack(files)
s_ppt<- writeRaster(s, filename="ppt_s", bylayer=TRUE, format="GTiff")


rastlist<- list.files(path = "C:/Users/sj20e051/Documents/biogeography/ecohealth/ppt/", pattern='*.tif$', all.files=TRUE, full.names=FALSE)
ppt<- raster::stack(paste0("C:/Users/sj20e051/Documents/biogeography/ecohealth/ppt/", rastlist))

lc_raster <- raster("C:/Users/sj20e051/Documents/biogeography/ecohealth/lc.tif")
topo_raster <- raster("C:/Users/sj20e051/Documents/biogeography/ecohealth/topo.tif")
lf_raster <- raster("C:/Users/sj20e051/Documents/biogeography/ecohealth/lc_modi.tif")

file.nc <- "gpw_v4_population_density_rev11_2pt5_min.nc"
importnetcdf <- raster(file.nc)

writeRaster(importnetcdf,filename="pop_den.tiff",format="GTiff",overwrite=TRUE,bylayer=TRUE)
pop_raster <-raster("C:/Users/sj20e051/Documents/biogeography/ecohealth/pop_den.tif")
extent(clim)


e = extent(res)

plot(res)

### crop to the extent
lc <- crop(lc_raster, e) 
alt_e <- crop(alt, e) 
lf <- crop(lf_raster, e) 
pop <- crop(pop_raster, e)
tmin <- crop(tmin_m, e)
tmax <- crop(tmax_m,e)
ppt <- crop(ppt_m,e)
clim1 <- stack(tmin,tmax)
clim2 <- addLayer(clim1, ppt)
clim <- addLayer(clim2,alt_e)

r12_lc<- resample(lc, clim,method= 'bilinear')
r12 <- addLayer(clim, lc)
r12_pop <- resample(pop, r12,method= 'bilinear')
r123 <- addLayer(r12, r12_pop)
r123_lf <- resample(lf, r123,method= 'bilinear')
bioclim.data <- addLayer(r123, r123_lf)
res(bioclim.data)
plot(bioclim.data)
#extract data
presence = p@coords
presvals <- extract(bioclim.data, presence, cellnumber=TRUE)
prevals_coords =cbind(presence,presvals)
absence = background@coords
absvals <- extract(bioclim.data, absence, cellnumber=TRUE)
absvals_coords =cbind(absence,absvals)
pb <- c(rep(1, nrow(prevals_coords)), rep(0, nrow(absvals_coords)))
pa <- data.frame(cbind(pb, rbind(prevals_coords, absvals_coords)))


ncelltot <- length(env)
pa <- data.frame(cbind(pb, rbind(prevals_coords, absvals_coords)))
coords <- data.frame(cbind(p@coords,background@coords))
plot(bioclim.data)
View(gr)
### extract data from the covariates##
grid$id <- 1:nrow(grid)
grid$Y <- grid$layer
grid$cellarea <- resolution*resolution
grid$tmin <- extract(bioclim.data$layer.1, coordinates(grid))
grid$tmax <- extract(bioclim.data$layer.2, coordinates(grid))
grid$ppt <- extract(bioclim.data$layer, coordinates(grid))
grid$lc <- extract(bioclim.data$lc, coordinates(grid))
grid$pop <- extract(bioclim.data$pop, coordinates(grid))
grid$lf <- extract(bioclim.data$lc_modi, coordinates(grid))
grid$alt <- extract(bioclim.data$alt, coordinates(grid))
gridmap <- raster::intersect(grid, res_mvd)
grid <- grid[grid$id %in% gridmap$id, ]
summary(grid)
##replacing NA in the covariates##
indNA_tmin <- which(is.na(grid$tmin))
grid$tmin[indNA_tmin] <- grid$tmin[indNA_tmin+1]
indNA_tmax <- which(is.na(grid$tmax))
grid$tmax[indNA_tmax] <- grid$tmax[indNA_tmax+1]
indNA_ppt <- which(is.na(grid$ppt))
grid$ppt[indNA_ppt] <- grid$ppt[indNA_ppt+1]
indNA_lf <- which(is.na(grid$lf))
grid$lf[indNA_lf] <- grid$lf[indNA_lf+1]
grid@data[is.na(grid@data)] <- 0
library(rgeos)
gridborder <- gUnaryUnion(grid)
tmap_mode("plot")
tm_shape(grid) +
  tm_polygons(col = c("Y", "tmin", "topo", "pop"), border.col = "transparent") +
  tm_shape(gridborder) + tm_borders() +
  tm_facets(ncol = 2) + tm_legend(legend.position = c("left", "bottom"))
library(INLA)

grid$id2 <- grid$id
View(grid@data)
formula <- Y ~ 1 + alt+ lf+ lc+ pop+
  f(id, model="rw2d", nrow = nrow, ncol = ncol) +
  f(id2, model="iid")
#tmin+ tmax+ ppt+ lf+ lc+ pop+
res <- inla(formula, family = "poisson", data = grid@data,
            E = cellarea, control.predictor = list(compute = TRUE), verbose = TRUE)
# Pull together coordinates and PA data into SpatialPointsDataFrame
dataframe = sp::SpatialPointsDataFrame(coords = p, data = data.frame(y = PA))

# Run the model.
model <- inlaSDM(dataframe, 
                 predictors, 
                 spatial = TRUE, 
                 cross_validation = FALSE,
                 meshvals = list(cutoff = 0.3, inner.max.edge = 1))
