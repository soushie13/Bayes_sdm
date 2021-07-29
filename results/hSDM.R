packages=c("hSDM","ggplot2","rasterVis","maptools","maps","dplyr","coda","doParallel","knitr","markdown","rgdal")
library(hSDM)
library(ggplot2)
library(rasterVis)
library(maptools)
library(maps)
library(dplyr)
library(coda)
library(doParallel)
library(knitr)
library(markdown)
library(raster)
library(rgeos)
library(rgdal)
## open xx occurrence data file##
getwd()
setwd("C:/")
dat <- read.csv(file = "xx.csv", header = TRUE)
obs.data <- dat[, c("x", "y")]
plot(obs.data)

##create the spatial extents from the reservoir extents for xx##
data(wrld_simpl)
res<- readOGR("REServoir.shp")
e <- extent(res)
rgeos::set_RGEOS_CheckValidity(2L)
res <- crop(wrld_simpl, e)
## add spatial buffers to the occurance data###
x1 = circles(obs.data, d=10000, lonlat=T)
x2 = circles(obs.data, d=50000, lonlat=T)
x3 = circles(obs.data, d=100000, lonlat=T)

x5 <- polygons(x1)
x10 <- polygons(x2)
x15 <- polygons(x3)

##generation of presence points at spatial buffer of 10km##
p = spsample(x5, 500, type='random', iter=1000)
plot(p)

##generation of pseudo absence points at spatial buffer of distribution of the reservoir##

# Randomly sample points (same number as our observed points)

background <- spsample(res,n= 1000,"random", iter= 5000) 
plot(background)

##climatic and elevation covariates##
tmin <- getData(name = "worldclim",var = "tmin", res = 2.5, path = "tmin/")
tmin_m <- mean(tmin)
tmax <- getData(name = "worldclim",var = "tmax", res = 2.5, path = "tmax/")
tmax_m <- mean(tmax)
ppt <- getData(name = "worldclim",var = "prec", res = 2.5, path = "ppt/")
ppt_m <- mean(ppt)
alt <- getData(name = "worldclim",var = "alt", res = 2.5, path = "")
## land and population raster
lc_raster <- raster("lc.tif")
topo_raster <- raster("topo.tif")
lf_raster <- raster("lc_modi.tif")
pop_raster <-raster("pop_den.tif")

### crop environmental predictors to the extent
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

## Extract environmental values and cell number for observations


pa$Presences <- pa$pb
pa$Trials <- c(10)
## omit rows with missing data (primarily ocean pixels)
pa.cc=na.omit(pa)


## Normalized continuous covariates
pa.norm <- pa.cc
Mean <- vector()
Sd <- vector()
for (i in c(5:11)) {
  m <- mean(pa.cc[,i],na.rm=TRUE)
  s <- sd(pa.cc[,i],na.rm=TRUE)
  Mean <- c(Mean,m)
  Sd <- c(Sd,s)
pa.norm[,i] <- (pa.cc[,i]-m)/s
}
## Data-frame with mean and sd for each variable
df.mean.sd <- as.data.frame(rbind(Mean,Sd))
names(df.mean.sd) <- names(pa.norm)[c(5:11)]

## Raster stack for predictions (with normalized covariates)
env <- bioclim.data
for (i in c(5:11)) {
  var.name <- names(pa.norm)[i] ## Variable name
  w <- which(names(env)==var.name) ## Position in the stack 
  m <- df.mean.sd[1,var.name] ## Mean
  s <- df.mean.sd[2,var.name] ## Sd
  orig <- values(subset(env,w)) ## Original values
  trans <- (orig-m)/s ## Transformed values
  env[[w]][] <- trans
}

## Select only grid cells with no NA
env.df.pred <- as.matrix(env)
w <- complete.cases(env.df.pred) ## Note: w will be used to obtain the cell identifier for predictions in iCAR model
env.df.pred.complete <- as.data.frame(env.df.pred[w,])

## Make a cluster for parallel MCMCs
nchains <- 2
ncores <- nchains ## One core for each MCMC chains
cores<-detectcores()
clust <- makeCluster(ncores)
registerDoParallel(clust)

## Starting values and random seed
seed <- 1234
set.seed(seed)
beta.start <- runif(nchains,-1,1)
gamma.start <- runif(nchains,-1,1)
Vrho.start <- runif(nchains,0,10)
seed.mcmc <- round(runif(nchains,0,1e6))
pa.norm$Trials <- c(1)
pa.norm$Presences <- pa.norm$pb
##===============================================
##
## 1. Model with environmental variables
##
##===============================================

## layer.1: mean minimum temperature
## layer.2: mean maximum temperature
## layer: mean precipitation 
## alt: elevation above sea level
## pop: population density
## lc: land cover
## lc_modi: land cover from modifications in human activity


## hSDM model using Zero Inflated Binomial (ZIB)
mod.ZIB.env <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.ZIB(presences=pa.norm$Presences,
                  trials=pa.norm$Trials,
                  suitability=~layer.1+layer.2+layer+lc+alt+lc_modi+pop_den,     
                  observability=~1,
                  data=pa.norm,
                  suitability.pred=env.df.pred.complete,
                  burnin=10000,
                  mcmc=10000, thin=10,
                  beta.start=beta.start[i],
                  gamma.start=gamma.start[i],
                  mubeta=0, Vbeta=1.0E6,
                  mugamma=0, Vgamma=1.0E6,
                  seed=seed.mcmc[i], verbose=1,
                  save.p=0)
  return(mod)
}
## Extract list of MCMCs from output
ZIB.env.mcmc <- mcmc.list(lapply(mod.ZIB.env,"[[","mcmc"))

## Outputs summary
ZIB.env.stat <- summary(ZIB.env.mcmc)$statistics
sink(file="results/ZIB.env.xx_mcmc_summary.txt")
summary(ZIB.env.mcmc)
cat(rep("\n",3))
gelman.diag(ZIB.env.mcmc)
#gelman.plot(ZIB.env.mcmc)
sink()
## Deviance
deviance.ZIB.env <- ZIB.env.stat["Deviance","Mean"]

## Detection probability
gamma.hat <- ZIB.env.stat["gamma.(Intercept)","Mean"]
delta.est <- inv.logit(gamma.hat) 
delta.est
## Plot trace and posterior distributions
pdf("results/ZIB.env.xx_mcmc_trace.pdf")
plot(ZIB.env.mcmc)
dev.off()
## Prediction on the landscape
prob.p.z <- subset(bioclim.data,1) ## create a raster for predictions
values(prob.p.z)[w] <- mod.ZIB.env[[1]]$prob.p.pred ## assign predicted values
values(prob.p.z)[!w] <- NA ## set NA where no environmental data
## Plot the predictions
pdf(file="results/ZIB.env_xx_predictions.pdf")
plot(prob.p.z)
#plot(pa.norm[pa.norm$pb==0,],pch=".",col=grey(0.5),add=TRUE)
#plot(pa.norm$pb,pch=3,add=TRUE)
points(obs.data$x,obs.data$y, pch=3, add=TRUE)
dev.off()
## Export the results as GeoTIFF
writeRaster(prob.p,filename="results/ZIB_env_pred_mod_xx.tif",overwrite=TRUE)

##===============================================
##
## 2. Binomial model  
##===============================================


## binomial model
## hSDM model using Binomial for perfect detection
mod.binomial <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.binomial(presences=pa.norm$Presences,
                  trials=pa.norm$Trials,
                  suitability=~layer.1+layer.2+layer+lc+alt+lc_modi+pop_den,
                  data=pa.norm,
                  suitability.pred=env.df.pred.complete,
                  burnin=5000,
                  mcmc=5000, thin=5,
                  beta.start=beta.start[i],
                  mubeta=0, Vbeta=1.0E6,
                  seed=seed.mcmc[i], verbose=1,
                  save.p=0)
  return(mod)
}

## Extract list of MCMCs from output
binomial.env.mcmc <- mcmc.list(lapply(mod.binomial,"[[","mcmc"))

sink(file="results/co_mcmc_summary.txt")
summary(binomial.env.mcmc)
sink()
## Outputs summary
bionomial.env.stat <- summary(binomial.env.mcmc)$statistics
sink(file="results/binomial_xx_mcmc_summary.txt")
summary(binomial.env.mcmc)
cat(rep("\n",3))
gelman.diag(binomial.env.mcmc)
sink()
## Deviance
deviance.bionomial.env <- bionomial.env.stat["Deviance","Mean"]

## Plot trace and posterior distributions
pdf("results/bionomial.env_xx_mcmc_trace.pdf")
plot(binomial.env.mcmc)
dev.off()

## Prediction on the landscape
prob.p.bi <- subset(bioclim.data,1) ## create a raster for predictions
values(prob.p.bi)[w] <- mod.binomial[[1]]$theta.pred ## assign predicted values
values(prob.p.bi)[!w] <- NA ## set NA where no environmental data
## Plot the predictions
pdf(file="results/binomial.env_xx_predictions.pdf")
plot(prob.p.bi)
plot(pa.norm[pa.norm$pb==0,],pch=".",col=grey(0.5),add=TRUE)
plot(pa.norm[pa.norm$pb>0,],pch=3,add=TRUE)
dev.off()

## Export the results as GeoTIFF
writeRaster(prob.p,filename="results/binomial_pred_mod_xx.tif",overwrite=TRUE)

##===============================================
##
## 3. Binomial iCAR model  
##===============================================


## Landscape and neighbors
ncells <- ncell(bioclim.data)
neighbors.mat <- adjacent(bioclim.data, cells=c(1:ncells), directions=8, pairs=TRUE, sorted=TRUE)
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
adj <- neighbors.mat[,2]
cells.pred <- which(w) ## Vector w indicates the cells with environmental information (without NA)

## binomial icar model
## hSDM model using Binomial icar for perfect detection
mod.binomial.icar <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.binomial.iCAR(presences=pa.norm$Presences,
                       trials=pa.norm$Trials,
                       suitability=~layer.1+layer.2+layer+lc+alt+lc_modi+pop_den,
                       data=pa.norm,
                       ## Spatial structure
                       spatial.entity=pa.norm$cells,
                       n.neighbors=n.neighbors,
                       neighbors=adj,
                       suitability.pred=env.df.pred.complete,
                       spatial.entity.pred=cells.pred,
                       burnin=5000,
                       mcmc=5000, thin=5,
                       beta.start=beta.start[i],
                       Vrho.start=Vrho.start[i],
                       ## Priors
                       priorVrho="Uniform",
                       mubeta=0, Vbeta=1.0E6,
                       Vrho.max=10,
                       seed=seed.mcmc[i], verbose=1,
                       save.p=0)
  return(mod)
}

## Extract list of MCMCs from output
binomial.icar.mcmc <- mcmc.list(lapply(mod.binomial.icar,"[[","mcmc"))
sink(file="results/binomial_icar_co_mcmc_summary.txt")
summary(binomial.icar.mcmc)
sink()
## Outputs summary
bionomial.icar.stat <- summary(binomial.icar.mcmc)$statistics
sink(file="results/binomial.icar_xx_mcmc_summary.txt")
summary(binomial.icar.mcmc)
cat(rep("\n",3))
gelman.diag(binomial.icar.mcmc)
sink()
## Deviance
deviance.bionomial.icar <- bionomial.icar.stat["Deviance","Mean"]

## Plot trace and posterior distributions
pdf("results/bionomial.icar_xx_mcmc_trace.pdf")
plot(binomial.icar.mcmc)
dev.off()
## Spatial random effects
rho <- subset(bioclim.data,1) ## create a raster
values(rho) <- mod.binomial.icar[[1]]$rho.pred
pdf(file="results/xx_binomial.iCAR_random_effects.pdf")
plot(rho)
dev.off()
## Prediction on the landscape
prob.p.b <- subset(bioclim.data,1) ## create a raster for predictions
values(prob.p.b)[w] <- mod.binomial.icar[[1]]$theta.pred ## assign predicted values
values(prob.p.b)[!w] <- NA ## set NA where no environmental data
## Plot the predictions
pdf(file="results/binomial.icar_xx_predictions.pdf")
plot(prob.p.b)
plot(pa.norm[pa.norm$pb==0,],pch=".",col=grey(0.5),add=TRUE)
plot(pa.norm[pa.norm$pb>0,],pch=3,add=TRUE)
dev.off()

## Export the results as GeoTIFF
writeRaster(prob.p.b,filename="results/binomial_icar_pred_mod_xx.tif",overwrite=TRUE)

##===============================================
##
## 4. ZIB iCAR model  
##===============================================

## ZIB.iCAR model
mod.ZIB.iCAR <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.ZIB.iCAR(
    ## Observations
    presences=pa.norm$Presences,
    trials=pa.norm$Trials,
    ## Habitat
    suitability=~layer.1+layer.2+layer+lc+alt+lc_modi+pop_den,
    observability=~1,
    ## Data-set
    data=pa.norm,
    ## Spatial structure
    spatial.entity=pa.norm$cells,
    n.neighbors=n.neighbors,
    neighbors=adj,
    ## Predictions
    suitability.pred=env.df.pred.complete,
    spatial.entity.pred=cells.pred,
    ## Chains
    burnin=5000, mcmc=5000, thin=5,
    ## Starting values
    beta.start=beta.start[i],
    gamma.start=gamma.start[i],
    Vrho.start=Vrho.start[i],
    ## Priors
    priorVrho="Uniform",
    #priorVrho=10,
    mubeta=0, Vbeta=1.0E6,
    mugamma=0, Vgamma=1.0E6,
    Vrho.max=10,
    ## Various
    seed=seed.mcmc[i], verbose=1,
    save.rho=0, save.p=0) ## Set save.p=1 to save predictive posterior for each spatial cell
  return(mod)
}

## Extract list of MCMCs from output
ZIB.iCAR.mcmc <- mcmc.list(lapply(mod.ZIB.iCAR,"[[","mcmc"))
sink(file="results/zib_icar_co_mcmc_summary.txt")
summary(ZIB.iCAR.mcmc)
sink()
## Outputs summary
ZIB.iCAR.stat <- summary(ZIB.iCAR.mcmc)$statistics

sink(file="results/ZIB_iCAR_xx_mcmc_summary.txt")
summary(ZIB.iCAR.mcmc)
cat(rep("\n",3))
gelman.diag(ZIB.iCAR.mcmc)
sink()
## Deviance
deviance.ZIB.iCAR <- ZIB.iCAR.stat["Deviance","Mean"]

## Plot trace and posterior distributions
pdf("results/ZIB_iCAR_xx_mcmc_trace.pdf")
plot(ZIB.iCAR.mcmc)
dev.off()

## Spatial random effects
rho <- subset(bioclim.data,1) ## create a raster
values(rho) <- mod.ZIB.iCAR[[1]]$rho.pred
pdf(file="results/ZIB_iCAR_xx_random_effects.pdf")
plot(rho)
dev.off()

## Prediction on the landscape
summary(mod.ZIB.iCAR[[1]]$prob.p.pred)
prob.p <- subset(bioclim.data,1) ## create a raster for predictions
values(prob.p)[w] <- mod.ZIB.iCAR[[1]]$prob.p.pred ## assign predicted values
values(prob.p)[!w] <- NA ## set NA where no environmental data
## Plot the predictions
pdf(file="results/ZIB.iCAR_xx_predictions.pdf")
plot(prob.p,zlim=c(0,1))
dev.off()
#= Summary plots


## Export the results as GeoTIFF
writeRaster(prob.p,filename="results/xx_pred_mod_ZIB_iCAR.tif",overwrite=TRUE)

##===============================================
##
## Model comparison based on deviance
##
##===============================================

## Null model
mod.ZIB.null <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.ZIB(presences=pa.norm$Presences,
                  trials=pa.norm$Trials,
                  suitability=~1,
                  observability=~1,
                  data=pa.norm,
                  suitability.pred=env.df.pred.complete,
                  burnin=5000,
                  mcmc=5000, thin=5,
                  beta.start=beta.start[i],
                  gamma.start=gamma.start[i],
                  mubeta=0, Vbeta=1.0E6,
                  mugamma=0, Vgamma=1.0E6,
                  seed=seed.mcmc[i], verbose=1,
                  save.p=0)
  return(mod)
}
## Stop cluster
stopCluster(clust)

## Extract list of MCMCs from output
ZIB.null.mcmc <- mcmc.list(lapply(mod.ZIB.null,"[[","mcmc"))

## Deviance
ZIB.null.stat <- summary(ZIB.null.mcmc)$statistics
deviance.null <- ZIB.null.stat["Deviance","Mean"]


prob.p.null <- subset(bioclim.data,1) ## create a raster for predictions
values(prob.p.null)[w] <- mod.ZIB.null[[1]]$prob.p.pred ## assign predicted values
values(prob.p.null)[!w] <- NA 

##= Table of deviance
dev.tab <- data.frame(Model=rep(NA,5),Deviance=rep(0,5),Perc=rep(0,5))
dev.tab$Model <- c("NULL","env","binomial","binomial.icar", "ZIB.icar")
dev.tab$Deviance <- c(deviance.null, deviance.ZIB.env,deviance.bionomial.env, deviance.bionomial.icar, deviance.ZIB.iCAR)
dev.tab$Perc <- round(100*(dev.tab$Deviance[1]-dev.tab$Deviance)/(dev.tab$Deviance[1]-dev.tab$Deviance[5]))
##= Export
sink(file="results/xx_deviance.txt")
dev.tab
sink()

##================================================================
##
## TSS (True Skill Statistics) and SDA (Species Distribution Area)
##
##================================================================

## Function to compute threshold dependent indexes
Index.fun <- function(obs,prob.p,thresh) {
  ## Transform probabilities into {0,1} given threshold
  pred <- ifelse(prob.p>=thresh,1,0)
  ## contingency table (pred/obs)
  n00 <- sum(pred==0 & obs==0,na.rm=TRUE)
  n11 <- sum(pred==1 & obs==1,na.rm=TRUE)
  n01 <- sum(pred==0 & obs==1,na.rm=TRUE)
  n10 <- sum(pred==1 & obs==0,na.rm=TRUE)
  ## Threshold  dependent indexes
  OA <- (n11+n00)/(n11+n10+n00+n01) ## Overall accuracy
  Sensitivity <- n11/(n11+n01)
  Specificity <- n00/(n00+n10)
  TSS <- Sensitivity+Specificity-1
  return(list(OA=OA,TSS=TSS,Sens=Sensitivity,Spe=Specificity))
}

## Suitable sites
pa$Suit <- 0
pa$Suit[pa$pb>0] <- 1
pa.obs.data <- pa[, c("x", "y")]
## Extract predicted probability of presence
pa$prob.p <- extract(prob.p,pa.obs.data)

## TSS as a function of the threshold
OA <- vector()
TSS <- vector()
thresh.seq <- seq(0,1,by=0.01)
for (i in 1:length(thresh.seq)) {
  Index <- Index.fun(pa$Suit,pa$prob.p,thresh.seq[i])
  OA[i] <- Index$OA
  TSS[i] <- Index$TSS
}

## Identifying the threshold maximizing TSS
maxTSS <- max(TSS,na.rm=TRUE)
maxTSS
w.th <- which(TSS==maxTSS)
thresh.maxTSS <- mean(thresh.seq[w.th],na.rm=TRUE)
OA.thresh.maxTSS <- Index.fun(pa$Suit,pa$prob.p,thresh.maxTSS)$OA 
tss.df <- data.frame(masTSS=maxTSS,OA=OA.thresh.maxTSS,prob=thresh.maxTSS)

## Plot evolution of TSS with threshold
pdf(file="results/xx_TSS.pdf")
plot(thresh.seq,TSS,type="l",xlab=c("probability threshold"),ylab="TSS")
abline(v=thresh.maxTSS)
dev.off()

## SDA based on maxTSS
SDA <- prob.p
SDA[SDA>=thresh.maxTSS] <- 1
SDA[SDA<thresh.maxTSS] <- 0
pdf(file="results/xx_SDA.pdf")
plot(SDA,legend=FALSE)
dev.off()

## Export the result
writeRaster(SDA,filename="results/xx_SDA.tif",overwrite=TRUE)
write.table(round(tss.df,2),file="results/xx_TSS.txt",row.names=FALSE,sep="\t")

## Estimating SDA area (in km2)
n.pix <- sum(values(SDA),na.rm=TRUE)
area.SDA <- n.pix*res(SDA)[1]*res(SDA)[2]/1.0e+6
area.SDA
predscale=scale_fill_gradientn(values=c(0,.5,1),colours=c('white','darkgreen','green'),na.value="transparent")
gplot(prob.p)+geom_raster(aes(fill=value)) +
  predscale+
  coord_equal()+
  ggcoast+gx+gy+ylab("Latitude")+xlab("Longitude")+
  labs(col = "p(presence)")+
  coord_equal()
ggsave("results/ggpt_co.pdf")


save(mod.binomial.icar,pa.norm, neighbors.mat, file= "results/binomial_coronaviridae.RData")
save(mod.ZIB.iCAR,pa.norm, neighbors.mat, file= "results/ZIB_coronaviridae.RData")
###combining models by weighted mean##
models <- stack(prob.p.null,prob.p.z,prob.p.bi,prob.p.b, prob.p)
names(models) <- c("null","ZIB","Binomial","Binomial iCAR", "ZIB.iCAR")
plot(models)
m <- mean(models)
plot(m, main='Mean')
dev <- sapply(list(ZIB.env.stat,bionomial.env.stat, bionomial.icar.stat, ZIB.env.stat), function(x) x["Deviance"])
weight <- dev.tab$Perc
m2 <- weighted.mean( models, weight)
plot(m2)
## Export the results as GeoTIFF
writeRaster(m2,filename="results/xx_weighted_model.tif",overwrite=TRUE)


##
## Final SDA map with Google Map background (with ggmap::get_map() function)
##
##==========================================================================

## Reproject in Lat/Long (epsg:4326)
GM.crs <- CRS("+init=epsg:4326")
## SDA
SDA.GM <- projectRaster(from=SDA,crs=GM.crs,method="ngb")
SDA.df <- as.data.frame(SDA.GM,xy=TRUE,na.rm=TRUE)
loc <- extent(SDA.GM)[c(1,3,2,4)]
names(SDA.df) <- c("x","y","pres")
sda.df <- SDA.df[SDA.df$pres==1,]
library(ggmap)
## Plot with ggplot2
bg <- ggmap(get_map(location=loc,zoom=4,maptype="terrain",source="stamen"))
g.sda <- geom_tile(mapping=aes(x,y,fill=pres),data=sda.df) 
 
gg.plot <- bg + 
  g.sda +  scale_fill_gradientn(colours = c("red"))## + g.abs ## If we want the absences to be plot

#+ coord_cartesian()  coord_equal() 
## Save as png image file
ggsave(filename="xx_SDA_ggmap.png",plot=gg.plot,device="png",path="results/",width=10,height=7,units="cm",dpi=300)
ggsave("results/xx_sda_pred.pdf")
plot(bg)
plot(SDA.GM, add = T, legend=FALSE, color="red")

t = tm_shape(prob.p) + 
  tm_raster(style = "fisher", title = "Hotspots for xx",
            palette = "YlOrRd",
            legend.hist = FALSE)+
  tm_legend(outside = TRUE)+
    tm_shape(res)+ tm_polygons(col="grey", border.col=NULL, alpha = 0.3)+
  tm_compass(type="arrow", position= c("LEFT", "BOTTOM"))+
    tm_scale_bar(position =  c("LEFT", "BOTTOM"))

  tmap_save(t, "results/Hotspots for xx.pdf") # height interpreted in inches
  