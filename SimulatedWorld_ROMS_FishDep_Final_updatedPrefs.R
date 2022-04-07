#Code updated 11/8/21 by Melissa Karp
#This code is the operating model for a HMS archetype representing a Pelagic mobile predator-like species
#It uses average spring conditions from downscales ROMS projections using Hadley 
#IMPORTANT: The ROMS-NEMUCSC projection data are available upon request from Mercedes Pozo Buil, or from NOAA’s ERDDAP data servers at the following URL: 
#https://oceanview.pfeg.noaa.gov/erddap/search/index.html?&searchFor=CCS+ROMS.


#We include a trophic interaction between predator (PMP) and prey, but note the estmation model will use chl-a as a proxy for prey information:
#prey: distribution and suitability driven by SST and zooplankton integrated across top 200m
#predator: distribution and abundance driven by SST, MLD, and Species A
#Simulates fishery dependent sampling bias by building a fishing location suitability raster using virtual species
#Biases Include:
# Random Sampling (control)
# Preferential Sampling - fishing suitability a function of target species habitat suitability in Time Y-1
# Distance from Port - fishing location suitability function of distance from closest port & target species habitat suitability 
# Bycatch avoidance - adds behavior of fishermen actively trying to avoid areas of high bycatch risk, while still targeting species of interest 
# Spatial closure - fishing locations limited by spatial closure 

#hadley directory
  dir <- "~/DisMAP project/Location, Location, Location/Location Workshop/ROMS/hadley_11_2_21" #this is the local drive where you have the ROMS data stored

#*******************Closed Area Simulation********************#
#Code courtesy of Owen Liu 
# Closed area designation can be a set of coordinates (a vector as c(xmin,xmax,ymin,ymax) in lat/lon), or a simple features (sf) polygon object
# The function takes the given area and converts to a raster on the same grid as the grid used for the species simulation
# The resulting raster is made up of zeroes (closed area) and ones, such that it could be multiplied by a habitat suitability or
# sampling raster, in order to exclude those areas from the sample
# default box is southern California, approximately where the current gillnet loggerhead closure is located
define_closed_area <- function(closed_area=c(-122,-117,31,35)){
  library(fasterize)
  library(raster)
  library(sf)
  library(dplyr)
  #template raster
  rst_files <- list.files(paste0(dir,'/sst_spring_avg'), full.names = TRUE, pattern=".grd")
  rst <- raster(rst_files[[1]])
  if(class(closed_area)=='numeric'){
    crds <- matrix(c(closed_area[1],closed_area[3],closed_area[1],closed_area[4],
                     closed_area[2],closed_area[4],closed_area[2],closed_area[3],
                     closed_area[1],closed_area[3]),ncol=2,byrow=T)
    pol <- st_polygon(list(crds)) %>% st_sfc(crs = 4326) %>% st_sf()
    rst_out <- fasterize(pol,rst)
  }
  # if(sf %in% class(closed_area)){
  #   rst_out <- closed_area %>% st_transform(4326) %>% fasterize(rst)
  # }
  rst_out <- rst_out==0
  rst_out[is.na(rst_out)] <- 1
  return(rst_out)
}
closed_area_1<-define_closed_area(closed_area=c(-126,-118,36,41))
closed_area_2<-define_closed_area(closed_area=c(-127,-118,35,42))
closed_area_3<-define_closed_area(closed_area=c(-129,-118,34,43))
#plot(closed_area_1*ROMS)#ROMS raster comes from the DistancetoPort.R code
#plot(closed_area_2*ROMS)
#plot(closed_area_3*ROMS)

###### *******Simulate Species Distribution, Abundance and Sampling********######

SimulateWorld_ROMS_FishDepFun_Final <- function(dir, nsamples){
  #dir is the local directory that points to where ROMS data is stored 
  #nsamples let's you choose the number of samples each year
  
  #required libraries
  library(virtualspecies)
  library(scales)
  library(raster)
  library(glmmfields)
  library(dplyr)
  library(splitstackshape)
  
  #----Create output file----
  #This will be the information passed to the estimation model
  #grid_had<-21912 #21912 non-NA grid cells in ROMS hadley
  #grid_gfdl<-22065
  output <- as.data.frame(matrix(NA, nrow=21912*121,ncol=25)) 
  colnames(output) <- c("lon","lat","year","pres","suitability_t", "spC_suitability", 
                        "spA_suitability", "random_sampled", "pref_sampled", 
                        "dist_sampled_npo", "dist_sampled_npn", "dist_sampled_mpo", "dist_sampled_mpn",
                        "dist_sampled_spo", "dist_sampled_spn", "dist_sampled_allo", "dist_sampled_alln",
                        "BY_sampled_2", "Closed_sampled_1", "Closed_sampled_2", "Closed_sampled_3", 
                        "sst","zoo_200", "mld", "chl_surface")
  
  #----Load in rasters----
  #These are the average spring conditions from the downscaled 'had' earth system model
  files_sst <- list.files(paste0(dir,'/sst_spring_avg'), full.names = TRUE, pattern=".grd") #should be 1452 files
  files_mld <- list.files(paste0(dir,'/ild_0.5C'), full.names = TRUE, pattern=".grd")
  files_zoo <- list.files(paste0(dir,'/zoo_200m'), full.names = TRUE, pattern=".grd")
  files_chl_surface <- list.files(paste0(dir,'/chl_surface'), full.names = TRUE, pattern=".grd")
  years <- seq(1980,2100,1)
  
  #---Load in Distance to port file -----
  dist_to_ports<-read.csv("~/DisMAP project/Location, Location, Location/Location Workshop/Dist_to_Ports_had.csv")
  
  #----Loop through each year----
  for (y in 1:121){
    print(paste0("Creating environmental simulation for Year ",years[y]))
    
    #Load in environmental rasters for a specific year
    sst <- raster(files_sst[y])
    mld <- raster(files_mld[y])
    zoo <- raster(files_zoo[y])
    chl_surface<-raster(files_chl_surface[y])
    chl_surface<-log(chl_surface)
    
    #Optional: plot environmental layers
    # par(mfrow=c(1,4))
    # plot(sst, main= 'SST')
    # plot(chl_surface, main = 'Chl-a Surface')
    # plot(mld, main = 'MLD')
    # plot(zoo, main = 'Zooplankton')
    
    #Environmental at time t-1 (previous year).
    #If year 1, then just use at time t.
    if (y>1){
      sst_t1 <- raster(files_sst[y-1])
      mld_t1<-raster(files_mld[y-1])
      zoo_t1<-raster(files_zoo[y-1])
      
    } else {
      sst_t1 <- sst
      mld_t1 <- mld
      zoo_t1 <- zoo
    }
    
    #---- ASSIGN RESPONSE CURVES ---- ####
    
    #----SPECIES A (prey): assign response curves----
    #Species A: likes high zooplankton and medium temps
    
    #Stack rasters
    spA_stack <- stack(sst, zoo)
    names(spA_stack) <- c('sst', 'zoo')
    
    #Assign preferences
    spA_parameters <- formatFunctions(sst = c(fun="dnorm",mean=14,sd=7),
                                      zoo = c(fun="logisticFun",alpha=-10,beta=45))
    spA_suitability <- generateSpFromFun(spA_stack,parameters=spA_parameters, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spA_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spA_parameters$sst$args[1], mean=spA_parameters$sst$args[1], sd=spA_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_zoo <- 1 / (1 + exp(((zoo@data@max) - spA_parameters$zoo$args[2])/spA_parameters$zoo$args[1]))
    ref_max <- ref_max_sst * ref_max_zoo #simple multiplication of layers.
    spA_suitability$suitab.raster <- (1/ref_max)*spA_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    # spA_suitability$suitab.raster
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    
    
    #----SPECIES B (Pelagic mobile predator): assign response curves----
    #Species B: likes to eat Species A, and warmer temperatues & shallow MLD
    
    #Stack rasters
    spB_stack <- stack(sst, mld, spA_suitability$suitab.raster)
    names(spB_stack) <- c('sst',"mld", "spA")
    
    #Assign preferences
    spB_parameters <- formatFunctions(sst = c(fun="dnorm",mean=17,sd=4),
                                      mld = c(fun="dnorm",mean=50,sd=30),
                                      spA = c(fun="logisticFun",alpha=-0.15,beta=0.4))
    spB_suitability <- generateSpFromFun(spB_stack,parameters=spB_parameters, rescale = FALSE,rescale.each.response = FALSE)
    #plot(spB_suitability$suitab.raster) #plot habitat suitability
    virtualspecies::plotResponse(spB_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spB_parameters$sst$args[1], mean=spB_parameters$sst$args[1], sd=spB_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_mld <- dnorm(spB_parameters$mld$args[1], mean=spB_parameters$mld$args[1], sd=spB_parameters$mld$args[2])
    ref_max_spA <- 1 / (1 + exp(((spA_suitability$suitab.raster@data@max) - spB_parameters$spA$args[2])/spB_parameters$spA$args[1]))
    ref_max <- ref_max_sst * ref_max_mld * ref_max_spA
    spB_suitability$suitab.raster <- (1/ref_max)*spB_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum temp is encountered
    # print(spB_suitability$suitab.raster)
    plot(spB_suitability$suitab.raster) #plot habitat suitability
    
    #-------Species B: Suitability for PREVIOUS YEAR (y-1)---------#
    #Sp A suitability y-1
    
    #Stack rasters
    spA_stack_t1 <- stack(sst_t1, zoo_t1)
    names(spA_stack_t1) <- c('sst_t1', 'zoo_t1')
    
    #Assign preferences
    spA_parameters_t1 <- formatFunctions(sst_t1 = c(fun="dnorm",mean=14,sd=7),
                                         zoo_t1 = c(fun="logisticFun",alpha=-10,beta=45))
    spA_suitability_t1 <- generateSpFromFun(spA_stack_t1,parameters=spA_parameters_t1, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
    # plot(spA_suitability_t1$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spA_suitability_t1) #plot response curves
    
    #manually rescale
    ref_max_sst_t1 <- dnorm(spA_parameters_t1$sst_t1$args[1], mean=spA_parameters_t1$sst_t1$args[1], sd=spA_parameters_t1$sst_t1$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_zoo_t1 <- 1 / (1 + exp(((zoo_t1@data@max) - spA_parameters_t1$zoo_t1$args[2])/spA_parameters_t1$zoo_t1$args[1]))
    ref_max_t1 <- ref_max_sst_t1 * ref_max_zoo_t1 #simple multiplication of layers.
    spA_suitability_t1$suitab.raster <- (1/ref_max_t1)*spA_suitability_t1$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    # spA_suitability_t1$suitab.raster
    #plot(spA_suitability_t1$suitab.raster) #plot habitat suitability
    
    #Sp B suitability in year y-1
    spB_stack_t1 <- stack(sst_t1, mld_t1, spA_suitability_t1$suitab.raster)
    names(spB_stack_t1) <- c('sst_t1',"mld_t1", "spA_t1")
    
    #Assign preferences
    spB_parameters_t1 <- formatFunctions(sst_t1 = c(fun="dnorm",mean=17,sd=4),
                                         mld_t1 = c(fun="dnorm",mean=50,sd=30),
                                         spA_t1 = c(fun="logisticFun",alpha=-0.15,beta=0.4))
    spB_suitability_t1 <- generateSpFromFun(spB_stack_t1,parameters=spB_parameters_t1, rescale = FALSE,rescale.each.response = FALSE)
    #plot(spB_suitability_t1$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(spB_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst_t1 <- dnorm(spB_parameters_t1$sst_t1$args[1], mean=spB_parameters_t1$sst_t1$args[1], sd=spB_parameters_t1$sst_t1$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_mld_t1 <- dnorm(spB_parameters_t1$mld_t1$args[1], mean=spB_parameters_t1$mld_t1$args[1], sd=spB_parameters_t1$mld_t1$args[2])
    ref_max_spA_t1 <- 1 / (1 + exp(((spA_suitability_t1$suitab.raster@data@max) - spB_parameters_t1$spA_t1$args[2])/spB_parameters_t1$spA_t1$args[1]))
    ref_max_t1 <- ref_max_sst_t1 * ref_max_mld_t1 * ref_max_spA_t1
    spB_suitability_t1$suitab.raster <- (1/ref_max_t1)*spB_suitability_t1$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum temp is encountered
    # print(spB_suitability_t1$suitab.raster)
    plot(spB_suitability_t1$suitab.raster) #plot habitat suitability
    
    #-----Species C: Turtle-like (BYCATCH SPECIES)------
    #Species C: current year 
    
    #Stack rasters
    spC_stack <- stack(sst, zoo)
    names(spC_stack) <- c('sst', 'zoo')
    
    #Assign preferences
    spC_parameters <- formatFunctions(sst = c(fun="dnorm",mean=25,sd=10),
                                         zoo = c(fun="logisticFun",alpha=-6,beta=50))
    spC_suitability <- generateSpFromFun(spC_stack,parameters=spC_parameters, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
    #plot(spC_suitability$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(spC_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spC_parameters$sst$args[1], mean=spC_parameters$sst$args[1], sd=spC_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_zoo <- 1 / (1 + exp(((zoo@data@max) - spC_parameters$zoo$args[2])/spC_parameters$zoo$args[1]))
    ref_max <- ref_max_sst * ref_max_zoo #simple multiplication of layers.
    spC_suitability$suitab.raster <- (1/ref_max)*spC_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    #plot(spC_suitability$suitab.raster) #plot habitat suitability
    
   
  #----CONVERT SP B SUITABILITY TO PRESENCE-ABSENCE----####
    
    #Use a specific function to convert suitability (0-1) to presence or absence (1 or 0)
    set.seed(y) #*********set seed to change every year, based on value of y. This will keep within year constant, but allow for noise between years
    #this will insure that the PA raster created for year Y and year Y-1 will be the same as we are using the pa.method = "probability"
    suitability_PA <- virtualspecies::convertToPA(spB_suitability, PA.method = "probability", beta = 0.4,
                                                  alpha = -0.07, species.prevalence = NULL, plot = FALSE)
    # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    # plot(suitability_PA$pa.raster)
    
###########-------BUILD FISHING SUITABILITY RASTERS FOR THE DIFFERENT SAMPLING SCENARIOS------######
    
  ##### 1. Preference for Target Species Habitat Suitability (fisher targeting behavior)
    Tar_stack <- stack(spB_suitability_t1$suitab.raster)
    names(Tar_stack) <- c('spB_t1')
    
    #1.3:Target_0.7
    Tar_params_3 <- formatFunctions(spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.7))
    Tar_suitability_3 <- generateSpFromFun(Tar_stack,parameters=Tar_params_3, rescale = FALSE,rescale.each.response = FALSE)
    # plot(Tar_suitability_3$suitab.raster) 
    # virtualspecies::plotResponse(Tar_suitability_3) #plot response curves
    #manually rescale
    Tar_reclassify_3<- Tar_suitability_3$suitab.raster/maxValue(Tar_suitability_3$suitab.raster)
    #plot(Tar_rescale)
    plot(Tar_reclassify_3)
    
    # 
  ##### 2. Distance to Port and Target Species Suitability considered
    #get "minimum distance from a port" for each cell
   
   #2.1: Northern ports only
    dist_to_ports$dist_min_NP <- apply(dist_to_ports[,c("dp4", "dp5")], 1, FUN=min)
    dist_min_raster_NP<- rasterFromXYZ(dist_to_ports[,c("lon","lat","dist_min_NP")])
    dist_min_km_NP<-dist_min_raster_NP/1000
    #plot(dist_min_km_NP)
    
    Dist_stack_NP <- stack(spB_suitability_t1$suitab.raster, dist_min_km_NP)
    names(Dist_stack_NP) <- c('spB_t1', "min_dist")
    spB<-rasterToPoints(spB_suitability_t1$suitab.raster)
    #Northern port - nearshore fishery only (50miles)
    Dist_params_NPn <- formatFunctions(min_dist = c(fun="logisticFun", alpha=50, beta=110),#100 kms (nearshore fishery active between 50-80 miles)
                                  spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.7))
    Dist_suitability_NPn <- generateSpFromFun(Dist_stack_NP,parameters=Dist_params_NPn, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Dist_suitability_NPn$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(Dist_suitability_NPn) #plot response curves
    Dist_reclassify_NPn<-Dist_suitability_NPn$suitab.raster/maxValue(Dist_suitability_NPn$suitab.raster)
    plot(Dist_reclassify_NPn)
 
    #Northern port - offshore fishery (300miles)
    Dist_params_NPo <- formatFunctions(min_dist = c(fun="logisticFun", alpha=50, beta=480),#offshore at arond 300miles
                                       spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.7))
    Dist_suitability_NPo <- generateSpFromFun(Dist_stack_NP,parameters=Dist_params_NPo, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Dist_suitability_NPo$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(Dist_suitability_NPo) #plot response curves
    Dist_reclassify_NPo<-Dist_suitability_NPo$suitab.raster/maxValue(Dist_suitability_NPo$suitab.raster)
    plot(Dist_reclassify_NPo)
    
   #2.2: Middle Ports Only 
    dist_to_ports$dist_min_MP <- apply(dist_to_ports[,c("dp2", "dp3")], 1, FUN=min)
    dist_min_raster_MP<- rasterFromXYZ(dist_to_ports[,c("lon","lat","dist_min_MP")])
    dist_min_km_MP<-dist_min_raster_MP/1000
    #plot(dist_min_km_MP)
   
    Dist_stack_MP <- stack(spB_suitability_t1$suitab.raster, dist_min_km_MP)
    names(Dist_stack_MP) <- c('spB_t1', "min_dist")
    
    #Middle ports - nearshore fishery only (50miles)
    Dist_params_MPn <- formatFunctions(min_dist = c(fun="logisticFun", alpha=50, beta=110),
                                       spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.7))
    Dist_suitability_MPn <- generateSpFromFun(Dist_stack_MP,parameters=Dist_params_MPn, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Dist_suitability_MPn$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(Dist_suitability_MPn) #plot response curves
    Dist_reclassify_MPn<-Dist_suitability_MPn$suitab.raster/maxValue(Dist_suitability_MPn$suitab.raster)
    plot(Dist_reclassify_MPn)
    
    #Middle ports - offshore fishery (300miles)
    Dist_params_MPo <- formatFunctions(min_dist = c(fun="logisticFun", alpha=50, beta=480),
                                       spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.7))
    Dist_suitability_MPo <- generateSpFromFun(Dist_stack_MP,parameters=Dist_params_MPo, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Dist_suitability_MPo$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(Dist_suitability) #plot response curves
    Dist_reclassify_MPo<-Dist_suitability_MPo$suitab.raster/maxValue(Dist_suitability_MPo$suitab.raster)
    plot(Dist_reclassify_MPo)
    
   #2.3: Southern Port only 
    dist_min_raster_SP<- rasterFromXYZ(dist_to_ports[,c("lon","lat","dp1")])
    dist_min_km_SP<-dist_min_raster_SP/1000
    #plot(dist_min_km_SP)
    
    Dist_stack_SP <- stack(spB_suitability_t1$suitab.raster, dist_min_km_SP)
    names(Dist_stack_SP) <- c('spB_t1', "min_dist")
    
    #Southern port - nearshore fishery only (50miles)
    Dist_params_SPn <- formatFunctions(min_dist = c(fun="logisticFun", alpha=50, beta=110),
                                       spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.7))
    Dist_suitability_SPn <- generateSpFromFun(Dist_stack_SP,parameters=Dist_params_SPn, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Dist_suitability_SPn$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(Dist_suitability_SPn) #plot response curves
    Dist_reclassify_SPn<-Dist_suitability_SPn$suitab.raster/maxValue(Dist_suitability_SPn$suitab.raster)
    plot(Dist_reclassify_SPn)
    
    #Southern port - offshore fishery (300miles)
    Dist_params_SPo <- formatFunctions(min_dist = c(fun="logisticFun", alpha=50, beta=480),
                                       spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.7))
    Dist_suitability_SPo <- generateSpFromFun(Dist_stack_SP,parameters=Dist_params_SPo, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Dist_suitability_SPo$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(Dist_suitability_SPo) #plot response curves
    Dist_reclassify_SPo<-Dist_suitability_SPo$suitab.raster/maxValue(Dist_suitability_SPo$suitab.raster)
    plot(Dist_reclassify_SPo)
    
   #2.4: All Ports/Coastwide
    dist_to_ports$dist_min_All <- apply(dist_to_ports[,c("dp1", "dp2", "dp3", "dp4", "dp5")], 1, FUN=min)
    dist_min_raster_All<- rasterFromXYZ(dist_to_ports[,c("lon","lat","dist_min_All")])
    dist_min_km_All<-dist_min_raster_All/1000
    #plot(dist_min_km_All)
  
    Dist_stack_All <- stack(spB_suitability_t1$suitab.raster, dist_min_km_All)
    names(Dist_stack_All) <- c('spB_t1', "min_dist")
    
   #All port - nearshore fishery only (50miles)
    Dist_params_Alln <- formatFunctions(min_dist = c(fun="logisticFun", alpha=50, beta=110),
                                       spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.7))
    Dist_suitability_Alln <- generateSpFromFun(Dist_stack_All,parameters=Dist_params_Alln, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Dist_suitability_Alln$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(Dist_suitability_Alln) #plot response curves
    Dist_reclassify_Alln<-Dist_suitability_Alln$suitab.raster/maxValue(Dist_suitability_Alln$suitab.raster)
    plot(Dist_reclassify_Alln)
    
   #All port - offshore fishery (300miles)
    Dist_params_Allo <- formatFunctions(min_dist = c(fun="logisticFun", alpha=50, beta=480),
                                       spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.7))
    Dist_suitability_Allo <- generateSpFromFun(Dist_stack_All,parameters=Dist_params_Allo, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Dist_suitability_Allo$suitab.raster) #plot habitat suitability
    #virtualspecies::plotResponse(Dist_suitability) #plot response curves
    Dist_reclassify_Allo<-Dist_suitability_Allo$suitab.raster/maxValue(Dist_suitability_Allo$suitab.raster)
    plot(Dist_reclassify_Allo)
    
  ##### 3. BYCATCH RISK and target spp suitability considered

    Bycatch_stack <- stack(spB_suitability_t1$suitab.raster, spC_suitability$suitab.raster)
    names(Bycatch_stack) <- c('spB_t1', "spC")

    #Bycatch_2 preferences
    Bycatch_params_2 <- formatFunctions(spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.7),
                                        spC = c(fun="logisticFun", alpha= 0.05, beta=0.5))
    BY_suitability_2 <- generateSpFromFun(Bycatch_stack,parameters=Bycatch_params_2, rescale = FALSE,rescale.each.response = FALSE)
    #plot(Bycatch_suitability_2$suitab.raster) 
    #virtualspecies::plotResponse(BY_suitability_2) #plot response curves
    By_reclassify_2<-BY_suitability_2$suitab.raster/maxValue(BY_suitability_2$suitab.raster)
    plot(By_reclassify_2)
    
  ##### 4. Static Closed Area ######
   
    #4.1: Smallest closed area
    Opt_closed_1<-Tar_suitability_3$suitab.raster*closed_area_1
    #plot(Opt_closed_1)
    closed_reclassify_1<-Opt_closed_1/maxValue(Opt_closed_1)
    plot(closed_reclassify_1)
    
    #4.2: Medium closed area
    Opt_closed_2<-Tar_suitability_3$suitab.raster*closed_area_2
    #plot(Opt_closed_2)
    closed_reclassify_2<-Opt_closed_2/maxValue(Opt_closed_2)
    plot(closed_reclassify_2)
    
    #4.3: Largest closed area
    Opt_closed_3<-Tar_suitability_3$suitab.raster*closed_area_3
    #plot(Opt_closed_3)
    closed_reclassify_3<-Opt_closed_3/maxValue(Opt_closed_3)
    plot(closed_reclassify_3)
    
  ###########-----SAMPLE PRESENCES AND ABSENCES-----##########
    
    #******Random Sampling of nsamples*******
    presence.points.random <- sampleOccurrences(suitability_PA,n = nsamples,type = "presence-absence",
                                                detection.probability = 1, error.probability=0, plot = TRUE,
                                                sample.prevalence = NULL)
    
    #convert to dataframe
    pres_df_random <- cbind(as.data.frame(presence.points.random$sample.points$x),as.data.frame(presence.points.random$sample.points$y))
    colnames(pres_df_random) <- c("x","y")
    pres_df_random$random_sampled <- 1
    
    #expand dataframe to include all possible locations
    df_full <- as.data.frame(rasterToPoints(sst)[,1:2]) #picking an example raster to extract lat and lon from
    df_full_2 <- left_join(df_full, pres_df_random, by=c('x','y'))
    df_full_2$random_sampled <- ifelse(is.na(df_full_2$random_sampled),0,df_full_2$random_sampled)
    
    #****** 1.1 Preferential Sampling of nsamples - based on Precence of target species 0.5**********
    presence.points.pref1<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                             detection.probability = 1, bias = "manual", weights = Tar_reclassify_3, plot = TRUE)
    #convert to dataframe
    pres_df_pref1 <- cbind(as.data.frame(presence.points.pref1$sample.points$x),as.data.frame(presence.points.pref1$sample.points$y))
    colnames(pres_df_pref1) <- c("x","y")
    pres_df_pref1$pref_sampled <- 1
    #add to dataframe
    df_full_3 <- left_join(df_full_2, pres_df_pref1, by=c('x','y'))
    df_full_3$pref_sampled <- ifelse(is.na(df_full_3$pref_sampled),0,df_full_3$pref_sampled)
    
    #******* 2.1 Sampling of nsamples based on Habitat Suitability & DISTANCE TO PORT NPoff********** 
    presence.points.dist_npo<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                detection.probability = 1, bias = "manual", weights= Dist_reclassify_NPo, plot = TRUE)
    #convert to dataframe
    pres_df_dist_npo <- cbind(as.data.frame(presence.points.dist_npo$sample.points$x),as.data.frame(presence.points.dist_npo$sample.points$y))
    #pres_df_dist <- cbind(as.data.frame(Dist_data$lons[Dist_data$cell_num %in% cells_dist]),as.data.frame(Dist_data$lats[Dist_data$cell_num %in% cells_dist]))
    colnames(pres_df_dist_npo) <- c("x","y")
    pres_df_dist_npo$dist_sampled_npo <- 1
    
    #add to dataframe
    df_full_4 <- left_join(df_full_3, pres_df_dist_npo, by=c('x','y'))
    df_full_4$dist_sampled_npo <- ifelse(is.na(df_full_4$dist_sampled_npo),0,df_full_4$dist_sampled_npo)
    
    #******* 2.2 Sampling of nsamples based on Habitat Suitability & DISTANCE TO PORT NPnear********** 
    presence.points.dist_npn<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                detection.probability = 1, bias = "manual", weights= Dist_reclassify_NPn, plot = TRUE)
    #convert to dataframe
    pres_df_dist_npn <- cbind(as.data.frame(presence.points.dist_npn$sample.points$x),as.data.frame(presence.points.dist_npn$sample.points$y))
    #pres_df_dist <- cbind(as.data.frame(Dist_data$lons[Dist_data$cell_num %in% cells_dist]),as.data.frame(Dist_data$lats[Dist_data$cell_num %in% cells_dist]))
    colnames(pres_df_dist_npn) <- c("x","y")
    pres_df_dist_npn$dist_sampled_npn <- 1
    
    #add to dataframe
    df_full_5 <- left_join(df_full_4, pres_df_dist_npn, by=c('x','y'))
    df_full_5$dist_sampled_npn <- ifelse(is.na(df_full_5$dist_sampled_npn),0,df_full_5$dist_sampled_npn)
    
    #******* 2.3 Sampling of nsamples based on Habitat Suitability & DISTANCE TO PORT MPoff********** 
    presence.points.dist_mpo<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                detection.probability = 1, bias = "manual", weights= Dist_reclassify_MPo, plot = TRUE)
    #convert to dataframe
    pres_df_dist_mpo <- cbind(as.data.frame(presence.points.dist_mpo$sample.points$x),as.data.frame(presence.points.dist_mpo$sample.points$y))
    #pres_df_dist <- cbind(as.data.frame(Dist_data$lons[Dist_data$cell_num %in% cells_dist]),as.data.frame(Dist_data$lats[Dist_data$cell_num %in% cells_dist]))
    colnames(pres_df_dist_mpo) <- c("x","y")
    pres_df_dist_mpo$dist_sampled_mpo <- 1
    
    #add to dataframe
    df_full_6 <- left_join(df_full_5, pres_df_dist_mpo, by=c('x','y'))
    df_full_6$dist_sampled_mpo <- ifelse(is.na(df_full_6$dist_sampled_mpo),0,df_full_6$dist_sampled_mpo)
    
    #******* 2.4 Sampling of nsamples based on Habitat Suitability & DISTANCE TO PORT MPnear********* 
    presence.points.dist_mpn<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                detection.probability = 1, bias = "manual", weights= Dist_reclassify_MPn, plot = TRUE)
    #convert to dataframe
    pres_df_dist_mpn <- cbind(as.data.frame(presence.points.dist_mpn$sample.points$x),as.data.frame(presence.points.dist_mpn$sample.points$y))
    #pres_df_dist <- cbind(as.data.frame(Dist_data$lons[Dist_data$cell_num %in% cells_dist]),as.data.frame(Dist_data$lats[Dist_data$cell_num %in% cells_dist]))
    colnames(pres_df_dist_mpn) <- c("x","y")
    pres_df_dist_mpn$dist_sampled_mpn <- 1
    
    #add to dataframe
    df_full_7 <- left_join(df_full_6, pres_df_dist_mpn, by=c('x','y'))
    df_full_7$dist_sampled_mpn <- ifelse(is.na(df_full_7$dist_sampled_mpn),0,df_full_7$dist_sampled_mpn)
    
    #******* 2.5 Sampling of nsamples based on Habitat Suitability & DISTANCE TO PORT SPoff********* 
    presence.points.dist_spo<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                detection.probability = 1, bias = "manual", weights= Dist_reclassify_SPo, plot = TRUE)
    #convert to dataframe
    pres_df_dist_spo <- cbind(as.data.frame(presence.points.dist_spo$sample.points$x),as.data.frame(presence.points.dist_spo$sample.points$y))
    #pres_df_dist <- cbind(as.data.frame(Dist_data$lons[Dist_data$cell_num %in% cells_dist]),as.data.frame(Dist_data$lats[Dist_data$cell_num %in% cells_dist]))
    colnames(pres_df_dist_spo) <- c("x","y")
    pres_df_dist_spo$dist_sampled_spo <- 1
    
    #add to dataframe
    df_full_8 <- left_join(df_full_7, pres_df_dist_spo, by=c('x','y'))
    df_full_8$dist_sampled_spo <- ifelse(is.na(df_full_8$dist_sampled_spo),0,df_full_8$dist_sampled_spo)
    
    #******* 2.6 Sampling of nsamples based on Habitat Suitability & DISTANCE TO PORT SPnear********* 
    presence.points.dist_spn<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                detection.probability = 1, bias = "manual", weights= Dist_reclassify_SPn, plot = TRUE)
    #convert to dataframe
    pres_df_dist_spn <- cbind(as.data.frame(presence.points.dist_spn$sample.points$x),as.data.frame(presence.points.dist_spn$sample.points$y))
    #pres_df_dist <- cbind(as.data.frame(Dist_data$lons[Dist_data$cell_num %in% cells_dist]),as.data.frame(Dist_data$lats[Dist_data$cell_num %in% cells_dist]))
    colnames(pres_df_dist_spn) <- c("x","y")
    pres_df_dist_spn$dist_sampled_spn <- 1
    
    #add to dataframe
    df_full_9 <- left_join(df_full_8, pres_df_dist_spn, by=c('x','y'))
    df_full_9$dist_sampled_spn <- ifelse(is.na(df_full_9$dist_sampled_spn),0,df_full_9$dist_sampled_spn)
    
    #******* 2.7 Sampling of nsamples based on Habitat Suitability & DISTANCE TO PORT Alloff********* 
    presence.points.dist_allo<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                 detection.probability = 1, bias = "manual", weights= Dist_reclassify_Allo, plot = TRUE)
    #convert to dataframe
    pres_df_dist_allo <- cbind(as.data.frame(presence.points.dist_allo$sample.points$x),as.data.frame(presence.points.dist_allo$sample.points$y))
    #pres_df_dist <- cbind(as.data.frame(Dist_data$lons[Dist_data$cell_num %in% cells_dist]),as.data.frame(Dist_data$lats[Dist_data$cell_num %in% cells_dist]))
    colnames(pres_df_dist_allo) <- c("x","y")
    pres_df_dist_allo$dist_sampled_allo <- 1
    
    #add to dataframe
    df_full_10 <- left_join(df_full_9, pres_df_dist_allo, by=c('x','y'))
    df_full_10$dist_sampled_allo <- ifelse(is.na(df_full_10$dist_sampled_allo),0,df_full_10$dist_sampled_allo)
    
    #******* 2.8 Sampling of nsamples based on Habitat Suitability & DISTANCE TO PORT Allnear********* 
    presence.points.dist_alln<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                 detection.probability = 1, bias = "manual", weights= Dist_reclassify_Alln, plot = TRUE)
    #convert to dataframe
    pres_df_dist_alln <- cbind(as.data.frame(presence.points.dist_alln$sample.points$x),as.data.frame(presence.points.dist_alln$sample.points$y))
    #pres_df_dist <- cbind(as.data.frame(Dist_data$lons[Dist_data$cell_num %in% cells_dist]),as.data.frame(Dist_data$lats[Dist_data$cell_num %in% cells_dist]))
    colnames(pres_df_dist_alln) <- c("x","y")
    pres_df_dist_alln$dist_sampled_alln <- 1
    
    #add to dataframe
    df_full_11 <- left_join(df_full_10, pres_df_dist_alln, by=c('x','y'))
    df_full_11$dist_sampled_alln <- ifelse(is.na(df_full_11$dist_sampled_alln),0,df_full_11$dist_sampled_alln)
    
    #****** 3.1 Closed Area Sampling - based on HABITAT SUITABILITY of target species & Small closed area**********
    presence.points.close1<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                              detection.probability = 1, bias = "manual", weights = closed_reclassify_1)
    close_df_pref_1 <- cbind(as.data.frame(presence.points.close1$sample.points$x),as.data.frame(presence.points.close1$sample.points$y))
    colnames(close_df_pref_1) <- c("x","y")
    close_df_pref_1$Closed_sampled_1 <- 1
    
    #add to dataframe
    df_full_12 <- left_join(df_full_11, close_df_pref_1, by=c('x','y'))
    df_full_12$Closed_sampled_1 <- ifelse(is.na(df_full_12$Closed_sampled_1),0,df_full_12$Closed_sampled_1)
    
    #****** 3.2 Closed Area Sampling - based on HABITAT SUITABILITY of target species & Medium closed area**********
    presence.points.close2<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                              detection.probability = 1, bias = "manual", weights = closed_reclassify_2)
    close_df_pref_2 <- cbind(as.data.frame(presence.points.close2$sample.points$x),as.data.frame(presence.points.close2$sample.points$y))
    colnames(close_df_pref_2) <- c("x","y")
    close_df_pref_2$Closed_sampled_2 <- 1
    
    #add to dataframe
    df_full_13 <- left_join(df_full_12, close_df_pref_2, by=c('x','y'))
    df_full_13$Closed_sampled_2 <- ifelse(is.na(df_full_13$Closed_sampled_2),0,df_full_13$Closed_sampled_2)
    
    #****** 3.3 Closed Area Sampling - based on HABITAT SUITABILITY of target species & Large closed area**********
    presence.points.close3<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                              detection.probability = 1, bias = "manual", weights = closed_reclassify_3)
    close_df_pref_3 <- cbind(as.data.frame(presence.points.close3$sample.points$x),as.data.frame(presence.points.close3$sample.points$y))
    colnames(close_df_pref_3) <- c("x","y")
    close_df_pref_3$Closed_sampled_3 <- 1
    #add to dataframe
    df_full_14 <- left_join(df_full_13, close_df_pref_3, by=c('x','y'))
    df_full_14$Closed_sampled_3 <- ifelse(is.na(df_full_14$Closed_sampled_3),0,df_full_14$Closed_sampled_3)
    
    #******* 4.1 Sampling of nsamples based on Habitat Suitability & BYCATCH RISK********** 
    
    presence.points.Bycatch2<-sampleOccurrences(suitability_PA, n=nsamples, type="presence-absence",
                                                detection.probability = 1, bias = "manual", weights= By_reclassify_2, plot = TRUE)
    #convert to dataframe
    pres_df_BY2 <- cbind(as.data.frame(presence.points.Bycatch2$sample.points$x),as.data.frame(presence.points.Bycatch2$sample.points$y))
    colnames(pres_df_BY2) <- c("x","y")
    pres_df_BY2$BY_sampled_2 <- 1
    #add to dataframe
    df_full_15 <- left_join(df_full_14, pres_df_BY2, by=c('x','y'))
    df_full_15$BY_sampled_2 <- ifelse(is.na(df_full_15$BY_sampled_2),0,df_full_15$BY_sampled_2)
    
    colnames(df_full_15)[1:2] <- c("lon","lat")
    
    #----EXTRACT DATA for each year----
    
    print("Extracting suitability")
    ei <- 21912*y #end location in output grid to index to
    se <- ei - (21912-1) #start location in output grid to index to
    output$lat[se:ei] <- rasterToPoints(sst)[,2]
    output$lon[se:ei] <- rasterToPoints(sst)[,1]
    output$year[se:ei] <- rep(years[y],21912)
    output$pres[se:ei] <- rasterToPoints(suitability_PA$pa.raster)[,3] 
    output$suitability_t[se:ei] <- rasterToPoints(spB_suitability$suitab.raster)[,3]  #extract points from suitability file
    output$spC_suitability[se:ei] <- rasterToPoints(spC_suitability$suitab.raster)[,3] 
    output$spA_suitability[se:ei] <- rasterToPoints(spA_suitability$suitab.raster)[,3] 
    output$sst[se:ei] <-  rasterToPoints(sst)[,3]   #extract points from suitability file
    output$zoo_200[se:ei] <-  rasterToPoints(zoo)[,3] 
    output$mld[se:ei] <-  rasterToPoints(mld)[,3] 
    output$chl_surface[se:ei] <-  rasterToPoints(chl_surface)[,3] 
    #temporary storage: make sure sample coordinates match output coordinates
    temp <-   left_join(output[se:ei,1:7], df_full_15[1:16], by=c('lon','lat')) 
    # temp_pref_1 <-   left_join(output[se:ei,1:6], df_full_19[,c(1,2,4)], by=c('lon','lat'))
    # temp_pref_2 <- left_join(output[se:ei,1:6], df_full_19[,c(1,2,5)], by=c('lon','lat'))
    # 
    # temp_Dist <-   left_join(output[se:ei,1:6], df_full_19[,c(1,2,10)], by=c('lon','lat'))
    # temp_BY<- left_join(output[se:ei,1:6], df_full_19[,c(1,2,12)], by=c('lon','lat'))
    # temp_ByDist<- left_join(output[se:ei,1:6], df_full_19[,c(1,2,14)], by=c('lon','lat'))
    # temp_Closed<- left_join(output[se:ei,1:6], df_full_19[,c(1,2,16)], by=c('lon','lat'))
    #assign to output
    output$random_sampled[se:ei] <- temp$random_sampled
    output$pref_sampled[se:ei] <- temp$pref_sampled
    output$dist_sampled_npo[se:ei] <- temp$dist_sampled_npo
    output$dist_sampled_npn[se:ei] <- temp$dist_sampled_npn
    output$dist_sampled_mpo[se:ei] <- temp$dist_sampled_mpo
    output$dist_sampled_mpn[se:ei] <- temp$dist_sampled_mpn
    output$dist_sampled_spo[se:ei] <- temp$dist_sampled_spo
    output$dist_sampled_spn[se:ei] <- temp$dist_sampled_spn
    output$dist_sampled_allo[se:ei] <- temp$dist_sampled_allo
    output$dist_sampled_alln[se:ei] <- temp$dist_sampled_alln
    output$BY_sampled_2[se:ei]<- temp$BY_sampled_2
    output$Closed_sampled_1[se:ei]<- temp$Closed_sampled_1
    output$Closed_sampled_2[se:ei]<- temp$Closed_sampled_2
    output$Closed_sampled_3[se:ei]<- temp$Closed_sampled_3
  }
  
  #Average monthly biomass available to CCS is: 1.18x10^5 ± (0.13x10^5 se) mt (from Desiree Tommasi)
  mean_spatial <- log((118000/21912)*5) #Kind of had to make these up to get biomass pop that matches 1.18, and with error low enough that the EM does ok.
  se_spatial <- log((13000)/10000)
  output$abundance <- ifelse(output$pres==1,rlnorm(nrow(output),mean_spatial, se_spatial)*output$suitability_t,0)
  
  # #--------PRINT and SAVE CSV------------------#####
  
  print('Saving csv to working directory')
  write.csv(output, 'FisheryDependent_OM_Simulation_11_8_21_had_newPrefs.csv',row.names = FALSE)
  return(output)
}

#test<-SimulateWorld_ROMS_FishDepFun_Final(dir=dir, nsamples = 100)
