## Calculating distance from port for each cell in ROMS raster layer 
library(raster)
library(rgdal)
setwd("~/DisMAP project/Location, Location, Location/Location Workshop/ROMS/hadley_11_2_21/sst_spring_avg")
#setwd("~/DisMAP project/Location, Location, Location/Location Workshop/ROMS/gfdl/sst_spring_avg")

#create output data frame for all the distance to port data
output <- as.data.frame(matrix(NA, nrow=33300,ncol=7)) #33300 cells for Had, 33666 for gfdl
colnames(output) <- c("lon", "lat", "dp1", "dp2", "dp3", "dp4", "dp5")

p1<- c(-117.1441, 32.6717)#San Diego bay, CA
p2<- c(-122.001620, 36.965719)#Santa Cruz, CA
p3<- c(-123.050618, 38.334302)#Bodega Bay, CA
p4<- c(-124.292000, 43.383975)#Garibaldi, OR
p5<- c(-124.114934, 46.911534)#Westport, WA

#read in ROMS sst data and convert to raster layer, and plot Port locations on Map
ROMS<-raster("sst_monthly_had_SpringMean_2100.grd")
plot(ROMS)
points(-117.1441, 32.6717, pch=0, cex=2)
points(-122.001620, 36.965719, pch=0, cex=2)
points(-123.050618, 38.334302, pch=0, cex=2)
points(-124.292000, 43.383975, pch=0, cex=2)
points(-124.114934, 46.911534, pch=0, cex=2)
 rNA <- sum(!is.na(ROMS))  
freq(ROMS, value=NA) # counts number of NA cells in raster
## Calculate distance of every non-NA cell to each Port
dp1<- distanceFromPoints(ROMS, p1) ## distances are in meters I think?
dp2<- distanceFromPoints(ROMS, p2)
dp3<- distanceFromPoints(ROMS, p3)
dp4<- distanceFromPoints(ROMS, p4)
dp5<- distanceFromPoints(ROMS, p5)

#Extract distance values for each cell and put in output file
output$lat <- rasterToPoints(dp1)[,2]
output$lon <- rasterToPoints(dp1)[,1]
output$dp1 <- rasterToPoints(dp1)[,3]
output$dp2 <- rasterToPoints(dp2)[,3]
output$dp3 <- rasterToPoints(dp3)[,3]
output$dp4 <- rasterToPoints(dp4)[,3]
output$dp5 <- rasterToPoints (dp5)[,3]
head(output)

write.csv(output, "~/DisMAP project/Location, Location, Location/Location Workshop/Dist_to_Ports_had.csv")

## Now have a distance dataframe where fishemen can only land fish in CA ports, so 
 # eliminate dp3 and dp4 (the OR and WA ports)
#Extract distance values for each cell and put in output file
output2 <- as.data.frame(matrix(NA, nrow=33300,ncol=5))
colnames(output2) <- c("lon", "lat", "dp1", "dp2", "dp3")

output2$lat <- rasterToPoints(dp1)[,2]
output2$lon <- rasterToPoints(dp1)[,1]
output2$dp1 <- rasterToPoints(dp1)[,3]
output2$dp2 <- rasterToPoints(dp2)[,3]
output2$dp3 <- rasterToPoints(dp3)[,3]
head(output2)
write.csv(output2, "~/DisMAP project/Location, Location, Location/Location Workshop/Dist_to_Ports_CAports.csv")


##Now Just WA and OR ports
#Extract distance values for each cell and put in output file
output3 <- as.data.frame(matrix(NA, nrow=33300,ncol=4))
colnames(output3) <- c("lon", "lat", "dp4", "dp5")

output3$lat <- rasterToPoints(dp4)[,2]
output3$lon <- rasterToPoints(dp4)[,1]
output3$dp4 <- rasterToPoints(dp4)[,3]
output3$dp5 <- rasterToPoints(dp5)[,3]

head(output3)
write.csv(output3, "~/DisMAP project/Location, Location, Location/Location Workshop/Dist_to_Ports_WAORports.csv")


#Assign preferences - unused code from simulation distance section 
# linear.function <- function(x, a, b)
# {
#   a*x + b
# }
# Dist_params <- formatFunctions(min_dist = c(fun="linear.function",a=-0.0007, b=1),
#                                spB_t1 = c(fun="logisticFun",alpha=-0.05,beta=0.5))