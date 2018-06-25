# Set directories
mainDir <- "S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS"
setwd(mainDir)

#=======================================================================#
# Inital preparation (Fusion and LAStools)
#=======================================================================#
# Merge LAZ files into a single point cloud
system(paste(file.path("C:","Fusion", "mergedata.exe"),
             "S:\COS\PyroGeog\amartinez\UAV_seedlings\Lidar\fullLAZ\filelist.txt",
             file.path(mainDir, "Seedling_merge.las"),
             sep = " "))

# Tile pointcloud into tiles of 2 million points 
system(paste(file.path("C:","Fusion", "lastile.exe"),
             "-i", file.path(mainDir, "Seedling_merge.las"),
             "-o", file.path(mainDir, "Tiles", "tile.las"),
             "-refine 20000000", "-cores 6",
             sep = " "))
system(paste("dir/s/b",
             file.path(mainDir, "Tiles", "*.las"),
             ">", file.path(mainDir, "TileList.txt"),
             sep = " "))


# Clip tiles to seedling area
system(paste(file.path("C:","Fusion", "clipdata.exe"),
             "/shape:0 /index",
             file.path(mainDir, "TileList.txt"),
             file.path(mainDir, "MoscowMtn_clip.las"),
             "2380665 1893999 2381231 1894722",
             sep=" "))

# View LAS information using LAStools and Fusion
system(paste(file.path("C:","LAStools", "bin", "lasinfo.exe"),
             "-i", file.path(mainDir, "LAS_norm.las"),
             "-last_only",
             "-cd",
             sep = " "))
system(paste(file.path("C:","Fusion", "catalog.exe"),
             "/density:0.125,5,10",
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))

#=======================================================================#
# Create ground model -- 0.25 ft (3 in) resolution
#=======================================================================#
# Load packages
library(lidR)
library(rLiDAR)
#library(raster)
#library(rgeos)
#library(rgdal)

# Import LAS
las <- readLAS("S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/MoscowMtn_clip.las")

# Identify ground points (takes approx. 20 hours to run on this dataset)
dem.las <- grid_terrain(las, res = 0.25, method = "knnidw", k = 10, p = 2,
                    model = gstat::vgm(0.59, "Sph", 874), keep_lowest = FALSE)
dem.raster <- as.raster(dem.las)
writeRaster(dem.raster, file.path(mainDir, "Ground", "pmf_dem.tif"))

#=======================================================================#
# Create canopy model -- 0.25 ft (3 in) resolution
#=======================================================================#
# Remove power lines
hist(las@data$Z)
las.pwr.rm <- lasfilter(las, Z < 3050)
writeLAS(las.pwr.rm, file.path(mainDir, "LAS_pwr_rm.las"))

# Create "spike free" canopy height model
system(paste(file.path("C:","LAStools", "bin", "las2dem.exe"),  # CHM with normalized LAS
             "-i", file.path(mainDir, "LAS_pwr_rm.las"),
             "-spike_free 0.45",
             "-step 0.25",
             "-o",  file.path(mainDir, "chm", "CHM_spike_free.asc"),
             sep=" "))
chm.sf <- raster(file.path(mainDir, "chm", "CHM_spike_free.asc"))

############Untested########
#=======================================================================#
# Locate trees
#=======================================================================#
# Identify tree tops
tree.top <- tree_detection(chm2, ws = 25, hmin = 1) #6 foot radius, 12 inch minimum seedling 
tree.top.pts <- rasterToPoints(tree.top, spatial = F)
shapefile(rasterToPoints(tree.det.sf, spatial = T), 
          filename = file.path(mainDir, "Trees", "tree_tops.shp"), overwrite = T)

# Classify trees in point cloud
lastrees_silva(las.pwr.rm, chm.sf, tree.top, max_cr_factor = 0.8)
writeLAS(las.pwr.rm, file.path(mainDir, "LAS_pwr_rm_class.las"))

# Canopy
canopy <- ForestCAS(chm.sf, tree.top, maxcrown = 6, exclusion = 0)




############################################################################
#                                  END                                     #
############################################################################
# Import canopy height model
chm <- raster(file.path(mainDir, "chm", "MoscowMtn_clip_CHM.asc"))
plot(chm)
summary(chm)

# Smooth canopy (set parameters and execute)
ws <- 5 # 3x3 window size
filter <- "Gaussian" # Gaussian filter type
sigma <- 0.5
sCHM <- CHMsmoothing(chm, filter, ws, sigma)
plot(sCHM)

# Individual tree detection list (set parameters and execute)
fws <- 5 # 3x3 fixed window size
minht <- 1 # 0.5 ft. minimum height above ground
loc.trees <- FindTreesCHM(sCHM, fws, minht)
loc <- subset(loc.trees, height < 7)
coords <- loc[,1:2]
data <- as.data.frame(loc[,3])
trees <- SpatialPointsDataFrame(coords, data = data, proj4string = CRS("+init=epsg:6453"))
plot(trees)
plot(loc.obs)
shapefile(trees, "Trees/trees.shp", overwrite = TRUE)


loc.obs <- readOGR("S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/GIS/SeedClip.shp")


###
# 2
system(paste(file.path("C:","Fusion", "groundfilter.exe"),
             "/gparam:0 /wparam:1 /smooth:3 /iterations:10",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts2.las"),
             0.25,
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))
system(paste(file.path("C:","Fusion", "gridsurfacecreate.exe"),
             "/smooth:5", 
             file.path(mainDir, "Ground", "MoscowMtn_DEM2.dtm"),
             "0.25 F F 2 11 2 0",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts2.las"),
             sep=" "))
system(paste(file.path("C:","Fusion", "dtm2ascii.exe"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM2.dtm"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM2.asc"),
             sep=" "))
#
# 3
system(paste(file.path("C:","Fusion", "groundfilter.exe"),
             "/gparam:-2 /wparam:2.5 /smooth:3 /iterations:10",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts3.las"),
             0.25,
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))
system(paste(file.path("C:","Fusion", "gridsurfacecreate.exe"),
             "/smooth:5", 
             file.path(mainDir, "Ground", "MoscowMtn_DEM3.dtm"),
             "0.25 F F 2 11 2 0",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts3.las"),
             sep=" "))
system(paste(file.path("C:","Fusion", "dtm2ascii.exe"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM3.dtm"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM3.asc"),
             sep=" "))
# best so far
# 4
system(paste(file.path("C:","Fusion", "groundfilter.exe"),
             "/gparam:-4 /wparam:4.5 /smooth:3 /iterations:10",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts4.las"),
             0.25,
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))
system(paste(file.path("C:","Fusion", "gridsurfacecreate.exe"),
             "/smooth:5", 
             file.path(mainDir, "Ground", "MoscowMtn_DEM4.dtm"),
             "0.25 F F 2 11 2 0",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts4.las"),
             sep=" "))
system(paste(file.path("C:","Fusion", "dtm2ascii.exe"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM4.dtm"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM4.asc"),
             sep=" "))
#
# 5
system(paste(file.path("C:","Fusion", "groundfilter.exe"),
             "/gparam:-2 /wparam:2.5 /median:5 /iterations:10",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts5.las"),
             0.25,
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))
system(paste(file.path("C:","Fusion", "gridsurfacecreate.exe"),
             "/smooth:5", 
             file.path(mainDir, "Ground", "MoscowMtn_DEM5.dtm"),
             "0.25 F F 2 11 2 0",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts5.las"),
             sep=" "))
system(paste(file.path("C:","Fusion", "dtm2ascii.exe"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM5.dtm"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM5.asc"),
             sep=" "))


system(paste(file.path("C:","Fusion", "lda2ascii.exe"),
             file.path(mainDir, "MoscowMtn_clip.las"),
             file.path(mainDir, "MoscowMtn_clip.txt"),
             "1",
             sep=" "))
# 
# 4
system(paste(file.path("C:","Fusion", "groundfilter.exe"),
             "/gparam:-4 /wparam:4.5 /smooth:3 /iterations:10",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts4.las"),
             0.25,
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))
system(paste(file.path("C:","Fusion", "gridsurfacecreate.exe"),
             "/smooth:5", 
             file.path(mainDir, "Ground", "MoscowMtn_DEM4.dtm"),
             "0.25 F F 2 11 2 0",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts4.las"),
             sep=" "))
system(paste(file.path("C:","Fusion", "dtm2ascii.exe"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM4.dtm"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM4.asc"),
             sep=" "))


####
library(lidR)
las <- readLAS("S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/MoscowMtn_clip.las")
ws <- c(1,2,4,8,16,32,64,128)
th <- seq(0.1, 2, length.out = length(ws))
las.ground <- lasground(las, MaxWinSize = 20, Slope = 0.3, InitDist = 0.15, MaxDist = 10, CellSize = 0.2)

####
las <- readLAS("S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/MoscowMtn_clip.las")
ws <- c(1,2,4,8,16,32,64,128)
th <- seq(0.1, 2, length.out = length(ws))
las.ground <- lasground(las, "pmf", ws, th)

###
las <- readLAS("S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/MoscowMtn_clip.las")
u <- util_makeZhangParam(b = 2, dh0 = 0.15, dhmax = 3, s = 0.35, max_ws = 15)
Sys.time()
lasground(las, "pmf", u$ws, u$th)
writeLAS(las, "S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/LAS_gr_class.las")
Sys.time()

dem <- grid_terrain(las, res = 0.25, method = "knnidw", k = 10, p = 2,
             model = gstat::vgm(0.59, "Sph", 874), keep_lowest = FALSE)
dem2 <- as.raster(dem)
writeRaster(dem2, "S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/Ground/pmf_dem.tif")

##############
############## Contuinuing on
##############
las.norm <- lasnormalize(las, dem, method = "knnidw", k = 10, p= 2, model = gstat::vgm(0.59, "Sph", 874), copy = T)
las.filt <- lasfilter(las.norm, Z < 150)
writeLAS(las.filt, "S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/LAS_norm.las")
chm <- grid_canopy(las.filt, res = 0.25, na.fill = "knnidw", k = 10, p = 2)
chm2 <- as.raster(chm)
writeRaster(chm2, "S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/chm/chm.tif", overwrite = T)

tree.det <- tree_detection(chm, ws = 25, 1) #18 inch radius, 12 inch minimum seedling 
writeRaster(tree.det, "S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/chm/tree_det.tif", overwrite = T)
lastrees_silva(las.filt, chm2, tree.det)
tree.met <- tree_metrics(las.filt, func = quantile(Z, probs = 0.9))
tree.met.seed <- subset(tree.met, V1 < 6)



system(paste(file.path("C:","LAStools", "bin", "las2dem.exe"),  # CHM with normalized LAS
      "-i", file.path(mainDir, "LAS_norm.las"),
      "-spike_free 0.5",
      "-step 0.25",
      "-o",  file.path(mainDir, "chm", "CHM_spike_free.asc"),
      sep=" "))
chm.sf <- raster(file.path(mainDir, "chm", "CHM_spike_free.asc"))
tree.det.sf <- tree_detection(chm.sf, ws = 17, 1) #18 inch radius, 12 inch minimum seedling 
summary(tree.det.sf)
tree.point.sf <- rasterToPoints(tree.det.sf, spatial = T)
shapefile(tree.point.sf, filename = file.path(mainDir, "Trees", "tree_sf"))

tree.met.sf <- tree_metrics(las.filt, func = quantile(Z, probs = 0.9))
tree.met.sf.seed <- subset(tree.met, V1 < 6)

###
###
###
paste(file.path("C:","LAStools", "bin", "las2dem.exe"),  # CHM with normalized LAS
             "-i", file.path(mainDir, "LAS_norm.las"),
             "-spike_free 0.5",
             #"-subcircle 0.15",
             "-step 0.25",
             "-o",  file.path(mainDir, "chm", "chm_sf6.asc"),
             sep=" ")

paste(file.path("C:","LAStools", "bin", "las2dem.exe"),   # CHM with original LAS
      "-i", file.path(mainDir, "MoscowMtn_clip.las"),
      "-spike_free 0.5",
      #"-subcircle 0.15",
      "-step 0.25",
      "-o",  file.path(mainDir, "chm", "chm_orig_sf.asc"),
      sep=" ")

paste(file.path("C:","LAStools", "bin", "lasview.exe"),
      "-i", file.path(mainDir, "chm", "chm_sf3.asc"),
      "-spike_free 0.45",
      sep = " ")

paste(file.path("C:","LAStools", "bin", "lasinfo.exe"),
      "-i", file.path(mainDir, "LAS_norm.las"),
      "-last_only",
      "-cd",
      sep = " ")

chm.sf <- raster(file.path(mainDir, "chm", "chm_sf6.asc"))
tree.det <- tree_detection(chm.sf, ws = 17, 1)
####

