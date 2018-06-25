# Set directories
mainDir <- "S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS"
setwd(mainDir)

#        #
# Fusion #
#        #

# ClipData
# Clip tiles to seedling area
system(paste(file.path("C:","Fusion", "clipdata.exe"),
             "/shape:0 /index",
             file.path(mainDir, "TileList_clip.txt"),
             file.path(mainDir, "MoscowMtn_clip.las"),
             "2380665 1893999 2381231 1894722",
             sep=" "))

# Catalog
# Produce desriptive report of lidar set
system(paste(file.path("C:","Fusion", "catalog.exe"),
             "/density:0.125,5,10",
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))

# GroundFilter
# Approximate the ground's surface (bare-earth points)
system(paste(file.path("C:","Fusion", "groundfilter.exe"),
             "/gparam:0 /wparam:1 /iterations:10",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts.las"),
             0.25,
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))

# GridSurfaceCreate
# Compute the elevation of each grid cell using the average elevation of all 
# points within the cell
system(paste(file.path("C:","Fusion", "gridsurfacecreate.exe"),
             "/smooth:5", 
             file.path(mainDir, "Ground", "MoscowMtn_DEM.dtm"),
             "0.25 F F 2 11 2 0",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts.las"),
             sep=" "))

# DTM2ASCII
# Convert the DEM from DTM into ASCII raster
system(paste(file.path("C:","Fusion", "dtm2ascii.exe"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM.dtm"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM.asc"),
             sep=" "))

# CanopyModel
# Create a canopy surface model using a LIDAR point cloud by CanopyModel assigning
# the elevation of the highest return within each grid cell to the grid cell center
system(paste(file.path("C:","Fusion", "canopymodel.exe"),
             "/ascii /outlier:-1,150",
             paste("/ground:",file.path(mainDir, "Ground", "MoscowMtn_DEM.dtm"), sep=""),
             file.path(mainDir, "chm", "MoscowMtn_clip_CHM.dtm"),
             "0.5 F F 2 11 2 0",
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))

# Load packages
library(rLiDAR)
library(raster)
library(rgeos)
library(rgdal)

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
chm <- grid_canopy(las.filt, res = 0.25, na.fill = "knnidw", k = 10, p = 2)
chm2 <- as.raster(chm)
writeRaster(chm2, "S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/chm/chm.tif", overwrite = T)

tree.det <- tree_detection(chm, ws = 17, 1) #18 inch radius, 12 inch minimum seedling 
writeRaster(tree.det, "S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS/chm/tree_det.tif", overwrite = T)

tree.seg <- lastrees_silva(las, chm2, tree.det, )