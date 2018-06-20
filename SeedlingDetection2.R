# Set directories
mainDir <- "S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/LAS"
setwd(mainDir)

# Fusion
# Merge tiles for seedling area
#system(paste(file.path("C:","Fusion", "MergeData.exe"),
#             "/index",
#             file.path(mainDir, "TileList_clip.txt"),
#             file.path(mainDir, "MoscowMtn_clip.las"),
#             sep=" "))

# Fusion
# Clip tiles to seedling area
system(paste(file.path("C:","Fusion", "clipdata.exe"),
             "/shape:0 /index",
             file.path(mainDir, "TileList_clip.txt"),
             file.path(mainDir, "MoscowMtn_clip2.las"),
             "[2380665 1893999 2381231 1894722]",
             sep=" "))

# Fusion
# Produce desriptive report of lidar set
system(paste(file.path("C:","Fusion", "catalog.exe"),
             "/coverage /density:0.25,5,10",
             paste("/projection:", file.path(mainDir, "MoscowMtn_clip.prj"), sep = ""),
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))

# Fusion
# Approximate the ground's surface (bare-earth points)
system(paste(file.path("C:","Fusion", "groundfilter.exe"),
             "/gparam:0 /wparam:1 /tolerance:1 /iterations:10",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts.las"),
             0.5,
             file.path(mainDir, "MoscowMtn_clip.las"),
             sep=" "))

# Fusion
# Compute the elevation of each grid cell using the average elevation of all 
# points within the cell
system(paste(file.path("C:","Fusion", "gridsurfacecreate.exe"),
             "/smooth:5", 
             file.path(mainDir, "Ground", "MoscowMtn_DEM.dtm"),
             "0.5 F F 2 11 2 0",
             file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts.las"),
             sep=" "))

# Fusion
# Convert the DEM ifrom DTM into ASCII raster
system(paste(file.path("C:","Fusion", "dtm2ascii.exe"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM.dtm"),
             file.path(mainDir, "Ground", "MoscowMtn_DEM.asc"),
             sep=" "))

# Fusion
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

# Import canopy height model
chm <- raster(file.path(mainDir, "chm", "MoscowMtn_clip_CHM.asc"))
plot(chm)
summary(chm)

# Smooth canopy (set parameters and execute)
ws <- 3 # 3x3 window size
filter <- "Gaussian" # Gaussian filter type
sigma <- 0.5
sCHM <- CHMsmoothing(chm, filter, ws, sigma)
plot(sCHM)

# Individual tree detection list (set parameters and execute)
fws <- 5 # 3x3 fixed window size
minht <- 0.5 # 0.5 ft. minimum height above ground
loc <- FindTreesCHM(sCHM, fws, minht)
loc.seedling5 <- subset(loc, height < 15)
coords <- loc.seedling5[,1:2]
data <- as.data.frame(loc.seedling5[,3])
trees5 <- SpatialPointsDataFrame(coords, data = data, proj4string = CRS("+init=epsg:6453"))
plot(trees5)

