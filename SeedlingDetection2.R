#=======================================================================#
#=======================================================================#
#             Using UAV lidar to identify conifer seedlings             #
#                         NRS 504 Final Project                         #
#                            Anthony Martinez                           #
#=======================================================================#
#=======================================================================#

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
library(rLiDAR)
library(lidR)
library(raster)
library(rgeos)
library(rgdal)
library(data.table)

# Import LAS
las <- lidR::readLAS(file.path(mainDir, "MoscowMtn_clip.las"))

# Remove flight lines
hist(las@data$Z) # identify cutoff threshold
las <- lasfilter(las, Z < 3050)

# Identify ground points (takes approx. 20 hours to run on this dataset)
ws <- c(1,2,4,8,16,32,64,128)
th <- seq(0.1, 2, length.out = length(ws))
las.ground <- lasground(las, "pmf", ws, th)
writeLAS(las.ground, file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts.las"))
dem.grid <- grid_terrain(las.ground, res = 0.25, method = "knnidw", k = 10, p = 2, 
                        keep_lowest = FALSE)
dem.raster <- as.raster(dem.grid)
writeRaster(dem.raster, file.path(mainDir, "Ground", "pmf_dem2.tif"))

#=======================================================================#
# Create canopy model -- 0.25 ft (3 in) resolution
#=======================================================================#
# Create "spike free" canopy height model (with LAStools)
system(paste(file.path("C:","LAStools", "bin", "las2dem.exe"),  # CHM with normalized LAS
             "-i", file.path(mainDir, "LAS_pwr_rm.las"),
             "-spike_free 0.45", # freeze distance: ~ 3x the average pulse spacing
             "-step 0.25",
             "-o",  file.path(mainDir, "chm", "CHM_spike_free.asc"),
             sep=" "))
chm.sf <- raster(file.path(mainDir, "chm", "CHM_spike_free.asc"))
chm.sf <- CHMsmoothing(chm.sf, filter = "Gaussian", ws = 3, sigma = 0.5)
writeRaster(chm.sf, file.path(mainDir, "chm", "CHM_spike_free_sm.asc"))
chm.sf <- raster(file.path(mainDir, "chm", "CHM_spike_free_sm.asc"))

#=======================================================================#
# Locate trees
#=======================================================================#
# Identify tree tops -- 6 foot diameter, 12 inch minimum seedling
tree.top.ras <- tree_detection(chm.sf, ws = 25, hmin = 1)
tree.top.pts <- rasterToPoints(tree.top.ras, spatial = T)
tree.top <- as.data.frame(rasterToPoints(tree.top.ras, spatial = F))
colnames(tree.top) <- c("x", "y", "id")
shapefile(tree.top.pts, filename = file.path(mainDir, "Trees", "tree_tops.shp"), 
          overwrite = T)

# Classify trees in point cloud
#lastrees_silva(las.pwr.rm, chm.sf, tree.top, max_cr_factor = 0.8)
#writeLAS(las.pwr.rm, file.path(mainDir, "LAS_pwr_rm_class.las"))

# Determine tree heights
chm.norm <- chm.sf - dem.raster
writeRaster(chm.norm, file.path(mainDir, "chm", "chm_norm.tif"), overwrite = T)
chm.norm <- raster(file.path(mainDir, "chm", "chm_norm.tif"))
tree.top.pts@data$ht <- extract(chm.norm, tree.top.pts)

############################################################################
#  Accuracy assessment
############################################################################
# Import sampled trees and tree-absences
tree.abs <- readOGR("S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/GIS/RandomPoints.shp")
tree.samp <- readOGR("S:/COS/PyroGeog/amartinez/UAV_seedlings/Lidar/GIS/SampledTrees.shp", stringsAsFactors = FALSE)
proj4string(tree.top.pts) <- proj4string(tree.samp)

SapID <- function(x, y, z){
  buff <- gBuffer(y, byid = T, width = z)
  a <- gContains(buff, tree.top.pts, byid = T)
  b <- as.data.frame(which(a, arr.ind = T))
  colnames(b) <- c("tree.top.id", "sample.no")
  b$sample.name <- y$OBJNAME[b[,2]]
  b$tree.top.ht <- x@data$ht[b$tree.top.id]
  b$sample.ht <- y@data$Ht__in_[b$sample.no]/12
  b$diff <- abs(b$tree.top.ht - b$sample.ht)
  b <- data.table(b)
  b <- b[ , .SD[which.min(diff)], by = tree.top.id]
  b <- b[ , .SD[which.min(diff)], by = sample.name]
  b$no <- z
  buff.abs <- gBuffer(tree.abs, byid = T, width = z)
  d <- gContains(buff.abs, tree.top.pts, byid = T)
  e <- as.data.frame(which(d, arr.ind = T))
  colnames(e) <- c("tree.top.id", "abs.sample.no")
  acc <- c(nrow(e), nrow(b))
  names(acc) <- c("False.Pos", "True.Pos")
  assign(paste0("tree.samp.", z), tree.samp[unique(b$sample.no),], envir = .GlobalEnv)
  assign(paste0("tree.top.", z), x[b$tree.top.id,], envir = .GlobalEnv)
  assign(paste0("tree.abs.", z), tree.abs[unique(e$abs.sample.no),], envir = .GlobalEnv)
  assign(paste0("tree.id.", z), data.frame(b), envir = .GlobalEnv)
  assign(paste0("accuracy.", z), acc)
  paste0("created objects: ", "tree.samp.", z, ", tree.top.", z, ", tree.abs.", z, " tree.id.", z, ", & accuracy.", z)
}

SapID(tree.top.pts, tree.samp, 1)
SapID(tree.top.pts, tree.samp, 2)
SapID(tree.top.pts, tree.samp, 3)
SapID(tree.top.pts, tree.samp, 4)
SapID(tree.top.pts, tree.samp, 5)
shapefile(tree.samp.5, file.path(mainDir, "Trees", "Tree_samp_5.shp"))


############################################################################
#  Visualization and figures
############################################################################
sap.ht <- tree.top.pts@data$ht[tree.top.pts@data$ht < 9]
hist(sap.ht, xlab = "Sapling height (ft.)", main = "")
summary(sap.ht)
tree.id.4
boxplot(tree.id.1$diff, tree.id.2$diff, tree.id.3$diff, tree.id.4$diff, tree.id.5$diff,
        xlab = "Buffer radius (ft.)", ylab = "Height difference (ft.)", names = 1:5)

a <- tree.id.5
a$no[which(a$tree.top.id %in% Reduce(intersect, list(tree.id.5$tree.top.id,tree.id.3$tree.top.id)))] <- 4

for(i in rev(1:5)){
  a$no[which(a$tree.top.id %in% 
               Reduce(intersect, list(tree.id.5$tree.top.id, 
                                      get(paste0("tree.id.", i))$tree.top.id)))] <- 6-i
}
as.numeric(rownames(tree.top.5@data))
tree.top.5@data$no <- a$no
shapefile(tree.top.5, file.path(mainDir, "Trees", "Tree_top_5.shp"), overwrite = T)

hist(las@data$Z, breaks = 50, col = "gray40", border = "gray40", yaxt='n', xaxt = 'n', ann = F)

############################################################################
#                                  END                                     #
############################################################################
#import
las.pwr.rm <- lidR::readLAS(file.path(mainDir, "LAS_pwr_rm_class.las"))
las.pwr.rm <- lidR::readLAS(file.path(mainDir, "LAS_pwr_rm_class.las"))
dem.raster <- raster(file.path(mainDir, "Ground", "pmf_dem.tif"))
las.ground <- lidR::readLAS(file.path(mainDir, "Ground", "MoscowMtn_clip_groundpts.las"), 
                            filter = "-drop_z_above 3050")
dem.raster <- raster(file.path(mainDir, "Ground", "pmf_dem.tif"))


x = tree.top.pts
y = tree.samp
z = 5

