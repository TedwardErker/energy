library(raster)
  library(doParallel)
  library(foreach)

  build.regions <- shapefile("../DD/building.regions.3m16angles.shp")

#  build.regions <- shapefile("/Users/erker/g/projects/energy/DD/building.regions.3m16angles_utcProj_shps/1.shp")

  dir.create("../DD/building.regions.3m16angles.sections")
  chunks <- 5000
  a.s <- split(1:length(build.regions), f = rep(1:chunks, each = ceiling(length(build.regions)/chunks)))

  build.regions.list <- lapply(a.s, function(i) build.regions[i,])

  cl <- makeCluster(20)
  registerDoParallel(cl)

  foreach(build.region = build.regions.list, .packages = c("sp","raster")) %dopar% {
      shapefile(build.region, paste0("../DD/building.regions.3m16angles.sections/building.regions_subset",build.region@data$BUILDIN[1]), overwrite = T)
      }
