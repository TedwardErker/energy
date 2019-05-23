library(dplyr)
library(rgeos)
library(rgdal)
 library(raster)
  library(doParallel)
  library(foreach)


    build.res.primary.parcels <- readOGR(dsn = "../DD/", layer = "buildings.w.parcels.energy.propinfo")
    colnames(build.res.primary.parcels@data) <- readRDS("../DD/colnames.buildings.w.parcels.energy.propinfo.rds")
    build.res.primary.parcels@data$Parcel <- readRDS("../DD/parcelID.buildings.w.parcels.energy.propinfo.rds")
    build.res.primary.parcels@data$BUILDINGINFO <- readRDS("../DD/BUILDINGINFO.buildings.w.parcels.energy.propinfo.rds")

  build.res.primary.parcels.reduced.data <- build.res.primary.parcels
  build.res.primary.parcels.reduced.data@data <-   build.res.primary.parcels.reduced.data@data %>% dplyr::select(BUILDINGFO)


    build.regions <- shapefile("../DD/building.regions.3m16angles.shp")

build.regions@data$region.area <- round(gArea(build.regions, byid = T),2)

cores <- 44

chunks <- 5000
a.s <- split(1:length(build.regions), f = rep(1:chunks, each = ceiling(length(build.regions)/chunks)))

build.regions.list <- lapply(a.s, function(i) build.regions[i,])

dir.create("../DD/building.around.building.3m16angles/")

cl <- makeCluster(cores)
registerDoParallel(cl)

out <- foreach(build.region = build.regions.list[4702], .packages = c("plyr","dplyr","magrittr","sp","raster","rgeos"), .combine = "rbind") %dopar% {

    ## in order to get build.regions.list[4702] to work...
     br <- gBuffer(build.region,width = .0001, byid = T)
     d <- build.region@data
     build.region <- as(br, "SpatialPolygonsDataFrame")
     build.region@data <- d

    o <- raster::intersect(build.region, build.res.primary.parcels.reduced.data)
    o@data$area <- round(gArea(o, byid = T),2)

    df <- left_join(build.region@data, o@data) %>%
        dplyr::mutate(area = ifelse(is.na(area), 0, area)) %>%
        ungroup() %>%
        group_by(BUILDIN, dstnc__, ovr_bld, directn, region.area) %>%
        dplyr::summarize(area = sum(area))

    saveRDS(df, paste0("../DD/building.around.building.3m16angles/",build.region@data$BUILDIN[1],".rds"))
df



}

closeAllConnections()
