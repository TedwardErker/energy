library(plyr)
library(ascii)
library(broom)
library(tidyr)
library(stringr)
library(raster)
library(rgeos)
library(rgdal)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(GGally)

1267.5 / 2.20462 / 1000 *12 /44

kWh2kgC <- function(kWh) {
    kWh * .1567988
}

.1*14.46

therm2kgC <- function(therms) {
    therms * 1.446
}

extract.polygons.parallel <- function(x, sp.polygons, cores) {

    n <- ceiling(length(sp.polygons) / cores)
    v <- 1:length(sp.polygons)
    l <- split(v, ceiling(seq_along(v)/n))

          cl <- makeCluster(cores)
          registerDoParallel(cl)

  out <- foreach(i = seq_along(l), .packages = "raster") %dopar% {
        sp.sub <- sp.polygons[l[[i]],]
        raster::extract(x, sp.sub)
    }
        closeAllConnections()

return(out)
  }

options(asciiType = "org")
ascii.nowarn.print <- function(x,...) {
                                        #op <- options(warn = -1)
                                        #      on.exit(options(op))

    suppressWarnings(print(ascii(x,...)))

}

makeRegionsAroundPolygon <- function(poly, buffer.widths, theta = c(-67.5,-22.5,22.5,67.5)) {
    bc <- gCentroid(poly)
    bcb <- sapply(buffer.widths, function(w) gBuffer(poly, width = w, byid = T))
    bcbd <- sapply(length(bcb):1, function(i) {
        if(i >1) {
            gDifference(bcb[[i]], bcb[[(i-1)]])
        } else {
            bcb[[i]]
        }
    })
    bcb <- do.call(bind, bcbd)

    theta.radians <- theta*pi/180
    r <- max(buffer.widths) + sqrt(gArea(poly)/pi) * 3

    y <- r*sin(theta.radians)
    x <- r*cos(theta.radians)

    x1 <- coordinates(bc)[,1 ] + x
    y1 <- coordinates(bc)[,2] + y

    x2 <- coordinates(bc)[,1 ] - x
    y2 <- coordinates(bc)[,2] - y

    line = SpatialLines(list(Lines(list(Line(cbind(c(x1,x2),c(y1,y2)))), ID="line")))

    line = SpatialLines(sapply(1:length(x1), function(i) list(Lines(Line(cbind(c(x1[i],x2[i]),c(y1[i],y2[i]))), ID=i))))
    proj4string(line) <- crs(poly)


    o <- lapply(1:length(bcb), function(i) {
        lpi <- gIntersection(bcb[i,], line)
        blpi <- gBuffer(lpi, width = 0.0001)  # create a very thin polygon buffer of the intersected line
        dpi <- gDifference(bcb[i,], blpi)              # split using gDifference
        disaggregate(dpi)
    })

    o <- do.call(bind, o)



    dist <- round(gDistance(poly, o, byid = T),-1)[,1]
    over.building <- 0 == round(gDistance(bc, o, byid = T),0)[,1]

    direction <- sapply(1:length(o), function(i) {
        diff.coords <- coordinates(bc) - coordinates(gCentroid(o[i,]))
        round(atan2(diff.coords[1], diff.coords[2]) * 180/pi, 0)
    })

    angles <- c(0, 45, 90, 135, 180, -180, -135, -90, -45)
    closest.angles <- sapply(direction, function(dir) which(abs(angles - dir) == min(abs(angles - dir))))
    direction <- angles[closest.angles]
    direction <- mapvalues(direction, from = angles, to = c("s","sw","w","nw","n","n","ne","e","se"))
    o <- SpatialPolygonsDataFrame(o, data = data.frame(BUILDINGFO = poly@data$BUILDINGFO,
                                                       distance.from.building = dist,
                                                       over.building = over.building,
                                                       direction = direction))
    o
}

poly <- readWKT("POLYGON ((78.6 41.7, 91.3 41.7, 91.3 26.6, 78.7 26.6, 78.6 41.7))")

  poly <- as(poly, "SpatialPolygonsDataFrame")
poly@data$BUILDINGFO = 0


poly <- build.res.primary.parcels[27940,]
buffer.widths <- seq(0,51,3)
out <- makeRegionsAroundPolygonManyAngles(poly, buffer.widths)

buffer.widths <- seq(0,51,3)
  theta = seq(0,157.5,22.5)

      makeRegionsAroundPolygonManyAngles <- function(poly, buffer.widths, theta = seq(0,157.5,22.5)) {

      bc <- gCentroid(poly)
      bcb <- sapply(buffer.widths, function(w) gBuffer(poly, width = w, byid = T))
      bcbd <- sapply(length(bcb):1, function(i) {
          if(i >1) {
              gDifference(bcb[[i]], bcb[[(i-1)]])
          } else {
              bcb[[i]]
          }
      })
      bcb <- do.call(bind, bcbd)

      theta.radians <- theta*pi/180
      r <- max(buffer.widths) + sqrt(gArea(poly)/pi) * 3

      y <- r*sin(theta.radians)
      x <- r*cos(theta.radians)

      x1 <- coordinates(bc)[,1 ] + x
      y1 <- coordinates(bc)[,2] + y

      x2 <- coordinates(bc)[,1 ] - x
      y2 <- coordinates(bc)[,2] - y

      line = SpatialLines(list(Lines(list(Line(cbind(c(x1,x2),c(y1,y2)))), ID="line")))

      line = SpatialLines(sapply(1:length(x1), function(i) list(Lines(Line(cbind(c(x1[i],x2[i]),c(y1[i],y2[i]))), ID=i))))
      proj4string(line) <- crs(poly)


      o <- lapply(1:length(bcb), function(i) {
          lpi <- gIntersection(bcb[i,], line)
          blpi <- gBuffer(lpi, width = 0.001)  # create a very thin polygon buffer of the intersected line
          dpi <- gDifference(bcb[i,], blpi)              # split using gDifference
          disaggregate(dpi)
      })

      o <- do.call(bind, o)

      dist <- gDistance(poly, o, byid = T)
          dist <- sapply(dist, function(d) buffer.widths[which(abs(buffer.widths - d) == min(abs(buffer.widths - d)))])

#This over.building code also seems to make errors because some of the regions are quite close.
      over.building <- 0 == round(gDistance(bc, o, byid = T),0)[,1]

      direction <- sapply(1:length(o), function(i) {
          diff.coords <- coordinates(bc) - coordinates(gCentroid(o[i,]))
          atan2(diff.coords[1], diff.coords[2]) * 180/pi
      })

  direction[direction < 0] <- 360 + direction[direction < 0]  # get rid of negative angles, but on 0-360 scale

  angles <- seq(11.25, 348.75, 22.5)

      closest.angles <- sapply(direction, function(dir) which(abs(angles - dir) == min(abs(angles - dir))))

          # shifting so that angles originate in north, grow in a clockwise direction
          direction <- angles[closest.angles] + 180
          direction[direction > 360] <-         direction[direction > 360] - 360

  direction <- as.character(direction)

      o <- SpatialPolygonsDataFrame(o, data = data.frame(BUILDINGFO = poly@data$BUILDINGFO,
                                                         distance.from.building = dist,
                                                         over.building = over.building,
                                                         direction = direction))
      o
  }

makeNESWRegionsAroundPolygon <- function(poly, buffer.widths, theta = c(45,-45)) {
    bc <- gCentroid(poly)
    bcb <- sapply(buffer.widths, function(w) gBuffer(poly, width = w, byid = T))
    bcbd <- sapply(length(bcb):1, function(i) {
        if(i >1) {
            gDifference(bcb[[i]], bcb[[(i-1)]])
        } else {
            bcb[[i]]
        }
    })
    bcb <- do.call(bind, bcbd)

    theta.radians <- theta*pi/180
    r <- max(buffer.widths) + sqrt(gArea(poly)/pi) * 3

    y <- r*sin(theta.radians)
    x <- r*cos(theta.radians)

    x1 <- coordinates(bc)[,1 ] + x
    y1 <- coordinates(bc)[,2] + y

    x2 <- coordinates(bc)[,1 ] - x
    y2 <- coordinates(bc)[,2] - y

    line = SpatialLines(list(Lines(list(Line(cbind(c(x1,x2),c(y1,y2)))), ID="line")))

    line = SpatialLines(sapply(1:length(x1), function(i) list(Lines(Line(cbind(c(x1[i],x2[i]),c(y1[i],y2[i]))), ID=i))))
    proj4string(line) <- crs(poly)


    o <- lapply(1:length(bcb), function(i) {
        lpi <- gIntersection(bcb[i,], line)
        blpi <- gBuffer(lpi, width = 0.0001)  # create a very thin polygon buffer of the intersected line
        dpi <- gDifference(bcb[i,], blpi)              # split using gDifference
        disaggregate(dpi)
    })

    o <- do.call(bind, o)


    dist <- gDistance(poly, o, byid = T)[,1]
    closest.dist <- sapply(dist, function(d) which(abs(buffer.widths - d) == min(abs(buffer.widths - d))))
    dist <- buffer.widths[closest.dist]
    over.building <- 0 == round(gDistance(bc, o, byid = T),0)[,1]

    direction <- sapply(1:length(o), function(i) {
        diff.coords <- coordinates(bc) - coordinates(gCentroid(o[i,]))
        atan2(diff.coords[1], diff.coords[2]) * 180/pi
    })

    angles <- c(0, 90, 180, -180, -90)
    closest.angles <- sapply(direction, function(dir) which(abs(angles - dir) == min(abs(angles - dir))))
    direction <- angles[closest.angles]
                                        #    direction <- mapvalues(direction, from = angles, to = c("s","sw","w","nw","n","n","ne","e","se"))
        direction <- mapvalues(direction, from = angles, to = c("s","w","n","n","e"))
    o <- SpatialPolygonsDataFrame(o, data = data.frame(BUILDINGFO = poly@data$BUILDINGFO,
                                                       distance.from.building = dist,
                                                       over.building = over.building,
                                                       direction = direction))
    o
}

makeRegionsAroundPolygonCentroid <- function(polygon, buffer.widths, theta) {

    bc <- gCentroid(polygon)
    bcb <- sapply(buffer.widths, function(w) gBuffer(bc, width = w, byid = T))
    bcbd <- sapply(length(bcb):1, function(i) {
        if(i >1) {
            gDifference(bcb[[i]], bcb[[(i-1)]])
        } else {
            bcb[[i]]
        }
    })
    bcb <- do.call(bind, bcbd)

    theta.radians <- theta*pi/180
    r <- max(widths)

    y <- r*sin(theta.radians)
    x <- r*cos(theta.radians)

    x1 <- coordinates(bc)[,1 ] + x
    y1 <- coordinates(bc)[,2] + y

    x2 <- coordinates(bc)[,1 ] - x
    y2 <- coordinates(bc)[,2] - y

    line = SpatialLines(list(Lines(list(Line(cbind(c(x1,x2),c(y1,y2)))), ID="line")))

    line = SpatialLines(sapply(1:length(x1), function(i) list(Lines(Line(cbind(c(x1[i],x2[i]),c(y1[i],y2[i]))), ID=i))))
    proj4string(line) <- crs(polygon)


    o <- lapply(1:length(bcb), function(i) {
        lpi <- gIntersection(bcb[i,], line)
        blpi <- gBuffer(lpi, width = 0.000001)  # create a very thin polygon buffer of the intersected line
        dpi <- gDifference(bcb[i,], blpi)              # split using gDifference
        disaggregate(dpi)
    })

    o <- do.call(bind, o)
    o
}

build.shapefile.path <- "../RD/dane_cty_building_footprint_2010/BuildingFootprint.shp"
parcels.shapefile.path <- "../RD/dane_cty_parcels_2014/Dane_Parcels_2014.shp"
energy.data.frame.path <- "../RD/energy.dataframe-2016-04-12-15-40-59.rds"
assessors.property.path <- "../RD/Assessor_Property_Information.csv"
mad.utc.path <- "../RD/madison_utc/3_ClassifiedUrbanArea.tif"

energy <- readRDS(energy.data.frame.path)

ass.prop <- read.csv(assessors.property.path, stringsAsFactors = F)

ass.prop.energy <- left_join(energy, ass.prop) %>%
  mutate(PropertyAd = toupper(Address))

parcels <- shapefile(parcels.shapefile.path)

parcels.energy <- subset(parcels, PropertyAd %in% ass.prop.energy$PropertyAd)

parcels.energy@data <- left_join(parcels.energy@data, ass.prop.energy, by = c("PropertyAd"))

parcels.energy <- subset(parcels.energy, complete.cases(kWh_High) & kWh_High >0)


parcels.energy@data <- parcels.energy@data %>%
          dplyr::select(Parcel,
                      Address,
                 kWh_High,
                 kWh_Mon_Avg_last12mo,
                 Cost_kWh_High,
                 Cost_kWh_Avg,
                 Therms_High,
                 Therms_Mon_Avg_last12mo,
                 Cost_Therms_High,
                 Cost_Therms_Avg,
                 Electric_for_Heating,
                 Current.Year.Land.Value,
                 Current.Year.Improvement.Value,
                 Lot.Size.Sq.Ft,
                 Water.Frontage,
                 Elementary.School,
                 Middle.School,
                 High.School,
                 Ward,
                 State.Assembly.District,
                 Home.Style,
                 Year.Built,
                 Stories,
                 Dwelling.Units,
                 Bedrooms,
                 Full.Baths,
                 Half.Baths,
                 Fireplaces,
                 Central.Air,
                 First.Floor.Living.Area,
                 Second.Floor.Living.Area,
                 Third.Floor.Living.Area,
                 Above.Third.Floor.Living.Area,
                 Total.Living.Area,
                 Finished.Attic,
                 Finished.Basement,
                 Total.Basement,
                 Crawl.Space,
                 Exterior.Wall.1,
                 Exterior.Wall.2,
                 Foundation,
                 Roof,
                 Roof.Replaced.Year,
                 Garage.1,
                 Garage.Stalls.1,
                 Garage2,
                 Garage.Stalls.2,
                 Driveway) %>%
          mutate(Electric_for_Heating = factor(Electric_for_Heating),
                 Water.Frontage = factor(Water.Frontage),
                 Elementary.School = factor(Elementary.School),
                 Middle.School = factor(Middle.School),
                 High.School = factor(High.School),
                 Ward = factor(Ward),
                 State.Assembly.District = factor(State.Assembly.District),
                 Home.Style = factor(Home.Style),
                 Central.Air = factor(Central.Air),
                 Exterior.Wall.1 = factor(Exterior.Wall.1),
                 Exterior.Wall.2 = factor(Exterior.Wall.2),
                 Foundation = factor(Foundation),
                 Roof = factor(Roof),
                 Garage.1 = factor(Garage.1),
                 Garage2 = factor(Garage2),
                 Driveway = factor(Driveway)) %>%
          mutate(Current.Year.Improvement.Value.th = as.numeric(str_replace(Current.Year.Improvement.Value,"\\$",""))/1000,
                 Current.Year.Land.Value.th = as.numeric(str_replace(Current.Year.Land.Value,"\\$",""))/1000) %>%
          dplyr::select(-Current.Year.Improvement.Value, -Current.Year.Land.Value)

writeOGR(parcels.energy, dsn = "../DD/", layer = "parcels.energy.propinfo", driver = "ESRI Shapefile", overwrite = T)
saveRDS(colnames(parcels.energy@data), "../DD/colnames.parcels.energy.propinfo.rds")
saveRDS(parcels.energy$Parcel, "../DD/parcelID.parcels.energy.propinfo.rds")

parcels.energy <- shapefile("../DD/parcels.energy.propinfo.shp")
parcels.energy$Parcel <- readRDS("../DD/parcelID.parcels.energy.propinfo.rds")
colnames(parcels.energy@data) <- readRDS("../DD/colnames.parcels.energy.propinfo.rds")
  parcels.energy <- spTransform(parcels.energy, crs("+init=epsg:32616"))

build <- shapefile(build.shapefile.path)
    build.res <- subset(build, BuildingUs == "Residential")
    build.res.primary <- subset(build.res, BuildingCl == "Primary")
  build.res.primary <- spTransform(build.res.primary, crs("+init=epsg:32616"))
writeOGR(build.res.primary, dsn = "../DD/", layer = "buildings_residential_primary", driver = "ESRI Shapefile")

build.res.primary <- shapefile("../DD/buildings_residential_primary.shp")

e <- new("Extent"
         , xmin = 294022.535226926
         , xmax = 317479.19202743
         , ymin = 4763205.41481759
         , ymax = 4781697.04891745
           )

  e.parcels.energy <- crop(parcels.energy, e)

  library(alphahull)
  library(igraph)

  g <- geom(e.parcels.energy)
  i <- sample(1:nrow(g), 20000)

  a <- ashape(unique(g[i,5:6]), alpha = 2000)

  ag = graph.edgelist(cbind(as.character(a$edges[, "ind1"]), as.character(a$edges[,
      "ind2"])), directed = FALSE)
  plot(ag)

  if (!is.connected(ag)) {
      stop("Graph not connected")
  }
  if (any(degree(ag) != 2)) {
      stop("Graph not circular")
  }
  if (clusters(ag)$no > 1) {
      stop("Graph composed of more than one circle")
  }

  cutg = ag - E(ag)[1]
  # find chain end points
  ends = names(which(degree(cutg) == 1))
  path = get.shortest.paths(cutg, ends[1], ends[2])[[1]]
  # this is an index into the points
  pathX = as.numeric(V(ag)[path[[1]]]$name)
  # join the ends
  pathX = c(pathX, pathX[1])
  # now show the alpha shape plot with our poly on top
  plot(a, lwd = 10, col = "gray")
  # get the points from the ashape object
  lines(a$x[pathX, ], lwd = 2)


  a.sp <- SpatialPolygonsDataFrame(SpatialPolygons(list(Polygons(list(Polygon(a$x[pathX,])),1))), data = data.frame(d = 1))

  a.sp.buf <- gBuffer(a.sp, width = 100)

build.res.primary.crp <- crop(build.res.primary, a.sp.buf)

writeOGR(build.res.primary.crp, dsn = "../DD/", layer = "buildings_residential_primary_crp2Mad", driver = "ESRI Shapefile")

build.res.primary.crp <- shapefile("../DD/buildings_residential_primary_crp2Mad.shp")

parallel.crop <- function(sp1,sp2, cores) {
    a <- 1:length(sp1)
    a.s <- split(a, rep(1:cores,each = ceiling(length(a)/cores)))

    cl <- makeCluster(cores)
    registerDoParallel(cl)

    out <- foreach(i = a.s, .packages = c("sp","raster")) %dopar% {
        o <- crop(sp1[i,], sp2)
    }
    closeAllConnections()
    return(out)
}

out <- parallel.crop(build.res.primary.crp, parcels.energy, 40)

build.res.primary.parcels <- do.call("rbind",out[sapply(out, FUN = function(o) !is.null(o))])

build.res.primary.parcels <- spTransform(build.res.primary.parcels, crs("+init=epsg:32616"))

                                        #remove the fragments of buildings
o.area <- data.frame(id = build.res.primary.crp@data$BUILDINGFO, o.area = gArea(build.res.primary.crp, byid = T))
n.area <- data.frame(id = build.res.primary.parcels@data$BUILDINGFO, n.area =gArea(build.res.primary.parcels, byid = T))

compare.area <- left_join(n.area, o.area) %>%
    mutate(area.ratio = n.area / o.area,
           keep = area.ratio > .6)

build.res.primary.parcels <- build.res.primary.parcels[compare.area$keep,]





df  <- over(build.res.primary.parcels, parcels.energy)
build.res.primary.parcels@data <- cbind(build.res.primary.parcels@data, df)


                                        # If there are two buildings in a parcel, remove the smaller of the two

parcels.w2buidings <- build.res.primary.parcels@data %>%
    group_by(Parcel) %>%
    summarize(n = n()) %>%
    filter(n > 1) %>%
    pull(Parcel)

build.to.remove <- build.res.primary.parcels@data %>%
    filter(Parcel %in% parcels.w2buidings) %>%
    group_by(Parcel) %>%
    dplyr::select(BUILDINGFO, Shape_area) %>%
    arrange(Parcel) %>%
    filter(Shape_area == min(Shape_area)) %>%
    pull(BUILDINGFO)


build.res.primary.parcels <- build.res.primary.parcels[ ! build.res.primary.parcels@data$BUILDINGFO %in% build.to.remove,]

                                        # If there are buildings that are in multiple parcels??

writeOGR(build.res.primary.parcels, dsn = "../DD/", layer = "buildings.w.parcels.energy.propinfo", driver = "ESRI Shapefile", overwrite = T)

saveRDS(colnames(build.res.primary.parcels@data), "../DD/colnames.buildings.w.parcels.energy.propinfo.rds")
saveRDS(build.res.primary.parcels@data$Parcel, "../DD/parcelID.buildings.w.parcels.energy.propinfo.rds")
saveRDS(build.res.primary.parcels@data$BUILDINGINFO, "../DD/BUILDINGINFO.buildings.w.parcels.energy.propinfo.rds")
saveRDS(build.res.primary.parcels@data, "../DD/data.buildings.w.parcels.energy.propinfo.rds")

build.res.primary.parcels <- readOGR(dsn = "../DD/", layer = "buildings.w.parcels.energy.propinfo")
colnames(build.res.primary.parcels@data) <- readRDS("../DD/colnames.buildings.w.parcels.energy.propinfo.rds")
build.res.primary.parcels@data$Parcel <- readRDS("../DD/parcelID.buildings.w.parcels.energy.propinfo.rds")
build.res.primary.parcels@data$BUILDINGINFO <- readRDS("../DD/BUILDINGINFO.buildings.w.parcels.energy.propinfo.rds")

build.res.primary.parcels <- bind(build.res.primary.parcels, gBuffer(build.res.primary.parcels[27940,], width = .001))

build.res.primary.parcels <- build.res.primary.parcels[-27940,]

cores <- 45
      widths <- seq(0,51,3)

  chunks <- 80
a <- 1:length(build.res.primary.parcels)
a.s <- split(a, rep(1:chunks,each = ceiling(length(a)/chunks)))

      dir.create("../DD/buildings.regions")

      cl <- makeCluster(cores)
      registerDoParallel(cl)


      foreach(i = a.s, .packages = c("plyr","sp","raster","rgeos"), .combine = "rbind") %dopar% {
          o <- lapply(i, function(j) {
              print(j)
              op <- makeRegionsAroundPolygonManyAngles(build.res.primary.parcels[j,], widths)
          })
          o <- do.call(bind, o)
          shapefile(o, paste0("../DD/buildings.regions/",i[1],".shp"), overwrite = T)
      }

      closeAllConnections()


    shp.files <- list.files("../DD/buildings.regions", pattern = ".*shp", full.names = T)

      cl <- makeCluster(cores)
      registerDoParallel(cl)

    out <- foreach(shp.file = shp.files, .packages = c("sp","raster","rgeos")) %dopar% {
        shapefile(shp.file)
        }

    out <- do.call(bind, out)

colnames(out@data) <- c("BUILDINGFO","distance.from.building","over.building","direction")

writeOGR(out, dsn = "../DD/", layer = "building.regions.3m16angles", driver = "ESRI Shapefile", overwrite = T)

#  unlink("../DD/buildings.regions.15mNESW", recursive = T)

out <- shapefile("../DD/building.regions.3m16angles.shp")
colnames(out@data) <- c("BUILDINGFO","distance.from.building","over.building","direction")

utc <- raster(mad.utc.path)

out2 <- spTransform(out, proj4string(utc))
shapefile(out2,"../DD/building.regions.3m16angles_utcProj.shp")

out2 <- shapefile("../DD/building.regions.3m16angles_utcProj.shp")
colnames(out2@data) <- c("BUILDINGFO","distance.from.building","over.building","direction")

#out2 <- shapefile("../DD/building.regions.3m16angles_utcProj_shps/1.shp")
#  colnames(out2@data) <- c("BUILDINGFO","distance.from.building","over.building","direction")
#out2 <- out2[1:1000,]

cores <- 42


  dir.create("../DD/extracted.cover3m16angles")

  chunks <- 100
  a.s <- split(1:length(out2), f = rep(1:chunks, each = ceiling(length(out)/chunks)))

  o <-      lapply(a.s, function(i) {
      res.set <- extract.polygons.parallel(utc, out2[i,], cores = cores)
      res.u <- unlist(res.set, recursive = F)
      saveRDS(res.u, paste0("../DD/extracted.cover3m16angles/eachRegionAroundBuilding_set",i[1],".rds"))
  })

res.u.list <- list.files("../DD/extracted.cover3m16angles", full.names = T)

                                        #put them in the correct order for reconstruction
i <- as.numeric(str_extract(res.u.list, "[0-9]+"))
res.u.list <- res.u.list[order(i)]


res.u <- lapply(res.u.list, readRDS)
res.u.u <- unlist(res.u, recursive = F)

sum.tree.per.region <- sapply(res.u.u, function(x) sum(x == 1))
n.pixels.per.region <- sapply(res.u.u, function(x) length(x))

  out2@data <- cbind(out2@data, sum.tree.per.region, n.pixels.per.region)

#fix distances  I should only have to do this once, I need to fix the code above, but I didn't re run because it takes 3 days.
buffer.widths <- seq(0,51,3)
out2@data$distance.from.building <- sapply(out2@data$distance.from.building, function(d) buffer.widths[which(abs(buffer.widths - d) == min(abs(buffer.widths - d)))])






            o <- out2@data %>% unite_("direction.dist.overbuilding", c("direction","distance.from.building", "over.building"))





            # some regions, say south 10m from a house are made or more than one polygon.  This combines them.
            o <- o %>% group_by(BUILDINGFO, direction.dist.overbuilding) %>%
                summarize(sum.tree.per.region = sum(sum.tree.per.region),
                          n.pixels.per.region = sum(n.pixels.per.region))

      saveRDS(o, "../DD/sum.tree_npixels_byregion_atbuildings_3m16angles.rds")


      o <- select(o,  -n.pixels.per.region)

      o <- spread(o, key = direction.dist.overbuilding, value = sum.tree.per.region)

build.res.primary.parcels <- shapefile("../DD/buildings.w.parcels.energy.propinfo", stringsAsFactors = F)

  colnames(build.res.primary.parcels@data) <- readRDS("../DD/colnames.buildings.w.parcels.energy.propinfo.rds")
build.res.primary.parcels@data$BUILDINGFO <- as.numeric(build.res.primary.parcels@data$BUILDINGFO)

build.res.primary.parcels.tree <- build.res.primary.parcels[build.res.primary.parcels@data$BUILDINGFO %in% o$BUILDINGFO, ]
build.res.primary.parcels.tree@data <- left_join(build.res.primary.parcels.tree@data, o)

colnames(build.res.primary.parcels.tree@data)[57:344] <- paste0("d",colnames(build.res.primary.parcels.tree@data)[57:344],"_tree")

shapefile(as(build.res.primary.parcels.tree, "SpatialPolygons"), "../DD/buildings.w.parcels.energy.propinfo_3m16angles", overwrite = T)
saveRDS(colnames(build.res.primary.parcels.tree@data), "../DD/colnames.buildings.w.parcels.energy.propinfo_3m16angles.rds")
saveRDS(build.res.primary.parcels.tree$Parcel, "../DD/Parcel.buildings.w.parcels.energy.propinfo_3m16angles.rds")
saveRDS(build.res.primary.parcels.tree@data, "../DD/data_buildings.w.parcels.energy.propinfo_3m16angles.rds")

res.u.list <- list.files("../DD/building.around.building.3m16angles", full.names = T)


res.u <- lapply(res.u.list, readRDS)

res.u.u <- bind_rows(res.u)

saveRDS(res.u.u, "../DD/building.around.building.3m16angles/buildings.in.regions.around.buildings.rds")

buildings.in.regions <- readRDS("../DD/building.around.building.3m16angles/buildings.in.regions.around.buildings.rds")

                                        # remove over building because it is meaningless

buildings.in.regions <- buildings.in.regions %>% dplyr::filter(ovr_bld == 0)

abuildings.in.regions <- buildings.in.regions %>% unite_("direction.dist.overbuilding", c("directn", "dstnc__", "ovr_bld"))

buildings.in.regions <- abuildings.in.regions %>% mutate(direction.dist.overbuilding = paste0(direction.dist.overbuilding, "_building"))

#remove trouble buildings
buildings.in.regions <- filter(buildings.in.regions, BUILDIN != 190246)

b <- buildings.in.regions %>% select(-region.area)

            # some regions, say south 10m from a house are made or more than one polygon.  This combines them.
            b <- b %>% group_by(BUILDIN, direction.dist.overbuilding) %>%
                summarize(area = sum(area))


bs <- spread(b, key = direction.dist.overbuilding, value = area)

buildings.propinfo.energy.tree <- shapefile("../DD/buildings.w.parcels.energy.propinfo_3m16angles")

colnames(buildings.propinfo.energy.tree@data) <- readRDS("../DD/colnames.buildings.w.parcels.energy.propinfo_3m16angles.rds")

buildings.propinfo.energy.tree@data <- readRDS("../DD/data_buildings.w.parcels.energy.propinfo_3m16angles.rds")

                                          #remove trouble buildings
  buildings.propinfo.energy.tree <- buildings.propinfo.energy.tree[which(buildings.propinfo.energy.tree@data$BUILDIN != 190246),]

d <- left_join(buildings.propinfo.energy.tree@data, bs, c("BUILDINGFO" = "BUILDIN"))
   buildings.propinfo.energy.tree@data <- d

buildings.propinfo.energy.tree.building <- buildings.propinfo.energy.tree
shapefile(buildings.propinfo.energy.tree.building, "../DD/buildings.w.parcels.energy.propinfo.tree.building.shp", overwrite = T)
saveRDS(colnames(buildings.propinfo.energy.tree.building@data), "../DD/buildings.w.parcels.energy.propinfo.tree.building.colnames.rds")
saveRDS(buildings.propinfo.energy.tree.building@data, "../DD/buildings.w.parcels.energy.propinfo.tree.building.data.rds")

library(plyr)
  library(ascii)
  library(broom)
  library(tidyr)
  library(stringr)
  library(raster)
  library(rgeos)
  library(rgdal)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(dplyr)
  library(lme4)
library(GGally)

library(plyr)
  library(ascii)
  library(broom)
  library(tidyr)
  library(stringr)
  library(raster)
  library(rgeos)
  library(rgdal)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(dplyr)
  library(lme4)
library(GGally)
