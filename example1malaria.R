########################
# MAP COUNTRY
########################

map <- ne_states(country = "Madagascar", returnclass = "sf")
map <- map[order(map$name), ]


########################
# POPULATION RASTER
########################

r <- terra::rast("data/mdg_pd_2012_1km.tif")
map$popraster <- terra::extract(r, map, sum, na.rm = TRUE)$mdg_pd_2012_1km


########################
# AREAL DATA ESTIMATION
########################

# positive, total, prevalence
da <- read.csv("data/m.csv", header = T)
da <- da[order(da$name), ]
da$pos <- round(da$X2012)
da$total <- round(map$popraster) # total population
da$prev <- da$pos / da$total # prevalence
da <- cbind(map, da)
mapview(da, zcol = "prev")

########################
# ELEVATION RASTER
########################

elv <- rast("data/dataexample1/elv.tif")
da$cov <- terra::extract(elv, map, mean, na.rm = TRUE)[[2]]


########################
# POINT DATA ESTIMATION
########################

point <- getPR(country = "Madagascar", species = "Pf")
point <- point[point$year_start == 2012,]
point <- data.frame(prev = point$pr, Longitude = point$longitude, Latitude = point$latitude, Positive = point$positive,
                    Examined = point$examined )
point <- na.omit(point)
point$pos <- round(point$Positive)
point$total <- round(point$Examined)
point$prev <- point$pos/point$total
point$cov <- terra::extract(elv, point[, c("Longitude", "Latitude")])

point <- st_as_sf(point, coords = c("Longitude", "Latitude"), crs = 4326)
dp <- point
dp <- st_filter(dp, map)
mapview(dp, zcol = "prev")




########################
# PROJECTION UTM
########################

sf_use_s2(FALSE)
crsproj <- 29738

dp <- st_transform(dp, crs = crsproj)
da <- st_transform(da, crs = crsproj)
boundaryregion <- st_as_sf(st_union(da))

########################
# POINT DATA PREDICTION
########################
# Continuous surface
bb <- unname(attributes(st_geometry(boundaryregion))$bbox)
x <- seq(bb[1] - 1, bb[3] + 1, length.out = 100)
y <- seq(bb[2] - 1, bb[4] + 1, length.out = 100)
coop <- expand.grid(x, y)
coop_sf <- sf::st_as_sf(coop, coords = c('Var1','Var2'), crs = crsproj)
dpcontsurface <- st_filter(coop_sf, boundaryregion)
ggplot() + geom_sf(data = dpcontsurface)

coordTransf <- st_coordinates(st_transform(dpcontsurface, 4326))
dpcontsurface$cov <- terra::extract(elv, coordTransf[, c("X", "Y")])[[1]]



########################
# FIT MODEL
########################

# point data for estimation
depoint <- dp
depoint$value <- depoint$pos
depoint$numtrials <- depoint$total

# areal data for estimation
dearea <- da
dearea$value <- dearea$pos
dearea$numtrials <- dearea$total



########################
# Mesh
mesh <- fnCreateMesh(depoint, boundaryregion, 100)
plot(mesh)
points(as.matrix(st_coordinates(depoint)[ , c(1, 2)]), col = 2)

# Fit model "binomial"
respre <- fnPredictMelding(depoint =  depoint, dearea = dearea,
                           dppoint = dpcontsurface, dparea = NULL,
                           boundaryregion = boundaryregion,  mesh = mesh, model = "binomial")

res <- respre[[3]]
summary(res)

# pred_mean
dppred <- respre[[1]]
dppred$x <- st_coordinates(dppred)[, 1]
dppred$y <- st_coordinates(dppred)[, 2]

########################
# Exceedance probabilities

index <- inla.stack.index(stack = respre[[4]], tag = "pred1")$data

dppred$excprob1 <- sapply(res$marginals.fitted.values[index],
                          FUN = function(marg){1-inla.pmarginal(q = 0.1, marginal = marg)})


###############
# PLOTS
###############

theme_set(theme_minimal(base_size = 20))
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
bound <- st_union(map)
bound <- st_transform(bound, 4326)

min_value <- min(depoint$prev, dearea$prev, dppred$pred_mean, 0)
max_value <- max(depoint$prev, dearea$prev, dppred$pred_mean)
(color_limits <- c(min_value, max_value))


fnPlot <- function(dppred, vble, titlelegend, color_limits){
  ggplot(dppred) + geom_sf() + geom_raster(data = as.data.frame(dppred, xy = TRUE), aes_string(x = "x", y = "y", fill = vble)) +
    scale_fill_viridis_c(option = "viridis", limits = color_limits) + labs(fill = titlelegend, x = "", y = "") +
    geom_sf(data = bound, aes(geometry = geometry), col = "black", linewidth = .1, fill = NA) + coord_sf(label_graticule = "----") +
    labs(color = titlelegend)
}

# Posterior mean
fnPlot(dppred, "pred_mean", "Pred prevalence", color_limits)
ggsave("plots/malaria-predmean.png")

# Exceedance probabilities
fnPlot(dppred, "excprob1", "P(prev > 0.1)", c(0, 1))
ggsave("plots/malaria-excprob.png")

# Point data of malaria 
ggplot(data = bound) + geom_sf() + geom_sf(data = depoint, aes(color = prev), size = 2) +
  scale_fill_viridis_c(option = "viridis", limits = c(0, 1)) +
  labs(color = "Prevalence") + coord_sf(label_graticule = "----")
ggsave("plots/malaria-pointdata.png")

# Areal data of malaria
ggplot(data = bound) + geom_sf(data = dearea, aes(fill = prev)) + labs(fill = "Prevalence") +
  labs(color = "Prevalence") + coord_sf(label_graticule = "----")
ggsave("plots/malaria-arealdata.png")

# Altitude as a raster
polygons <- as.polygons(elv, dissolve = TRUE)
sf_object <- st_as_sf(polygons)
ggplot(bound) + geom_sf(data = sf_object, aes(color = sf_object$MDG_elv_msk)) + labs(color = "Altitude") +
  labs(color = "Altitude") + coord_sf(label_graticule = "----")
ggsave("plots/malaria-altitude.png")

