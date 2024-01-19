##########################
# POINT DATA PM2.5
##########################

crsproj <- "+proj=utm +zone=29 +nord +units=km"
depoint <- st_read('data/depoint_cov.shp')  %>%  st_transform(crsproj)

##########################
# AREAL DATA PM2.5
##########################

dearea <- st_read('data/dearea_cov.shp') %>%  st_transform(crsproj)

####################################################
# MAP
####################################################

map <- ne_countries(country = "United Kingdom", scale = "large", returnclass = "sf")
map <- st_transform(map, crsproj)
boundaryregion <- st_as_sf(st_union(map))

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

########################
# FIT MODEL
########################

# point data for estimation
depoint <- depoint
depoint$value <- depoint$value
depoint$cov <- rep(1, nrow(depoint))
depoint$numtrials <- rep(1, nrow(depoint))

# areal data for estimation
dearea <- dearea
dearea$value <- dearea$value
dearea$cov <- rep(1, nrow(dearea))
dearea$numtrials <- rep(1, nrow(dearea))

# point data for prediction
dpcontsurface$cov <- rep(1, nrow(dpcontsurface))

########################
# Mesh
mesh <- fnCreateMesh(depoint, boundaryregion, 25)
plot(mesh)
points(as.matrix(st_coordinates(depoint)[ , c(1, 2)]), col = 2)

# Fit model "normal"
respre <- fnPredictMelding(depoint =  depoint, dearea = dearea,
                           dppoint = dpcontsurface, dparea = NULL,
                           boundaryregion = boundaryregion,  mesh = mesh, model = "gaussian")

res <- respre[[3]]
summary(res)



###########################
# PLOTS
###########################


theme_set(theme_minimal(base_size = 20))
theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
# Define color scale limits
dearea_clipped <- st_intersection(dearea, boundaryregion)
min_value <- min(depoint$value, dearea_clipped$value)
max_value <- max(depoint$value, dearea_clipped$value)
color_limits <- c(min_value, max_value)

# point data
ggplot(data = map) + geom_sf() +
  geom_sf(data = depoint, aes(geometry = geometry, color = value), size = 3) +
  scale_color_viridis_c(option = "viridis", limits = color_limits) +  labs(color = expression(PM[2.5])) +
  geom_sf(data = boundaryregion, fill = NA, col = "black", linewidth = 0.5) + coord_sf(label_graticule = "----", crs = crsproj, xlim = c(500, 1230.3))
ggsave("plots/airpollution-pointdata.png")

# areal data
ggplot(data = map) + geom_sf() +
  geom_sf(data = dearea_clipped, aes(fill = value)) +
  scale_fill_viridis_c(option = "viridis", limits = color_limits) +  labs(fill = expression(PM[2.5])) +
  geom_sf(data = map, fill = NA, col = "black", linewidth = 0.5)  + coord_sf(label_graticule = "----", xlim = c(500, 1230.3))
ggsave("plots/airpollution-arealdata.png")


# pred_mean
dppred <- respre[[1]]
dppred$x <- st_coordinates(dppred)[, 1]
dppred$y <- st_coordinates(dppred)[, 2]


fnPlot <- function(dppred, vble, titlelegend){
ggplot(boundaryregion) +
  geom_tile(data = dppred, aes_string(x = "x", y = "y", fill = vble)) +
  scale_fill_viridis_c(option = "viridis", limits = color_limits) + labs(fill = titlelegend, x = "", y = "")  +
    geom_sf(data = map, fill = NA, col = "black", linewidth = 0.5) + coord_sf(label_graticule = "----", xlim = c(500, 1230.3))
}

fnPlot(dppred, "pred_mean", expression(paste("Pred ", PM[2.5])))
ggsave("plots/airpollution-predmean.png")
fnPlot(dppred, "pred_ll", "2.5 perc.")
ggsave("plots/airpollution-predll.png")
fnPlot(dppred, "pred_ul", "97.5 perc.")
ggsave("plots/airpollution-predul.png")

