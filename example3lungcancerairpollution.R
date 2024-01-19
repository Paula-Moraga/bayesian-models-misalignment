#################################
# Get map of Alabama
#################################

map <- st_read("data/dataexample1/map.shp")
map <- st_transform(map, 4326)


####################################################
# Read lung cancer data counties Alabama
####################################################

d <- read_excel("data/dataexample1/d1.xlsx")

#################################
# Read population of Alabama (raster)
#################################

r <- rast("data/dataexample1/r.tif")

# Sum population in counties of Alabama map
map$popraster <- terra::extract(r, map, sum, na.rm = TRUE)$usa_pd_2002_1km


################################
# Merge map and areal data, and calculate E and SIR
#################################

map <- merge(map, d[, c("area", "Count")], by.x = "NAME", by.y = "area")
map$Y <- map$Count
map$E <- round(map$popraster * sum(map$Y)/sum(map$popraster))
map$SIR <- map$Y/map$E


####################################################
# Read PM2.5 (data.frame)
####################################################

dpol <- st_read("data/dataexample1/covariate.shp")

##################
# Calculate PM25 weighted by population within counties
##################

dpol$pm25Xpop <- as.vector(dpol$cov * terra::extract(r, st_coordinates(dpol)))[[1]]
dpol$pop <- as.vector(terra::extract(r, st_coordinates(dpol)))[[1]]
inter <- st_intersects(map, dpol)
map$pm25xpop <- sapply(inter, FUN = function(x){sum(dpol[x, ]$pm25Xpop, na.rm = TRUE)})
map$sumpop <- sapply(inter, FUN = function(x){sum(dpol[x, ]$pop, na.rm = TRUE)})
map$pm25weightedpop <- map$pm25xpop/map$sumpop


######################
# Models
######################


fnFitModel <- function(typemodel){
  map1 <- as(map, "Spatial")
  nb <- poly2nb(map1)
  nb2INLA("map1.adj", nb)
  g <- inla.read.graph(filename = "map1.adj")
  map1$re_u <- 1:nrow(map1)
  map1$re_v <- 1:nrow(map1)
  map1$b <- 1:nrow(map1)
  
  if(typemodel == "BYM"){formula <- Y ~ pm25weightedpop +
    f(re_u, model = "besag", graph = g, scale.model = TRUE) + f(re_v, model = "iid")}
  
  if(typemodel == "iid"){formula <- Y ~ pm25weightedpop +f(re_v, model = "iid")}
  
  res <- inla(formula, family = "poisson", data = map1@data, E = E, control.predictor = list(compute = TRUE),
              control.compute = list(dic = T, waic = T, cpo = T, return.marginals.predictor = TRUE), verbose = FALSE)
  return(res)
}

# Fit models
fitBYM <- fnFitModel("BYM")
fitiid <- fnFitModel("iid")

# Compare model based on DIC and WAIC
data.frame(model = c("BYM", "iid"),
DIC = c(fitBYM$dic$dic, fitiid$dic$dic),
WAIC = c(fitBYM$waic$waic,fitiid$waic$waic))

summary(fitBYM)
  

###########################
# PLOTS
###########################

theme_set(theme_minimal(base_size = 20))
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")

# Observed cases Y, SIR, RR, iid random effect, spatial random effect

map$RR <- fitBYM$summary.fitted.values[, "mean"]
map$reiid <- fitBYM$summary.random$re_v$mean
map$rebesag <- fitBYM$summary.random$re_u$mean
map$exc <- sapply(fitBYM$marginals.fitted.values,
                  FUN = function(marg){1 - inla.pmarginal(q = 1.2, marginal = marg)})

# PM2.5 weighted by population
ggplot(data = map) + geom_sf(aes(fill = pm25weightedpop)) + labs(fill = "WNO3-") + coord_sf(label_graticule = "----")
ggsave("plots/lc-pm25weightedpop.png")
ggplot(data = map) + geom_sf(aes(fill = Y)) +  labs(fill = "Y") + coord_sf(label_graticule = "----")
ggsave("plots/lc-Y.png")
ggplot(data = map) + geom_sf(aes(fill = E)) +  labs(fill = "E") + coord_sf(label_graticule = "----")
ggsave("plots/lc-E.png")
ggplot(data = map) + geom_sf(aes(fill = SIR)) +  labs(fill = "SIR") + coord_sf(label_graticule = "----")
ggsave("plots/lc-SIR.png")
ggplot(data = map) + geom_sf(aes(fill = RR)) +  labs(fill = "RR") + coord_sf(label_graticule = "----")
ggsave("plots/lc-RR.png")
ggplot(data = map) + geom_sf(aes(fill = reiid)) +  labs(fill = "unstructured") + coord_sf(label_graticule = "----")
ggsave("plots/lc-reiid.png")
ggplot(data = map) + geom_sf(aes(fill = rebesag)) +  labs(fill = "structured") + coord_sf(label_graticule = "----")
ggsave("plots/lc-rebesag.png")
ggplot(data = map) + geom_sf(aes(fill = exc)) +  labs(fill = "P(RR > 1.2)") + coord_sf(label_graticule = "----")
ggsave("plots/lc-exc.png")

# Plot PM25 data.frame of locations with ggplot2
dpol$lon <- st_coordinates(dpol)[, 1]
dpol$lat <- st_coordinates(dpol)[, 2]
ggplot(data = map) + geom_sf(col = NA, fill = NA) +
  geom_point(data = dpol, aes(x = lon, y = lat, color = cov)) + labs(color = "", x = "", y = "") +
  labs(color = "WNO3-") + coord_sf(label_graticule = "----")
ggsave("plots/lc-pm25.png")

