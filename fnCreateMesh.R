fnCreateMesh <- function(de1, boundaryregion, maxedgedenominator) {
  # maxedgedenominator = 25 (air pollution UK), 100 (malaria)
  
  location <- NULL
  if(!is.null(de1)){
    #de1 <- st_transform(de1, crs = 4326)
    location <- as.matrix(st_coordinates(de1)[, c(1, 2)])
  }
  
  maxedge <-  fnMaxEdgeMeshFromBoundary(boundaryregion)
  
  bdsp  <- as(boundaryregion, 'Spatial')
  mesh <- inla.mesh.2d(loc = location, boundary = bdsp, max.edge = c(maxedge/10, maxedge), cutoff = maxedge/maxedgedenominator)
  
  return(mesh)
}

fnMaxEdgeMeshFromBoundary <- function(bd){
  maxedge <- 0.33 * max(attributes(st_geometry(bd))$bbox[3] - attributes(st_geometry(bd))$bbox[1],
                        attributes(st_geometry(bd))$bbox[4] - attributes(st_geometry(bd))$bbox[2])
  return(maxedge)
}




