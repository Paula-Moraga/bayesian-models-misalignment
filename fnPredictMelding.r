fnPredictMelding <- function(depoint = NULL, dearea = NULL, dppoint = NULL, dparea = NULL, boundaryregion,
                             mesh = NULL, priorspdesigma = NULL, priorspderange = NULL, model = "gaussian"){
  
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # Check inputs
  fnCheckInputsMelding(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  de1ToF <- !is.null(de1)
  de2ToF <- !is.null(de2)
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  if(model == "gaussian"){ # air pollution UK
    maxedgedenominator <- 25
  }else{ # malaria Madagascar
    maxedgedenominator <- 100    
  }

  # Mesh
  if(is.null(mesh)){
    mesh <- fnCreateMesh(de1, boundaryregion, maxedgedenominator)
  }
  
  # Create spde and index
  if(!is.null(priorspdesigma) & !is.null(priorspderange)){
    fnCheckPrior(priorspdesigma, priorspderange)
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = priorspderange, prior.sigma = priorspdesigma)
  }else{
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
  }
  indexs <- inla.spde.make.index("s", spde$n.spde)
  
  
  # Projection matrices for points (estimation point and prediction point)
  if(de1ToF){Ae1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(de1)[ , c(1,2)]))}
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dp1)[ , c(1,2)]))}
  
  # Create projection matrix A for areas (estimation area and prediction area)
  if(de2ToF){Ae2 <- fnProjectionMatrixArea(de2, mesh)}
  if(dp2ToF){Ap2 <- fnProjectionMatrixArea(dp2, mesh)}
  
  # Create stk.full
  stk.e1 = NULL
  stk.e2 = NULL
  stk.p1 = NULL
  stk.p2 = NULL
  
  # estimation point (de1), estimation area (de2), prediction point (dp1), prediction area (dp2)
  if(de1ToF){stk.e1 <- inla.stack(tag = "est1", data = list(y = de1$value, numtrials = de1$numtrials), A = list(1, Ae1), effects = list(data.frame(b0 = rep(1, nrow(de1)), cov = de1$cov), s = indexs$s))}
  if(de2ToF){stk.e2 <- inla.stack(tag = "est2", data = list(y = de2$value, numtrials = de2$numtrials), A = list(1, Ae2), effects = list(data.frame(b0 = rep(1, nrow(de2)), cov = de2$cov), s = indexs$s))}
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA, numtrials = NA), A = list(1, Ap1), effects = list(data.frame(b0 = rep(1, nrow(dp1)), cov = dp1$cov), s = indexs$s))}
  if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA, numtrials = NA), A = list(1, Ap2), effects = list(data.frame(b0 = rep(1, nrow(dp2)), cov = dp2$cov), s = indexs$s))}
  # construct stack full with the data we have # stk.full <- inla.stack(stk.e1, stk.p1)
  stk.full <- do.call(inla.stack, list(stk.e1, stk.e2, stk.p1, stk.p2)[c(de1ToF, de2ToF, dp1ToF, dp2ToF)])
  
  
  
  # Specify formula melding
  formula <- y ~ 0 + b0 + f(s, model = spde)
  
  # Call inla()
  if(model == "gaussian"){
    res <- inla(formula, family = "gaussian", data = inla.stack.data(stk.full),
                control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))
  }
  if(model == "binomial"){
    formula <- y ~ 0 + b0 + cov + f(s, model = spde)
  res <- inla(formula, family = "binomial", Ntrials = numtrials, data = inla.stack.data(stk.full),
              control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)),
              control.compute = list(return.marginals.predictor = TRUE))
  }
  
  # Retrieve predictions points
  # Predictions points
  if(dp1ToF){dp1 <- fnRetrievePredictions(stk.full, res, "pred1", dp1)}
  # Predictions areas
  if(dp2ToF){dp2 <- fnRetrievePredictions(stk.full, res, "pred2", dp2)}
  
  return(list(dp1,dp2,res,stk.full))
}

#####################################################################
# END MAIN FUNCTION
#####################################################################



#####################################################################
# INI AUXILIARY FUNCTIONS
#####################################################################




# Create projection matrix A for areas
fnProjectionMatrixArea <- function(de2, mesh){
  
  meshcoo <- data.frame(X = mesh$loc[, 1], Y = mesh$loc[, 2])
  
  meshin <- meshcoo %>% st_as_sf(coords = c("X", "Y"), dim = "XY") %>%
    st_set_crs(crs(de2)) %>% st_cast("MULTIPOINT")
  
  # find points in mesh n area.sf
  locin <- st_join(meshin, de2, left = F)
  
  block <- rep(0, nrow(locin))
  
  for(i in 1:nrow(de2)) {
    block[as.vector(which(!is.na(st_join(locin, de2[i,], left = T)$value.y)))] <- i
  }
  
  A <- inla.spde.make.A(mesh = mesh,
                        loc = as.matrix(st_coordinates(locin[, 1]))[, c(1, 2)],
                        block = block, block.rescale = "sum")
  return(A)
}


# Retrieve predictions
fnRetrievePredictions <- function(stack, res, tag, dataset){
  index <- inla.stack.index(stack = stack, tag = tag)$data
  dataset$pred_mean <- res$summary.fitted.values[index, "mean"]
  dataset$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
  dataset$pred_ul <- res$summary.fitted.values[index, "0.975quant"]
  return(dataset)
}

# Prior Check
fnCheckPrior <- function(prior.range, prior.sigma){
  if(length(prior.range) != 2 ||sum(prior.range >0) != 2 || prior.range[2] > 1){
    stop("The parameters of prior.sigma are not valid, the length 2 vector contains spatial range (>0) and 
         prior probability of the deviation (0~1)")
  }
  if(length(prior.sigma) != 2 || sum(prior.sigma >0) != 2 || prior.sigma[2] > 1){
    stop("The parameters of prior.sigma are not valid, the length 2 vector contains marginal standard deviation (>0) and 
         prior probability of the deviation (0~1)")
  }
}