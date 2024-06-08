#' @name ESM_Bp.Sampling
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @title Ensemble of Small Models: Sampling background points using 4 different methods.
#' @description This function generates background following the 4 different methods described by Steen et al (2024).
#' Two are generated in the environmental space whereas the remaining is in the geographic space (see details).
#'  
#' @param env a \code{SpatRaster} of at least one layer. if \emph{method = "strat.geo"}, a minimum of 2 layers are needed.
#' @param n.points \code{integer}. The number of background to be selected. \emph{Note that this number can change depending on the technique}
#' @param method \code{character}. one of: "rand.geo", "strat.geo", "rand.env" or "strat.env". \emph{see Details}.
#' @param digit.val.env \code{integer}. The number of digit to keep to remove too similar environmental values.
#' Only needed when method = "rand.env". \emph{Default: 1}.
#' @param aggr.fact.geo \code{integer}. The aggregating factor to generate the checkerboard. Only needed 
#' when method = "strat.geo". \emph{see Details}. \emph{Default: 5}.
#' @param n.strat.env \code{integer}. The number of classes to create for each environmental layer. Only needed
#' when method = "strat.env". \emph{see Details}. \emph{Default: 3}.
#' @param To.plot \code{logical}. Should the background point plotted on a map (in black).
#' @param xy.pres A two-column \code{matrix} or \code{data.frame} containing the coordinates of the species occurrences.
#' Optional and only used when To.plot = TRUE. Plot the occurrences on the map (in aquamarine). \emph{Default: NULL}.
#' 
#' @details
#' The "rand.geo" method corresponds to a random selection of background points in the geographic space. 
#' This is the most used techniques in SDM studies. This selection can also be stratified (method ="strat.geo"). To 
#' realise this stratification, a checkerboard is created with a pixel of a certain resolution. The size of these
#' pixels can be modified with the argument 'aggr.fact.geo'. For each pixel of the checkerboard the same number of 
#' background points are randomly selected to reach the value of \emph{n.points}. Selecting background points can also 
#' be performed in the environmental space (see Steen et al, 2024). The first method is the full random background
#' point selection (method = "rand.env"). To do so, a PCA is first performed and the two fist axes are kept to 
#' reflect the environmental space. Then, a grid of 100*100 pixel is created. To reduce too similar points, 
#' the PCA scores of each observation are rounded at a certain digit (argument digit.val.env) and only 
#' one observation among the similar ones is kept. Several points are then randomly selected in each pixel of this grid.
#' By doing this, the entire environmental is captured avoiding the over-abundance of some common environment. The
#' second method is the stratified selection of background point in the environmental spaces (method = "strat.env").
#' To perform this selection, each environmental layer is converted into classes of \emph{n.strat.env} (Default: 3)
#' by dividing the range of the environmental predictor by \emph{n.strat.env} of the same size. For each pixel 
#' constituting the grid, we then combined all the predictor together, creating a combination of classes. For each 
#' of the different possible combination of classes available, the same number of background point is selected  to 
#' reach the value of \emph{n.points}.
#'  
#' @return a \code{matrix} containing the selected background points. The two first columns correspond to the coordinates.
#' The other columns correspond to the environmental values of these points and the last column correspond to the geographic 
#' or environmental classes (named "BigClass) to which the observation belongs (only for method = "strat.geo" or "strat.env")
#'  
#' @examples 
#' library(terra)
#' library(ecospat)
#' env <- terra::rast(system.file("extdata","ecospat.testEnv.tif",package="ecospat"))
#' # Selection full random in the geographic space
#' Bp <- ESM_Bp.Sampling(env = env,
#'                       n.points = 1000,
#'                       method = "rand.geo",
#'                       To.plot = FALSE)
#'                        
#' # Selection stratified  in the geographic space                     
#' Bp <- ESM_Bp.Sampling(env = env,
#'                       n.points = 1000,
#'                       method = "strat.geo",
#'                       aggr.fact.geo = 2,
#'                       To.plot = FALSE)
#'                        
#' # Selection full random in the environmental space                     
#' Bp <- ESM_Bp.Sampling(env = env,
#'                       n.points = 1000,
#'                       method = "rand.env",
#'                       digit.val.env = 2,
#'                       To.plot = FALSE)       
#'                                         
#' # Selection stratified in the environmental space                     
#' Bp <- ESM_Bp.Sampling(env = env,
#'                       n.points = 1000,
#'                       method = "strat.env",
#'                       n.strat.env = 3,
#'                       To.plot = FALSE)                          
#' @references 
#' Steen,B., Broennimann, O., Maiorano, L., Guisan, . 2024. How sensitive are species distribution models to different background point 
#' selection strategies? A test with species at various equilibrium levels. 
#' \emph{Ecological Modelling}. \bold{493}, 110754. \doi{10.1016/j.ecolmodel.2024.110754}.
#' 
#' @export

ESM_Bp.Sampling <- function(env,
                            n.points = 10000,
                            method = "rand.geo",
                            digit.val.env = 1,
                            aggr.fact.geo = 5,
                            n.strat.env = 3,
                            To.plot = FALSE,
                            xy.pres = NULL){
  
  if(!inherits(env,"SpatRaster")){
    stop("env must be a SpatRaster from terra package.")
  }
  if(!is.numeric(n.points) | n.points%%1 != 0 |  n.points <= 0){
    stop("n.points must be a positive integer")
    
  }
  if(!(method %in% c("rand.geo", "rand.env", "strat.geo", "strat.env")) | length(method)>1){
    stop("method must be one of: 'rand.geo', 'rand.env', 'strat.geo' or, 'strat.env'")
  }
  
  ## Start the process ----
  
  if(method == "rand.geo"){
    coord.Env <- terra::crds(env)
    if(n.points > nrow(coord.Env)){
      cat("\nn.points is greater than the number of pixels with values. All pixels are taken")
      bp <- coord.Env
    }else{
      bp <- coord.Env[sample(1:nrow(coord.Env), size = n.points, replace = FALSE),]
    }
    bp <-  cbind.data.frame(bp,terra::extract(env,bp))
    
    
  }else if(method == "rand.env"){
    
    if(terra::nlyr(env)<2){
      stop("when method == rand.env or strat.env, a minimum of 2 environnmental layers is needed in env.")
    }
    if(digit.val.env %%1 !=0 | digit.val.env < 0){
      stop("digit.val.env must be an integer greater than 0")
    }
    
    env.dat <- terra::as.data.frame(env, xy = T)
    env.dat <- stats::na.omit(env.dat)
    env.pca <- ade4::dudi.pca(env.dat[,-c(1:2)], scale = TRUE, scannf = FALSE, nf = 2)
    Total.variance <- 100 * (cumsum(env.pca$eig/sum(env.pca$eig)))[2]
    cat(paste("\nThe two first PCA axes explained", round(Total.variance,1),"% of the total variance"))
    env.score.round <- round(env.pca$li,digit.val.env)
    # density_2d <- ks::kde(env.score.round,compute.cont=TRUE)
    # aa <- grDevices::contourLines(density_2d$eval.points[[1]], density_2d$eval.points[[2]], density_2d$estimate, level = 0.001)
    # bb <- lapply(aa, function(x){do.call(cbind,x)[,c(2:3)]})
    # cc <- terra::vect(bb,type = "polygons",crs="")
    # cc<-terra::aggregate(cc)
    grid.env <- terra::rast(terra::ext(apply(env.score.round,2, range)), ncols = 100, nrows = 100,crs="") #vals = 1:(100*100)
    # ee <- terra::mask(vide,cc)
    
    env.score.round.filt <- unique(env.score.round)
    
    cell.pos <- terra::cellFromXY(grid.env, env.score.round.filt)
    cell.env <- unique(cell.pos)
    if(n.points<length(cell.env)){
      cat("\nThere are more environmental classes than n.points. Thus,the number of background points will be equal to the number of classes")
      n.ObsPerClass = 1
    }else{
      n.ObsPerClass <- ceiling(n.points/length(cell.env))
    }
    
    ToKeep <- c()
    ToPrint <- TRUE
    
    ##Sampling
    for(i in 1:length(cell.env)){
      pointVal <- which(cell.pos == cell.env[i])
      if(length(pointVal) < n.ObsPerClass){
        if(ToPrint){
          cat("\nSome Classes have less observation than the number of observation per class. Thus, all the observations for these classes will be sampled. 
              Note that the number of background points sampled will be less than asked.")
        }
        ToPrint <- FALSE
        ToKeep <- c(ToKeep,pointVal)
      }else{
        if(length(pointVal)==1){
          ToKeep <- c(ToKeep,pointVal)
        }else{
          ToKeep <- c(ToKeep,sample(pointVal,n.ObsPerClass, replace = FALSE))
          
        }
      }
    }
    
    
    bp <- env.dat[rownames(env.score.round.filt)[ToKeep],]
    
    
  }else if(method == "strat.geo"){
    if(aggr.fact.geo %%1 !=0 | aggr.fact.geo <= 1){
      stop("aggr.fact.geo must be an integer greater than 1")
    }
    coord.Env <- terra::crds(env)
    Checkboard <- terra::aggregate(terra::subset(env,1), fact = aggr.fact.geo, na.rm=T)
    Checkboard.pos <- terra::cellFromXY(Checkboard, coord.Env)
    Checkboard.pos.class <- unique(Checkboard.pos)
    
    if(n.points<length(Checkboard.pos.class)){
      cat("There are more geographic classes than n.points. Thus,the number of background points will be equal to the number of classes")
      n.ObsPerClass = 1
    }else{
      n.ObsPerClass <- ceiling(n.points/length(Checkboard.pos.class))
    }
    
    ToKeep <- c()
    ToPrint <- TRUE
    
    ##Sampling
    for(i in 1:length(Checkboard.pos.class)){
      pointVal <- which(Checkboard.pos == Checkboard.pos.class[i])
      if(length(pointVal) < n.ObsPerClass){
        if(ToPrint){
          cat("Some Classes have less observation than the number of observation per class. Thus, all the observations for these classes will be  sampled. 
              Note that the number of background points sampled will be less than  asked.")
        }
        ToPrint <- FALSE
        ToKeep <- c(ToKeep,pointVal)
      }else{
        if(length(pointVal)==1){
          ToKeep <- c(ToKeep,pointVal)
        }
        ToKeep <- c(ToKeep,sample(pointVal,n.ObsPerClass, replace = FALSE))
      }
    }
    
    bp <- terra::extract(env,coord.Env[ToKeep,])
    
    bp <- cbind.data.frame(coord.Env[ToKeep,],
                           bp, 
                           BigClass = Checkboard.pos[ToKeep])
    
  }else{
    
    ##Strat in the environment ----
    
    env.dat <- terra::as.data.frame(env, xy = T)
    env.dat <- stats::na.omit(env.dat)
    
    ## Create classes per predictor
    env.class <- apply(as.data.frame(env.dat[,-c(1:2)]), 2, 
                       function(x,n.strat.env){
                         classes <- seq(min(x),max(x),
                                        length.out = n.strat.env+1)
                         classes[1] = classes[1]-1
                         classes[n.strat.env+1] = classes[n.strat.env+1]+1 ##to avoid NAs
                         return(as.numeric(cut(x,classes)))
                       }, 
                       n.strat.env = n.strat.env)
    ## Generate the classes combining all the predictors
    Class.Tot <- as.factor(apply(env.class, 1, paste, collapse = "" ))
    
    Cat.Class.Tot <- levels(Class.Tot)
    
    if(n.points<length(Cat.Class.Tot)){
      cat("\nThere are more environmental classes than n.points. Thus,the number of background points will be equal to the number of classes")
      n.ObsPerClass = 1
    }else{
      n.ObsPerClass <- ceiling(n.points/length(Cat.Class.Tot))
    }
    
    ToKeep <- c()
    ToPrint <- TRUE
    
    ##Sampling
    for(i in 1:length(Cat.Class.Tot)){
      pointVal <- which(Class.Tot == Cat.Class.Tot[i])
      if(length(pointVal) < n.ObsPerClass){
        if(ToPrint){
          cat("\nSome Classes have less observation than the number of observation per class. Thus, all the observations for these classes will be  sampled. 
              Note that the number of background points sampled will be less than you asked.")
        }
        ToPrint <- FALSE
        ToKeep <- c(ToKeep,pointVal)
      }else{
        if(length(pointVal)==1){
          ToKeep <- c(ToKeep,pointVal)
        }
        ToKeep <- c(ToKeep,sample(pointVal,n.ObsPerClass, replace = FALSE))
      }
    }
    
    bp <- cbind.data.frame(env.dat[ToKeep,],
                           BigClass = Class.Tot[ToKeep])
    
    
  }
  if(To.plot){
    terra::plot(terra::subset(env,1), col= "grey80", legend = F)
    terra::points(bp[,1:2], pch = 20)
    if(!is.null(xy.pres)){
      if(ncol(xy.pres) !=2){
        stop("\nxy.pres should only contain to columns")
      }
      terra::points(xy.pres, pch = 19, col = "aquamarine3")
    }
  }
  return(bp)
}



####

#' @name ESM_Range.Shift
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @title Ensemble of Small Models: Range Shift between projections
#' @description This function compare ESM or SDM binary projections between two time periods by comparing the change 
#' in suitable/Unsuitable pixels. It also returns percentage of loss and gained
#' 
#' @param proj.curr \code{SpatRaster}. One or more ESM/SDM binary projections at current time. 
#' \emph{Note that if the number of layer is greater than one for proj.curr, the number of layer 
#' in proj.fut must be the same as comparisons will be between between map i of proj.curr and map i of proj.fut.}
#' @param proj.fut \code{SpatRaster}. One or more ESM/SDM binary projections at future time.
#' 
#' @details
#' The comparisons betwween SDM/ESM binary predictions are made following this formula: 
#' \bold{RangeShift = proj.curr + 2*proj.fut}. 
#' Therefore, when the value is equal to 0. It correspond to a pixel considered as "absent"/"unsuitable" in both 
#' time period. 
#' When it is 1, pixels considered as "presence" in proj.curr  but not in proj.fut \bold{(pixel lost)}.
#' When it is 2, pixels considered as "presence" in proj.fut  but not in proj.curr \bold{(pixel gained)}.
#' When it is 3, pixels are considered as "present" in both time periods \bold{(pixel Stable)}.
#'   
#' 
#' 
#' @return \code{list} containing:
#' \itemize{
#' \item{RangeShift}: \code{SpatRaster}. The range comparisons between current and future (\emph{see Details for the values}).
#' \item{RangeShift.table}: \code{data.frame}. Summary of the comparisons between maps. containing 9 columns:
#' \itemize{
#' \item{layer}. The layer names of proj.fut used to make the comparisons.
#' \item{Pixel.Absence.Stable}. The number of pixels that are considered as absent for both time period.
#' \item{Pixel.Presence.Current.Only}. The number of pixels considered as "presence" in proj.curr but not in proj.fut.
#' \item{Pixel.Presence.Future.Only}. The number of pixels considered as "presence" in proj.fut but not in proj.curr.
#' \item{Pixel.Presence.Stable}. The number of pixels considered as "presence" both in proj.curr and proj.fut.
#' \item{Current.Range.size}. The number of pixels considered as "presence" in proj.curr.
#' \item{Future.Range.size}. The number of pixels considered as "presence" in proj.fut assuming that species is
#' not limited by their dispersal capacities.
#' \item{Perc.Loss}. The percentage of pixel considered as "presence" in proj.curr and becoming "absent"
#' in proj.fut compared to the current range size. 
#' \item{Perc.Gain}. The percentage of pixel considered as "absent" in proj.curr and becoming 
#' "present" in proj.fut compared to the current range size. 
#' }
#' }
#' 
#' @seealso [ESM_Modeling], [ESM_Ensemble.Modeling],  [ESM_Projection], [ESM_Ensemble.Projection], [ESM_Threshold]
#' @export

ESM_Range.Shift <- function(proj.curr,
                            proj.fut){
  
  if(!(inherits(proj.curr,"SpatRaster") & inherits(proj.fut,"SpatRaster"))){
    stop("proj.curr and proj.fut must be SpatRaster a accepted")
  }
  
  if(terra::nlyr(proj.curr) > 1 & terra::nlyr(proj.curr) != terra::nlyr(proj.fut)){
    stop("When proj.curr has more than one layer the number of layer in proj.fut must match terra::nlyr(proj.curr)")
  }
  shift <- proj.curr + 2 * proj.fut
  results <- terra::freq(shift, wide=T)
  n.column <- ncol(results)
  while(ncol(results) < 5){
    
    results <- cbind(results,0)
  } 
  colnames(results)[(n.column+1):5] = setdiff(c("layer",0:3), colnames(results))
  results <- results[,c(1,order(colnames(results)[2:5])+1)]
  colnames(results) <-c("layer","Pixel.Absence.Stable","Pixel.Presence.Current.only",
                        "Pixel.Presence.Future.Only","Pixel.Presence.Stable")
  results$layer <- names(proj.fut)
  results$Current.Range.size <- results$Pixel.Presence.Current.only+ results$Pixel.Presence.Stable
  results$Future.Range.size <- results$Pixel.Presence.Future.Only+ results$Pixel.Presence.Stable
  results$Perc.Loss <- round(100*results$Pixel.Presence.Current.only/
                               (results$Pixel.Presence.Current.only+results$Pixel.Presence.Stable),
                             2)
  
  results$Perc.Gain <- round(100*results$Pixel.Presence.Future.Only/
                               (results$Pixel.Presence.Current.only+results$Pixel.Presence.Stable),
                             2)  
  
  return(list(RangeShift = shift,
              RangeShift.table = results))
}