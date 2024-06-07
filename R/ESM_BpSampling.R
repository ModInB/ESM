#' @name ESM_Bp.Sampling
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @title Ensemble of Small Models: Sampling background points using 4 different methods.
#' @references 
#' Steen,B., Broennimann, O., Maiorano, L., Guisan, . 2024. How sensitive are species distribution models to different background point 
#' selection strategies? A test with species at various equilibrium levels. 
#' \emph{Ecological Modelling}. \bold{493}, 110754. \doi{10.1016/j.ecolmodel.2024.110754}.
#' 
#' @export

ESM_Bp.Sampling <- function(env,
                            n.points = 10000,
                            method = "rand.geo",
                            aggr.fact.geo = 5,
                            n.strat.env = 3,
                            round.val.env = 2,
                            To.plot = FALSE,
                            xy.pres = NULL){
  
  if(!inherits(env,"SpatRaster")){
    stop("env should be a SpatRaster from terra package.")
  }
  if(!is.numeric(n.points) | n.points%%1 != 0 |  n.points <= 0){
    stop("n.points should be a positive integer")
    
  }
  if(!(method %in% c("rand.geo", "rand.env", "strat.geo", "strat.env")) | length(method)>1){
    stop("method should be one of: 'rand.geo', 'rand.env', 'strat.geo' or, 'strat.env'")
  }
  
  ## Start the process ----
  
  if(method == "rand.geo"){
    coord.Env <- terra::crds(env)
    if(n.points > n.cells){
      cat("\nn.points is greater than the number of pixels with values. All pixels are taken")
      test.bp <- coord.Env
    }else{
      test.bp <- coord.Env[sample(1:nrow(coord.Env), size = n.points, replace = FALSE),]
    }
    test.bp <-  cbind.data.frame(test.bp,
                                 terra::extract(env,test.bp))
    
    
  }else if(method == "rand.env"){
    
    if(terra::nlyr(env)<2){
      stop("when method == rand.env or strat.env, a minimum of 2 environnmental layers is needed in env.")
    }
    
    test <- terra::as.data.frame(env, xy = T)
    test <- stats::na.omit(test)
    test.env.pca <- ade4::dudi.pca(test[,-c(1:2)], scale = TRUE, scannf = FALSE, nf = 2)
    Total.variance <- 100 * (cumsum(test.env.pca$eig/sum(test.env.pca$eig)))[2]
    cat(paste("\nThe two first PCA axes explained", round(Total.variance,1),"% of the total variance"))
    env.score.round <- round(test.env.pca$li,round.val.env)
    # density_2d <- ks::kde(env.score.round,compute.cont=TRUE)
    # aa <- grDevices::contourLines(density_2d$eval.points[[1]], density_2d$eval.points[[2]], density_2d$estimate, level = 0.001)
    # bb <- lapply(aa, function(x){do.call(cbind,x)[,c(2:3)]})
    # cc <- terra::vect(bb,type = "polygons",crs="")
    # cc<-terra::aggregate(cc)
    vide <- terra::rast(terra::ext(apply(env.score.round,2, range)), ncols = 100, nrows = 100,crs="") #vals = 1:(100*100)
    # ee <- terra::mask(vide,cc)
    ee <- vide
    env.score.round.filt <- unique(env.score.round)

    cell.pos <- terra::cellFromXY(ee, env.score.round.filt)
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
          cat("\nSome Classes have less observation than the number of observation per class. Thus, all the observations for these classes will be  sampled. 
              Note that the number of background points sampled will be less than you asked.")
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
  
    
    test.bp <- test[rownames(env.score.round.filt)[ToKeep],]
      

  }else if(method == "strat.geo"){
    
    coord.Env <- terra::crds(env)
    Checkboard <- terra::aggregate(terra::subset(env,1), fact = aggr.fact, na.rm=T)
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
    
    test.bp <- terra::extract(env,coord.Env[ToKeep,])
    
    test.bp <- cbind.data.frame(coord.Env[ToKeep,],
                                test.bp, 
                                BigClass = Checkboard.pos[ToKeep])
    
  }else{
    
    ##Strat in the environment ----
    
    test <- terra::as.data.frame(env, xy = T)
    test <- stats::na.omit(test)
    
    ## Create classes per predictor
    Test.class <- apply(as.data.frame(test[,-c(1:2)]), 2, 
                        function(x,n.strat.env){
                          classes <- seq(min(x),max(x),
                                         length.out = n.strat.env+1)
                          classes[1] = classes[1]-1
                          classes[n.strat.env+1] = classes[n.strat.env+1]+1 ##to avoid NAs
                          return(as.numeric(cut(x,classes)))
                          }, 
                        n.strat.env = n.strat.env)
    ## Generate the classes combining all the predictors
    Class.Tot <- as.factor(apply(Test.class, 1, paste, collapse = "" ))
      
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
    
    test.bp <- cbind.data.frame(test[ToKeep,],
                                BigClass = Class.Tot[ToKeep])
    
    
  }
  if(To.plot){
    terra::plot(terra::subset(env,1), col= "grey80", legend = F)
    terra::points(test.bp[,1:2], pch = 20)
    if(!is.null(xy.pres)){
      if(ncol(xy.pres) !=2){
        stop("\nxy.pres should only contain to columns")
      }
      terra::points(xy.pres, pch = 19, col = "aquamarine3")
    }
  }
  return(test.bp)
}
