#' @name Bp_Sampling
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @title Ensemble of Small Models: Sampling background points using 4 different methods.
#' @description This function generates background following the 4 different methods described by Steen et al (2024).
#' Two are generated in the environmental space whereas the remaining is in the geographic space (see details).
#' 
#' @param env a \code{SpatRaster} of at least one layer. if \emph{method = "strat.geo"}, a minimum of 2 layers are needed.
#' @param n.points \code{integer}. The number of background to be selected. \emph{Note that this number can change depending on the technique}
#' @param method \code{character}. one of: "rand.geo", "strat.geo", "rand.env" or "strat.env". \emph{see Details}. \emph{Default: "rand.geo"}.
#' @param digit.val.env \code{integer}. The number of digit to keep to remove too similar environmental values.
#' Only needed when method = "rand.env". \emph{Default: 1}.
#' @param res.grid.env \code{numeric}. Number of rows and columns to generate the grid (see Details). Only needed
#' when method = "rand.env". \emph{see Details}. \emph{Default: 100}.
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
#' reflect the environmental space. Then, a grid of res.grid.env*res.grid.env (default 100\*100) pixel is created. Note that 
#' increasing this value will strongly increase computation time. To reduce too similar points, 
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
#' \donttest{
#' library(terra)
#' env <- terra::unwrap(ESM_Env)
#' # Selection full random in the geographic space
#' Bp <- Bp_Sampling(env = env,
#'                   n.points = 1000,
#'                   method = "rand.geo",
#'                   To.plot = FALSE)
#'                        
#' # Selection stratified  in the geographic space                     
#' Bp <- Bp_Sampling(env = env,
#'                   n.points = 1000,
#'                   method = "strat.geo",
#'                   aggr.fact.geo = 2,
#'                   To.plot = FALSE)
#'                        
#' # Selection full random in the environmental space                     
#' Bp <- Bp_Sampling(env = env,
#'                   n.points = 1000,
#'                   method = "rand.env",
#'                   digit.val.env = 2,
#'                   res.grid.env = 100,
#'                   To.plot = FALSE)       
#'                                         
#' # Selection stratified in the environmental space                     
#' Bp <- Bp_Sampling(env = env,
#'                   n.points = 1000,
#'                   method = "strat.env",
#'                   n.strat.env = 3,
#'                   To.plot = FALSE)       
#' }               
#' @references 
#' Steen,B., Broennimann, O., Maiorano, L., Guisan, . 2024. How sensitive are species distribution models to different background point 
#' selection strategies? A test with species at various equilibrium levels. 
#' \emph{Ecological Modelling}. \bold{493}, 110754. \doi{10.1016/j.ecolmodel.2024.110754}.
#' @seealso \code{\link{ESM_Modeling}}
#' @export

Bp_Sampling <- function(env,
                        n.points = 10000,
                        method = "rand.geo",
                        digit.val.env = 1,
                        res.grid.env = 100,
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
      stop("when method == 'rand.env' or 'strat.env', a minimum of 2 environnmental layers is needed in env.")
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
    grid.env <- terra::rast(terra::ext(apply(env.score.round,2, range)), ncols = res.grid.env, nrows = res.grid.env,crs="") #vals = 1:(100*100)
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
#' The comparisons between SDM/ESM binary predictions are made following this formula: 
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
#' @seealso \code{\link{ESM_Projection}}, \code{\link{ESM_Ensemble.Projection}}, \code{\link{ESM_Threshold}}, \code{\link{ESM_Binarize}}
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
  if(length(names(proj.fut)) == length(unique(names(proj.fut)))){
    names(shift) = names(proj.fut)
    }
  results <- terra::freq(shift, wide=T)
  n.column <- ncol(results)
  while(ncol(results) < 5){
    
    results <- cbind(results,0)
    colnames(results)[ncol(results)] = setdiff(c("layer",0:3), colnames(results))[1]
  } 
  
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


#' @name ESM_Binarize
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @title Binarize probability values using a threshold
#' @description This function binarizes probability values based on a specific threshold
#' 
#' @param proj a \code{SpatRaster}, \code{data.frame}, \code{matrix}, or \code{numeric} containing the data to binarize. 
#' @param thr \code{numeric}. threshold to binarize the probabilities. \bold{must be a single value}. 
#' 
#' @details
#' \describe{
#' Probabilities strictly below the threshold will be assigned to 0 (
#' environmentally unsuitable for the species) while probabilities greater or 
#' equal to the threshold will be 1 (environementally suitable).
#' }
#' 
#' @return 
#' An object of the same class as proj: \code{SpatRaster}, \code{data.frame}, \code{matrix}, or \code{numeric}.
#' 
#' @seealso \code{\link{ESM_Projection}}, \code{\link{ESM_Ensemble.Projection}}, \code{\link{ESM_Threshold}} ,\code{\link{ESM_Range.Shift}}
#' 
#' @examples 
#' ## Generate a vector to binarize
#' proj <- seq(0,1000, by = 100)
#' # Binraization of the vector
#' ESM_Binarize(proj = proj, thr = 400)
#' 
#' @export

ESM_Binarize <- function(proj,
                         thr){
  
  if(length(thr)>1 | !is.numeric(thr)){
    stop("thr must contain a signle numeric")
  }
  if(inherits(proj,"SpatRaster")){
    
    rclmat <- matrix(c(0, thr, 0,thr, Inf, 1),
                     ncol=3, 
                     byrow=TRUE)
    
    proj.bin <- terra::classify(proj, rclmat, include.lowest=TRUE, right = FALSE)
    
  }else if(is.data.frame(proj) | is.matrix(proj) | is.numeric(proj)){
    
    proj.bin <- sapply(proj, FUN = function(x,thr){as.integer(x>=thr)}, thr = thr)
    
  }else{
    stop("proj must be either a SpatRaster, a data.frame, a matrix, or a numeric")
  }
  return(proj.bin)
}

#' @name ESM_Generate.ODMAP
#' @title Generates and fills ODMAP table
#' @author Flavien Collart \email{flaviencollart@hotmail.com}
#' @description
#' The function generates an ODMAP table to report your Modelling procedure 
#' using the outputs of \code{\link{ESM_Modeling}} and asks you some
#' questions to help you fill this table.
#' @param ESM.Mod The object returned by \code{\link{ESM_Modeling}}. \emph{Note that the object can also be NULL}
#' @param ESM.ensembleMod The object returned by \code{\link{ESM_Ensemble.Modeling}}. \emph{Note that the object can also be NULL}
#' @param ask.to.fill \code{logical}. If TRUE, teh functions will ask you 20 questions to fill some parts of the table.
#' @details
#' ODMAP (Overview, Data, Model, Assessment, Prediction) is a standard protocol
#' to report your study on SDMs proposed by Zurell et al (2020).This function
#' helps you to fill a part of this ODMAP table. However, some parts cannot be 
#' filled automatically. We've added a column in this table to help you to fill
#' the other lines on your own. 
#' 
#' @return a \code{data.frame} containing the 4 columns and 84 lines of the ODMAP table.
#' some values can be filled depending on the arguments. The fifth column is an help
#' to fill the different lines.
#' 
#' \code{\link{ESM_Ensemble.Modeling}} but also with some questions that will be
#' asked.
#' 
#' @seealso \code{\link{ESM_Ensemble.Modeling}}, \code{\link{ESM_Modeling}}
#' @examples
#' # A simple example where no values are filled in the table
#'  ODMAP_Table <- ESM_Generate.ODMAP(ESM.Mod = NULL,
#'                                    ESM.ensembleMod = NULL,
#'                                    ask.to.fill = FALSE)
#' # To see another example, see in ?ESM_Modeling
#'  
#' @references 
#' Zurell, D., Franklin, J., König, C., Bouchet, P.J., Dormann, C.F., Elith, J., 
#' Fandos, G., Feng, X., Guillera-Arroita, G., Guisan, A., Lahoz-Monfort, J.J., 
#' Leitão, P.J., Park, D.S., Peterson, A.T., Rapacciuolo, G., Schmatz, D.R., 
#' Schröder, B., Serra-Diaz, J.M., Thuiller, W., Yates, K.L., Zimmermann, N.E. 
#' and Merow, C. (2020), A standard protocol for reporting species distribution 
#' models. \emph{Ecography}, \bold{43}, 1261-1277. \doi{10.1111/ecog.04960}
#' @importFrom utils packageVersion
#' @export
ESM_Generate.ODMAP <- function(ESM.Mod = NULL,
                               ESM.ensembleMod = NULL,
                               ask.to.fill = TRUE){
  
  ## Creation of an empty table ####
  ODMAP <- as.data.frame(matrix("",ncol = 5, nrow = 84))
  colnames(ODMAP) = c("section","subsection","element","Value", "Information")
  ODMAP$section = c(rep("Overview",25),
                    rep("Data",33),
                    rep("Model",14),
                    rep("Assessment", 5),
                    rep("Prediction", 7))
  
  ODMAP$subsection = c(rep("Autorship",4),
                       rep("Model objective",2),
                       "Focal Taxon",
                       "Location",
                       rep("Scale of Analysis",5),
                       rep("Biodiversity data",2),
                       "Predictors","Hypotheses",
                       "Assumptions",
                       rep("Algorithms",3),
                       "Workflow",
                       rep("Software",3),
                       rep("Biodiversity data",12),
                       rep("Data partitioning",3),
                       rep("Predictor variables",10),
                       rep("Transfer data",8),
                       "Variable pre-selection",
                       "Multicollinearity",
                       rep("Model settings",2),
                       rep("Model estimates",3),
                       rep("Model selection",3),
                       rep("Analysis and Correction of non-independence",3),
                       "Threshold selection",
                       rep("Performance statistics",3),
                       rep("Plausibility check",2),
                       rep("Prediction output",2),
                       rep("Uncertainty quantification",5))
  
  ODMAP$element = c("Study title","Author names","Contact ","Study link",
                    "Model objective","Target output","Focal Taxon","Location",
                    "Spatial extent","Spatial resolution","Temporal extent",
                    "Temporal resolution","Boundary","Observation type",
                    "Response data type","Predictor types","Hypotheses",
                    "Model assumptions","Modelling techniques","Model complexity",
                    "Model averaging","Model workflow","Software","Code availability",
                    "Data availability","Taxon names","Taxonomic reference system",
                    "Ecological level","Data sources","Sampling design","Sample size",
                    "Clipping","Scaling","Cleaning","Absence data","Background data",
                    "Errors and biases","Training data","Validation data",
                    "Test data","Predictor variables","Data sources","Spatial extent",
                    "Spatial resolution","Coordinate reference system","Temporal extent",
                    "Temporal resolution","Data processing","Errors and biases",
                    "Dimension reduction","Data sources","Spatial extent",
                    "Spatial resolution","Temporal extent","Temporal resolution",
                    "Models and scenarios","Data processing","Quantification of Novelty",
                    "Variable pre-selection","Multicollinearity","Model settings (fitting)",
                    "Model settings (extrapolation)","Coefficients","Parameter uncertainty",
                    "Variable importance","Model selection","Model averaging","Model ensembles",
                    "Spatial autocorrelation","Temporal autocorrelation","Nested data",
                    "Threshold selection","Performance on training data",
                    "Performance on validation data","Performance on test data","Response shapes",
                    "Expert judgement","Prediction unit","Post-processing",
                    "Algorithmic uncertainty","Input data uncertainty",
                    "Parameter uncertainty","Scenario uncertainty","Novel environments")
  
  ## Add Information to certain complicated rows to fill
  ODMAP$Information = c(rep("",3),
                        "Link to study (e.g., DOI, web address)",
                        rep("",8),
                        "Natural, political and/or rectangle",
                        "citizen science, field survey, GPS tracking, range map and/or standardised montoring data",
                        "","climatic, edaphic, habitat and/or topographic",
                        "Hypotheses about species-environment relationships",
                        "Critical model assumptions (see table 2 in Zurell et al 2020. 10.1111/ecog.04960)",
                        "","Justification of the model complexity choice","","Complete conceptual description of the different modelling steps (fitting, evaluation, prediction)",
                        "", "Here add if possible a github link with you code", 
                        "Here add if possible the link to your data","",
                        "","", "Communties, individuals, OTU, populations and/or species",
                        "Details on species data source: e.g., URL/DOI, accession date, database version",
                        "Sampling design: spatial design: e.g. random, uniform, stratified), temporal design, nestedness",
                        "Country/region mask, if applicable","Details on scaling, if applicable: e.g., rasterisation of polygon maps, spatial and temporal thinning, measures to address spatial uncertainties",
                        "Details on data cleaning/filtering steps, if applicable: e.g., taxonomically, outlier presence/treatment",
                        "Details on absence collection if applicable", "Details on background data derivation, if applicable. See ESM::Bp.Sampling",
                        "Details on potential errors and biases in data, if applicable: e.g., detection probability, misidentification potential, geo-referencing errors, sampling bias",
                        rep("",3),
                        "Here write the list of available predictors  (before the selection of predictor)",
                        "Here you can write DOI or refs for your predictors",
                        rep("",4),"Temporal resolution of raw data, if applicable",
                        "Details on data processing and on spatial, temporal and thematic scaling: e.g. upscaling/downscaling, transformations, normalisations, thematic aggregations (e.g. of land cover classes), measures to address spatial uncertainties",
                        "Details on measurements errors and bias, when known",
                        "Details on dimension reduction of variable set, if applicable - if model-based, this should be contained in Model section (element: Details on pre-selection of variables)",
                        "Details on data sources: e.g., URL/DOI, accession date, database version",
                        "Spatial extent of transfer data", 
                        "Spatial resolution of transfer data", 
                        "Temporal extent of transfer data", 
                        "Temporal resolution of transfer data", 
                        "Models and scenarios used", 
                        "Details on data processing and scaling (see section Predictor variables)", 
                        "Quantification of novel environmental conditions and novel environmental combinations: e.g., distance to training data.  Here, if temporal projection or projection onto a new geographical extent, we recommend to perform ecospat::ecospat.climan or ecospat::ecopat.mess", 
                        "Details on pre-selection of variables, if applicable", "Methods for identifying and dealing with multicollinearity (Dormann, et al. 2013) or justification if multicollinearity is not explicitly dealt with",
                        "","", "Here, you could use ESM::ESM_VariableContributions" , 
                        "Details on quantification of parameter uncertainty, e.g. resampling", 
                        "Here, you could use ESM::ESM_VariableContributions" , 
                        "Model selection strategy: e.g. information-theoretic approach for variable selection, shrinkage and regularization", 
                        "Method for model averaging: e.g. derivation of weights", 
                        "Ensemble method: e.g. initial conditions (input data), model classes, model parameters, boundary conditions", 
                        "Method for addressing spatial autocorrelation in residuals", 
                        "Method for addressing temporal autocorrelation in residuals", 
                        "Method to account for nested data: e.g., fixed and random effects", 
                        "Here, you could use ESM::ESM_Threshold",
                        "","" , "", "Response plots, e.g. partial response plots, evaluation strips, inflated response plots. Check ESM::ESM_Response.Plot", 
                        "Expert judgements, e.g. map display","", "Post-processing, e.g. clipping, reprojection",
                        "Algorithmic uncertainity, if applicable", 
                        "Uncertainty in input data, if applicable", 
                        "Effect of parameter uncertainty, error propagation, if applicable", 
                        "Uncertainty in scenarios (e.g. climate models, land use models, storylines)", 
                        "Visualization/treatment of novel environments: e.g., masking). Here, if temporal projection or projection onto a new geographical extent, we recommend to perform ecospat::ecospat.climan or ecospat::ecopat.mess")
  colnames(ODMAP)[5] = "Information to correctly fill the 'Value' column"
  
  ## Start pre-fill the table based on the provided object ####
  ODMAP$Value[23] = paste(R.version.string, "With ESM package version", packageVersion("ESM"))
  if(!is.null(ESM.Mod)){
    ODMAP$Value[19] = paste(ESM.Mod$model.info$models, collapse = ", ")
    ODMAP$Value[26] = ESM.Mod$data$sp.name
    N.obs <- sum(ESM.Mod$data$resp)
    prevalence <- round(100*N.obs/length(ESM.Mod$data$resp),2)
    ODMAP$Value[31] = paste(N.obs,"occurrences were used in the models with a prevalence of",
                            prevalence,"%. This prevalence was afterwards set to", 100*ESM.Mod$model.info$prevalence,
                            "% in the models.")
    
    ODMAP$Value[41] = paste0(colnames(ESM.Mod$data$env.var),collapse = ", ")
    if(ESM.Mod$data$env.info$type == "SpatRaster"){
      ODMAP$Value[c(9,43)] = paste(names(as.vector(ESM.Mod$data$env.info$extent)),as.vector(ESM.Mod$data$env.info$extent), sep =" = ",collapse = ", ")
      ODMAP$Value[c(10,44)] = paste(c("xres = ", "yres = "),ESM.Mod$data$env.info$res, collapse = ", ")
      EPSG <- terra::crs(ESM.Mod$data$env.info$proj,describe=TRUE)[,c("authority","code")]
      proj4String <- terra::crs(ESM.Mod$data$env.info$proj,proj=TRUE)
      ODMAP$Value[45] = paste(paste0(EPSG,collapse = ":"),
                              "; PROJ-string = ",proj4String)
    }
    ODMAP$Value[60] = .printModelParameters(model.options = ESM.Mod$model.info$models.options,models = ESM.Mod$model.info$models)
    cv.method <- ESM.Mod$cv.method
    if(cv.method=="block"){
      cv.method <- "-fold"
    }
    ODMAP$Value[38] = paste("All bivariate models were evaluated using", ncol(ESM.Mod$cv.split.table)-1,ESM.Mod$cv.method,
                            "cross-validation.")
    if(cv.method == "split-sampling"){
      cv.ratio = round(sum(ESM.Mod$cv.split.table[,1])/nrow(ESM.Mod$cv.split.table),2)
      ODMAP$Value[38] = paste(ODMAP$Value[38],100*cv.ratio, "% were employed to calibrate the model and the remaining to evaluate it.")
    }
    if(ESM.Mod$model.info$pooling){
      ODMAP$Value[38] = paste(ODMAP$Value[38], "The pooling method as described in Collart & Guisan (2023. Ecol Inform) was afterwards applied to evaluate each model.")
    }
    ODMAP$Value[39] = ODMAP$Value[38]
    
    ODMAP$Value[74] = ODMAP$Value[73] = paste0(colnames(ESM.Mod$biva.calibration[[1]]), collapse = ", ")
    
  }
  if(!is.null(ESM.ensembleMod)){
    ODMAP$Value[67] = "weighted mean"
    if(nrow(ESM.ensembleMod$EF.algo$weights.algo)>1){
      ODMAP$Value[c(21,68)] = paste("The ensemble models were generated using a weighted mean, weighting each bivariate model by its", ESM.ensembleMod$weighting.score, "value for each modelling algorithm and removing all models having a performance <",ESM.ensembleMod$threshold,".")
      
    }else{
      ODMAP$Value[c(21,68)] = paste("The ensemble models were generated using a weighted mean, weighting each bivariate model by its", ESM.ensembleMod$weighting.score, "value for each modelling algorithmand removing all models having a performance <",ESM.ensembleMod$threshold,". An ensemble was then
                              realized between the algorithms by applying a weighted mean, weighting each ensemble by its", ESM.ensembleMod$weighting.score,"value.")
    }
  }
  
  ## Ask information to fill the cells ####
  if(ask.to.fill){
    cat("Starting the filing of ODMAP (~20 questions), please enter your answers without quote. You can simply let the value empty by pressing enter. If you want to stop the process but keep your data that you've already filled, write 'stop' (without the quote) during a question.")
    ## Overview ----
    ### Study ----
    x <- readline(prompt = "Enter your study title:")
    if(x == "stop"){
      return(ODMAP)
    }
    
    ODMAP$Value[1] = x
    
    ### Authors ----
    x <- readline(prompt = "Enter the authors:")
    if(x == "stop"){
      return(ODMAP)
    }
    
    ODMAP$Value[2] = x
    
    ### Contacts ----
    x <- readline(prompt = "Enter your email address:")
    if(x == "stop"){
      return(ODMAP)
    }
    
    ODMAP$Value[3] = x
    
    ### Objectives ----
    x <- readline(prompt = "What is your model objetive? write:\n\t1 for Inference and explanation\n\t2 for Mapping and Interpolation\n\t3 for Forecast and Transfer")
    if(x == "stop"){
      return(ODMAP)
    }
    x <- as.numeric(x)
    Objective <- x
    if(is.na(x) | x>3|x<=0){
      cat("\n Not a correct answer. This question will be passed.")
      Objective <- 1
    }else{
      ODMAP$Value[5] = c("Inference and explanation","Mapping and Interpolation","Forecast and Transfer")[x]
    }
    
    if(Objective>1){
      x <- readline(prompt = "What is your main target output? write:\n\t1 for suitable vs. unsuitable habitat\n\t2 continuous habitat suitability index")
      if(x == "stop"){
        return(ODMAP)
      }
      x <- as.numeric(x)
      if(is.na(x) | x>2|x<=0){
        cat("\n Not a correct answer. This question will be passed.")
      }else{
        ODMAP$Value[6] = c("suitable vs. unsuitable habitat","continuous habitat suitability index")[x]
      }
      
      
    }
    
    ## Focal Taxon
    x <- readline(prompt = "What is your focal taxon? (If you want to stop the process, write stop)")
    if(x == "stop"){
      return(ODMAP)
    }
    ODMAP$Value[7] = x
    
    ## Study area
    x <- readline(prompt = "What is your study area? (e.g., Switzerland, Europe)")
    if(x == "stop"){
      return(ODMAP)
    }
    ODMAP$Value[8] = x
    
    ##Temporal extent
    x <- readline(prompt = "What is your temporal extent of your observation data? (e.g, 1981-2010)")
    if(x == "stop"){
      return(ODMAP)
    }
    ODMAP$Value[11] = x
    
    ## Boundary
    x <- readline(prompt = "What is the boundary type of your study area? write:\n\t1 for Natural \n\t2 for Political \n\t3 for Rectangle")
    if(x == "stop"){
      return(ODMAP)
    }
    x <- as.numeric(x)
    if(is.na(x) | x>3|x<=0){
      cat("\n Not a correct answer. This question will be passed.")
    }else{
      ODMAP$Value[13] = c("Natural","Political","Rectangle")[x]
    }
    
    ## Observation data
    x <- readline(prompt = "What is the type of your observation data? Allow several answers (e.g, write 13 to select citizen science data and GPS tracking)\n\t1 citizen science \n\t2 field survey \n\t3 GPS tracking \n\t4 range map \n\t5 standardised monitoring data")
    if(x == "stop"){
      return(ODMAP)
    }
    if(is.na(x)){
      cat("\n Not a correct answer. This question will be passed.")
    }else{
      ToKeep <- sapply(as.character(1:5),grep, x= x)==1
      ToKeep <- !is.na(ToKeep)
      if(sum(ToKeep)==0){
        cat("\n Not a correct answer. This question will be passed.")
        
      }else{
        ODMAP$Value[14] = paste(c("citizen science", "field survey", "GPS tracking", "range map", "standardised monitoring data")[ToKeep],collapse = ", ")
      }
      
    }
    
    ## Response variables
    x <- readline(prompt = "What is the type of your response data?\n\t1 for presence/absence \n\t2 for presence-only")
    if(x == "stop"){
      return(ODMAP)
    }
    x <- as.numeric(x)
    if(is.na(x) | x>2 | x<0){
      cat("\n Not a correct answer. This question will be passed.")
    }else{
      ODMAP$Value[15]  = c("presence/absence","presence-only")[x]
    }
    
    ##Predictor type
    
    x <- readline(prompt = "What are the categories of your predictors? Allow several answers (e.g, write 14 to select climate and topographic predictors)\n\t1 climate \n\t2 edaphic \n\t3 habitat (i.e, land-cover) \n\t4 topographic")
    if(x == "stop"){
      return(ODMAP)
    }
    if(is.na(x)){
      cat("\n Not a correct answer. This question will be passed.")
    }else{
      ToKeep <- sapply(as.character(1:4),grep, x = x)==1
      ToKeep <- !is.na(ToKeep)
      if(sum(ToKeep)==0){
        cat("\n Not a correct answer. This question will be passed.")
        
      }else{
        ODMAP$Value[16] = paste(c("climate", "edaphic", "habitat", "topographic")[ToKeep],collapse = ", ")
      }
      
    }
    
    
    x <- readline(prompt = "What is the ecological level of your observation data? Allow several answers (e.g., write 23 for Individuals and OTUs) \n\t1 Communities \n\t2 Individual \n\t3 OTUs \n\t4 Population \n\t5 Species")
    if(x == "stop"){
      return(ODMAP)
    }
    if(is.na(x)){
      cat("\n Not a correct answer. This question will be passed.")
    }else{
      ToKeep <- sapply(as.character(1:5),grep, x = x)==1
      ToKeep <- !is.na(ToKeep)
      if(sum(ToKeep)==0){
        cat("\n Not a correct answer. This question will be passed.")
        
      }else{
        ODMAP$Value[28] = paste(c("Communties", "Individuals", "OTU", "Populations", "species")[ToKeep],collapse = ", ")
      }
      
    }
    
    x <- readline(prompt = "Did you clip your observation data to a single region? \n\t1 for Yes \n\t2 for False \n\t stop to stop")
    if(x == "stop"){
      return(ODMAP)
    }
    x <- as.numeric(x)
    if(is.na(x) | x>2 | x<0){
      cat("\n Not a correct answer. This question will be passed.")
    }else{
      ODMAP$Value[32]  = c("The observation data were clipped","No clipping was made")[x]
    }
    
    
    x <- readline(prompt = "What was the temporal extent of your predictor data for the modeling? (e.g., 1981-2010)")
    if(x == "stop"){
      return(ODMAP)
    }
    ODMAP$Value[46] = x
    if(Objective == 3){
      x <- readline(prompt = "If you projected your models, what was the other temporal extent? (e.g, 2041-2070)")
      if(x == "stop"){
        return(ODMAP)
      }
      ODMAP$Value[54] = x
      
      x <- readline(prompt = "What were the scenarios? (e.g, UKESM ssp 5-8.5)")
      if(x == "stop"){
        return(ODMAP)
      }
      ODMAP$Value[56] = x
    }
    
  }
  cat("\n\n Your ODMAP have been generated with your information. Please try to fill the other cells. We added help to fill the different sections that could not be filled.")
  return(ODMAP)
}

.printModelParameters <- function(model.options, models){
  ToReturn <- NULL
  for(i in 1:length(models)){
    if(models[i] == "GLM"){
      if(model.options$GLM$test=="none"){
        sentence <- paste("GLM was performed under the",model.options$GLM$family$family, "family with",model.options$GLM$type,
                          "terms")
      }else{
        sentence <- paste("GLMs were performed under the",model.options$GLM$family$family, "family with",model.options$GLM$type,
                          "terms, selecting the best structure with a stepwise AIC")
      }
      ToReturn <- c(ToReturn, 
                    sentence)
    }else if(models[i] == "ANN"){
      sentence <- paste0("ANNs were performed with ", model.options$ANN$size, " units in the hidden layers, a weight decay of ", model.options$ANN$decay,
                         ", an initial random weight of ", model.options$ANN$rang, ", and a maximum of ", model.options$ANN$maxit,
                         " iterations, all other parameters was set as default")
      ToReturn <- c(ToReturn, 
                    sentence)
    }else if(models[i] == "CTA"){
      # xval = 5, minbucket = 5, minsplit = 5, cp = 0, maxdepth =25
      sentence <- paste0("CTAs were performed with a minimum of ", model.options$CTA$control$minsplit,
                         " observations existing in a node in order for a split to be attempted, a minimum of ",
                         model.options$CTA$control$minbucket, " observations in a terminal node, a complexity of ",
                         model.options$CTA$control$cp,", ", model.options$CTA$control$maxcompete,
                         " competitor splits retained in the output, ", model.options$CTA$control$maxsurrogate,
                         " surrogate splits retained in the output, a maximal depth of ", model.options$CTA$control$maxdepth,
                         ", and ",model.options$CTA$control$xval," internal cross-validations, all other parameters were set as default")
      ToReturn <- c(ToReturn, 
                    sentence)
      
    }else if(models[i] == "GAM"){
      sentence <- paste0("GAMs were performed under a ", model.options$GAM$family," distribution with ",
                         model.options$GAM$smooth.k, " dimension(s) used to represent the smooth term, ", model.options$GaM$smooth.bs,
                         " as a penalized smoothing basis, ", paste(model.options$GaM$optimizer,collapse = " and "), 
                         " as the optimization method(s) to use to optimize the smoothing parameter estimation criterion, all other parameters were set as default.")
      ToReturn <- c(ToReturn, 
                    sentence)
    }else if(models[i] == "GBM"){
      sentence <- paste0("GBMs were performed under a ", model.options$GBM$distribution," distribution with ",
                         model.options$GBM$n.trees, " trees, a maximum interaction depth of ", model.options$GBM$interaction.depth,
                         ", a minimum of ", model.options$GBM$n.minobsinnode," observation in the termal nodes of trees, a learning rate of ",
                         model.options$GBM$shrinkage,", a fraction of ", model.options$GBM$bag.fraction,
                         " of the training set observations randomly selected to propose the next tree in the expansion, a first fraction of ",
                         model.options$GBM$train.fraction, " to fit the model and ",model.options$GBM$cv.folds," internal cross-vaidations, all other parameters was set as default.")
      ToReturn <- c(ToReturn, 
                    sentence)
    }else{
      ToReturn <- c(ToReturn, 
                    "MAXNET were performed with all default parameters.")
    }
  }
  ToReturn <- paste(ToReturn,collapse = "; ")
  return(ToReturn)
}

