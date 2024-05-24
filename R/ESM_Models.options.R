
#' @export
ESM_Models.Options <- function(ANN = NULL,
                               CTA = NULL,
                               GLM = NULL,
                               GBM = NULL
                               ){
  ## Default options
  opt <- list(ANN = list(size = 8,
                         decay = 0.001,
                         rang = 0.1,
                         maxit = 200
                         ),
              CTA = list(na.action= rpart::na.rpart,
                         method = "class",
                         x = FALSE,
                         y = TRUE,
                         model = FALSE,
                         control = list(xval = 5, 
                                        minbucket = 5,
                                        minsplit = 5,
                                        cp = 0, 
                                        maxdepth = 25,
                                        maxcompete = 4,
                                        maxsurrogate = 5,
                                        usesurrogate = 2,
                                        surrogatestyle = 0)
                         ),
              GLM = list(type = 'quadratic',
                         myFormula = NULL,
                         test = 'none',
                         family = binomial(link = 'logit')
                         ),
              GBM = list(distribution = 'bernoulli',
                         n.trees = 1000,
                         interaction.depth = 4,
                         n.minobsinnode = 5,
                         shrinkage = 0.005,
                         bag.fraction = 0.5,
                         train.fraction = 1,
                         cv.folds = 3,
                         verbose = FALSE,
                         n.cores = NULL
                          )
              )
  
  ### Adapted Code from biomod2 4.2.4
  if (!is.null(ANN)){
    if(!is.null(ANN$size)){
      if(!is.integer(ANN$size)){
        stop("ANN$size should be an integer")
      }
      opt$ANN$size = ANN$size
    }
    if(!is.null(ANN$decay)){
      if(!is.numeric(ANN$decay)){
        stop("ANN$decay should be a numeric")
      }
      opt$ANN$decay = ANN$decay 
    }
    if(!is.null(ANN$rang)){
      if(!is.numeric(ANN$rang)){
        stop("ANN$rang should be a numeric")
      }
      opt$ANN$rang = ANN$rang 
    }
    if(!is.null(ANN$maxit)){
      if(!is.integer(ANN$maxit)){
        stop("ANN$maxit should be a numeric")
      }
      opt$ANN$maxit = ANN$maxit 
    }
  }
  if (!is.null(CTA)){
    if(!is.null(CTA$na.action)){
      opt$CTA$na.action = CTA$na.actio
    }
    if(!is.null(CTA$method)){
      if(CTA$method != "class"){
        stop("CTA$method should be 'class'.")
      }
    }
    if(!is.null(CTA$x)){
      if(!is.logical(CTA$x)){
        stop("CTA$x should be logical")
      }
      opt$CTA$x = CTA$x 
    }
    if(!is.null(CTA$y)){
      if(!is.logical(CTA$y)){
        stop("CTA$y should be logical")
      }
      opt$CTA$y = CTA$y 
    }
    if(!is.null(CTA$model)){
      if(!is.logical(CTA$model)){
        stop("CTA$model should be logical")
      }
      opt$CTA$model = CTA$model 
    }
    if(!is.null(CTA$control)){
      if(!is.list(CTA$control) | length(CTA$control) != 9 ){
        stop("CTA$control should be a list with 9 elements. Please use the function rpart::rpart.control to generate this object")
      }
      cat("\n CTA$control has been changed, Please note this parameter is not checked")
      opt$CTA$control = CTA$control
    }
  }
  if (!is.null(GLM)) {
    if (!is.null(GLM$type)) {
      if(!(any(GLM$type %in% c("linear","quadratic","polynomial")))){
        stop("GLM$type should be either 'linear', 'quadratic' or 'polynomial'")
      }
      opt$GLM$type <- GLM$type
    }
    if (!is.null(GLM$myFormula)) {
      opt$GLM$myFormula <- GLM$myFormula
    }
    if (!is.null(GLM$test)) {
      if(!(any(GLM$test %in% c("none","AIC")))){
        stop("GLM$type should be either 'none' or 'AIC'")
      }
      opt$GLM$test <- GLM$test
    }
  }
  if (!is.null(GBM)) {
    if (!is.null(GBM$n.trees)) {
      if(!is.integer(GBM$n.trees)){
        stop("GBM$n.trees should be an integer")
      }
      opt$GBM$n.trees <- GBM$n.trees
    }
    if (!is.null(GBM$interaction.depth)) {
      if(!is.integer(GBM$interaction.depth)){
        stop("GBM$interaction.depth should be an integer")
      }
      opt$GBM$interaction.depth <- GBM$interaction.depth
    }
    if (!is.null(GBM$n.minobsinnode)) {
      if(!is.integer(GBM$n.minobsinnode)){
        stop("GBM$n.minobsinnode should be an integer")
      }
      opt$GBM$n.minobsinnode <- GBM$n.minobsinnode
    }
    if (!is.null(GBM$shrinkage)) {
      if(!is.numeric(GBM$shrinkage)){
        stop("GBM$shrinkage should be a numeric")
      }
      opt$GBM$shrinkage <- GBM$shrinkage
    }
    if (!is.null(GBM$bag.fraction)) {
      if(!is.numeric(GBM$bag.fraction | GBM$bag.fraction>1 | GBM$bag.fraction<0)){
        stop("GBM$bag.fraction should be a numeric and comprised between 0 and 1")
      }
      opt$GBM$bag.fraction <- GBM$bag.fraction
    }
    if (!is.null(GBM$train.fraction)) {
      if(!is.numeric(GBM$train.fraction | GBM$train.fraction>1 | GBM$train.fraction<0)){
        stop("GBM$train.fraction should be a numeric and comprised between 0 and 1")
      }
      opt$GBM$train.fraction <- GBM$train.fraction
    }
    if (!is.null(GBM$cv.folds)) {
      if(!is.integer(GBM$cv.folds)){
        stop("GBM$cv.folds should be an integer")
      }
      opt$GBM$cv.folds <- GBM$cv.folds
    }
    if (!is.null(GBM$verbose)) {
      if(!is.logical(GBM$verbose)){
        stop("GBM$v should be logical")
      }
      opt$GBM$verbose <- GBM$verbose
    }
    if (!is.null(GBM$n.cores)) {
      if(!is.integer(GBM$n.cores)){
        stop("GBM$n.cores should be an integer")
      }
      opt$GBM$n.cores <- GBM$n.cores
    }
  }
  ####
  
  return(opt)
  
}




