#' @name ESM_Models.Options
#' @author Flavien Collart from the code of biomod2::BIOMOD_ModelingOptions
#' 
#' @title Model parameters for ESMs
#' 
#' @description Generate a list of model parameters
#' 
#' @param ANN a \code{list}. with the objects:
#' \itemize{
#' \item{size}: an \code{integer}. The numer of units in the hidden layers. \emph{Default}: 8.
#' \item{decay}: a \code{numeric}. The weight decay. \emph{Default}: 0.001.
#' \item{rang}: a \code{numeric}. Initial random weights. \emph{Default}: 0.1.
#' \item{maxit}: an \code{integer}.The maximal number of iterations. \emph{Default}: 200.
#' }
#' 
#' @param CTA a \code{list}. with the objects:
#' \itemize{
#' \item{na.action}: The default action to remove observations with NA in the response variable. \emph{Default}: \code{\link[rpart]{na.rpart}}
#' \item{method}: a \code{character}. Only "class" is available
#' \item{x}: a \code{logical}. Keep a copy of the x matrix. \emph{Default}: FALSE.
#' \item{y}: a \code{logical}. Keep a copy of the y matrix. \emph{Default}: TRUE.
#' \item{model}: a \code{logical}. Keep a copy of the model frame. \emph{Default}: FALSE.
#' \item{control}: a \code{list}. a list of options obtained with \code{\link[rpart]{rpart.control}}. \emph{Default}: $xval = 5, minbucket = 5, minsplit = 5, cp = 0, maxdepth =25. for the other parameters see \code{\link[rpart]{rpart.control}} 
#' }
#' 
#' @param GAM a \code{list}. with the objects:
#' \itemize{
#' \item{smooth.k}: \code{integer}. The dimension of the basis used to represent the smooth term. \emph{Default}: 4.
#' \item{smooth.bs}: \code{character}. The (penalized) smoothing basis to use.  \emph{Default}: 'tp'. see \code{\link[mgcv]{smooth.terms}}.
#' \item{method}: \code{character}. The smoothing parameter estimation method.  \emph{Default}: 'GCV.Cp'.
#' \item{optimizer}: \code{character}. An array specifying the numerical optimization method to use to optimize the smoothing parameter estimation criterion. \emph{Default}: c('outer','newton').
#' \item{family}: error distribution family. Only \code{stats}[binomial(link = 'logit')] was tested thus we discourage to change this parameter to another distribution family.
#' \item{control}: a \code{list}. a list of options obtained with \code{\link[mgcv]{gam.control}}.
#' \item{scale}: \code{numeric}. scale parameter. \emph{Default}: 0. If negative, it means that the scale parameter is unknown. 0 means 
#' that the scale is equal to 1 for a binomial distribution.
#' \item{select}: \code{logical}. should gam add an extra penalty to each term so that it can be penalized to zero? Use \emph{gamma} to increase level of penalization. \emph{Default}: FALSE.
#' \item{gamma}: \code{numerical}. Increasing this number produce smoother models (but could also produce overfitting issues). \emph{Default}: 1.
#' \item{knots}: \code{list}. knot values to be used for basis construction. \emph{Default}: \code{NULL}.
#' \item{H}: Fixed quadratic penalty on the parameters. \emph{Default}: \code{NULL}.
#' }
#' 
#' @param GLM a \code{list}. with the objects:
#' \itemize{
#' \item{type}: \code{character}. Either "linear", "quadratic" or "cubic". \emph{Default}: "quadratic".
#' \item{test}: \code{character}. Either "AIC" or "none". Perform or not a step AIC selection for the formula. \emph{Default}: "none".
#' \item{family}: error distribution family. Only \code{stats}[binomial(link = 'logit')] was tested thus we discourage to change this parameter to another distribution family.
#' }
#' 
#' @param GBM a \code{list}. with the objects:
#' \itemize{
#' \item{n.trees}: \code{integer}. Number of trees to fit. \emph{Default}: 1000.
#' \item{interaction.depth}: \code{integer}. The maximum depth of each tree \code{numeric}. Number of trees to fit. \emph{Default}: 4.
#' \item{n.minobsinnode}: \code{integer}. The minimum number of observations in the terminal nodes of the trees. \emph{Default}: 5.
#' \item{shrinkage}: \code{numeric}. The learning rate   \emph{Default}: 0.005.
#' \item{bag.fraction}: \code{numeric}. The fraction of the training set observations randomly selected to propose the next tree in the expansion. \emph{Default}: 0.5.
#' \item{train.fraction}: \code{numeric}. The first fraction of observations to fit the model. \emph{Default}: 1.
#' \item{cv.folds}: \code{integer}. Number of cross-validations to perform. \emph{Default}: 3.
#' \item{verbose}: \code{logical}. print or not the progress. \emph{Default}: FALSE.
#' \item{n.cores}: \code{logical}. Number of CPU cores to use. \emph{Default}: \code{NULL}.
#' }
#' 
#' @details For the arguments of each modeling technique, please refer to the manual of \code{\link[nnet]{nnet}}, \code{\link[rpart]{rpart}}, \code{\link[mgcv]{gam}}, \code{\link[stats]{glm}}, and \code{\link[gbm]{gbm}}.
#' @return a \code{list} of parameters for ESM.
#' @examples 
#' ## Perform a GLM with step AIC to select the best structure 
#' ## and allows linear, quadratic and cubic terms
#' models.options = ESM_Models.Options(GLM=list(test="AIC",type="cubic"))
#' ## Perform a GLM with linear and quadratic terms 
#' ## but does not perform a selection of the structure
#' models.options = ESM_Models.Options(GLM=list(test="none",type="quadratic"))
#' 
#' @seealso \code{\link{ESM_Modeling}}, \code{\link{ESM_Ensemble.Modeling}}, \code{\link{ESM_Projection}}, \code{\link{ESM_Ensemble.Projection}}
#' 
#' @export

ESM_Models.Options <- function(ANN = NULL,
                               CTA = NULL,
                               GAM = NULL,
                               GLM = NULL,
                               GBM = NULL){
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
              GAM=list(smooth.k = 4,
                       smooth.bs = "tp",
                       method = "GCV.Cp",
                       optimizer = c('outer','newton'),
                       family = stats::binomial(link = 'logit'),
                       control = mgcv::gam.control(),
                       scale = 0,
                       select = FALSE,
                       gamma = 1,
                       knots = NULL,
                       H = NULL
                       ),
              GLM = list(type = 'quadratic',
                         test = 'none',
                         family = stats::binomial(link = 'logit')
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
  
  ## ANN ----
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
  ## CTA----
  if (!is.null(CTA)){
    if(!is.null(CTA$na.action)){
      opt$CTA$na.action = CTA$na.action
    }
    if(!is.null(CTA$method)){
      if(CTA$method != "class"){
        stop("CTA$method should be 'class'.")
      }
      opt$CTA$method = CTA$method
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
  
  ## GAM ----
  if(!is.null(GAM)){
    if(!is.null(GAM$smooth.k)){
      if(!is.numeric(GAM$smooth.k)){
        stop("GAM$smooth.k should be an integer")
      }
      opt$GAM$smooth.k = GAM$smooth.k
    }
    if(!is.null(GAM$smooth.bs)){
      if(!is.character(GAM$smooth.bs) | length(GAM$smooth.bs) > 1 | 
         (!any(GAM$smooth.bs %in% c('tp', 'ts', 'ds', 'cr', 'cs',
                                    'cc', 'sos','bs','ps','cp','re',
                                    'mrf','gp','so','sf','sw', 'te',
                                    'ti','t2','ad','sz', 'fs')
               )
                                                                     )
         ){
        stop("GAM$smooth.bs should be a character and must be either 'tp', 'ts', 'ds', 'cr', 'cs', 
             'cc', 'sos','bs','ps','cp','re','mrf','gp','so','sf','sw', 'te','ti','t2','ad','sz' or 'fs'. ")
      }
      opt$GAM$smooth.bs = GAM$smooth.bs
    }
    if(!is.null(GAM$method)){
      if(!is.character(GAM$method)| length(GAM$method) > 1 | 
         (!any(GAM$method %in% c('GCV.Cp','GACV.Cp','NCV','QNCV','REML','P-REML','ML','P-ML')))
         ){
        stop("GAM$method must be a character and must be either 'GCV.Cp','GACV.Cp','NCV','QNCV','REML','P-REML','ML', or'P-ML'")
      }
      opt$GAM$method = GAM$method
    }
    if(!is.null(GAM$optimizer)){
      if(!all(GAM$optimizer %in% c('outer','newton','bfgs','optim','nlm','efs'))){
        stop("GAM$optimizer must be a character array from c('outer','newton','bfgs','optim','nlm','efs') ")
      }
      opt$GAM$optimizer = GAM$optimizer
    }
    if (!is.null(GAM$family)) {
      opt$GAM$family <- GAM$family
    }
    if(!is.null(GAM$control)){
      if(!is.list(GAM$control)){
        stop("GAM$control must be a list. Please use the function mgcv::gam.control to generate this object")
      }
      cat("\n gam$control has been changed, Please note this parameter is not checked")
      opt$GAM$control = GAM$control
    }
    if(!is.null(GAM$scale)){
      if(!is.numeric(GAM$scale)){
        stop("GAM$scale should be an integer")
      }
      opt$GAM$scale = GAM$scale
    }
    if(!is.null(GAM$select)){
      if(!is.logical(GAM$select)){
        stop("CTA$select must be logical")
      }
      opt$GAM$select = GAM$select 
    }
    if(!is.null(GAM$gamma)){
      if(!is.numeric(GAM$gamma)){
        stop("GAM$gamma must be a numeric")
      }
      opt$GAM$gamma = GAM$gamma
    }
    if(!is.null(GAM$knots)){
      if(!is.list(GAM$knots)){
        stop("GAM$knots must be a list")
      }
      opt$GAM$knots = GAM$knots
    }
    if(!is.null(GAM$H)){
      opt$GAM$H = GAM$H
    }
    
  }
  
  ## GLM ----
  if (!is.null(GLM)) {
    if (!is.null(GLM$type)) {
      if(!(any(GLM$type %in% c("linear","quadratic","cubic")))){
        stop("GLM$type should be either 'linear', 'quadratic' or 'cubic'")
      }
      opt$GLM$type <- GLM$type
    }
    if (!is.null(GLM$family)) {
      opt$GLM$family <- GLM$family
    }
    if (!is.null(GLM$test)) {
      if(!(any(GLM$test %in% c("none","AIC")))){
        stop("GLM$type should be either 'none' or 'AIC'")
      }
      opt$GLM$test <- GLM$test
    }
  }
  ## GBM ----
  if (!is.null(GBM)) {
    if (!is.null(GBM$distribution)) {
      opt$GBM$distribution <- GBM$distribution
    }
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




