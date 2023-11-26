ESM_Models.options <- function(GLM = NULL,
                               # GAM = NULL,
                               GBM = NULL
                               ){
  ## Default options
  opt <- list(GLM = list(type = 'quadratic',
                         myFormula = NULL,
                         test = 'none',
                         family = binomial(link = 'logit')),
              GBM = list(distribution = 'bernoulli',
                         n.trees = 1000,
                         interaction.depth = 4,
                         n.minobsinnode = 5,
                         shrinkage = 0.005,
                         bag.fraction = 0.5,
                         train.fraction = 1,
                         cv.folds = 3,
                         verbose = FALSE,
                         n.cores = 1
                          ))
  
  fam_GLM = c("binomial", "gaussian", "Gamma", "inverse.gaussian", 
                        "poisson", "quasi", "quasibinomial", "quasipoisson")
  if (!is.null(GLM)) {
    if (!is.null(GLM$type)) {
      opt$GLM$type <- GLM$type
    }
    if (!is.null(GLM$myFormula)) {
      opt$GLM$myFormula <- GLM$myFormula
    }
    if (!is.null(GLM$test)) {
      opt$GLM$test <- GLM$test
    }
    if (!is.null(GLM$family)) {
      fam.test <- TRUE
      if (inherits(GLM$family, "family")) {
        opt$GLM$family <- GLM$family
      }
      else if (is.character(GLM$family)) {
        if (!unlist(strsplit(GLM$family, "[/(]"))[1] %in% 
            fam_GLM) {
          fam.test <- FALSE
        }
        if (grepl(")", GLM$family)) {
          opt$GLM$family <- eval(parse(text = GLM$family))
        }
        else {
          opt$GLM$family <- eval(parse(text = paste0(GLM$family, 
                                                     "()")))
        }
      }
      else {
        fam.test <- FALSE
      }
      if (!fam.test) {
        cat("\n!!! invalid GLM$family given -> binomial(link = 'logit') was automatically set up !!!")
        opt$GLM$family <- binomial(link = "logit")
      }
    }
  }
  if (!is.null(GBM)) {
    if (!is.null(GBM$distribution)) {
      opt$GBM$distribution <- GBM$distribution
    }
    if (!is.null(GBM$n.trees)) {
      opt$GBM$n.trees <- GBM$n.trees
    }
    if (!is.null(GBM$interaction.depth)) {
      opt$GBM$interaction.depth <- GBM$interaction.depth
    }
    if (!is.null(GBM$n.minobsinnode)) {
      opt$GBM$n.minobsinnode <- GBM$n.minobsinnode
    }
    if (!is.null(GBM$shrinkage)) {
      opt$GBM$shrinkage <- GBM$shrinkage
    }
    if (!is.null(GBM$bag.fraction)) {
      opt$GBM$bag.fraction <- GBM$bag.fraction
    }
    if (!is.null(GBM$train.fraction)) {
      opt$GBM$train.fraction <- GBM$train.fraction
    }
    if (!is.null(GBM$cv.folds)) {
      opt$GBM$cv.folds <- GBM$cv.folds
    }
    if (!is.null(GBM$verbose)) {
      opt$GBM$verbose <- GBM$verbose
    }
    if (!is.null(GBM$n.cores)) {
      opt$GBM$n.cores <- GBM$n.cores
    }
    else {
      opt$GBM$n.cores <- 1
    }
  }
  
  ####
  
  
  return(opt)
  
}




