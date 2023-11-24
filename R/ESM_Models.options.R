ESM_Models.options <- function(GLM = NULL
                               # GAM = NULL,
                               # GBM = NULL
                               ){
  opt <- new("ESM.models.options")
  fam_GLM = c("binomial", "gaussian", "Gamma", "inverse.gaussian", 
                        "poisson", "quasi", "quasibinomial", "quasipoisson")
  if (!is.null(GLM)) {
    if (!is.null(GLM$type)) {
      opt@GLM$type <- GLM$type
    }
    if (!is.null(GLM$interaction.level)) {
      opt@GLM$interaction.level <- GLM$interaction.level
    }
    if (!is.null(GLM$myFormula)) {
      opt@GLM$myFormula <- GLM$myFormula
    }
    if (!is.null(GLM$test)) {
      opt@GLM$test <- GLM$test
    }
    if (!is.null(GLM$family)) {
      fam.test <- TRUE
      if (inherits(GLM$family, "family")) {
        opt@GLM$family <- GLM$family
      }
      else if (is.character(GLM$family)) {
        if (!unlist(strsplit(GLM$family, "[/(]"))[1] %in% 
            fam_GLM) {
          fam.test <- FALSE
        }
        if (grepl(")", GLM$family)) {
          opt@GLM$family <- eval(parse(text = GLM$family))
        }
        else {
          opt@GLM$family <- eval(parse(text = paste0(GLM$family, 
                                                     "()")))
        }
      }
      else {
        fam.test <- FALSE
      }
      if (!fam.test) {
        cat("\n!!! invalid GLM$family given -> binomial(link = 'logit') was automatically set up !!!")
        opt@GLM$family <- binomial(link = "logit")
      }
    }
  }
  
}


setClass("ESM.models.options",
         representation(GLM = "list"
                        # GBM = "list",
                        # GAM = "list"
         ),
         prototype(GLM = list(type = 'quadratic',
                              interaction.level = 0,
                              myFormula = NULL,
                              test = 'none',
                              family = binomial(link = 'logit'))
                   # GBM = list(distribution = 'bernoulli',
                   #            n.trees = 2500,
                   #            interaction.depth = 7,
                   #            n.minobsinnode = 5,
                   #            shrinkage = 0.001,
                   #            bag.fraction = 0.5,
                   #            train.fraction = 1,
                   #            cv.folds = 3,
                   #            keep.data = FALSE,
                   #            verbose = FALSE,
                   #            # class.stratify.cv = 'bernoulli',
                   #            perf.method = 'cv',
                   #            n.cores = 1),
                   # GAM = list(algo = "GAM_mgcv",
                   #            type = "s_smoother",
                   #            k = NULL,
                   #            interaction.level = 0,
                   #            myFormula = NULL,
                   #            family = binomial(link = 'logit'),
                   #            control = list(epsilon = 1e-06, trace = FALSE, maxit = 100),
                   #            method = "GCV.Cp",
                   #            optimizer = c("outer", "newton"),
                   #            select = FALSE,
                   #            knots = NULL,
                   #            paraPen = NULL)
         ),
         validity = function(object) {
           # test <- TRUE
           # 
           # ## GLM ------------------------------------------------------------
           # test <- .fun_testIfIn(test, "GLM$type", object@GLM$type,
           #                       c("simple", "quadratic", "polynomial", "user.defined"))
           # test <- .fun_testIfPosInt(test, "GLM$interaction.level", object@GLM$interaction.level)
           # 
           # if (!is.null(object@GLM$myFormula) && 
           #     !inherits(object@GLM$myFormula, "formula")) {
           #   cat("\nGLM$myFormula must be NULL or a formula object")
           #   test <- FALSE
           # }
           # 
           # test <- .fun_testIfIn(test, "GLM$test", object@GLM$test, c("AIC", "BIC", "none"))
           # fam <- "none"
           # if (!inherits(object@GLM$family, "family")) {
           #   cat("\nGLM$family must be a valid family object")
           #   test <- FALSE
           # }
           # if (!is.list(object@GLM$control)) {
           #   cat("\nGLM$control must be a list like that returned by glm.control")
           #   test <- FALSE
           # }
           # 
           # ## GBM ------------------------------------------------------------
           # test <- .fun_testIfIn(test, "GBM$distribution", object@GBM$distribution, 
           #                       c("bernoulli", "huberized", "multinomial", "adaboost"))
           # # test <- .fun_testIfPosInt(test, "GBM$n.trees", object@GBM$n.trees)
           # if (!is.numeric(object@GBM$n.trees)) {
           #   cat("\nGBM$n.trees must be a integer")
           #   test <- FALSE
           # } else {
           #   if (object@GBM$n.trees < 0 | floor(object@GBM$n.trees) != object@GBM$n.trees) {
           #     cat("\nGBM$n.trees must be a positive integer")
           #     test <- FALSE
           #   }
           # }
           # test <- .fun_testIfPosInt(test, "GBM$interaction.depth", object@GBM$interaction.depth)
           # test <- .fun_testIfPosInt(test, "GBM$n.minobsinnode", object@GBM$n.minobsinnode)
           # test <- .fun_testIfPosNum(test, "GBM$shrinkage", object@GBM$shrinkage)
           # test <- .fun_testIf01(test, "GBM$bag.fraction", object@GBM$bag.fraction)
           # test <- .fun_testIf01(test, "GBM$train.fraction", object@GBM$train.fraction)
           # test <- .fun_testIfPosInt(test, "GBM$cv.folds", object@GBM$cv.folds)
           # if (!is.logical(object@GBM$keep.data)) {
           #   cat("\nGBM$keep.data must be a logical")
           #   test <- FALSE
           # }
           # if (!is.logical(object@GBM$verbose)) {
           #   cat("\nGBM$verbose must be a logical")
           #   test <- FALSE
           # }
           # # test <- .fun_testIfIn(test, "GBM$class.stratify.cv", object@GBM$class.stratify.cv, c('bernoulli', 'multinomial'))
           # test <- .fun_testIfIn(test, "GBM$perf.method", object@GBM$perf.method, 
           #                       c('OOB', 'test', 'cv'))
           # 
           # ## GAM ------------------------------------------------------------
           # test <- .fun_testIfIn(test, "GAM$algo", object@GAM$algo,
           #                       c('GAM_mgcv', 'GAM_gam', 'BAM_mgcv'))
           # test <- .fun_testIfIn(test, "GAM$type", object@GAM$type,
           #                       c('s_smoother', 's', 'lo', 'te'))
           # if (!is.null(object@GAM$k)) {
           #   if (!is.numeric(object@GAM$k)) {
           #     cat("\nGAM$k must be a integer")
           #     test <- FALSE
           #   } else {
           #     if (object@GAM$k < -1 | object@GAM$k %% 1 != 0) {
           #       cat("\nGAM$k must be > -1")
           #       test <- FALSE
           #     }
           #   }
           # }
           # test <- .fun_testIfPosInt(test, "GAM$interaction.level", object@GAM$interaction.level)
           # if (!is.null(object@GAM$myFormula) &&
           #     (!inherits(object@GAM$myFormula, "formula"))) {
           #   cat("\nGAM$myFormula must be NULL or a formula object")
           #   test <- FALSE
           # }
           # 
           # if (!inherits(object@GAM$family, "family")) {
           #   cat("\nGAM$family must be a valid family object")
           #   test <- FALSE
           # }
           # if (!is.list(object@GAM$control)) {
           #   cat("\nGAM$control must be a list like that returned by gam.control")
           #   test <- FALSE
           # }
           # test <- .fun_testIfIn(test, "GAM$method", object@GAM$method,
           #                       c("GCV.Cp", "GACV.Cp", "REML", "P-REML", "ML", "P-ML"))
           # if (any(!object@GAM$optimizer %in% 
           #         c("perf", "outer", "newton", "bfgs", "optim", "nlm", "nlm.fd"))) {
           #   cat("\nGAM$optimizer bad definition (see ?mgcv::gam)")
           #   test <- FALSE
           # }
           # if (!is.logical(object@GAM$select)) {
           #   cat("\nGAM$select must be a logical")
           #   test <- FALSE
           # }
           
         }
)


