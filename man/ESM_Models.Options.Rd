\name{ESM_Models.Options}
\alias{ESM_Models.Options}

\title{
Model Options for ESMs
}
\description{
Generate a list of model parameters
}
\usage{
ESM_Models.Options(GLM = NULL, 
                   GBM = NULL)
}
\arguments{
  \item{GLM}{
\code{list}. with the objects:
\itemize{
\item{type}: \code{character}. Either "linear", "quadratic" or "polynomial". \emph{Default}: "quadratic".
\item{myFormula}: \code{formula}. Custom formula. The response variable should be called resp. Example: as.formula("resp~predictor") \emph{Default}: \code{NULL}.
\item{test}: \code{character}. Either "AIC" or "none". Perform or not a step AIC selection for the formula. \emph{Default}: "none".
}
}
  \item{GBM}{
\code{list}. with the objects (see \link[gbm]{gbm}):
\itemize{
\item{n.trees}: \code{integer}. Number of trees to fit. \emph{Default}: 1000.
\item{interaction.depth}: \code{integer}. The maximum depth of each tree \code{numeric}. Number of trees to fit. \emph{Default}: 4.
\item{n.minobsinnode}: \code{integer}. The minimum number of observations in the terminal nodes of the trees. \emph{Default}: 5.
\item{shrinkage}: \code{numeric}. The learning rate   \emph{Default}: 0.005.
\item{bag.fraction}: \code{numeric}. The fraction of the training set observations randomly selected to propose the next tree in the expansion. \emph{Default}: 0.5.
\item{train.fraction}: \code{numeric}. The first fraction of observations to fit the model. \emph{Default}: 1.
\item{cv.folds}: \code{integer}. Number of cross-validations to perform. \emph{Default}: 3.
\item{verbose}: \code{logical}. print or not the progress. \emph{Default}: FALSE.
\item{n.cores}: \code{logical}. Number of CPU cores to use. \emph{Default}: \code{NULL}.
}
}
}
\details{
For the arguments of each modeling technique, please refer to the manual of \code{\link[stats]{glm}} and \code{\link[gbm]{gbm}}.
}
\value{
a \code{list} of parameters for ESM
}

\author{
Flavien Collart from the code of biomod2::BIOMOD_ModelingOptions (version 4.2-4)
}
\examples{
models.options = ESM_Models.Options(GLM=list(test="AIC",
                                             type="polynomial"))
}
