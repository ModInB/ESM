\name{ESM_Pooling.Evaluation}
\alias{ESM_Pooling.Evaluation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ensemble of Small Models: Evaluation via the Pooling procedure
}
\description{
This function evaluates the Ensemble of Small Models by pooling the different runs of the cross validation as in Collart et al (2021) and recommended for rare species by Collart & Guisan (2023).
}
\usage{
ESM_Pooling.Evaluation(ESM.Mod, 
                       ESM.ensembleMod, 
                       EachSmallModels = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ESM.Mod}{
The object returned by \code{ESM_Modeling()}
}
  \item{ESM.ensembleMod}{
The object returned by \code{ESM_Ensemble.Modeling()}
}
  \item{EachSmallModels}{
\code{logical}. Should the individual bivariate models be evaluated by the pooling procedure?
}
}
\details{
Because a minimum sample size is needed to evaluate models (see Collart & Guisan, 2023), this function uses the approach from Collart et al.(2021), which consists to pool the suitability values of the hold-out data (evaluation dataset) across replicates. As the same data point (presence or absence or background point) is presumably sampled in several replicates, the suitability values for each data point is consequently averaged across replicates where they were sampled. This procedure generates a series of independent suitability values with a size approximately equal (as some data points may not have been sampled by chance in any of the \emph{n} replicates) to that of the number of data point.
}
\value{
a \code{list} containing:
\itemize{
  \item{ESM.evaluations}{a \code{matrix} with the evaluation scores for the ESMs based on the different modelling algorithms and based on the consensus across the modelling algorithms (called here "ensemble")}
  \item{ESM.fit}{a \code{matrix} of predicted values resulting from the pooling procedure and used to compute the evaluation scores. The column \emph{resp} is where the species occurs or not}
  \item{ESM.evaluations.bivariate.models}{a \code{list} containing a matrix of evaluation scores for each bivariate models (generated only if EachSmallModels = T)}
  \item{ESM.fit.bivariate.models}{a \code{list} containing a matrix of of predicted values resulting from the pooling procedure for each bivariate models (generated only if EachSmallModels = T)}
  }
}
\references{
Collart, F., & Guisan, A. (2023). Small to train, small to test: Dealing with low sample size in model evaluation. \emph{Ecological Informatics}. \bold{75}, 102106. \doi{10.1016/j.ecoinf.2023.102106}

Collart, F., Hedenas, L., Broennimann, O., Guisan, A. and Vanderpoorten, A. 2021. Intraspecific differentiation: Implications for niche and distribution modelling. \emph{Journal of Biogeography}. \bold{48}, 415-426. \doi{10.1111/jbi.14009}
}
\author{
Flavien Collart \email{flaviencollart@hotmail.com} with contributions of Olivier Broennimann \email{olivier.broennimann@unil.ch}
}

\seealso{
\code{\link[ESM]{ESM_Ensemble.Modeling}}, \code{\link[ESM]{ESM_Modeling}}
}