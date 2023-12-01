\name{ESM_Variable.Contributions}
\alias{ESM_Variable.Contributions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variable contribution in ESMs
}
\description{
Calculates the variable contribution of each variable and method in an ESM model
}
\usage{
ESM_Variable.Contributions(ESM.Mod, 
                           ESM.ensembleMod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ESM.Mod}{
The object returned by \code{ESM_Modeling()}
}
  \item{ESM.ensembleMod}{
The object returned by \code{ESM_Ensemble.Modeling()}
}
}
\details{
Calculates the ratio between sum of weights of bivariate models where a focal variable was used and sum of weights of bivariate models where the focal variable was not used. The ratio is corrected for the number of models with or without the focal variable. This ratio gives an indication on the proportional contribution of the variable in the final ensemble model. A value of higher than 1 indicate that the focal variable has a higher contribution than average.
In the case of multiple methods (e.g., GLM, GAM...), the contributions are counted per method. For ensemble model, the contributions are then weighted means (based on the weighting score as chosen in ecospat.ESM.EnsembleModeling()) of single methods
    
}

\value{
Returns a dataframe with contribution values (i.e. proportional contribution) by variable and model

}
\author{
Olivier Broennimann \email{Olivier.Broennimann@unil.ch} with contributions of Heidi Mod \email{Heidi.mod1@gmail.com} and Daniel Scherrer \email{daniel.j.a.scherrer@gmail.com} and the adaptations of Flavien Collart \email{flaviencollart@hotmail.com}
}

\seealso{
\code{\link[ESM]{ESM_Modeling}}
}
