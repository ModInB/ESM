% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ESM_Tools.R
\name{ESM_Generate.ODMAP}
\alias{ESM_Generate.ODMAP}
\title{Generates and fills ODMAP table}
\usage{
ESM_Generate.ODMAP(ESM.Mod = NULL, ESM.ensembleMod = NULL, ask.to.fill = TRUE)
}
\arguments{
\item{ESM.Mod}{The object returned by \code{\link{ESM_Modeling}}. \emph{Note that the object can also be NULL}}

\item{ESM.ensembleMod}{The object returned by \code{\link{ESM_Ensemble.Modeling}}. \emph{Note that the object can also be NULL}}

\item{ask.to.fill}{\code{logical}. If TRUE, teh functions will ask you 20 questions to fill some parts of the table.}
}
\value{
a \code{data.frame} containing the 4 columns and 84 lines of the ODMAP table.
some values can be filled depending on the arguments. The fifth column is an help
to fill the different lines.

\code{\link{ESM_Ensemble.Modeling}} but also with some questions that will be
asked.
}
\description{
The function generates an ODMAP table to report your Modelling procedure 
using the outputs of \code{\link{ESM_Modeling}} and asks you some
questions to help you fill this table.
}
\details{
ODMAP (Overview, Data, Model, Assessment, Prediction) is a standard protocol
to report your study on SDMs proposed by Zurell et al (2020).This function
helps you to fill a part of this ODMAP table. However, some parts cannot be 
filled automatically. We've added a column in this table to help you to fill
the other lines on your own.
}
\examples{
# A simple example where no values are filled in the table
 ODMAP_Table <- ESM_Generate.ODMAP(ESM.Mod = NULL,
                                   ESM.ensembleMod = NULL,
                                   ask.to.fill = FALSE)
# To see another example, see in ?ESM_Modeling
 
}
\references{
Zurell, D., Franklin, J., König, C., Bouchet, P.J., Dormann, C.F., Elith, J., 
Fandos, G., Feng, X., Guillera-Arroita, G., Guisan, A., Lahoz-Monfort, J.J., 
Leitão, P.J., Park, D.S., Peterson, A.T., Rapacciuolo, G., Schmatz, D.R., 
Schröder, B., Serra-Diaz, J.M., Thuiller, W., Yates, K.L., Zimmermann, N.E. 
and Merow, C. (2020), A standard protocol for reporting species distribution 
models. \emph{Ecography}, \bold{43}, 1261-1277. \doi{10.1111/ecog.04960}
}
\seealso{
\code{\link{ESM_Ensemble.Modeling}}, \code{\link{ESM_Modeling}}
}
\author{
Flavien Collart \email{flaviencollart@hotmail.com}
}
