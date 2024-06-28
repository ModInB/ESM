### Dataset documentation

#' Species and environmental data to perform ESM with the ESM package
#' 
#' Contains sampling of 3 bryophytes in 413 plots in the Western Swiss Alps
#'
#' @format A data.frame with 413 plots and 10 variables:
#' \describe{
#' \item{x}{X coordinate in the CRS EPSG:2056}
#' \item{y}{y coordinate in the CRS EPSG:2056}
#' \item{Tayloria_serrata}{Presence/Absence data of \emph{Tayloria serrata}}
#' \item{Thuidium_assimile}{Presence/Absence data of \emph{Thuidium assimile}}
#' \item{Tortella_pseudofragilis}{Presence/Absence data of \emph{Tortella pseudofragilis}}
#' \item{ddeg0}{Growing degree-days above 0C}
#' \item{mind68}{moisture index for month June to August}
#' \item{srad68}{solar radiation for month June to August}
#' \item{slope25}{average of slopes at 25m resolution}
#' \item{topos25}{average of topographic positions at 25m resolution}
#' }
#' 
#' @author Flavien Collart
#' @examples 
#' data(ESM_species.env)
"ESM_species.env"


#' Environmental SpatRester to perform ESM with the ESM package
#' 
#' Contains stacj of 5 predictors in the Western Swiss Alps
#'
#' @format A SpatRaster with 5 layers:
#' \describe{
#' \item{ddeg0}{Growing degree-days above 0C}
#' \item{mind68}{moisture index for month June to August}
#' \item{srad68}{solar radiation for month June to August}
#' \item{slope25}{average of slopes at 25m resolution}
#' \item{topos25}{average of topographic positions at 25m resolution}
#' }
#' @author Flavien Collart and Olivier Broennimann
#' @examples 
#' library(terra)
#' data(ESM_Env)
#' ESM_Env <- terra::unwrap(ESM_Env)
"ESM_Env"