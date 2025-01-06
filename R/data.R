### Dataset documentation

#' Species and environmental data to perform ESM with the ESM package
#' 
#' Contains sampling of 3 bryophytes in 413 plots in the Western Swiss Alps
#'
#' @format A data.frame with 413 plots and 10 variables:
#' \describe{
#' \item{x}{x coordinate in the CRS EPSG:2056}
#' \item{y}{y coordinate in the CRS EPSG:2056}
#' \item{Campylophyllum_halleri}{Presence/Absence data of \emph{Campylophyllum halleri}}
#' \item{Dicranella_subulata}{Presence/Absence data of \emph{Dicranella subulata}}
#' \item{Solenostoma_gracillimum}{Presence/Absence data of \emph{Solenostoma gracillimum}}
#' \item{ch_edaphic_eivdescombes_pixel_r}{}
#' \item{sradY}{Annual solar radiation}
#' \item{bio3_tiso_8110_LV95}{}
#' \item{bio15_ps_8110_LV95}{}
#' }
#' 
#' @author Flavien Collart
#' @references
#' Collart, F., Kiebacher, T., Quetsch, M., Broennimann, O., Guisan, A. & Vanderpoorten, A. 2024. To what extent can we predict 
#' variation of bryophyte and tracheophyte community composition at fine spatial scale along an elevation gradient?
#' \emph{STOTEN}, 926:171741. \doi{10.1016/j.scitotenv.2024.171741}.
#' @examples 
#' data(ESM_species.env)
"ESM_Species.Env"


#' Environmental SpatRaster to perform ESM with the ESM package
#' 
#' Contains stack of 4 predictors in the Western Swiss Alps
#'
#' @format A SpatRaster with 5 layers:
#' \describe{
#' \item{ch_edaphic_eivdescombes_pixel_r}{}
#' \item{sradY}{Annual solar radiation}
#' \item{bio3_tiso_8110_LV95}{}
#' \item{bio15_ps_8110_LV95}{}
#' }
#' @author Flavien Collart
#' @references
#' Collart, F., Kiebacher, T., Quetsch, M., Broennimann, O., Guisan, A. & Vanderpoorten, A. 2024. To what extent can we predict 
#' variation of bryophyte and tracheophyte community composition at fine spatial scale along an elevation gradient?
#' \emph{STOTEN}, 926:171741. \doi{10.1016/j.scitotenv.2024.171741}.
#' @examples 
#' library(terra)
#' data(ESM_Env)
#' ESM_Env <- terra::unwrap(ESM_Env)
#' 
"ESM_Env"