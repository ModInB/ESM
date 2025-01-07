### Dataset documentation

#' Species and environmental data to perform ESM with the ESM package
#' 
#' Contains sampling of 3 bryophytes in 413 plots in the Western Swiss Alps with
#' environmental values from Collart et al (2024)
#'
#' @format A data.frame with 413 plots and 10 variables:
#' \describe{
#' \item{x}{x coordinate in the CRS EPSG:2056}
#' \item{y}{y coordinate in the CRS EPSG:2056}
#' \item{Campylophyllum_halleri}{Presence/Absence data of \emph{Campylophyllum halleri}}
#' \item{Dicranella_subulata}{Presence/Absence data of \emph{Dicranella subulata}}
#' \item{Solenostoma_gracillimum}{Presence/Absence data of \emph{Solenostoma gracillimum}}
#' \item{ch_edaphic_eivdescombes_pixel_r}{Ecological indicator values reflecting 
#' Soil pH.}
#' \item{sradY}{Annual solar radiation. Sum of monthly solar radiations}
#' \item{bio3_tiso_8110_LV95}{Isothermalithy. It is equal to: Mean Diurnal 
#' Range/Temperature Annual Range * 100. The values have been averaged between 
#' 1981 and 2010}
#' \item{bio15_ps_8110_LV95}{Precipitation Seasonality. Computed using the 
#' coefficient of variation among monthly precipitation and averaged between
#' 1981-2010.}
#' }
#' 
#' @author Flavien Collart
#' @references
#' Collart, F., Kiebacher, T., Quetsch, M., Broennimann, O., Guisan, A. & Vanderpoorten, A. 2024. To what extent can we predict 
#' variation of bryophyte and tracheophyte community composition at fine spatial scale along an elevation gradient?
#' \emph{STOTEN}, 926:171741. \doi{10.1016/j.scitotenv.2024.171741}.
#' 
#' Külling, N., Adde, A., Fopp, F. et al. 2024. SWECO25: a cross-thematic raster 
#' database for ecological research in Switzerland. \emph{Sci Data}, 11:21. 
#' \doi{10.1038/s41597-023-02899-1}.
#' 
#' @examples 
#' data(ESM_species.env)
"ESM_Species.Env"


#' Environmental SpatRaster to perform ESM with the ESM package
#' 
#' Contains stack of 4 predictors at 25m resolution in the Western Swiss Alps.
#'
#' @format A SpatRaster with 5 layers:
#' \describe{
#' \item{ch_edaphic_eivdescombes_pixel_r}{Ecological indicator values reflecting 
#' Soil pH.}
#' \item{sradY}{Annual solar radiation. Sum of monthly solar radiations}
#' \item{bio3_tiso_8110_LV95}{Isothermalithy. It is equal to: Mean Diurnal 
#' Range/Temperature Annual Range * 100. The values have been averaged between 
#' 1981 and 2010}
#' \item{bio15_ps_8110_LV95}{Precipitation Seasonality. Computed using the 
#' coefficient of variation among monthly precipitation and averaged between
#' 1981-2010.}
#' }
#' @author Flavien Collart
#' @references
#' Collart, F., Kiebacher, T., Quetsch, M., Broennimann, O., Guisan, A. & Vanderpoorten, A. 2024. To what extent can we predict 
#' variation of bryophyte and tracheophyte community composition at fine spatial scale along an elevation gradient?
#' \emph{STOTEN}, 926:171741. \doi{10.1016/j.scitotenv.2024.171741}.
#' 
#' Külling, N., Adde, A., Fopp, F. et al. 2024. SWECO25: a cross-thematic raster 
#' database for ecological research in Switzerland. \emph{Sci Data}, 11:21. 
#' \doi{10.1038/s41597-023-02899-1}.
#' 
#' @examples 
#' library(terra)
#' data(ESM_Env)
#' ESM_Env <- terra::unwrap(ESM_Env)
#' 
"ESM_Env"