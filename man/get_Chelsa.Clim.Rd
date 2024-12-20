% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{get_Chelsa.Clim}
\alias{get_Chelsa.Clim}
\title{Download CHELSA data}
\usage{
get_Chelsa.Clim(
  var.names,
  time = "1981-2010",
  gcm = NULL,
  ssp = NULL,
  path = getwd(),
  extent = NULL,
  mask = FALSE,
  compress = TRUE,
  timeout = 300L
)
}
\arguments{
\item{var.names}{\code{character}. Vector containing either 'all' to download all the variables (works only at present time), 'all_fut'
for all variables that are also available for future time periods, 'bioclim' to download the 19 bioclimatic variables,
or the variable names to download. See details for more information.}

\item{time}{\code{character}. The time period(s) wanted. Must be: '1981-2010','2011-2040','2041-2070',and/or '2071-2100'.
\emph{Default: '1981-2010'}}

\item{gcm}{\code{character}. A vector containing the Global Circulation Models for which you want to download the datasets. 
Only needed when future time periods are provided. Must be: 'gfdl-esm4','ipsl-cm6a-lr','mpi-esm1-2-hr','mri-esm2-0' and/or,'ukesm1-0-ll'.
\emph{Default: NULL}.}

\item{ssp}{\code{character}. A vector containing the Shared Socio-economic Pathways for which you want to download the datasets. 
Only needed when future time periods are provided. Must be: '126', '370' and/or '585'.\emph{Default: NULL}.}

\item{path}{\code{character}.The path to store the climatic grids.}

\item{extent}{The extent to which you want to crop your climatic datasets. Extent must be either a 'SpatRaster', 'SpatVector', or'SpatExtent'.}

\item{mask}{\code{logical}. Do you want to mask your maps? \emph{Note that when TRUE, extent must be a 'SpatRaster' or a 'SpatVector'.}}

\item{compress}{\code{logical}. When extent is provided, do you want to compress the climatic grids? If TRUE, the compression
will be: "COMPRESS=DEFLATE", "PREDICTOR=2" and "ZLEVEL=6". \emph{Default: TRUE}.}

\item{timeout}{\code{integer}. The maximum downloading time allowed for each map in seconds. \emph{Default: 300L}.}
}
\value{
\code{character}. A vector containing the path to the downloaded files.
}
\description{
Download and crop climatic variables at different time periods from CHELSA v.2.1.
}
\details{
\describe{
get_chesla downloads and crops climatic variables at different time periods from CHELSA v.2.1.

The following list contains the different variables associated with their time period : 
\itemize{

\item{Present}: pet_perman_max, pet_perman_mean, pet_perman_min, pet_perman_range,
sfcWind_max, sfcWind_mean, sfcWind_min","sfcWind_range, vpd_max, vpd_mean, vpd_min, vpd_range,
hurs_max, hurs_mean, hurs_min, hurs_range, clt_max, clt_mean, clt_min, clt_range, cmi_max, cmi_mean,
cmi_min, cmi_range, rsds_1981-2010_mean, rsds_1981-2010_min, rsds_1981-2010_max, rsds_1981-2010_range.

\item{Present and Future}: the 19 bioclim variables (e.g. "bio1"), ai, fcf, fgd, gdd0, gdd10,
gddd5, gddlgd0, gddlgd5, gddlgd10, gdgfgd0, gdgfgd10, gdgfgd5, gsl, gsp, gst,
kg from 0 to 5 (e.g. "kg0"), lgd, ngd0, ngd10, ngd5, npp, scd, swd, swe.
}
The description of the different variables is available at the point 7 (p.11) of the following link :
\href{https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf}{CHELSA_tech_specification_V2.pdf}

The parameter "var.names" allows you to select the variables taken into account :
\itemize{

\item{all}: All available variables are selected. \emph{Note: Only works for present time}
\item{all_fut}: All the variables available in present AND future time are selected.
\item{bioclim}: the 19 bioclimatic variables are selected.
}
The variables available in the future are to be combined with a future scenario by choosing a global circular model (gcm) and a shared-socioeconomic pathway (ssp) :
\itemize{

\item{ssp126}: ssp1 (Sustainability - Taking the green road) combined with a radiative forcing of 2.6 W/m².
\item{ssp370}: ssp3 (Regional rivalry - A rocky road) combined with a radiative forcing of 7 W/m².
\item{ssp585}: ssp5 (Fossil-fueled development - Taking the highway) with a radiative forcing of 8.5 W/m². 
}
The complete description of the ssp is available in O'Neill et al. 2017. These ssp are combined with
Global circulation Models (GCM) to assess future climate.

If the files are to heavy and/or your connection is to slow, you may want to modified the timeout parameter to a greater period of time.
Note that you always have to end your number with "L" to make sure your time is an integer (e.g. 300L, 500L).

To crop a map to a different size than the original one, it is possible to add an extent and a mask parameter. If the extent of 
the map is large, you may want to compress it, reducing storage space without significantly increase the reading time in R.
}
}
\examples{
 \donttest{
# get_Chelsa.Clim("bio1")
}
}
\references{
Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, P., Kessler, M. (2017). Climatologies at high resolution for the Earth land surface areas. \emph{Scientific Data}. \bold{4}, 170122. \doi{10.1038/sdata.2017.122}.

Karger D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E, Linder, H.P., Kessler, M. (2018): Data from: Climatologies at high resolution for the earth’s land surface areas. \emph{EnviDat}. \doi{10.16904/envidat.228.v2.1}.

B. C. O’Neill, E. Kriegler, K. L. Ebi, E. Kemp-Benedict, K. Riahi, D. S. Rothman, B. J. van Ruijven, D. P. van Vuuren, J. Birkmann, K. Kok, M. Levy, W. Solecki. (2017). The roads ahead: Narratives for shared socioeconomic pathways describing world futures in the 21st century.
\emph{Global Environmental Change}. \bold{42} 169-180. \doi{10.1016/j.gloenvcha.2015.01.004}
}
\author{
Flavien Collart \email{flaviencollart@hotmail.com} with contributions of Adèle Hotermans
}
