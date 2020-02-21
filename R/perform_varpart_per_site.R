#' Perform Variation Partioning per Site
#'
#' @param data merged table with Metadata and OTU table
#' @param site the site you want to perform the analysis on
#' @param numberOTU number where the OTU table starts
#' @param names_spatial names for the Longitude and Lattitude data
#' @param ... more arguments for perform varpart, it need more than 2
#'
#' @return
#' @export
#'
#' @examples

perform_varpart_per_site <- function(data, site, numberOTU, names_spatial, ...){

  DT <- setDT(copy(data))
  set_site <- DT[Location %like% site]
  spatial <- set_site[, names_spatial, with=FALSE]
  sp.pcnm <- as.data.frame(scores(pcnm(dist(spatial))))
  OTU <- set_site[1:nrow(set_site), numberOTU:ncol(set_site)]
  var_dat <- varpart(OTU, sp.pcnm,...,  data=set_site)

  return(var_dat)

}
