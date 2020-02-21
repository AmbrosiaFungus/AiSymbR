#' Perform a Capscale on the data
#'
#' @param idata
#' @param norm_dat
#' @param dist_mat
#'
#' @return list
#' @export
#'
#' @examples

perform_cap_scale <- function(idata, norm_dat, dist_mat){

  # set up full and null models for ordistep
  cap1 <- capscale(norm_dat ~ ., data=idata, dist=as.character(dist_mat))
  cap0 <- capscale(norm_dat ~ 1, data=idata, dist=as.character(dist_mat))

  #perform forward and backward selection of explanatory variables

  env <- ordistep(cap0, scope=formula(cap1))

  anova_res <- env$anova

  return(list(env, cap0, cap1, anova))

}
