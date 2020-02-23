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

perform_cap_scale <- function(idata, norm_dat, env_list){

  env_data <- idata[, env_list, with = FALSE]
  # set up full and null models for ordistep
  cap1 <- capscale(norm_dat ~ ., data = env_data, dist = "bray")
  cap0 <- capscale(norm_dat ~ 1, data = env_data, dist = "bray")

  #perform forward and backward selection of explanatory variables

  env <- ordistep(cap0, scope=formula(cap1))

  anova_res <- env$anova

  return(list(env, cap0, cap1, anova))

}
