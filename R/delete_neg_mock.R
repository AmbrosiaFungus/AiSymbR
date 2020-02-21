#' Delete Mock-Community and Negatives
#'
#' @param df
#' @param names_mock
#' @param name_neg
#'
#' @return
#' @export
#'
#' @examples
delete_neg_mock <- function(df, names_mock, name_neg){

  df_clean <- as.data.table(df)

  result <- df_clean[!Sample_ID %like% get(names_mock) & !Sample_ID %like% name_neg]

  return(result)

}
