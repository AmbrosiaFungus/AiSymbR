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

delete_neg_mock <- function(df, name_mock, name_neg){

  DT <- as.data.table(df)

  result <- DT[!eval(Sample_ID) %like% name_mock & !eval(Sample_ID) %like% name_neg]

  return(result)

}

#test <- delete_neg_mock(df=bacteria, name_mock = "Mock", name_neg = "Neg" )
