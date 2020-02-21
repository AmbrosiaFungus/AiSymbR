#' Read data of Metadata and OTU
#'
#' @param path_df1 MetaData
#' @param path_df2 OTU table
#' @param number 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 

read_data <- function(path_df1, path_df2, number, ...){ 
  #df1 is your metadata and df2 is your distance matrix
  #number ist the number till which the distance matrix goes
  #in arguments I would like to put the parameters if it is ta
  df1 <- read.csv(file = as.character(path_df1), sep="\t")
  df2 <- read.csv(file = as.character(path_df2), sep=",", header = T)
  
  names(df1)[1] <- "SampleID" #name the first column as SampleID
  names(df2)[1] <- "SampleID" #same for the other df 
  
  tab <- df1[match(df2$SampleID, as.character(df1$SampleID)),]
  summary(df2$SampleID==tab$SampleID)
  
  #merge the two Dataframes
  
  new_df <- merge(tab,df2, by.x="SampleID", by.y = "SampleID")
  scaledOTU <- scale(new_df[1:nrow(new_df), as.numeric(number):ncol(new_df)], scale=T) #scales the matrix
  summary.new_df <- summary(scaledOTU, display=NULL)
  mean.new_df <- mean(scaledOTU)
  
  #return(list(new_df, scaledOTU, summary.new_df, mean.new_df))
  return(new_df)
  
  
  }

#test if its working
p1 <- "/home/robert/Projects/microbiome-a.incompertus/data/Metadata_with_Host.tsv"
p2 <- "/home/robert/Projects/microbiome-a.incompertus/data/OTU_table_97-fungi.csv"

t <- read_data(path_df1 = p1,
          path_df2 = p2,
          number = 23)



