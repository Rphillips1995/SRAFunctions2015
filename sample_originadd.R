#' Decreases number of NAs in the sample_origin by inputing sample_origins for those 
#' cell_lines that do not have a sample_origin
#' 
#' This function works by first matching the cell_lines from the reference dataframe to the cell lines
#' in metadata. The function then merges the dataframes so that there is a column with the missing 
#' tissue information. Finally, using an ifelse statement, the function adds sample_origins for those 
#' that match and do not have a sample_origin. The dataframe is then saved into the users working 
#' directory, where they can directly import it into their global environment. 
#'
#'
#' @return metadata An updated metadata dataframe.
#'

#' @examples
#' sample_originadd()

sample_originadd <- function(){
  
  #Read in the metadata as a different object so that we can make changes to it and 
  #override later. 
  df<-metadata
  #Must be converted to character so that when we run the ifelse function it will input 
  #the sample_origin and not the number of the factor in the list
  df$sample_origin <- as.character(df$sample_origin)
  #Manipulating the cell_lines so that they are uniform and easier to match with. 
  df2<-within(df,cell_line <- gsub('-','',with(df,cell_line)))
  df2<-within(df2,cell_line <- gsub(' ','',with(df2,cell_line)))
  df2<-within(df2,cell_line <-toupper(with(df2,cell_line)))
  #Giving the rows an ID will help us in reordering the rows into their original format. 
  df2$id<-1:nrow(df2)
  
  #This creates a dataframe of the cell lines in the reference table that match the cell lines    in md2
  match <- subset(ref_cell_lines,subset=(ref_cell_lines$cell_line %in% df2$cell_line)==TRUE)
  names(match) [1] <- 'cell_line'
  
  #The tables must be merged so that we can pull out the tissue . 
  table <- merge(df2,match,all=TRUE)
  table <- within(table,cell_line<-gsub('-','',with(table,cell_line)))
  table <- within(table,cell_line<-gsub(' ','',with(table,cell_line)))
  table$match<-table$cell_line %in% ref_cell_lines$cell_line
  
  table<-within(table, {
    sample_origin<-ifelse(table$match=='FALSE',
                          table$sample_origin,
                          toupper(table$tissue))
  })
  {
  table <- table[order(table$id),]
  table <- table[,c(2:70,1,71:77)]
  metadata <- table
  }
  save(metadata , file='metadata.Rda')
}



