#' This function builds a reference dataframe for cell_lines. 
#' 
#' This function uses the XML package to read the tables that exist on different URLS.
#' To build the reference dataframe of cell lines, we used tables from the Sigma
#' Aldrich table that contained over 900 different cell lines and their corresponding tissue
#' and species. The function saves the reference dataframe into the working directory. From the 
#' working directory you can load the reference dataframe into the global environment. 
#'
#'
#' @return ref_cell_lines A dataframe containing 913 different cell lines and corresponding tissues
#'
#' @keywords keywords
#' 
#' @examples
#' ref_cls()

ref_cls <- function(){
  #All urls will come from the sigma aldrich website. 
  #Getting a reference of cancer cell lines from sigma aldrich table. 
  url<-'http://www.sigmaaldrich.com/europe/life-science-offers/cell-cycle/sigma-ecacc-cell/cancer-cell-lines.html#BRC'
  cancer_cell_lines <- readHTMLTable(url, header = TRUE, 
                                     as.data.frame = TRUE, 
                                     which = 1, 
                                     stringsAsFactors = FALSE)
  
  #cancer_cell_lines clean up 
  cancer_cell_lines<-cancer_cell_lines[,-1] 
  cancer_cell_lines<-cancer_cell_lines[-c(1:18),] 
  cancer_cell_lines<-subset(cancer_cell_lines,subset=(V2!='NA')) 
  cancer_cell_lines<-within(cancer_cell_lines,V3<-toupper(with(cancer_cell_lines,V3)))
  cancer_cell_lines<-within(cancer_cell_lines,V4<-toupper(with(cancer_cell_lines,V4)))
  
  #There are rows where the tissue and Species were switched around. 
  cancer_cell_lines<-within(cancer_cell_lines,{
    V3<-ifelse(V3=='HUMAN',
               toupper(cancer_cell_lines$V4),
               cancer_cell_lines$V3)
  })
  
  
  
  #Getting a reference of cardio vascular cell lines 
  url2<-"http://www.sigmaaldrich.com/europe/life-science-offers/cell-cycle/sigma-ecacc-cell/cardiovascular-disease.html"
  cardio_cell_lines <- readHTMLTable(url2, header = TRUE, 
                                     as.data.frame = TRUE, 
                                     which = 1, 
                                     stringsAsFactors = FALSE)
  
  #Cleaning up cardio_cell_lines 
  cardio_cell_lines<- cardio_cell_lines[-c(1:9),]
  cardio_cell_lines<- cardio_cell_lines[,-1]
  cardio_cell_lines<- subset(cardio_cell_lines,subset=(V2!='NA'))
  
  #Getting a reference of diabetes and respiratory cell lines
  url3<-"http://www.sigmaaldrich.com/europe/life-science-offers/cell-cycle/sigma-ecacc-cell/diabetes-respiratory.html"
  diaresp_cell_lines <- readHTMLTable(url3, header = TRUE, 
                                      as.data.frame = TRUE, 
                                      which = 1, 
                                      stringsAsFactors = FALSE)
  
  #Cleaning up diaresp_cell_lines 
  diaresp_cell_lines <- diaresp_cell_lines[-c(1:9),]
  diaresp_cell_lines <- diaresp_cell_lines[,-1]
  diaresp_cell_lines <- subset(diaresp_cell_lines,subset=(V2!='NA'))
  #The names and order of the columns must be uniform throughout. 
  diaresp_cell_lines <- diaresp_cell_lines[,c(1,3,2)] 
  names(diaresp_cell_lines) [2] <- "V3"
  names(diaresp_cell_lines) [3] <- "V4"
  
  #Getting reference for Musculoskeletal cell lines 
  url4<-"http://www.sigmaaldrich.com/europe/life-science-offers/cell-cycle/sigma-ecacc-cell/musculoskeletal.html" 
  msstem_cell_lines <- readHTMLTable(url4, header = TRUE, 
                                     as.data.frame = TRUE, 
                                     which = 1, 
                                     stringsAsFactors = FALSE)
  
  #Cleaning up msstem_cell_lines 
  msstem_cell_lines <- msstem_cell_lines[-c(1:9),]
  msstem_cell_lines <- msstem_cell_lines[,-1]
  msstem_cell_lines <- subset(msstem_cell_lines,subset=(V2!="NA"))
  #The names and order of the columns must be uniform throughout. 
  msstem_cell_lines <- msstem_cell_lines[,c(1,3,2)]
  names(msstem_cell_lines) [2] <- "V3"
  names(msstem_cell_lines) [3] <- "V4"
  save(msstem_cell_lines,file='msmsstem_cell_lines.Rda')
  
  #Now we can create our final reference table 
  ref_cell_lines<-rbind(cancer_cell_lines,cardio_cell_lines,diaresp_cell_lines,msstem_cell_lines)
  #Finally we manipulate the cell lines just as we did for the metadata to make the cell lines more uniform and easier    
  #to search with. 
  ref_cell_lines <- within(ref_cell_lines,V2<-gsub('-','',with(ref_cell_lines,V2)))
  ref_cell_lines <- within(ref_cell_lines,V2<-gsub(' ','',with(ref_cell_lines,V2)))
  ref_cell_lines <- within(ref_cell_lines,V2<- toupper(with(ref_cell_lines,V2)))
  #There are duplicated cell_line in the table. We need to get rid of those. 
  refdups <- duplicated(ref_cell_lines[,1])
  ref_cell_lines <- ref_cell_lines[!refdups,]
  #Renaming 
  names(ref_cell_lines) [1] <- 'cell_line'
  names(ref_cell_lines) [2] <- 'tissue'
  names(ref_cell_lines) [3] <- 'species'
  { 
  #Saving the file
  save(ref_cell_lines,file='ref_cell_lines.Rda')
  }
}




