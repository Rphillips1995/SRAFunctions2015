#' Gets all accession numbers based on a search of cell line. 
#' 
#' This function first builds a dataframe that contains information on whether or not the cell lines matched one of the 
#' cell lines in the reference table. The function then extracts all accession numbers specific to the cell line searched. 
#' 
#' 
#' @param ell_line_name: Cell line for which you wish to find the sample_origin. THE FUNCTION DOES NOT SUPPORT - OR SPACES
#' BETWEEN THE LETTERS. WRITE THE CELL LINE AS ONE WORD WITH NO SPACES.  
#'
#' @return All: A list, separated by type of accession number, of all the accession numbers.  
#'
#' @examples
#' accession('MCF7')

accession <- function(cell_line_name){
  
  #Read in the md dataframe
  md <- read.table('all_illumina_sra_for_human.txt', 
                   sep = "\t",
                   header = TRUE, 
                   stringsAsFactors = FALSE, 
                   quote = "",
                   comment.char = "")
  
  #Clean up the data
  md2<-within(md,cell_line <- gsub('-','',with(md,cell_line)))
  md2<-within(md2,cell_line <- gsub(' ','',with(md2,cell_line)))
  md2<-within(md2,cell_line <-toupper(with(md2,cell_line)))
  
  #Getting a reference of cancer cell lines 
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
  
  #For some reason there are rows where V3 and V4 got switched
  
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
  msstem_cell_lines <- msstem_cell_lines[,c(1,3,2)]
  names(msstem_cell_lines) [2] <- "V3"
  names(msstem_cell_lines) [3] <- "V4"
  
  #Now we can create our final reference table 
  ref_cell_lines<-rbind(cancer_cell_lines,cardio_cell_lines,diaresp_cell_lines,msstem_cell_lines)
  
  #Finally we manipulate the cell lines just as we did for the metadata to make the cell lines more uniform and easier    to search with. 
  ref_cell_lines <- within(ref_cell_lines,V2<-gsub('-','',with(ref_cell_lines,V2)))
  ref_cell_lines <- within(ref_cell_lines,V2<-gsub(' ','',with(ref_cell_lines,V2)))
  ref_cell_lines <- within(ref_cell_lines,V2<- toupper(with(ref_cell_lines,V2)))
  
  #This creates a dataframe of the cell lines in the reference table that match the cell lines    in md2
  match <- subset(ref_cell_lines,subset=(ref_cell_lines$V2 %in% md2$cell_line)==TRUE)
  names(match) [1] <- 'cell_line'
  
  md.table <- merge(md2,match,all=TRUE)
  md.table <- within(md.table,cell_line<-gsub('-','',with(md.table,cell_line)))
  md.table <- within(md.table,cell_line<-gsub(' ','',with(md.table,cell_line)))
  md.table$match<-md.table$cell_line %in% ref_cell_lines$V2
  
  md.table<-within(md.table, {
    sample_origin<-ifelse(md.table$match=='FALSE',
                          md$sample_origin,
                          toupper(md.table$V3))
  })
  #View(md.table)
  {
  Subset <- as.data.frame(subset(md.table,
                                 subset=(cell_line==cell_line_name)))
  }
  runaccession <- unique(Subset[,2])
  sampleaccession <- unique(Subset[,3])
  expaccession <- unique(Subset[,4])
  studyaccession <- unique(Subset[,5])
  subaccession <- unique(Subset[,6])
  {
  All<-list(runaccession,sampleaccession,expaccession,studyaccession,subaccession)
  }
  {
  names(All) <- c('run_accession',
                  'sample_accession',
                  'experiment_accession',
                  'study_accession',
                  'submission_accession')
  }
  {
  return(All)
  }
  
}
