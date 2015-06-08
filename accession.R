#' Returns all accession numbers based on a search of cell line. 
#' 
#' This function creates a list of all accession numbers specific to the cell_line. This allows the user to further search 
#' studies that used this cell line. 
#' 
#' 
#' @param cell_line_name: Cell line for which you wish to find the accession number. THE FUNCTION DOES NOT SUPPORT - OR SPACES
#' BETWEEN THE LETTERS. WRITE THE CELL LINE AS ONE WORD WITH NO SPACES.  
#'
#' @return All: A list, separated by type of accession number, of all the accession numbers.  
#'
#' @examples
#' accession('MCF7')

accession <- function(cell_line_name){
  
  Subset <- as.data.frame(subset(metadata,
                                 subset=(cell_line==cell_line_name)))
  {
  runaccession <- unique(Subset[,1])
  sampleaccession <- unique(Subset[,2])
  expaccession <- unique(Subset[,3])
  studyaccession <- unique(Subset[,4])
  subaccession <- unique(Subset[,5])
  }
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
