#'  Returnsthe sample origin that is unique/specific to that cell line. 
#' 
#' This function searches for the sample_origin based on a search of cell line. 
#'
#' @param cell_line_name: Cell line for which you wish to find the sample_origin. THE FUNCTION DOES NOT SUPPORT - OR SPACES
#' BETWEEN THE LETTERS. WRITE THE CELL LINE AS ONE WORD WITH NO SPACES.  
#'
#' @return origin:  A character of the sample origin that is unique/specific to that cell line.  
#'
#' @examples
#' predictsampleorigin('AGS')
#' predictsampleorigin('MCF7')

predictsampleorigin <- function(cell_line_name){
  
  Subset <- as.data.frame(subset(metadata,
                                 subset=(cell_line==cell_line_name)))
  {  
  origin <- unique(with(Subset,sample_origin))
  }
  {
  return(origin)
  }
  
}
