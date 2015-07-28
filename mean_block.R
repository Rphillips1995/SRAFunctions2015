#' This function returns the mean coverage of genes for one or multiple BigWig files. 
#' 
#' The function first uses the extract block function to pull out the region of interest from 
#' the BigWig file. The mean of the region is then taken, and put into an empty vector. The
#' vector is then returned, so that the user may plot using whatever method preferred. 
#'
#' @param file : A single or vector of BigWig files. 
#' @param path :  The path to where the BigWig file is stored
#' @param chr : Chromosome for which the gene of interested is located. 
#' @param start : The numeric value for where on the chromosome the gene begins. 
#' @parm end : The numeric value for where on the chromosome the gene ends. 
#'
#' @return A plot in which the x axis represents the file, and the y axis represents the mean 
#' coverage of the gene. 
#'
#' 
#' @export
#' 
#' @examples
#' mean_block(file = 'NA20543_male_TSI_UNIGE_1-1-1.bw', 
#'            path = '/home/other/person/',
#'            chr = 6, 
#'            start = 30687978, 
#'            end = 30693203)



mean_block  <- function(file,path,chr,start,end) {
  extract.block <- function(files, chr, start, end, verbose = TRUE){
    rl <- IRanges::RangesList(IRanges::IRanges(start=start, end=end))
    names(rl) <- chr
    rles <- lapply(files, function(xx) {
      import(xx, as = "Rle", format = "bw", selection = BigWigSelection(rl))
    })
    megaMatrix <- do.call(cbind, lapply(rles, function(xx)      as.numeric(xx[[chr]][start:end])))
    megaMatrix
  }
  vec <- rep(NA,length(file))
  if(length(file)>=1)
  {
    for(i in 1:length(file))
    {
      print(paste0(path,file[i]))
      block <- extract.block(files = paste0(path, file[i]),
                             chr = chr,
                             start = start,
                             end = end,
                             verbose = TRUE)
      x <- mean(block)
      vec[i] <- x
    }
   return(vec) 
  }
}

    
