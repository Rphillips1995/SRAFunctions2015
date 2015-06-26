#' This function normalizes gene expression count data for BigWig files. 
#' 
#' The function first loads the extract block function needed to extract the region of interest from the bigwig
#' file. The function then ensures that the number of files is equal to the number of objects in the counts vector.
#' A for loop is then applied for each file that is inputted. The for loop first extracts the region of interest 
#' from the BigWig file. The data is normalized by dividing the count data by the number of total mapped reads, and 
#' then multiplying by 40000000. The mean is then taken to output the normalized value. The names of the files, and 
#' their corresponding normalized value is then outputted as a dataframe.
#'
#' @param file One or multiple BigWig files. MUST BE BIGWIG FILES.
#' @param chr Chromosome that the gene of interest is located on 
#' @param start The numeric value for where on the chromosome the gene begins
#' @param end The numeric value for where on the chromosome the gene ends 
#' @param counts This input can be a vector or a signle number. The numbers must be representative of the 
#' number of total mapped reads for that specific file. If a vector, the order of the numbers must be in the 
#' same order of the files inputted. 
#'
#' @return normalized_df A dataframe consisting of a column with the file, and corresponding normalized value. 
#'
#' @keywords keywords
#'
#' @export
#' 
#' @examples
#' normalized(file = "SRP040547_SRS582644_SRX501587_SRR1205725-1-1.bw",
#' chr = "chrX",
#' start = 73040486,
#' end = 73072588,
#' counts = 10973101)
#' 
#' normalized(file= vectoroffiles,
#' chr= "chrX",
#' start = 73040486,
#' end = 73072588,
#' counts= vectoroftotalmappedreads)



normalized <- function(file,chr,start,end,counts) {
  {
  extract.block <- function(files, chr, start, end, verbose = TRUE){
    rl <- IRanges::RangesList(IRanges::IRanges(start=start, end=end))
    names(rl) <- chr
    rles <- lapply(files, function(xx) {
      import(xx, as = "Rle", format = "bw", selection = BigWigSelection(rl))
    })
megaMatrix <- do.call(cbind, lapply(rles, function(xx)      as.numeric(xx[[chr]][start:end])))
megaMatrix
  }
}
if(class(file)=="character" & length(file)!=length(counts)){
  stop("The number of files must equal the length of the counts vector")
} 
normalized_values <- rep(NA,length(file))
if(length(file)>=1)
{
  for(i in 1:length(file))
  {
    total <- counts[i]
    block <- extract.block(files = file[i],
                           chr = chr,
                           start = start,
                           end = end,
                           verbose = TRUE)
    x <- mean(block/(total))* 40000000
    normalized_values[i] <- x
  }
  normalized_df  <- data.frame(file,normalized_values)
  return(normalized_df)
}

  }


