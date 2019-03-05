#' Function to get reverse complement
#' @param DNAstr String to get reverse complement of
#' @return Reverse compliment of the string
#' @examples
#' revcomp("ATCG")
#' @export
revcomp <- function(DNAstr) {
	step1 <- chartr("ACGT","TGCA",DNAstr)
	step2 <- unlist(strsplit(step1, split=""))
	step3 <- rev(step2)
	step4 <- paste(step3, collapse="")
	return(step4)
}

#' Function to reverse sequence
#' 
#' Used to correctly plot right flank in subsequence heatmaps
#' 
#' @param string Sequence to be reversed
#' @return Reversed string
#' @export
reverse_chars <- function(string){
	string_split = strsplit(as.character(string), split = "")
	reversed_split = string_split[[1]][nchar(string):1]
	paste(reversed_split, collapse="")
}
