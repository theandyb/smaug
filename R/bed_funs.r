#' BED to GRanges
#'
#' This function loads a BED-like file and stores it as a GRanges object.
#' The tab-delimited file must be ordered as 'chr', 'start', 'end', 'id', 'score', 'strand'.
#' The minimal BED file must have the 'chr', 'start', 'end' columns.
#' Any columns after the strand column are ignored.
#'
#' This is a forked version of the function written by Dave Tang
#' (https://github.com/davetang/bedr), with added support for bed files with headers
#'
#' @param file Location of your file
#' @param hd logical; indicates whether or not bed file has a header
#' @keywords BED GRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom utils read.table
#' @export
#' @examples
#' \dontrun{bed_to_granges('my_bed_file.bed')}
bed_to_granges <- function(file, hd=FALSE){
  df <- read.table(file,
                   header=hd,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  if(length(df)==3){
    gr <- with(df, GenomicRanges::GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GenomicRanges::GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GenomicRanges::GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GenomicRanges::GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}


#' Get recombination rate at each site
#' @param sites A dataframe with 10 columns (CHR, POS, Sequence, Dummy Variables for the 6 mutation types, and read depth)
#' @param rcrfile Path to file containing recombination rates
#' @importFrom bedr bedr.sort.region bedr.join.region
#' @importFrom dplyr select mutate arrange
#' @import magrittr
#' @return Vector containing recombination rates in order matching input dataset
#' @export
rcrCol <- function(sites, rcrfile){
  rcr <- read.table(rcrfile, header=T, stringsAsFactors=F)
  names(rcr) <- c("CHROM", "s", "e", "rate")
  rcr$CHROM <- as.character(rcr$CHROM)
  rcr <- bedr.sort.region(rcr,
                          check.chr=FALSE,
                          check.merge=FALSE,
                          check.zero.based=FALSE,
                          check.valid=FALSE,
                          verbose=FALSE)
  
  site_ranges <- sites %>%
    mutate(CHR=paste0("chr", CHR), START=POS-1, END=POS) %>%
    dplyr::select(CHR, START, END)
  
  site_ranges <- bedr.sort.region(site_ranges,
                                  check.chr=FALSE,
                                  check.merge=FALSE,
                                  check.zero.based=FALSE,
                                  check.valid=FALSE,
                                  verbose=FALSE)
  tmp <- bedr.join.region(site_ranges, rcr,
                          check.chr=FALSE,
                          check.merge=FALSE,
                          check.zero.based=FALSE,
                          check.sort=FALSE,
                          check.valid=FALSE,
                          verbose=FALSE) %>%
    dplyr::select(CHR, POS=END, rate) %>%
    mutate(CHR=as.numeric(gsub("chr", "", CHR)),
           rate=as.numeric(rate)) %>%
    arrange(CHR, POS)
  
  return(tmp$rate)
}

#' Get replication timing rate for each site
#' @param sites A dataframe with 10 columns (CHR, POS, Sequence, Dummy Variables for the 6 mutation types, and read depth)
#' @param repfile Path to file containing replication timing
#' @import magrittr
#' @return Vector with replication timing to be appended as column to input data.frame
#' @importFrom bedr bedr.sort.region bedr.join.region
#' @importFrom dplyr select mutate arrange filter
#' @importFrom utils read.table
#' @export
repCol <- function(sites, repfile){
  reptime <- read.table(repfile, header=F, stringsAsFactors=F, sep="\t")
  names(reptime) <- c("CHR", "END", "TIME")
  reptime <- reptime[!duplicated(reptime[,1:2]),] %>%
    filter(CHR<=22) %>%
    arrange(CHR, END) %>%
    mutate(CHR=paste0("chr", CHR), START=imputeStart(END)) %>%
    mutate(START=ifelse(START>END, 0, START)) %>%
    dplyr::select(CHROM=CHR, s=START, e=END, TIME)
  
  reptime <- bedr.sort.region(reptime,
                              check.chr=FALSE,
                              check.merge=FALSE,
                              check.zero.based=FALSE,
                              check.valid=FALSE,
                              verbose=FALSE)
  
  site_ranges <- sites %>%
    mutate(CHR=paste0("chr", CHR), START=POS-1, END=POS) %>%
    dplyr::select(CHR, START, END)
  
  site_ranges <- bedr.sort.region(site_ranges,
                                  check.chr=FALSE,
                                  check.merge=FALSE,
                                  check.zero.based=FALSE,
                                  check.valid=FALSE,
                                  verbose=FALSE)
  tmp <- bedr.join.region(site_ranges, reptime,
                          check.chr=FALSE,
                          check.merge=FALSE,
                          check.zero.based=FALSE,
                          check.sort=FALSE,
                          check.valid=FALSE,
                          verbose=FALSE) %>%
    dplyr::select(CHR, POS=END, TIME) %>%
    mutate(CHR=as.numeric(gsub("chr", "", CHR)),
           TIME=as.numeric(TIME)) %>%
    arrange(CHR, POS)
  
  return(tmp$TIME)
}

#' Check if site in bed file; returns 0 or 1
#' @param sites A dataframe with 10 columns (CHR, POS, Sequence, Dummy Variables for the 6 mutation types, and read depth)
#' @param bedfile Path to bedfile
#' @return Vector with binary indicators to be appended as column to input data.frame
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @export
binaryCol <- function(sites, bedfile){
  feat_ranges <- bed_to_granges(bedfile, hd=F)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$CHR),
                                        ranges=IRanges(start=sites$POS, end=sites$POS))
  out <- findOverlaps(site_ranges, feat_ranges, type="within", select="first")
  out[is.na(out)] <- 0
  return(as.integer(as.logical(out)))
}

#' Get GC content in 10kb window
#' @param sites A dataframe with 10 columns (CHR, POS, Sequence, Dummy Variables for the 6 mutation types, and read depth)
#' @param gcfile Path to file containing percent gc content
#' @return Vector with gc content to be appended as column to input data.frame
#' @importFrom utils read.table
#' @import magrittr
#' @importFrom dplyr select mutate arrange
#' @export
gcCol <- function(sites, gcfile){
  
  gcbins <- read.table(gcfile, header=F, stringsAsFactors=F)[,1:4]
  names(gcbins) <- c("CHR", "start", "end", "prop_GC")
  gcbins$CHR <- as.character(gcbins$CHR)
  site_tmp <- sites %>%
    dplyr::select(CHR, POS) %>%
    mutate(start=floor(POS/10000)*10000, end=ceiling(POS/10000)*10000)
  
  out1 <- merge(site_tmp, gcbins, by=c("CHR", "start")) %>%
    arrange(CHR, POS)
  return(out1$prop_GC)
}

#' Function determines start of interval from value in previous row
#' @param ends vector containing positions
#' @return Vector of starting positions
#' @export
imputeStart <- function(ends){
  starts<-c(0, ends[-(length(ends))])
  return(starts)
}