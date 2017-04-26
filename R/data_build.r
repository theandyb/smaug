##############################################################################
# Function gets basic summary data from input directory
#' getData
#' prepare input data for downstream processing
#'
#' @title getData
#' @param summfile input file containing list of sites
#' @param singfile file containing sites annotated with sample ID
#' @param bindir directory containing the motif counts
#' @return a list of four components: sites, bins, mct, & aggseq
#' @export
#' @examples
#' getdata("sites.txt", "sites_id.txt", "/path/to/motif_counts")
#' @author Jed Carlson
##############################################################################
getData <- function(summfile, singfile, bindir){

  if(!file.exists(summfile)){
    msg <- paste0("Summary file ", summfile, " does not exist.\n")
  	stop(msg)
  }

  # Read in per-site summary data
  cat("Reading summary file:", summfile, "...\n")
  sites <- read.table(summfile, header=T, stringsAsFactors=F)
  sites$Category2 <- sites$Category
  sites$Category <- gsub("cpg_", "", sites$Category)
  sites$BIN <- ceiling(sites$POS/binw)
  sites$MASK <- binaryCol(sites,
    paste0(parentdir, "/reference_data/testmask2.bed"))

  if(!missing(singfile)){
    cat("Annotating with sample ID...\n")
    inds <- read.table(singfile, header=T, stringsAsFactors=F)
    names(inds) <- c("CHR", "POS", "S", "ALT", "ID")
    sites <- merge(sites, inds, by=c("CHR", "POS", "ALT"))
  } else {
    cat("No ID file specified! Downstream scripts may fail.\n")
  }

  # summarize motif counts genome-wide
  cat("Counting motifs...\n")
  # Read in motif counts per chromosome
  p1 <- "motifs_full.txt"
  bins <- get_bins(bindir, p1)
  mct <- get_mct(bins)

  cat("Calculating K-mer rates...\n")
  aggseq <- get_aggseq(sites, mct)

  out <- list()
  out$sites <- sites
  out$bins <- bins
  out$mct <- mct
  out$aggseq <- aggseq
  return(out)
}

##############################################################################
# Function gets basic summary data from input directory
#' get_bins
#' Helper function for getData
#'
#' @title get_bins
#' @param bindir directory containing the motif counts (assumed 1 per chromosome)
#' @param pattern string or regular expression to match input files in bindir
#' @return concatenates files into a single data frame
#' @export
#' @examples
#' get_bins("/path/to/motif_counts", "motifs.txt")
#' @author Jed Carlson
##############################################################################
get_bins <- function(bindir, pattern){
  files <- list.files(path=bindir, pattern=pattern, full.names=T)
  if(length(files)!=22){
    msg <- paste0("Per-chromosome motif counts not initiated. \n")
    stop(msg)
  }

  cat("Reading bin files from:", bindir, "...\n")
  out <- do.call(`rbind`,
    lapply(files, read.table, header=T, stringsAsFactors=F))

  return(out)
}

##############################################################################
# Function gets basic summary data from input directory
#' get_mct
#' Helper function for getData
#'
#' @title get_mct
#' @param bins
#' @param pattern string or regular expression to match input files in bindir
#' @return counts total number of motifs in bin directory
#' @export
#' @examples
#' get_bins(bins)
#' @author Jed Carlson
##############################################################################
get_mct <- function(bins){
  out <- bins %>%
  	# dplyr::select(CHR, BIN, Sequence=MOTIF, nMotifs=COUNT) %>%
  	# dplyr::select(CHR, Motif, nMotifs=COUNT) %>%
  	group_by(Motif) %>%
  	summarise(nMotifs=sum(nMotifs))
  return(out)
}

##############################################################################
# Function gets basic summary data from input directory
#' get_aggseq
#' Helper function for getData
#'
#' @title get_aggseq
#' @param data directory containing the motif counts (assumed 1 per chromosome)
#' @param mct string or regular expression to match input files in bindir
#' @return calculates K-mer mutation rates from input data
#' @export
#' @examples
#' get_aggseq(sites, mct)
#' @author Jed Carlson
##############################################################################
get_aggseq <- function(data, mct){
  out <- data %>%
  	group_by(Motif, Category2, BIN) %>%
  	summarise(n=n()) %>%
  	summarise(nERVs=sum(n))
  out <- merge(out, mct, by="Motif")
  out$rel_prop <- out$nERVs/out$nMotifs
  return(out)
}
