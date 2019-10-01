#' Function gets basic summary data from input directory
#' @title getData
#' @param summfile input file containing list of sites
#' @param singfile file containing sites annotated with sample ID
#' @param bindir directory containing the motif counts
#' @param maskfile Path to genome mask file (BED format)
#' @param binw Binwidth
#' @return a list of four components: sites, bins, mct, & aggseq
#' @importFrom readr read_tsv
#' @export
#' @examples
#' \dontrun{("sites.txt", "sites_id.txt", "/path/to/motif_counts")}
#' @author Jed Carlson
getData <- function(summfile, singfile, bindir, maskfile, binw){

  if(!file.exists(summfile)){
    msg <- paste0("Summary file ", summfile, " does not exist.\n")
  	stop(msg)
  }

  # Read in per-site summary data
  cat("Reading summary file:", summfile, "...\n")
  sites <- read_tsv(summfile)
  sites$Category2 <- sites$Category
  sites$Category <- gsub("cpg_", "", sites$Category)
  sites$BIN <- ceiling(sites$POS/binw)
  sites$MASK <- binaryCol(sites, maskfile)

  if(!missing(singfile)){
    cat("Annotating with sample ID...\n")
    inds <- read_tsv(singfile)
    names(inds) <- c("CHR", "POS", "S", "ALT", "ID")
    sites <- merge(sites, inds, by=c("CHR", "POS", "ALT"))
    rm(inds)
  } else {
    cat("No ID file specified! Downstream scripts may fail.\n")
  }

  # summarize motif counts genome-wide
  cat("Counting motifs...\n")
  # Read in motif counts per chromosome
  p1 <- "motifs_full.txt"
  bins <- get_bins(bindir, p1)
  mct <- get_mct(bins)

  out <- list()
  out$sites <- sites
  out$bins <- bins
  out$mct <- mct

  return(out)
}

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
#' \dontrun{get_bins("/path/to/motif_counts", "motifs.txt")}
#' @author Jed Carlson
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

# Function gets basic summary data from input directory
#' get_mct
#' Helper function for getData
#'
#' @title get_mct
#' @param bins bins object from getData
#' @return counts total number of motifs in bin directory
#' @import magrittr
#' @importFrom dplyr summarise group_by
#' @export
#' @examples
#' \dontrun{get_bins(bins)}
#' @author Jed Carlson
get_mct <- function(bins){
  out <- bins %>%
  	group_by(Motif) %>%
  	summarise(nMotifs=sum(nMotifs))
  return(out)
}

# Function gets basic summary data from input directory
#' get_aggseq
#' Helper function for getData
#'
#' @title get_aggseq
#' @param data directory containing the motif counts (assumed 1 per chromosome)
#' @param mct string or regular expression to match input files in bindir
#' @return calculates K-mer mutation rates from input data
#' @import magrittr
#' @importFrom dplyr summarise group_by
#' @export
#' @examples
#' \dontrun{get_aggseq(sites, mct)}
#' @author Jed Carlson
get_aggseq <- function(data, mct){
  out <- data %>%
  	group_by(Motif, Category2, BIN) %>%
  	summarise(n=n()) %>%
  	summarise(nERVs=sum(n))
  out <- merge(out, mct, by="Motif")
  out$rel_prop <- out$nERVs/out$nMotifs
  return(out)
}
