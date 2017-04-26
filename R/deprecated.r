##############################################################################
# Function checks if elements in a exist in b
# Output is binary vector of length same as b
##############################################################################
toBin <- function(a,b){
	is.element(b,a)
}


##############################################################################
# Get standard error estimate (DEPRECATED)
##############################################################################
std <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

##############################################################################
# Function to read file from disk (DEPRECATED)
##############################################################################
make.data <- function(filename, chunksize, skiprows,...){
	conn<-NULL
	function(reset=FALSE){
		if(reset){
			if(!is.null(conn)) close(conn)
			conn<<-file(filename,open="r")
		} else{
			rval<-read.table(conn, nrows=chunksize, skip=skiprows,...)
			if ((nrow(rval)==0)) {
				close(conn)
				conn<<-NULL
				rval<-NULL
			}
			return(rval)
		}
	}
}

##############################################################################
# QQ plot in ggplot2 with qqline (DEPRECATED)
##############################################################################
ggQQ <- function (vec) {
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]

  d <- data.frame(resids = vec)

  ggplot(d, aes(sample = resids)) +
	stat_qq() +
	geom_abline(slope = slope, intercept = int)

}

##############################################################################
# Calculate standard error of auc curve
##############################################################################
seauc <- function(auc, n){
  q1<-auc/(2-auc)
  q2<-2*auc^2/(1+auc)
  num<-auc*(1-auc)+(n-1)*(q1-auc^2)+(n-1)*(q2-auc^2)
  denom<-n^2
  sqrt(num/denom)
}


##############################################################################
# Append columns to windowed count data for all motif lengths (DEPRECATED)
##############################################################################
getSubMotifs <- function(data, nts, b3){

	# nts <- ifelse(grepl("^AT", cat1), "A", "C")
	outdat <- data

	# Loop currently just runs for 7->5bp motif aggregation;
	# can run over 7->5->3 by setting last index to :1
	for(j in ((nbp-1)/2-1):2){

		# Specify iteration motif length
		mlength <- (j+1)*2+1

		# Define rule for substring evaluation
		# griddef <- paste(c(rep("bases", j), "nts", rep("bases", j)), collapse=",")

		griddef <- paste(c("bases", "bases", "nts", "b3", "bases"), collapse=",")

		# Evaluate substring rule and get vector of submotifs
		tris <- apply(eval(parse(text=paste("expand.grid(",griddef,")"))),
			1, paste, collapse="")

		# Loop through each substring and append column of
		# aggregated counts
		for(k in tris){
			# Generate regex string; j is fixed per iteration
			# (e.g., looking for internal 3-mers or 5-mers)
			# so we search for all 3-mers or 5-mers by allowing
			# any base preceding or following the internal motif
			# regtri <- paste0("^", "[A-Z]{", j, "}", i, "[A-Z]{", j, "}")
			regtri <- paste0("^[A-Z]", k, "[A-Z]")

			# Extract sequences matching this submotif
			z <- names(data)[grepl(regtri, names(data))]

			# Ensure motif match vector includes only sequences
			# corresponding to the appropriate motif length
			z <- z[nchar(head(gsub("_[A-Z]*", "", z)))==mlength]

			# Create column and append to df
			tripct <- data %>%
				dplyr::mutate_(.dots=setNames(paste(z, collapse="+"), k)) %>%
				dplyr::select_(.dots=k)
			outdat <- cbind(outdat, tripct)
		}
	}

	return(outdat)
}

##############################################################################
#' commandArgs parsing (DEPRECATED)
#' return a named list of command line arguments
#'
#' Usage:
#' call the R script thus
#'   ./myfile.R --args myarg=something
#' or
#'   R CMD BATCH --args myarg=something myfile.R
#'
#' Then in R do
#'   myargs <- getArgs()
#' and myargs will be a named list
#' > str(myargs)
#' List of 2
#' $ file : chr "myfile.R"
#' $ myarg: chr "something"
#'
#' @title getArgs
#' @param verbose print verbage to screen
#' @param defaults a named list of defaults, optional
#' @return a named list
#' @author Chris Wallace
##############################################################################
getArgs <- function(verbose=FALSE, defaults=NULL) {
	myargs <- gsub("^--","",commandArgs(TRUE))
	setopts <- !grepl("=",myargs)
	if(any(setopts))
	myargs[setopts] <- paste(myargs[setopts],"=notset",sep="")
	myargs.list <- strsplit(myargs,"=")
	myargs <- lapply(myargs.list,"[[",2 )
	names(myargs) <- lapply(myargs.list, "[[", 1)

	## logicals
	if(any(setopts))
	myargs[setopts] <- TRUE

	## defaults
	if(!is.null(defaults)) {
		defs.needed <- setdiff(names(defaults), names(myargs))
		if(length(defs.needed)) {
		  myargs[ defs.needed ] <- defaults[ defs.needed ]
		}
	}

	## verbage
	if(verbose) {
		cat("read",length(myargs),"named args:\n")
		print(myargs)
	}
	myargs
}
