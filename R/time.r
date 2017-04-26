##############################################################################
# Function to install and load packages
#' timefun
#' sources script with runtime info
#'
#' @title timefun
#' @param script path to script you wish to source
#' @export
#' @author Jed Carlson
##############################################################################
timefun <- function(script){
  ptm <- proc.time()

  source(script)

  tottime <- (proc.time()-ptm)[3]
  cat("Done (", tottime, "s)\n")
}
