##############################################################################
# Function to install and load packages
#' usePackage
#' Convenience function for installing and loading R packages from CRAN
#'
#' @title usePackage
#' @param package package to be installed
#' @export
#' @examples
#' usePackage("tidyverse")
#' @author Jed Carlson
##############################################################################
usePackage <- function(package) {
  if (!is.element(package, installed.packages()[,1])){
    install.packages(package, dep = TRUE)
		require(package, character.only = TRUE)
	} else {
		require(package, character.only = TRUE)
	}
}
