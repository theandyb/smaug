##############################################################################
# Negbin model processing functions
##############################################################################

#' Run negbin model for each formula specified in formlist
#' @param formlist List of formula objects
#' @param data Data on which you wish to run the model(s)
#' @return List of model objects returned by glm.nb
#' @importFrom MASS glm.nb
#' @export
runMod <- function(formlist, data){
	out <- lapply(formlist, function(x) glm.nb(x, data))
	names(out) <- names(formlist)
	return(out)
}

#' Build list of fitted values (with CHR/BIN names) for each model
#' @param modlist list containing model objects
#' @param data Data to be fitted
#' @importFrom stats fitted.values
#' @return list of fitted values
#' @export
getFits <- function(modlist, data){
	out <- lapply(modlist,
		function(x){
			y <- fitted.values(x)
			names(y) <- paste0(data$CHR, ".", data$BIN)
			y
		})
	names(out) <- names(modlist)
	return(out)
}

#' Build list of dataframes for each model
#' @param fitlist list containing fitted values
#' @param data Data to be fitted
#' @return list of fitted values
#' @export
buildDF <- function(fitlist, data){
	out <- lapply(fitlist,
		function(x){
			data.frame(CHR, Category2=cat1, BIN,
				exp=x,
				obs=data$obs,
				# res=names(fitlist),
				stringsAsFactors=F)
		}
	)
	names(out) <- names(fitlist)
	return(out)
}

#' Compute standard error for correlations
#' @param corval correlation
#' @param ct count?
#' @return list of fitted values
#' @export
corSE <- function(corval, ct){
	sqrt((1-corval^2)/(ct-2))
}
