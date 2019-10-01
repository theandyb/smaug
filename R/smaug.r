#' Helper library for smaug pipeline
#' 
#' This package contains functions used in Carlson, et al (2018)
#' for the evaluation of germline mutation rate heterogeneity across the 
#' genome using extremely rare variants.
#' 
"_PACKAGE"

# Mute warnings about internal objects
#' @importFrom utils globalVariables
globalVariables(c("Motif", "Category2", "BIN",
                       "cat1", "CHR", "POS",
                       "sig1", "sig2", "sig3",
                       "subtype", "value", "nMotifs","END", "ERV_rel_rate",
                  "ID", "START", "Sequence", "TIME", "Type", "n",
                  "cluster", "coef", "f", "nbp", "prob", "prop_diff4", 
                  "rate", "sig", "v2","v2a", "v3", "v4", "xhi", 
                  "xlo", "yhi", "ylo"))