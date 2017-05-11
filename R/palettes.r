##############################################################################
# Setup: define color palettes
##############################################################################
library(RColorBrewer)

#' @export
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#e7baea")
#' @export
iwhPalette <- c("#cd5431", "#a14ad9", "#67b03f", "#604dad", "#c79931",
  "#cc4498", "#4d9f83", "#b54f50", "#5d8cb6", "#7f7d48", "#a67abe", "#965571")
#' @export
myPaletteCat <- colorRampPalette(brewer.pal(12, "Paired"))
#' @export
myPaletteCatN <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
#' @export
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
#' @export
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
#' @export
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
#' @export
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
#' @export
myPaletteO <- colorRampPalette(rev(brewer.pal(11, "RdBu")), space="Lab")
#' @export
myPaletteBrBG <- colorRampPalette(rev(brewer.pal(11, "BrBG")), space="Lab")

#' @export
rb <- c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
#' @export
g <- myPaletteG(6)[1:3]
#' @export
rbg <- c(myPaletteB(6)[1:3],myPaletteR(6)[1:3], myPaletteG(6)[1:3])

#' @export
orderedcats <- c("AT_CG", "AT_GC", "AT_TA",
  "GC_AT", "GC_CG", "GC_TA",
  "cpg_GC_AT", "cpg_GC_CG", "cpg_GC_TA")
#Define colors if using this ordering
# cols <- myPaletteCat(12)[c(8,10,12,2,4,6,1,3,5)]

#' @export
orderedcats1 <- c("AT_GC", "AT_CG", "AT_TA",
  "GC_AT", "GC_TA", "GC_CG",
  "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")
#' @export
orderedcats2 <- c("A>G", "A>C", "A>T",
  "C>T (non-CpG)", "C>A (non-CpG)", "C>G (non-CpG)",
  "CpG>TpG", "CpG>ApG", "CpG>GpG")
#' @export
gp_cols <- myPaletteCat(12)[
  c(10,8,12,
    2,4,6,
    1,3,5)] #<- colors if using this ordering

# cols returns following color list:
# "#6A3D9A" "#FF7F00" "#B15928" # A>N
# "#1F78B4" "#33A02C" "#E31A1C" # C>N (non-CpG)
# "#A6CEE3" "#B2DF8A" "#FB9A99" # CpG>NpG
