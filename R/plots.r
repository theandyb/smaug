#' Function for plotting K-mer heatmaps
#' @param dat Data to be plotted
#' @param adj How many adjacent bases?
#' @return ggplot2 object with heatmap
#' @import ggplot2
#' @import magrittr
#' @importFrom dplyr mutate
#' @export
rrheat2 <- function(dat, adj){

  plotdat <- dat %>%
    mutate(v2=substr(Motif,1,adj),
      v2a=factor(as.character(lapply(as.vector(v2), reverse_chars))),
      v3=substr(Motif, adj+2, adj*2+1),
      v4=ERV_rel_rate,
      Category=Type,
      v5=factor(gsub("_", ">", Type)))

  nbox <- length(unique(plotdat$v2a))
  nint <- nbox/(4^(adj-1))
  xhi <- rep(1:(4^(adj-1)), 4^(adj-1))*nint+0.5
  xlo <- xhi-nint
  yhi <- rep(1:(4^(adj-1)),each=4^(adj-1))*nint+0.5
  ylo <- yhi-nint
  f <- data.frame(xlo,xhi,ylo,yhi)

  levs_a <- as.character(lapply(as.vector(levels(plotdat$v2a)),
		reverse_chars))

	p <- ggplot()+
	geom_tile(data=plotdat, aes(x=v3, y=v2a, fill=v4))+
	geom_rect(data=f, size=0.6, colour="grey70",
		aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
	scale_fill_gradientn("Relative Rate",
    colours=myPaletteO(11),
		trans="log10",
    breaks=0.0003*2^(0:9),
    labels=0.0003*2^(0:9),
		limits=c(0.0002, 0.2))+
	xlab("3' flank")+
	ylab("5' flank")+
  theme_classic()+
	theme(
		# legend.position="none",
		legend.title = element_text(size=18),
		legend.text = element_text(size=16),
		legend.key.height = unit(1.5, "cm"),
		axis.ticks.x = element_blank(),
		axis.ticks.y = element_blank(),
		axis.text.y = element_blank(),
		axis.text.x = element_blank(),
    axis.title.y = element_blank(),
		axis.title.x = element_blank())

	return(p)
}

#' Plot relative differences between ERVs and AV rates
#' @param dat Data to be plotted
#' @return ggplot2 object with heatmap
#' @import ggplot2
#' @export
rrheat3 <- function(dat){
	p<-ggplot()+
	geom_tile(data=dat, aes(x=v3, y=v2a, fill=prop_diff4))+
	geom_rect(data=f, size=0.6, colour="grey70",
		aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
	scale_fill_gradientn("Rp/Rs\n",
		colours=myPaletteBrBG(nbp),
		trans="log",
		breaks=c(0.5, 1, 2),
		labels=c("<0.5", "1", ">2"),
		limits=c(0.5, 2.2))+
		xlab("3' flank")+
		ylab("5' flank")+
	theme(
		legend.position="none",
		axis.ticks.x = element_blank(),
		axis.ticks.y = element_blank(),
		axis.text.y = element_blank(),
		axis.text.x = element_blank(),
    axis.title.y = element_blank(),
		axis.title.x = element_blank())

	return(p)
}

#' Display plot with inset
#' @param main ggplot2 object to be displayed
#' @param inset How many adjacent bases?
#' @param loc Where to place inset
#' @import ggplot2
#' @import magrittr
#' @importFrom dplyr mutate
#' @export
insetPlot <- function(main, inset, loc) {
	 print(main)
	 theme_set(theme_bw(base_size = 4))
	 print(inset, vp = loc)
	 theme_set(theme_bw())
 }

#' Multiple plot function
#'
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' @param ... ggplot objects passed in ..., or to plotlist (as list of ggplot objects)
#' @param plotlist List of plots to use (in lieu of above ggplot objects)
#' @param cols Number of columns in layout
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
#' @importFrom grid viewport pushViewport grid.newpage grid.layout
#' @export
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {

   # Make a list from the ... arguments and plotlist
   plots <- c(list(...), plotlist)

   numPlots = length(plots)

   # If layout is NULL, then use 'cols' to determine layout
   if (is.null(layout)) {
     # Make the panel
     # ncol: Number of columns of plots
     # nrow: Number of rows needed, calculated from # of cols
     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
   }

  if (numPlots==1) {
     print(plots[[1]])

   } else {
     # Set up the page
     grid.newpage()
     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

     # Make each plot, in the correct location
     for (i in 1:numPlots) {
       # Get the i,j matrix positions of the regions that contain this subplot
       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                       layout.pos.col = matchidx$col))
     }
   }
 }
