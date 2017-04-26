###############################################################################
# get signature loadings
#' @export
###############################################################################
get_loads <- function(widedat, nmfdat){
  sigloads <- data.frame(subtype=names(widedat),
    sig1=coef(nmfdat)[1,],
    sig2=coef(nmfdat)[2,],
    sig3=coef(nmfdat)[3,]) %>%
    mutate(sig1=sig1/sum(sig1), sig2=sig2/sum(sig2), sig3=sig3/sum(sig3)) %>%
    gather(subtype, value)

  names(sigloads) <- c("subtype", "sig", "value")
  sigloads <- sigloads %>%
    mutate(Category=substr(subtype, 1, 5),
      Sequence=substr(subtype, 7, 14))

  return(sigloads)
}

###############################################################################
# plot signature loadings
#' @export
###############################################################################
plot_loads <- function(sigloads){
  p <- ggplot(sigloads, aes(x=Sequence, y=value))+
    geom_bar(stat="identity")+
    facet_grid(sig~Category, scales="free_x")+
    # geom_label_repel(data=sigloads[sigloads$sig3>0.005,], aes(x=sig2, y=sig3, label=subtype))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90, hjust=1))
  return(p)
}

###############################################################################
# plot signature contribution across individuals
#' @export
###############################################################################
plot_ind_sigs <- function(sigdat){
  p <- ggplot(sigdat,
      aes(x=ID, y=prob, group=cluster, colour=cluster)) +
    geom_line()+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0), limits=c(0,1))+
    facet_wrap(~Study, scales="free_x", nrow=1)+
    ylab("signature contribution")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
      legend.position="bottom")

  return(p)
}
