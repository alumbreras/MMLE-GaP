
#' @title  Plot four matrices
#' @description Plot four matrices in the same row using ggplot
#' @details This method is useful to compare, for instance, the different W
#' given by different choices of $K$ (number of latent dimensions) in NMF
#' @param a set of matrices
plot_multiple_W <- function(..., cols=4){
  matrices <- list(...)
  
  # Convert the matrix to dataframes and cretae the plots with ggplot
  plots <- list()
  for (m in 1:length(matrices)){
    df <- melt(matrices[[m]], varnames = c('i', 'j'))
    
    p <- ggplot(df, aes(x=as.factor(j), y=as.factor(i))) +
      geom_raster(aes(fill = value)) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      scale_x_discrete("i", expand=c(0,0), breaks=seq(1, 1000, by=ifelse(max(df$j)>10, 20, 1)))+
      scale_y_discrete("j", expand=c(0,0), limits = rev(levels(as.factor(df$i)))) +
      geom_hline(yintercept=0.5+1:1000, colour='grey', size=0.2) +
      geom_vline(xintercept=0.5+1:1000, colour='grey', size=0.2) +
      theme_bw()  +
      theme(aspect.ratio = 1,
            legend.position="bottom")
    
    plots[[m]] <- p
  }
  
  # Use this function to allow different color scales for each matrix
  multiplot(plots, cols=cols) 
}