#' Extract spawning stock biomass (ssb) from gmacs run
#'
#' females
#'
#' @param M list object created by read_admb function
#' @return dataframe of observed and predicted survey
#' @export
#' 
.get_ssb_df <- function(M) {
    n <- length(M)
    mdf <- NULL
    for (i in 1:n)
    {
        A <- M[[i]]
        df <- data.frame(Model=names(M)[i], A$SSB)
       names(df) <- c("Model","Year","SSB","SE","lb","ub")
        mdf     <- rbind(mdf, df)
    }
    mdf$Model <- as.factor(df$Model)
    return(mdf)
}


#' Plot predicted spawning stock biomass (ssb)
#'
#'
#' @param M List object(s) created by read_admb function
#' @param xlab the x-label of the figure
#' @param ylab the y-label of the figure
#' @param ylim is the upper limit of the figure
#' @param alpha the opacity of the ribbon
#' @return Plot of model estimates of fit to survey 
#' @export
#' 
plot_ssb <- function(M, xlab = "Year", ylab = "Spawning biomass ", xlim=c(1980,2020), ylim = NULL, alpha = 0.1)
{
    xlab <- paste0("\n", xlab)
    ylab <- paste0(ylab, "\n")
    
    mdf <- .get_ssb_df(M)
    
    p <- ggplot(mdf) + labs(x = xlab, y = ylab) + xlim(xlim)
    
    if (is.null(ylim))
    {
        p <- p + expand_limits(y = 0)
    } else {
        p <- p + ylim(ylim[1], ylim[2])        
    }

    
    if (length(M) == 1)
    {
        p <- p + geom_line(aes(x = Year, y = SSB) )              + geom_ribbon(aes(x=Year, ymin=lb, ymax=ub), alpha = alpha)
    } else {
        p <- p + geom_line(aes(x = Year, y = SSB))  + geom_ribbon(aes(x=Year, ymin=lb, ymax=ub, fill=Model), alpha = alpha)
    }
    
    if(!.OVERLAY) 
        p <- p + facet_grid( .~Model) + guides(colour=FALSE)
        #p <- p + facet_wrap(~Model) + guides(colour=FALSE)
        #p <- p + ylim(0,3e6)        
    print(p + .THEME)
}
