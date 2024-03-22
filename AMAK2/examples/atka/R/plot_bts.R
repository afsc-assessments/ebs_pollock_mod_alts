#' Extract spawning stock biomass (ssb) from gmacs run
#'
#' females
#'
#' @param M list object created by read_admb function
#' @return dataframe of observed and predicted survey
#' @export
#' 
.get_bts_df <- function(M) {
    n <- length(M)
    mdf <- NULL
    for (i in 1:n)
    {
        A <- M[[i]]
        df <- data.frame(year = as.numeric(A$Obs_Survey_1[,1]))
        df$Model <- names(M)[i]
        df$pre  <- as.numeric(A$Obs_Survey_1[,3])
        df$ub  <- df$lb   <- df$obs  <- as.numeric(A$Obs_Survey_1[,2])
        for (j in 1:length(df$year))
        {
          if (!is.na(df$ub[j]))
          {
            #print(df$ub[j])
            #print(A$Obs_Survey_1[j,4] )
            df$lb[j]   <- df$obs[j] - 1.96*as.numeric(A$Obs_Survey_1[j,4] )
            df$ub[j]   <- df$obs[j] + 1.96*as.numeric(A$Obs_Survey_1[j,4] )
          }
        }
        mdf     <- rbind(mdf, df)
    }
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
plot_bts <- function(M, xlab = "Year", ylab = "Bottom trawl survey biomass index", xlim=c(1990,2018), ylim = NULL, alpha = 0.1)
{
    xlab <- paste0("\n", xlab)
    ylab <- paste0(ylab, "\n")
    
    mdf <- .get_bts_df(M[2])
    
    p <- ggplot(mdf) + labs(x = xlab, y = ylab) + xlim(xlim)
    
    if (is.null(ylim))
    {
        p <- p + expand_limits(y = 0)
    } else {
        p <- p + ylim(ylim[1], ylim[2])        
    }

    
    if (length(M) == 1)
    {
        mdf2 <- mdf[!is.na(mdf$obs),]
        p <- p + geom_line(aes(x = year, y = pre)) + geom_point(aes(x=year, y=obs))  #+
        geom_errorbar(data=mdf2,aes(x = year, ymax = ub, ymin = lb)) 
            #geom_ribbon(aes(x = year, ymax = ub, ymin = lb), alpha = alpha)
    } else {
        p <- p + geom_line(aes(x = year, y = pre, col = Model),size=1.2) + geom_point(aes(x=year, y=obs)) + 
            geom_errorbar(aes(x = year, ymax = ub, ymin = lb))
    }
    
    #if(!.OVERLAY) 
        #p <- p + facet_wrap(~Model) + guides(colour=FALSE)
        p <- p + facet_grid( Model~.) + guides(colour=FALSE)
        p <- p + ylim(0,3e6)        
    print(p + .THEME)
}
