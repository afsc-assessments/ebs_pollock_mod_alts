#' Extract spawning stock biomass (ssb) from gmacs run
#'
#' Spawning biomass may be defined as all males or some combination of males and
#' females
#'
#' @param M list object created by read_admb function
#' @return dataframe of spawning biomass
#' @export
#' 
idx=1;which="ind";
.get_mnage_df <- function(M,idx,which="ind")
{
    n <- length(M)
    mdf <- NULL
    i=1
    for (i in 1:n)
    {
        A <- M[[i]]
        if (which=="ind") {
          type="Survey" 
          fleet <- get("Index_names",A)[idx]
        } else {
          type="Fsh"
          fleet <- get("Fshry_names",A)[idx]
        }
        df <- data.frame(get(paste0("EffN_",type,"_",idx),A))
        names(df) <- c("year","Eff_N1","Eff_N2","obs","pre","sda","lb","ub")
        head(df)
        #df <- data.frame(year = A$EffN_fsh)
        df$Model <- names(M)[i]
        df$type  <- type
        df$fleet <- fleet
        mdf      <- rbind(mdf, df)
    }
    return(mdf)
}


#' Plot predicted mean age by gear type
#'
#' @param M List object(s) created by read_admb function
#' @param xlab the x-label of the figure
#' @param ylab the y-label of the figure
#' @param ylim is the upper limit of the figure
#' @return Plot of model estimates of mean age against observed (and implied confidence bounds) 
#' @export
#' 
#xlab = "Year"; title=NULL; type="ind"; idx=1; ylab = "Mean age"; ylim = NULL 
plot_mnage <- function(M, xlab = "Year", title=NULL, type="ind", idx=1, ylab = "Mean age", ylim = NULL )
{
    xlab <- paste0("\n", xlab)
    ylab <- paste0(ylab, "\n")
    
    mdf <- .get_mnage_df(M,idx=idx,which=type)
    
    p <- ggplot(mdf) + labs(x = xlab, y = ylab)
    
    if (is.null(title))
      title=mdf$fleet[1]

    if (is.null(ylim))
    {
        p <- p + expand_limits(y = 0)
    } else {
        p <- p + ylim(ylim[1], ylim[2])        
    }
    
    if (length(M) == 1)
    {
        p <- p + geom_line(aes(x = year, y = pre)) +geom_point(aes(x=year, y=obs)) + 
            geom_errorbar(aes(x = year, ymax = ub, ymin = lb),size=.5)
    } else {
        p <- p + geom_line(aes(x = year, y = pre, col = Model),size=1.2) + geom_point(aes(x=year, y=obs)) + 
            geom_errorbar(aes(x = year, ymax = ub, ymin = lb),size=.5)
    }
    p <- p + facet_grid(type ~ .) #+ guides(colour=FALSE)
    return(p + .THEME + ggtitle(title))
}
#plot_mnage(M,type='fsh',idx=1)
