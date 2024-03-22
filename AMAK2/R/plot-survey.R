#' Extract spawning stock biomass (ssb) from gmacs run
#'
#' Spawning biomass may be defined as all males or some combination of males and
#' females
#'
#' @param M list object created by read_admb function
#' @return dataframe of spawning biomass
#' @export
#' 
.get_srv_df <- function(M,idx=1)
{
  n <- length(M)
  mdf <- NULL
  i=1
  for (j in 1:length(idx))
  {

    for (i in 1:n)
    {
      A <- M[[i]]
      df <- data.frame(get(paste0("Obs_Survey_",idx),A) )
	  names(df) <- c("year","obs","pre","se","res1","res2")
      #1964 NA  5883.83 NA NA NA
      #df <- data.frame(year = A$yr_bts)
      df <- df %>% mutate(
            year = as.numeric(year),
            obs  = as.numeric(obs),
            pre  = as.numeric(pre),
            se   = as.numeric(se),
            res1 = as.numeric(res1),
            res2 = as.numeric(res2),
            year = as.numeric(year),
            lb   = as.numeric(ifelse(se=="NA",NA,obs/exp(2.*sqrt(log(1+se^2/obs^2))))),
            ub   = as.numeric(ifelse(se=="NA",NA,obs*exp(2.*sqrt(log(1+se^2/obs^2)))))
      )
      df$Model <- names(M)[i]
      mdf     <- rbind(mdf, df)
    }
  }
  return(mdf)
}


#' Plot predicted survey against observations
#'
#' Spawning biomass may be defined as all males or some combination of males and
#' females
#'
#' @param M List object(s) created by read_admb function
#' @param xlab the x-label of the figure
#' @param ylab the y-label of the figure
#' @param xlim is the year range
#' @param ylim is the upper limit of the figure
#' @return Plot of model estimates of spawning stock biomass 
#' @export
#' 
plot_survey <- function(M, which_survey=1, xlab = "Year", ylab=NULL, xlim= NULL, ylim = NULL, color="purple",overlay=FALSE)
{
    xlab <- paste0("\n", xlab)
    
    mdf <- .get_srv_df(M,idx=which_survey)
    
    p <- ggplot(mdf) + labs(x = xlab, y = ylab)
    
    if (is.null(ylab))
        ylab <-paste(M[[1]]$Index_names[which_survey],"Index")
    else
      ylab <- paste0(ylab, "\n")
    if (!is.null(xlim))
        p <- p + xlim(xlim[1], xlim[2])        
    if (is.null(ylim))
    {
        p <- p + expand_limits(y = 0)
    } else {
        p <- p + ylim(ylim[1], ylim[2])        
    }
    if (length(M) == 1)
    {
        p <- p + geom_line(aes(x = year, y = pre)) +geom_point(aes(x=year, y=obs),size=2,color="purple") + 
             geom_errorbar(aes(x = year, ymax = ub, ymin = lb),width=0.5) + ylab(ylab)
    } else {
        dw <- 0.5 
        p <- p + geom_line(aes(x = year, y = pre, col = Model),size=0.8,position=position_dodge(width=dw)) + geom_point(aes(x=year, y=obs,col=Model,fill=Model),position = position_dodge(width=dw)) + 
           geom_pointrange(aes(year, obs, ymax = ub, ymin = lb, color = Model,fill=Model), shape = 1, linetype = "solid", position = position_dodge(width = dw))
            #geom_errorbar(aes(x = year, ymax = ub, ymin = lb),width=0.5,position=position_dodge(width=0.9))
    }
    
    if(!overlay) 
      p <- p + guides(colour=FALSE)
    return(p + .THEME)
}
#.THEME=theme_few()
#plot_survey(M,which_survey=c(2),xlim=c(1990,2020))
#plot_survey(M,which_survey=c(5),xlim=c(1990,2020))
#plot_survey(M,which_survey=c(1,2))

