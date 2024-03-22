#' Extract selectivity from a model run
#'
#' @param M list object created by read_admb function
#' @return dataframe of spawning biomass
#' 
#which="ind";idx=1;styr=1980;
.get_sel <- function(M,which,idx,styr) {
    n <- length(M)
    mdf <- NULL
    i=1
    for (i in 1:n)
    {
      A   <- M[[i]]
      if (which=="ind") {
        type="Index" 
        fleet <- get("Index_names",A)[idx]
      } else {
        type="Fishery"
      fleet <- get("Fshry_names",A)[idx]
      }
	    sel <- get(paste0("sel_",which,"_",idx),A)[,-c(1,2)]
	    df  <- data.frame(year=A$Yr,sel=sel); 
	    lage <- length(sel[1,])
	    names(df) <- c("Year",1:lage)
	    df        <- gather(df,age,sel,2:(lage+1),-Year) %>% filter(Year>=styr) %>% mutate(age=as.numeric(age)) #+ arrange(age,yr)
      df$Model  <- names(M)[i]
      df$type   <- type
      df$fleet  <- fleet
      mdf       <- rbind(mdf, df)
    }
    return(mdf)
}
#' Plot selectivity over time                      
#'
#'
#' @param M List object(s) created by read_admb function
#' @param xlab the x-label of the figure
#' @param ylab the y-label of the figure
#' @param ylim is the upper limit of the figure
#' @param alpha the opacity of the ribbon
#' @return Plot of model estimates of spawning stock biomass 
#' @export
#' 
#styr=1964 ;styr=NULL; type="Fishery"; idx=1; alpha=0.2;scale=3.8;fill="purple"
plot_sel <- function(M,styr=NULL, type="Fishery", idx=1, alpha=0.2,scale=3.8,fill="purple")
{
	if (type=='Fishery') typshrt="fsh" else typshrt="ind"
	sdf <- .get_sel(M,which=typshrt,idx=idx,styr=styr)
	if (is.null(styr)) styr=min(sdf$Year)
	lage <- max(sdf$age)
  p1  <- ggplot(sdf,aes(x=age,y=as.factor(Year),height = sel)) + 
        geom_density_ridges(stat = "identity",scale = scale, alpha = alpha,fill=fill,color="black") + 
  	    .THEME +
         xlim(c(1,lage))+ ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$Year))))
	return(p1)
}