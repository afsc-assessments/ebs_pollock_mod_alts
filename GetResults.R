getwd()
library(r4ss)
library(tidyverse)
library(ggridges)
library(GOApollock)
library(ebswp)
library(here)
here()
#---Main pollock model run-----
here()
setwd(here("pm"))
mod_names <- c("Pollock_model")   
mod_dir <- c("m7")   
nextyr=2024
M<- get_results(rundir = "2023")
pm1<-(M[[1]])
#names(pm1)
pm_obj <- function(pm_run=pm1, yrmin=1964,yrmax=2024) {
  sel_pm <- as_tibble(pm_run$sel_fsh) |> rowwise() |> 
    mutate(across(everything(),~ ./max(c_across(everything())))) |>  ungroup()
  sel_pm <- cbind(1964:2023,sel_pm)
  names(sel_pm) = c("Year",1:15)
  sel_pm <- sel_pm |> filter(Year>=yrmin,Year<=yrmax) |> mutate(source="Pollock_model")
  ts <- rbind(
    data.frame(
      Year  = pm_run$SSB[,1],
      source= "pm",
      type  = "SSB",
      se    = pm_run$SSB[,3],
      value = pm_run$SSB[,2]),
  data.frame(
    Year =pm_run$R[,1],
    source= "pm",
    type  = "Recruits",
    se   =pm_run$R[,3],
    value=pm_run$R[,2])
)
return(list(sel=sel_pm,ts=ts))
}
pm_run<-pm_obj()
glimpse(pm_run)

#---SS run-----
# get FMSY
r1 <- SS_output(dir = here("ss","base")); 
SS_plots(r1)
#ss_run=r1; yrmin=1964;yrmax=2024; atype="Asel"; Fleet=NULL 
#atype=c("Asel2","Asel") 
r1$derived_quants

SS_obj <- function(ss_run=r1, yrmin=1964,yrmax=2023, atype=c("Asel")) {
  sel_ss <- ss_run$ageselex |> filter(Yr>=yrmin,Yr<=yrmax, Factor %in% atype) |> 
    mutate(fleet=ss_run$FleetNames[Fleet]) |> select(fleet,Year=Yr, c(9:23)) |> 
    mutate(source="SS") |> filter(fleet=="Fishery") |>  select(-fleet)

  #ts <- ss_run$sprseries |> filter(Era=="TIME",Yr>=yrmin,Yr<=yrmax ) |> 
  #  select(Year=Yr,source="SS",SSB,MnAge=MnAge_Smry,MnAgeCatch=MnAge_Catch,Recruits,Catch=Dead_Catch_B) 
tmp <- rbind(as_tibble(r1$derived_quants |> dplyr::filter(grepl("SSB_", Label))) |> 
  mutate(value=Value/2,se=StdDev,type="SSB",Year=as.numeric(substr(Label,5,9)),source="SS" ) |> 
  filter(!is.na(Year)) |> select(Year,source,type,se,value) ,
            as_tibble(r1$derived_quants |> dplyr::filter(grepl("Recr_", Label))) |> 
  return(list(sel=sel_ss,ts=tmp))
}
ss_run<-SS_obj()
ss_run$ts

#--Cole's in results from Cole's run----
#ssb <- filter(asdrep, name=='Espawnbio')
#unique(asdrep$name)
#ggplot(ssb, aes(year, est, ymin=lwr, ymax=upr)) + geom_line() +
#  geom_ribbon(alpha=.5) + labs(y='SSB (Mt)')
#rec <- filter(asdrep, name=='recruit')
#ggplot(rec, aes(year, est, ymin=lwr, ymax=upr)) + geom_pointrange() +
#  labs(y='Recruits (billions)')
#gp_sd=asdrep; gp_rep=arep; yrmin=1964;yrmax=2024 
gp_obj <- function(path=here("colepoll","2023"), gp_sd=asdrep, gp_rep=arep, yrmin=1964,yrmax=2024) {
  gp_rep <- read_pk_rep(path=path,model='gp', version='gp', styr=1964, endyr=2024)
  gp_sd  <- read_pk_std(path=path,model='gp', version='gp', styr=1964, endyr=2024)
  ssb <- filter(gp_sd, name=='Espawnbio')
  rec <- filter(gp_sd, name=='recruit')
  
  mysel <-  gp_rep$Fishery_selectivity
  mysel
  dim(mysel)
  new_columns <- matrix(rep(mysel[, 10], 5), nrow = nrow(mysel))
  
  # Combine the original matrix with the new columns
  mysel <- cbind(mysel, new_columns)
  
  sel_am <- as_tibble(cbind(1964:2023, as.data.frame(mysel)))
  names(sel_am) = c("Year",1:15)
  sel_am <- sel_am |> filter(Year>=yrmin,Year<=yrmax) |> mutate(source="GOA_model")
  ts <- rbind(
    data.frame(
      Year=ssb$year,
      source=ssb$version,
      type="SSB",
      se = ssb$se,
      value=ssb$est*500),
    data.frame(
      Year=rec$year,
      source=rec$version,
      type="Recruits",
      se = rec$se,
      value=rec$est*1000)
  )
  return(list(sel=sel_am,ts=ts))
}
gp_run<-gp_obj()

#--Read in results from SAM run----
load("SAM/poll23/baserun/model.RData")
SAM_obj <- function(myfit=fit, yrmin=1964,yrmax=2024) {
  ssb <-  data.frame(stockassessment::ssbtable(myfit))
  ssb$Year<- as.numeric(row.names(ssb))
  ssb <- ssb |> mutate(value=Estimate/2,type="SSB",se=(log(High)-log(Estimate)),source="SAM") 
   |>  select(Year,source,type,se,value)
  rec <-  data.frame(rectable(myfit))
  rec$Year<- as.numeric(row.names(rec))
  rec <- rec |> mutate(value=Estimate,type="Recruits",se=(log(High)-log(Estimate)),source="SAM") |>  select(Year,source,type,se,value)
  sel <-  data.frame(stockassessment::faytable(myfit))
  sel <- as.data.frame(t(apply(sel, 1, function(x) x / max(x))))
  sel <- cbind(as.numeric(row.names(sel)),sel)
  sel$source <- "SAM"
  names(sel)<-c("Year",1:15,"source")
  return(list(sel=sel,ts=rbind(rec,ssb)))
}
sam_run<- SAM_obj()

#--Read in results from AMAK2 run----
#setwd("~/_mymods/afsc-assessments/ebs_pollock_mod_alts/amak2/runs")
#jjmR::runJJM(models="amak",output="results",exec="../../src/amak")
#am1 <- jjmR::runit("amak",  pdf=TRUE,portrait=F,est=TRUE,exec="amak")
#am2 <- jjmR::readJJM(model="amak")
#names(am1[[1]][[1]][[5]][[1]])
AMAK_obj <- function(am_run=am2, yrmin=1964,yrmax=2024) {
  setwd(here("amak2","runs","base"))
  am_run <- PBSadmb::readRep("For_R_1"))
  setwd(here())
 # sel_am <- as_tibble(am_run[[1]]$ebswp$output$Stock_1$sel_fsh_1 )
  sel_am <- as_tibble(am_run$sel_fsh_1) 
  tmp <- (sel_am[,3:17]) |> rowwise() |> 
    mutate(across(everything(),~ ./max(c_across(everything())))) |>  ungroup()
  sel_am[,3:17]=tmp
  names(sel_am) = c("stock","Year",1:15)
  sel_am <- sel_am |> filter(Year>=yrmin,Year<=yrmax) |> select(-stock) |> 
    mutate(source="amak")
  sel_am
  ts <- rbind(
    data.frame(
      #Year  = am_run[[1]]$ebswp$output$Stock_1$SSB[,1],
      Year  = am_run$SSB[,1],
      source= "amak",
      type  = "SSB",
      #se    = am_run[[1]]$ebswp$output$Stock_1$SSB[,3]),
      #value = am_run[[1]]$ebswp$output$Stock_1$SSB[,2]),
      se    = am_run$SSB[,3],
      value = am_run$SSB[,2]),
  data.frame(
    Year =am_run$R[,1],
    #Year =am_run[[1]]$ebswp$output$Stock_1$R[,1],
    source= "amak",
    type  = "Recruits",
    #se   =am_run[[1]]$ebswp$output$Stock_1$R[,3])
    #value=am_run[[1]]$ebswp$output$Stock_1$R[,2])
    se   =am_run$R[,3],
    value=am_run$R[,2])
)
return(list(sel=sel_am,ts=ts))
}


#---Combine models---------
ss_run<-SS_obj()
gp_run<-gp_obj()
am_run<-AMAK_obj()
am_run$sel
sam_run$<-SAM_obj()
names(sam_run)
pm_run$sel
names(sam_run$sel)
all_sel <- rbind(sam_run$sel,pm_run$sel,am_run$sel,ss_run$sel,gp_run$sel)
dim(all_sel)
Plot_Sel_age()
Plot_Sel(lage=15)
all_ts <- rbind(sam_run$ts,pm_run$ts,am_run$ts,ss_run$ts,gp_run$ts)
Plot_SSB()
glimpse(pm_run$ts)
glimpse(am_run$ts)
glimpse(gp_run$ts)
glimpse(ss_run$ts)
Plot_SSB <- function(df=all_ts) {
  p1  <- df |> filter(Year<2024,Year>1963) |> ggplot(aes( x=Year,
        y= value,color=source)) + 
    geom_line(stat='identity') +
    geom_point(stat='identity') +
    ggthemes::theme_few() + ylab("SSB") + xlab("Year") +
    ylim(0,NA) +
    facet_grid(type~.,scales="free_y")  
  return(p1)
} 
Plot_SSB()
p1
sel=all_sel;fage=1;lage=15;nages=15
Plot_Sel_age <- function(sel=all_sel,fage=1,lage=15,nages=15) {
  sdf       <-  pivot_longer(sel,names_to="age",values_to="sel",cols=2:(nages+1)) %>% filter(Year>=1991,Year<2024) %>% mutate(age=as.numeric(age)) |> 
    filter(age>2, age<11)
  p1  <- sdf |> ggplot(aes( x=Year,
        y= sel,color=source)) + 
    geom_line(stat='identity') +
    geom_point(stat='identity') +
    ggthemes::theme_few() + ylab("Selectivity") + xlab("Year") +
    facet_wrap(.~age) ; 
  p1
}
Plot_Sel <- function(sel=all_sel,fage=1,lage=10,nages=15) {
  sdf <-  pivot_longer(all_sel,names_to="age",values_to="sel",cols=2:(nages+1)) |>  filter(Year>=1991,Year<2024) |>  mutate(age=as.numeric(age)) #|>   filter(age>2, age<11)
  p1  <- sdf |> ggplot(
    aes( x=age, y=as.factor(Year),height = sel)) + 
    geom_density_ridges(stat = "identity",scale = 2.8, 
                        alpha = .3, fill="goldenrod",color="tan") + 
    ggthemes::theme_few() + ylab("Year") + xlab("Age (years)") +
    scale_x_continuous(limits=c(fage,lage),breaks=fage:lage) +
    scale_y_discrete(limits=rev(levels(as.factor(sdf$Year)))) +
    facet_grid(.~source)
  return(p1)
}
Plot_Sel_age()

#--Other functions--------------
plot_sel <- function(Year=M$Yr,sel=M$sel_fsh, styr=1977, fage=NULL, lage=NULL, alpha=0.2,scale=3.8,fill="purple")
{
  df        <- data.frame(Year=Year,sel=sel );
  if (is.null(fage)) fage      <- 1
  if (is.null(lage)) lage      <- length(sel[1,])
  df <- df |> select(1:(lage-fage+2))
  names(df) <- c("Year",fage:lage)
  nages     <- length(fage:lage)
  names(df)
  sdf       <-  pivot_longer(df,names_to="age",values_to="sel",cols=2:(nages+1)) %>% filter(Year>=styr) %>% mutate(age=as.numeric(age)) #+ arrange(age,yr)
  p1  <- ggplot(sdf,aes(x=age,y=as.factor(Year),height = sel)) + geom_density_ridges(stat = "identity",scale = scale, alpha = alpha,
                                                                                     fill=fill,color="black") + ggthemes::theme_few() +
    ylab("Year") + xlab("Age (years)") +
    scale_x_continuous(limits=c(fage,lage),breaks=fage:lage) +
    scale_y_discrete(limits=rev(levels(as.factor(sdf$Year))))
  return(p1)
}