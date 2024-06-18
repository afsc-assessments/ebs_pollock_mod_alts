#length/age model process

#Generic growth data to create length at age 
max_age <- 15
Linf <- 100
K <- .3
a0 <- 0
cv <- .2
alpha <- 0.001
beta <- 3
rzero <- 10000
steep <- 1

#Set up length/age composition data ranges
comp_ages <- 0:max_age
comp_lengths <- 0:ceiling((1+3*cv)*Linf)

#Set up number of growth morphes to track
n_growth_morphs <- 7
length_props <- 1+(-3:3)*cv
rec_props <- dnorm(length_props,1,cv)
rec_props <- rec_props/sum(rec_props)

#Vector to specify the ages tracked in the model
#For simplicity start with assumption that ages in 
#comps and model are the same
ages <- comp_ages
n_ages <- length(ages)
#Function to get expected length at age from growth params
AtoL <- function(a,Linf,K){
  L <- Linf*(1-exp(-K*a))+0.001
}

#Function to get expected weight at length from growth params
LtoW <- function(L,alpha,beta){
  W <- alpha*L^beta
}

logistic_mature <- function(x, inflection_point, slope){
  prop <- 1 / (1 + exp(-1 * slope * (x - inflection_point)))
}

get_recruit <- function(ssb,rzero,steep,phi_0){
  recruits <- (0.8*rzero*steep*ssb)/(0.2*phi_0*rzero*(1.0-steep)+ssb*(steep-0.2))
}
#Vector of model years 
start_year <- 1
end_year <- 10
years <- start_year:end_year
n_years <- length(years)

#Set up indexing vectors for all tracking across years, ages, and lengths

#For simplicity start with the assumption that these are not time varying 
#Assume that an equal number of lengths are tracked
#for each age. This simplifies setting things up but I 
#don't think will be stricktly necessary

year <- sort(rep(years,n_ages*n_growth_morphs))

#Assume that an equal number of lengths are tracked
#for each age. This simplifies setting things up but I 
#don't think will be strickly necessary
ages_x_lengths <- sort(rep(ages,n_growth_morphs))

age <- rep(ages_x_lengths,n_years)

#For this example I'm setting the lengths tracked for each 
#age based on a VB growth curve and 3 standard deviations
#either side of the mean in 1 sd steps. This should mean that
#almost all individuals (99.7%) are inside these bounds 
#assuming a normal distribution.
lengths_x_ages <- AtoL(ages_x_lengths,
                       rep(Linf,length(ages_x_lengths)),
                       rep(K,length(ages_x_lengths))
                       )*rep(length_props,length(ages))

length <- rep(lengths_x_ages,n_years)

#For this example I'm setting the weights as fixed weight at length
weight_x_length <- LtoW(length,
                       rep(alpha,length(length)),
                       rep(beta,length(length)))

weight <- weight_x_length

#For this example set a single natural mortality for all years, ages, and lengths
nat_mort <- rep(0.2,n_years*n_ages*n_growth_morphs)

#For this example set constant 50% female
prop_female <- rep(0.5,n_years*n_ages*n_growth_morphs)

#For this example set proportion mature to logistic maturity by age
age_50_mat <- 6
slope <- 2
mature_x_age<-logistic_mature(ages,age_50_mat,slope)
prop_Mature <- rep(sort(rep(mature_x_age,n_growth_morphs)),n_years)

#Data frame to hold abundance data specific to year, age, and length this could 
#probably just be a vector of abundance if the indexes are used.
# Abundance <- data.frame(abundance=rep(0,n_years*n_ages*n_growth_morphs),
#                         year=year,
#                         age=age,
#                         length=length,
#                         weight=weight)

abundance <- rep(0,n_years*n_ages*n_growth_morphs)
unfished_abundance <- rep(0,n_years*n_ages*n_growth_morphs)

init_NAAL <- rep(0,n_years*n_ages*n_growth_morphs)
init_unfished_NAAL <- rep(0,n_years*n_ages*n_growth_morphs)
#Initialize ssb
ssb <- rep(NA,n_years)
unfished_ssb <- rep(NA,n_years)

#setup transition pairs
growth_source <- NULL
growth_sink <- NULL
for(i in seq_along(years)){
  for(j in rev(seq_along(ages))){
    for(k in 1:n_growth_morphs){
      if(i==1){
        #In first year no transfer so just leave in same year 
        #an initial population function with fill these.
        growth_source<-c(growth_source,which(year==years[i] & 
                                             age==ages[j])[k])
        
        growth_sink<-c(growth_sink,which(year==years[i] & 
                                         age==ages[j])[k])
      }else if(j==1){
        #The first age class will be filled using a recruitment function.
        growth_source<-c(growth_source,which(year==years[i] & 
                                               age==ages[j])[k])
        
        growth_sink<-c(growth_sink,which(year==years[i] & 
                                           age==ages[j])[k])
      }else{
        #For all other years this will define the proportional transfer
        #between age/length bins each year. Currently simplified to 100% 
        #along a range of growth morphs. Abundance will be adjusted by 
        #natural and fishing mortality impacts.
        growth_source<-c(growth_source,which(year==years[i-1] & 
                                               age==ages[j-1])[k])
        
        growth_sink<-c(growth_sink,which(year==years[i] & 
                                           age==ages[j])[k])
      }
    }
  }
}

#Specify the transfer of abundance between year/age/length bins
transition_proportion <- rep(1,length(growth_source))   

#Setup fleet and survey structure
n_fleets <- 2

if(n_fleets>0){
  fleets <- 1:n_fleets
}else{
  fleets <- NULL
}

n_surveys <- 1

if(n_surveys>0){
  surveys <- 1:n_surveys
}else{
  surveys <- NULL
}

fleet <- sort(rep(fleets,n_years*n_ages*n_growth_morphs))
fleet_year <- rep(year,n_fleets)
fleet_age <- rep(age,n_fleets)
fleet_length <- rep(length,n_fleets)

fish_mort <- rep(0,n_years*n_fleets*n_ages*n_growth_morphs)

total_mort <- rep(0,n_years*n_ages*n_growth_morphs)

catch_abun <- rep(0,n_years*n_fleets*n_ages*n_growth_morphs)

survey <- sort(rep(surveys,n_years*n_ages*n_growth_morphs))
survey_year <- rep(year,n_surveys)
survey_age <- rep(age,n_surveys)
survey_length <- rep(length,n_surveys)

survey_q <- rep(0,n_years*n_surveys*n_ages*n_growth_morphs)

survey_abun <- rep(0,n_years*n_surveys*n_ages*n_growth_morphs)



for(i in seq_along(Growth_transition[,1]))
{
  sink <- growth_sink[i]
  source <- growth_source[i]
  trans_prop <- transition_proportion[i]
  catch_sink <- which(catch_year==year[sink] &
                      catch_age==age[sink] &
                      catch_length==length[sink])
  catch_source <- which(catch_year==year[source] &
                        catch_age==age[source] &
                        catch_length==length[source])
  survey_sink <- which(survey_year==year[sink] &
                       survey_age==age[sink] &
                       survey_length==length[sink])
  survey_source <- which(survey_year==year[source] &
                         survey_age==age[source] &
                         survey_length==length[source])
  
  if(year[sink]==start_year){
    abundance[sink] <- init_NAAL[which(init_age==age[sink] & init_length==length[sink])]
    unfished_abundance[sink] <- init_unfished_NAAL[which(init_age==age[sink] & init_length==length[sink])]
  }else if(age[sink]==0){
    if(is.na(ssb[which(years==year[sink])])){
      ssb[which(years==year[sink])] <- sum(abundance[which(year==year[sink])]*
                                           weight[which(year==year[sink])]*
                                           prop_mature[which(year==year[sink])]*
                                           prop_female[which(year==year[sink])])
      
      unfished_ssb[which(years==year[sink])] <- sum(unfished_abundance[which(year==year[sink])]*
                                                    weight[which(year==year[sink])]*
                                                    prop_mature[which(year==year[sink])]*
                                                    prop_female[which(year==year[sink])])
    }
    abundance[sink] <- get_recruit[ssb[which(years==year[sink])],rzero,steep,unfished_ssb[1]/rzero]
    unfished_abundance[sink] <- get_recruit[unfished_ssb[which(years==year[sink])],rzero,steep,unfished_ssb[1]/rzero]
  }else if(age[sink]==max_age){
    
    total_mort[sink] <- sum(nat_mort[sink],fish_mort[catch_sink])
    
    catch[catch_sink]<-catch[catch_sink] + (fish_mort[catch_sink]/total_mort[sink])*abundance[sink]*(1-exp(-total_mort[sink]))
    
    survey_intercept[survey_sink]<-survey_intercept[survey_sink] + abundance[sink]*(1-exp(-survey_obs_frac[survey_sink])) 
    
    abundance[sink] <- abundance[sink]*exp(-total_mort[sink])
    
    total_mort[source] <- sum(nat_mort[source],fish_mort[catch_source])
    
    abundance[sink] <- abundance[sink] + trans_prop*abundance[source]*exp(-total_mort[source])

    catch[which(catch_year==year[sink] &
                  catch_age==age[sink] &
                  catch_length==length[sink])]<-catch[which(catch_year==year[sink] &
                                                                          catch_age==age[sink] &
                                                                          catch_length==length[sink])] + (fish_mort[which(catch_year==year[source] &
                                                                                                                                        catch_age==age[source] &
                                                                                                                                        catch_length==length[source])]/total_mort[source])*trans_prop*abundance[source]*(1-exp(-total_mort[source])) 
    
  }else{
    total_mort[source] <- sum(nat_mort[source],fish_mort[which(catch_year==year[source] &
                                                                     catch_age==age[source] &
                                                                     catch_length==length[source])])
    
    abundance[sink] <- trans_prop*abundance[source]*exp(-total_mort[source])
    
    catch[which(catch_year==year[sink] &
                  catch_age==age[sink] &
                  catch_length==length[sink])]<-catch[which(catch_year==year[sink] &
                                                                          catch_age==age[sink] &
                                                                          catch_length==length[sink])] + (fish_mort[which(catch_year==year[source] &
                                                                                                                                      catch_age==age[source] &
                                                                                                                                      catch_length==length[source])]/total_mort[source])*trans_prop*abundance[source]*(1-exp(-total_mort[source]))    
  }
}



#TODO: Still need to write up logic for interpolating between abundance values
#to calculate length comps.
#
for(i in seq_along(fleets)){
  
  for(j in seq_along(years)){
    
    for(k in seq_along(ages)){
      
      for(j in seq_along(comp_lengths)){
        
        temp_proportion <- (comp_lengths[j] - model_L_lower[j])/(model_L_upper[j] - model_L_lower[j])
        
        pred_comp[j] <- temp_proportion
      }
    } 
  }
}

for(i in seq_along(pred_obs)){
  
  
  
}

#
#


