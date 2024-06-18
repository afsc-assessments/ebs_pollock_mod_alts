library(dplyr) # sry Nathan! :D
library(tidyr)

# get unique ID (integer value) for each data set
get_id <- function(myaux){
  myaux <- myaux %>% 
    dplyr::group_by(obs_type, nll_type, fleet) %>% 
    dplyr::mutate(id = dplyr::cur_group_id(), .before = obs_type) %>% 
    dplyr::ungroup()
  return(myaux)
}


#Add a likelihood calculation index to the obs data frame
get_likelihood_index <- function(myaux) {
  # myaux <- input$obsdf
  myaux$likelihood_index <- NA
  like_index <- 1
  # Loop over all observations to identify multinomial
  # observations with linked likelihoods
  for(i in seq_along(myaux$obs_type)){
    if(i==1){
      #Set the first observation at index 1
      myaux$likelihood_index[i] <- like_index 
    }else if(myaux$obs_type[i]!=myaux$obs_type[i-1]){
      #If observation type changes then increment
      like_index <- like_index + 1
      myaux$likelihood_index[i] <- like_index 
    }else if(myaux$obs_type[i]<2){
      #All catches and indices are independent values so always increment
      like_index <- like_index + 1
      myaux$likelihood_index[i] <- like_index 
    }else if(myaux$year[i]!=myaux$year[i-1]){
      #Increment if the year of composition changes
      like_index <- like_index + 1
      myaux$likelihood_index[i] <- like_index 
    }else if(myaux$fleet[i]!=myaux$fleet[i-1]){
      #Increment if the fleet of composition changes
      like_index <- like_index + 1
      myaux$likelihood_index[i] <- like_index 
    }else{
      #If year and fleet don't change then comps will be combined for a 
      #single multinomial likelihood calculation
      myaux$likelihood_index[i] <- like_index 
    }
  }
  return(myaux)
}

# get prediction dataframe for observations (expanding to the full range of
# observations)
get_pred <- function(myaux, myinput) {
  # myaux <- dat$aux; myinput <- input
  lkup <- myaux %>% dplyr::distinct(id, obs_type, nll_type, fleet)
  myaux <- myaux %>% dplyr::filter(year %in% myinput$year)
  
  pred <- dplyr::bind_rows(
    # catch and index data
    myaux %>% 
      dplyr::filter(obs_type %in% c(0, 1)) %>% 
      dplyr::select(id, year, obs, obserror, age, len, likelihood_index) %>% 
      tidyr::complete(id, year = myinput$year),
    # age comp data
    myaux %>% 
      dplyr::filter(obs_type %in% c(2)) %>% 
      dplyr::select(id, year, obs, obserror, age, len, likelihood_index) %>% 
      tidyr::complete(id, year = myinput$year, nesting(age)), 
    # len comp data
    myaux %>% 
      dplyr::filter(obs_type %in% c(3)) %>% 
      dplyr::select(id, year, obs, obserror, age, len, likelihood_index) %>% 
      tidyr::complete(id, year = myinput$year, nesting(len))
  )
  pred <- pred %>% 
    dplyr::left_join(lkup) %>% 
    dplyr::relocate(names(lkup), .before = year) 
  
  return(pred)
}
