
#################################
###      weibull_wrapper     ####
#################################

weibull_wrapper <- function(x){
  
  taxon <- x$Species[1]
  observations <- x$DOY
  sample_size <- nrow(x)
  
  frst_hrb <- min(observations)
  lst_hrb <- max(observations)
  
  onset_est <- weib.limit(observations, upper=FALSE, alpha = .2)
  end_est   <- weib.limit(observations, upper=TRUE , alpha = .2)
  
  a <- bind_rows(
    data.frame(t(data.frame(onset_est))),
    data.frame(t(data.frame(end_est)))
  ) %>% 
    mutate(across(.cols = everything(), round, 0)) %>% 
    rownames_to_column('event') %>% 
    cbind(., taxon, sample_size, frst_hrb, lst_hrb) %>% 
    mutate(est_duration = round(estimate - lag(estimate, default = first(estimate)), 0)) %>% 
    mutate(upper.ci = if_else(upper.ci >= 305, 305, upper.ci))
  
  return(a)
}

#################################
####       herb_wrapper      ####
#################################

herb_wrapper <- function(x){
  
  taxon <- x$Species[1]
  observations <- x$DOY
  sample_size <- nrow(x)
  
  frst_hrb <- min(observations)
  lst_hrb <- max(observations)
  
  onset_est <- frst_hrb + ((lst_hrb - frst_hrb) * 0.75)
  end_est <- lst_hrb + ((lst_hrb - frst_hrb) * 1.25)
  
  a <- bind_cols(
    data.frame('onset_est' = t(data.frame(onset_est))),
    data.frame('end_est' = t(data.frame(end_est)))
  ) %>% 
    mutate(across(.cols = everything(), round, 0)) %>% 
    cbind(., taxon, sample_size, frst_hrb, lst_hrb) %>% 
    pivot_longer(!taxon:lst_hrb, names_to = 'event', values_to = 'estimate') %>% 
    mutate(est_duration = round(estimate - lag(estimate, default = first(estimate)), 0)) %>% 
    
    # DO SOMETHING HERE - SAY IF ESTIMATE ABOVE > 3 WEEKS THEN USE 3 WEEKS BEFORE AND AFTER
    # KNOWN EVENTS TO CONSTRAN THE ESTIMATE WITHIN REASON
  
  return(a)
}


##################################
####       daily_flowers      ####
##################################

daily_flowers <- function(x){
  
  # these create a dataframe of plant taxa which are flowering 
  # each day of a season, from the earliest flowering taxon to the last 
  # flowering taxon within the list
  
  results <- data.frame(
    seq(from = min(x$estimate, na.rm = T), to = max(x$estimate, na.rm = T), by = 1),
    x$taxon[1]
  )
  colnames(results) <- c('DOY', 'Taxon')
  return(results)
}