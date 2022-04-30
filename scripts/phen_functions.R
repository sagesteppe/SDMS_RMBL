
#################################
###      weibull_wrapper     ####
#################################

weibull_wrapper <- function(x){
  
  taxon <- x$Species[1]
  observations <- x$DOY
  sample_size <- nrow(x)
  
  frst_hrb <- min(observations)
  lst_hrb <- max(observations)
  
  onset_est <- weib_percentile_ci(observations = observations, iterations = 10, conf = 0.90,
                              percentile = 0.01, bootstraps = 250, parallelize = "multicore", ncpus = 16)
  tenth_est <- weib_percentile_ci(observations = observations, iterations = 10, conf = 0.90,
                              percentile = 0.1, bootstraps = 250, parallelize = "multicore", ncpus = 16)
  fifty_est <- weib_percentile_ci(observations = observations, iterations = 10, conf = 0.90,
                              percentile = 0.5, bootstraps = 250, parallelize = "multicore", ncpus = 16)
  ninety_est <- weib_percentile_ci(observations = observations, iterations = 10, conf = 0.90,
                                 percentile = 0.9, bootstraps = 250, parallelize = "multicore", ncpus = 16)
  end_est <- weib_percentile_ci(observations = observations, iterations = 10, conf = 0.90,
                            percentile = 0.99, bootstraps = 250, parallelize = "multicore", ncpus = 16)
  
  a <- bind_cols(
    data.frame(t(data.frame(onset_est))),
    data.frame(t(data.frame(tenth_est))),
    data.frame(t(data.frame(fifty_est))),
    data.frame(t(data.frame(ninety_est))),
    data.frame(t(data.frame(end_est)))
  ) %>% 
    purrr::set_names(c("onset", "tenth", "fifty", "ninety", "end")) %>% 
    mutate(across(.cols = everything(), round, 0)) %>% 
    rownames_to_column('event') %>% 
    pivot_longer(!event, names_to = 'metric', values_to = 'DOY') %>% 
    cbind(., taxon, sample_size, frst_hrb, lst_hrb) %>% 
    mutate(est_duration = round(
      filter(., event == 'estimate' & metric == "end")$DOY -
        filter(., event == 'estimate' & metric == "onset")$DOY )
      ) %>% 
    mutate(DOY = if_else(event == 'high_ci' & DOY >= 305, 305, DOY))
  
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
  
  #onset_est <- frst_hrb + ((lst_hrb - frst_hrb) * 0.80)
  #end_est <- lst_hrb + ((lst_hrb - frst_hrb) * 1.20)
  
  
  onset_est <- if_else(
    (lst_hrb - frst_hrb) * 0.80 >= 14, 
    frst_hrb - 14, 
    frst_hrb + ((lst_hrb - frst_hrb) * 0.80)
    )
  
  end_est <- if_else(
    (lst_hrb - frst_hrb) * 1.20 >= 14, 
    lst_hrb + 14, 
    lst_hrb + ((lst_hrb - frst_hrb) * 1.20)
    )
  
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
####       lwr_snw_bnd        ####
##################################

lwr_snw_bnd <- function(snow_data, site){
  
  # this serves to reduce the periodicity by which the weibull distribution
  # generates biologically impossible estimates of Phenological events based on 
  # restriction of onset events via Snow cover. (See
  # Iler et al. 2021 for details ).
  
  # INPUTS, snow_data = a list of files paths to .nc files of snow cover see
  # https://doi.org/10.6084/m9.figshare.5902381.v4 
  # site = a study location as an sf object. In this case a single site, 
  # this untested with multipoint/line/poly objects. 
  
  lwr_snw_bnd_res <- vector(mode = "list", length = length(snow_data))
  for (i in 1:length(snow_data)){
    
    x_year_snow <- brick(snow_data[i])
    lwr_snw_bnd_res[i] <- data.frame('Snow' = t(extract(x_year_snow, site))) %>% 
      rownames_to_column('Date') %>% 
      mutate(Date = str_remove(Date, 'X'),
             DOY = yday(as.Date(Date, format = "%Y.%m.%d")) +1,
             Snow_days = cumsum(Snow)) %>% 
      arrange(DOY) %>% 
      filter(Snow == 0, DOY <= 212) %>% 
      slice_max(Snow_days, with_ties = F)
  }
  
  lwr_snw_bnd_res <- data.frame('Date' = unlist(lwr_snw_bnd_res)) %>% 
    mutate(DOY = yday(as.Date(Date, format = "%Y.%m.%d")),
           third_quart = round(quantile(DOY, 0.75), 0)) 
  
  return(lwr_snw_bnd_res)
}


##################################
####       uppr_snw_bnd       ####
##################################

uppr_snw_bnd <- function(snow_data, site){
  
  # this serves to reduce the periodicity by which the weibull distribution
  # generates biologically impossible estimates of Phenological events based on 
  # restriction of onset events via Snow cover. (See
  # Iler et al. 2021 for details ).
  
  # INPUTS, snow_data = a list of files paths to .nc files of snow cover see
  # https://doi.org/10.6084/m9.figshare.5902381.v4 
  # site = a study location as an sf object. In this case a single site, 
  # this untested with multipoint/line/poly objects. 
  
  upr_snw_bnd_res <- vector(mode = "list", length = length(snow_data))
  for (i in 1:length(snow_data)){
    
    x_year_snow <- brick(snow_data[i])
    upr_snw_bnd_res[i] <- data.frame('Snow' = t(extract(x_year_snow, site))) %>% 
      rownames_to_column('Date') %>% 
      mutate(Date = str_remove(Date, 'X'),
             DOY = yday(as.Date(Date, format = "%Y.%m.%d")) +1,
             Snow_days = RcppRoll::roll_mean(Snow, n = 4, na.rm=TRUE, align="right", fill = NA)) %>% 
      group_by(Snow_days) %>% 
      arrange(DOY) %>% 
      filter(Snow == 1, DOY >= 212, Snow_days == 0.75) %>% 
      slice_min(Snow_days, with_ties = F)
  }
  
  upr_snw_bnd_res <- data.frame('Date' = unlist(upr_snw_bnd_res)) %>% 
    mutate(DOY = yday(as.Date(Date, format = "%Y.%m.%d")),
           first_quart = round(quantile(DOY, 0.25), 0)) 
  
  return(upr_snw_bnd_res)
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