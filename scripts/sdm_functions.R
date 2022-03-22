
########################
# ENSURE MULTIPOLYGONS #
########################

ensure_multipolygons <- function(X) { 
  # holla @ Josh Brian on GISStackExchange for this solution. 
  tmp1 <- tempfile(fileext = ".gpkg")
  tmp2 <- tempfile(fileext = ".gpkg")
  st_write(X, tmp1)
  ogr2ogr(tmp1, tmp2, f = "GPKG", nlt = "MULTIPOLYGON")
  Y <- st_read(tmp2)
  st_sf(st_drop_geometry(X), geom = st_geometry(Y))
}


##############
# PSEUDO REP #
##############

pseudo_rep_funct <- function(x){
  
  set1 <- x # need to feed in two sets of coordinates. 
  set2 <- x
  binomial <- x %>% pull(binomial)
  binomial <- binomial[1] # append to out
  nearest_neigh <- nngeo::st_nn(set1, set2, 
                                k = 2, # self and nearest neighbor 
                                progress = F,
                                maxdist = 128) # only search up to the distance of a cell.  
  nearest_neigh <- lapply(nearest_neigh, as.data.frame) # reformat data here....
  nearest_neigh <-Map(cbind, nearest_neigh , paste0("i_", 
                                                    1:length(nearest_neigh)))
  column_names <- c("neighbor", "individual")
  nearest_neigh <- lapply(nearest_neigh, setNames, column_names)
  
  nearest_neigh <- data.table::rbindlist(nearest_neigh) %>% 
    mutate(individual = str_remove(individual, "i_")) %>% 
    mutate(individual = as.numeric(individual)) %>% 
    group_by(individual) %>% 
    mutate(neigh_id = seq_along(individual)) %>% 
    filter(neigh_id != 1)
  
  if(nrow(nearest_neigh) == 0)
  {
    return(x)
  }
  else
  {
    dist_remove <- as.matrix(nearest_neigh[,1:2]) # we need to recover each pair of points within this distance
    dist_remove <- dist_remove[!duplicated(apply(dist_remove, 1, sort), 
                                           MARGIN = 2),]
    dist_remove <- as.data.frame(dist_remove)
    
    if(ncol(dist_remove) == 2){
      dist_remove <- dist_remove %>% pull(individual)
    }
    else{
      dist_remove <- dist_remove[2,]
    }
    
    focal <- x %>% 
      rowid_to_column() %>% 
      filter(!rowid %in% dist_remove) %>% 
      filter(!rowid %in% dist_remove) %>% 
      dplyr::select(-rowid)
  }
  
  return(focal)
} 

####################
# MULTI VARIOGRAM  #
####################

multi_variogram <- function(x){
  
  taxon <- x %>% na.exclude() %>% pull(binomial)
  taxon <- taxon[1]
  
  variogram = gstat::variogram(records_cell ~ 1, data = x)
  vario_fit <- gstat::fit.variogram(variogram, gstat::vgm(c("Exp", "Mat", "Sph")), fit.kappa = TRUE)
  vario_fit <- cbind(taxon, vario_fit)
  vario_fit <- vario_fit %>% 
    dplyr::select(!starts_with(c('ang','anis')))
  
  return(vario_fit)
}


################
# MORANS INDEX #
################

Morans_index_func <- function(x){
  
  y <- x %>% pull(records_cell)
  taxon <- x %>% na.exclude() %>% pull(binomial)
  taxon <- taxon[1]
  
  Moran <- moran.test(y, listw, randomisation = T) 
  mor.i <- Moran[["estimate"]][["Moran I statistic"]] 
  mor.exp. <- Moran[["estimate"]][["Expectation"]]
  mor.var <- Moran[["estimate"]][["Variance"]]
  mor.p <- Moran[["p.value"]] 
  
  out <-  data.frame(taxon, mor.i, mor.exp., mor.var, mor.p)
}




###############
# SA REDUCER  #
###############

SA_reducer <- function(x){
  set1 <- x
  set2 <- x
  binomial <- x %>% pull(binomial)
  binomial <- binomial[1]
  nearest_neigh <- nngeo::st_nn(set1, set2, # find the nearest neighbors with the range of auto-correlation.
                                k = 10, # 10 neighbors. 
                                maxdist = set1$range[1]) # find 10 closest records within range of SA. 
  nearest_neigh <- lapply(nearest_neigh, as.data.frame) # reformat data here....
  nearest_neigh <-Map(cbind, nearest_neigh , paste0("i_", 1:length(nearest_neigh)))
  column_names <- c("neighbor", "individual")
  nearest_neigh <- lapply(nearest_neigh, setNames, column_names)
  
  nearest_neigh <- data.table::rbindlist(nearest_neigh) %>% 
    mutate(individual = str_remove(individual, "i_")) %>% 
    mutate(individual = as.numeric(individual)) %>% 
    group_by(individual) %>% 
    mutate(neigh_id = seq_along(individual)) %>% 
    filter(neigh_id != 1)
  
  y <- set1 %>% 
    rowid_to_column() %>% 
    dplyr::select(rowid) # has all original observations. 
  
  outside_range <- nearest_neigh %>% # return the records outside the original range of spatial auto-correlation
    ungroup() %>% 
    distinct(individual) %>% 
    anti_join(y, ., by = c('rowid' = "individual")) %>% 
    rename('individual' = rowid) %>% 
    dplyr::select(individual)
  
  focal <- left_join(nearest_neigh, y, by = c("individual" =  "rowid")) %>% 
    st_as_sf() %>% 
    group_by(individual) %>% 
    arrange(individual)
  neighbors <- left_join(nearest_neigh, y, by = c("neighbor" =  "rowid")) %>% 
    st_as_sf() %>% 
    group_by(individual) %>% 
    arrange(individual)
  
  focal$dist <- st_distance(focal, neighbors, by_element = T)
  
  top_neighbors <- nearest_neigh %>% # remove the most connected neighbors, i.e. those with the most connections within the range.
    ungroup() %>% 
    add_count(neighbor) %>% 
    distinct(neighbor, .keep_all = T) %>% 
    slice_max(n, prop = 0.025, with_ties = F) %>% 
    pull(neighbor)
  focal <- focal %>% 
    filter(!individual %in% top_neighbors) %>% 
    filter(!neighbor %in% top_neighbors) 
  
  dist_remove <- focal %>% # now remove from a focal records any neighbors within 270 meters - ala spatial thinning. 
    filter(dist <= units::as_units(1000, "meter")) %>% 
    ungroup() %>% 
    arrange(individual) %>% 
    st_drop_geometry() %>% 
    dplyr::select(individual, neighbor) 
  dist_remove <- as.matrix(dist_remove) # we need to recover each pair of points within this distance
  dist_remove <- dist_remove[!duplicated(apply(dist_remove, 1, sort), MARGIN = 2),]
  dist_remove <- as.data.frame(dist_remove) %>% pull(individual) # now we only remove one member of the pair. 
  
  focal <- focal %>% 
    filter(!individual %in% dist_remove) %>% 
    filter(!neighbor %in% dist_remove)
  
  # now remove the 1% of the remaining most spatially auto-correlated points.
  testy <- as.data.frame(focal) %>%  
    ungroup() %>% 
    add_count(neighbor) %>% 
    dplyr::distinct(neighbor, .keep_all = T) %>% 
    slice_min(n, prop = 0.99, with_ties = F) %>% 
    st_as_sf() 
  
  testy_follow <- testy %>% # reducing the records from neighbors to the remaining individuals
    distinct(individual, .keep_all = T) %>% 
    dplyr::select(individual)
  
  out <- rbind(testy_follow, outside_range)# add the individuals with no SA to those pruned with SA. 
  out <- cbind(binomial, out)
  
}



###################
# TRUE ABSENCE ML #
###################

true_absence_ML <- function(x){ # for collecting true absence records from BLM land. 
  
  taxon <- x %>% distinct(binomial) %>%  pull(binomial)
  taxon <- taxon[1]
  TA_req <- round(x$no_record[1] * prop_blm, 0) # how many true absences needed?
  presence_PK <- x %>% pull(PlotKey) # which plots can absences not occur in because a presence is there?
  
  AIM_absence <- AIM_points %>% # remove plots with an occurrence of the taxon. 
    filter(!PlotKey %in% presence_PK)  %>% # make sure the plot does not have a presence record
    mutate('binomial' = taxon) %>% 
    mutate(no_record = TA_req)
    
  AIM_absence <- AIM_absence[sample(1:nrow(AIM_absence), size =  TA_req, replace = F),]
  
  out <- rbind(x, AIM_absence)
  return(out)
}


####################
# TRUE ABSENCE REG #
####################

true_absence_REG <- function(x){ # for collecting true absence records from BLM land. 
  
  taxon <- x %>% distinct(binomial) %>% pull(binomial)
  taxon <- taxon[1]
  TA_req <- round(1000 * prop_blm, 0) # For regression without ML we need 1k records * prop of blm land.  
  presence_PK <- x %>% pull(PlotKey) 
  
  AIM_absence <- AIM_points %>% 
    filter(!PlotKey %in% presence_PK) %>% 
    mutate('binomial' = taxon) %>% 
    mutate(no_record = TA_req) 
  
  AIM_absence <- AIM_absence[sample(1:nrow(AIM_absence), size =  TA_req, replace = F),]
  
  out <- rbind(x, AIM_absence) 
  return(out) 
}



###########################################
## BINOMIAL LOGISTIC REGRESSION ABSENCES ##
###########################################
BinLogReg_abs <- function(x, y){ # for collecting true absence records from BLM land. 
  
  #input 
  # x = dataset containing new presences
  # y = dataset containing the TRUE absences for modelling and ensembling
  
  taxon <- x %>% distinct(binomial) %>% pull(binomial)
  taxon <- taxon[1]
  presence_PK <- x %>% pull(PlotKey) %>% na.omit() 
  absences_PK <- y %>% pull(PlotKey) %>% na.omit()
  Pks <- c(presence_PK, absences_PK)
  samp_req <- nrow(x)
  
  AIM_to_samp <- AIM %>% # remove plots which have already been used
    filter(!PlotKey %in% Pks)  
  
  # verify that no occurrence of the taxon is in the plots to sample
  AIM_removals <- AIM_to_samp %>% 
    filter(binomial == taxon) %>% 
    distinct(PlotKey) %>% 
    pull(PlotKey) 
  # now remve the plots which could have another presence record in them
  AIM_to_samp <- AIM_to_samp %>% 
    filter(!PlotKey %in% AIM_removals)  %>% 
    mutate(binomial = taxon,
           occurrence = 0) 
    
  
  new_absence <- AIM_to_samp[sample(1:nrow(AIM_to_samp), size =  samp_req, replace = F),]
  
  out <- rbind(x, new_absence)
  return(out)
}
#####################
# RANDOM PA SPATIAL #
#####################

random_PA_spatial <- function(x){
  
  taxon <- x %>% pull(binomial)
  taxon <- taxon[1]
  no_ab <- x %>% filter(occurrence == 0) %>% nrow() 
  PA_req <- 1200 - no_ab # we will make 1200, and then remove those with high similarity to the target points. 
  
  x_buf <- st_buffer(x, 9000)
  blm_cp_sUB <- st_erase(x_buf, blm_cp_s)

  pseudo_abs <- st_sample(blm_cp_sUB, size = PA_req, type = 'random', by_polygon = F) %>% 
    st_as_sf() %>% 
    mutate(occurrence = 0) %>% 
    mutate(binomial = taxon) %>% 
    mutate(PlotKey = NA) %>% 
    rename('geometry' = x) %>% 
    mutate(no_record = nrow(.))
  
  x <- st_centroid(x)
  
  out <- rbind(x, pseudo_abs)
  
}


###############
# MEAN VALUES #
###############

mean_values <- function(x){
  
  binomial <- x %>% pull(binomial)
  binomial <- binomial[1]
  out <- raster::extract(WPDPV_rs, x, fun = mean, method = "simple")
  out1 <- as.data.frame(t(apply(out, 2, mean)))
  out2 <- as.data.frame(t(apply(out, 2, sd)))
  out3 <- as.data.frame(t(apply(out, 2, se)))
  out1$var <- "mean"  
  out2$var <- "sd" 
  out3$var <- "se"
  
  out1 <- cbind(binomial, out1)
  out2 <- cbind(binomial, out2)
  out3 <- cbind(binomial, out3)
  
  rbind(out1, out2, out3)
  
}




##############
# percentile #
##############
se <- function(x) sqrt(var(x)/length(x)) 

percentile <- function(x){
  
  my_quantile <- function(x, probs) {
    tibble(x = quantile(x, probs, na.rm = T))
  }
  #Input a dataframe of scaled and centered variables. 
  
  x1 <- st_drop_geometry(x)
  scaled_matrix <- as.matrix(x1[,5:20])
  maxs <- apply(scaled_matrix, 2, max)
  mins <- apply(scaled_matrix, 2, min)
  scaled_matrix <- scale(scaled_matrix, center = mins, scale = maxs - mins)
  
  # split into two matrices   - PRESENCE records generate target values
  #                           - BLM Absences are withheld from analysis
  #                           - The random absences are evaluated against PRESENCE
  
  
  # Presence Matrix 
  presence_matrix <- scaled_matrix[which(x$occurrence == 1), ] 
  percentiles <- as_tibble(presence_matrix) %>% 
    summarise(across(.cols = everything(), ~my_quantile(.x, probs = c(0.4, 0.6)))) %>% 
    as.matrix()
  
  lower <- as.numeric(percentiles[1,])
  higher <- as.character(percentiles[2,])
  
  # Absence matrix
  absence_matrix <- matrix(scaled_matrix[which(x$occurrence == 0 & is.na(x$PlotKey) ), ]
                           , byrow = T, ncol = 16)
  
  m1  <- absence_matrix
  m1[] <- c(0, 1)[(m1 > lower & m1 < higher) +1]
  m1 <- matrix(m1, byrow = T, ncol = 16)

  selections <- data.frame(twenty = rowSums(m1, na.rm = T))
  absence_matrix <- data.frame(absence_matrix)
  
  absence_matrix <- cbind(absence_matrix, selections)  %>% 
    mutate(ID  = row_number()) %>% 
    drop_na(1:8) %>% 
    slice_min(twenty, prop = 0.8, with_ties = F)
  
  # filter to the appropriate records
  
  base <- nrow(x) - nrow(m1)
  absence_positions <- absence_matrix %>% 
    mutate(ID  = ID + base) %>% 
    pull(ID)
  
  sub_sequence <- c(seq(from = 1, to = base, by = 1), absence_positions)
  
  results <- x[sub_sequence,]

 return(results)
  
}


#############
# ENVI DIST #
#############

envi_dist <- function(x){
  
  focal_taxon <- x %>% 
    pull(binomial) 
  focal_taxon <- focal_taxon[1] 
  
  focal_tax_mean <- mean_abiotic_rescaled %>% 
    filter(binomial == focal_taxon) %>% 
    filter(var == 'mean') %>% 
    dplyr::select(-var, -binomial, -X) 
  
  no_presence <- x %>% 
    pull(no.record)
  no_presence <- no_presence[1]
  true_absence <- nrow(x) - no_presence
  PA_req <- round(no_presence - true_absence, 0)
  
  dists <- apply(data_subset, 1, function(row) sum((row - focal_tax_mean) ** 2)) 
  dist.order <- order(dists)
  
  cell_nums <- round(seq(from = 10000, to = nrow(data_subset), length.out = PA_req),0)
  cell_values <- dist.order[cell_nums]
  raster_values <- data_subset[rownames(data_subset) %in% cell_values, ]
}


##############
# LINEAR SDM #
##############

Linear_SDM <- function(x){
  
  binomial <-  x %>% pull(binomial)
  binomial <- binomial[1]
  taxon <- x %>% 
    mutate(occurrence = case_when(occurrence == 2 | occurrence == 1 ~ 1,
                                  occurrence == 0 ~ 0)) %>% 
    dplyr::select(occurrence) %>% 
    as(., "Spatial")
  taxon = spTransform(taxon,geo_proj)
  
  
  sdm_data_obj <- sdmData(formula=occurrence~., 
                          train = taxon, 
                          predictors = WPDPV2)
  sdm_model <- sdm(formula = occurrence~., 
                   data = sdm_data_obj, 
                   methods = c('glm','gam'), replication='sub', test.percent=30, n=3)
  
  fname <- paste0(here(), '/results/maps/', binomial, "_glm_", Sys.time(),'.tif')
  fname <- gsub(' ', '_', fname)
  sdm_ensemble_prediction <- sdm::ensemble(sdm_model, WPDPV2, 
                                           setting = list(method = "weighted", stat = 'tss', opt = 2), 
                                           filename = fname)
  
  fname <- paste0(here(), '/results/stats/', binomial, "_glm_", Sys.time(),'.csv')
  fname <- gsub(' ', '_', fname)
  evaluation <- getEvaluation(sdm_model, stat=c('TSS','Kappa','AUC'), wtest=c('training','test'), opt = 1)
  write.csv(evaluation, file = fname)
  
  # Predictor variables importance
}



############################
###    Variance Calcs    ###
############################

variance_calcs  <- function(x){
  
  my_quantile <- function(x, probs) {
    tibble(x = quantile(x, probs, na.rm = T))
  }
  se <- function(x) sqrt(var(x)/length(x)) 
  not_all_na <- function(x) any(!is.na(x))
  
  #Input a dataframe of scaled and centered variables. 
  
  binomial <- x$binomial[1]
  x1 <- sf::st_drop_geometry(x)
  scaled_matrix <- as.matrix(x1[,5:20])
  maxs <- apply(scaled_matrix, 2, max)
  mins <- apply(scaled_matrix, 2, min)
  scaled_matrix <- scale(scaled_matrix, center = mins, scale = maxs - mins)
  
  # Select only presence records and then calculate 
  # measures of dispersion for them
  presence_matrix <- scaled_matrix[which(x$occurrence == 1), ] 
  
  standard_err <- as_tibble(presence_matrix) %>% 
    summarise(across(.cols = everything(), ~ se(.x))) %>% 
    as.matrix()

  variance <- as_tibble(presence_matrix) %>% 
    summarise(across(.cols = everything(), ~ var(.x))) %>% 
    as.matrix()
  
  percentiles <- as_tibble(presence_matrix) %>% 
    summarise(across(.cols = everything(), 
                     ~my_quantile(.x, probs = c(0.25, 0.75)))) %>% 
    mutate(across(.cols = everything(), 
                  ~ lag(.x))) %>% 
    as.matrix() 
  
  percentiles <- percentiles[2,]

  # prettify results
  results <- rbind(standard_err, variance, percentiles)
  length_col <- as.data.frame(results) %>% select(where(not_all_na))
  length_col <- ncol(length_col)
  results <- rowSums(results, na.rm = T)/length_col
  results <- cbind(binomial, results)
  results <- cbind(t(data.frame('SE','VAR','IQR')), results)
  row.names(results) <- NULL
  colnames(results) <- c('Dispersion', 'binomial', 'Value')
  results <- as.data.frame(results) %>% mutate(across(.cols = 3:ncol(results), as.numeric))
  
  return(results)
}



#######################
###    total_area   ###
#######################


total_area <- function(x){
  
  binomial <- x$binomial[1]
  presences <- x %>% 
    mutate(occurrence = ifelse(occurrence == 2, 0, occurrence)) %>% 
    filter(occurrence == 1) %>% 
    sf::st_union()
  
  hull <- sf::st_convex_hull(presences)
  area <- as.numeric(sf::st_area(hull))
  results <- data.frame(cbind(binomial, area))
  
  return(results)
}


############################
####    gbif_collector  ####
############################
 
# wrote one and what happened?


gbif_collector <- function(dataset){
  
  downloaded <- occ_data(
    scientificName = dataset, 
    hasGeospatialIssue = F,
    geometry = wkt_boundary,
    kingdomKey = 6,
    limit=1000
    )
  
  data <- downloaded[['data']][,c(1:8, 14:17,26,32,38,44,51,55,62,70)]
  data$endOfRecords <- downloaded[['meta']]$endOfRecords
  
  return(data)
  Sys.sleep(1)
}
