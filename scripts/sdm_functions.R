
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
  TA_req <- round(x$no.record[1] * prop_blm, 0) # how many true absences needed?
  presence_PK <- x %>% pull(PlotKey) # which plots can absences not occur in because a presence is there?
  
  AIM_absence <- AIM_points %>% # remove plots with an occurrence of the taxon. 
    filter(!PlotKey %in% presence_PK) %>% # make sure the plot does not have a presence record
    dplyr::sample_n(TA_req, replace = F) %>%   # sample this many plots
    mutate('binomial' = taxon) %>% 
    mutate(no.record = TA_req)
  
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
    dplyr::sample_n(TA_req, replace = F) %>% 
    mutate('binomial' = taxon) %>% 
    mutate(no.record = TA_req) 
  
  out <- rbind(x, AIM_absence) 
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
  
  pseudo_abs <- st_sample(blm_cp_sUB, size = PA_req, type = 'random') %>% 
    st_as_sf() %>% 
    mutate(occurrence = 0) %>% 
    mutate(binomial = taxon) %>% 
    mutate(PlotKey = NA) %>% 
    rename('geometry' = x) %>% 
    mutate(no.record = nrow(.))
  
  out <- rbind(x, pseudo_abs)
  
}


###############
# MEAN VALUES #
###############

mean_values <- function(x){
  
  binomial <- x %>% pull(binomial)
  binomial <- binomial[1]
  out <- raster::extract(bioclim.ml.scaled, x, fun = mean, method = "simple")
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
  
  binomial <- x@data['binomial'][1]
  
  sdm_data_obj <- sdmData(formula=occurrence~., 
                          train = x, 
                          predictors = bioclim)
  sdm_data_obj_1 <- sdm(formula = occurrence~., 
                        data = sdm_data_obj, 
                        methods = c('glm','gam'), replication='sub', test.percent=90, n=2)
  ensemb <- ensemble(sdm_data_obj_1, bioclim, 
                     setting = list(method = "weighted", stat = 'tss', opt = 2), 
                     filename = paste0(binomial, "_","lm", Sys.time(),'.tif'))
  
  # Predictor variables importance
  
  PrVar1 <- cbind("GLM_1","training",sdm_data_obj_1@models[["occurrence"]][["glm"]][["1"]]@varImportance[["training"]]@varImportance)
  #PrVar2 <- cbind("GLM_2","training",sdm_data_obj_1@models[["occurrence"]][["glm"]][["2"]]@varImportance[["training"]]@varImportance)
  PrVar3 <- cbind("GAM_1","training",sdm_data_obj_1@models[["occurrence"]][["gam"]][["3"]]@varImportance[["training"]]@varImportance)
  #PrVar4 <- cbind("GAM_2","training",sdm_data_obj_1@models[["occurrence"]][["gam"]][["4"]]@varImportance[["training"]]@varImportance)
  
  PrVar5 <- cbind("GLM_1","test", sdm_data_obj_1@models[["occurrence"]][["glm"]][["1"]]@varImportance[["test.dep"]]@varImportance)
  #PrVar6 <- cbind("GLM_2","test", sdm_data_obj_1@models[["occurrence"]][["glm"]][["2"]]@varImportance[["test.dep"]]@varImportance)
  PrVar7 <- cbind("GAM_3","test", sdm_data_obj_1@models[["occurrence"]][["gam"]][["3"]]@varImportance[["test.dep"]]@varImportance)
  #PrVar8 <- cbind("GAM_4","test", sdm_data_obj_1@models[["occurrence"]][["gam"]][["4"]]@varImportance[["test.dep"]]@varImportance)
  
  Mt_train1 <- cbind("GLM_1","training",sdm_data_obj_1@models[["occurrence"]][["glm"]][["1"]]@evaluation[["training"]]@threshold_based)
  Mt_train2 <- cbind("GAM_1","training",sdm_data_obj_1@models[["occurrence"]][["gam"]][["2"]]@evaluation[["training"]]@threshold_based)
  Mt_test1  <- cbind("GLM_1", "test",sdm_data_obj_1@models[["occurrence"]][["glm"]][["1"]]@evaluation[["test.dep"]]@threshold_based)
  Mt_test2  <- cbind("GAM_1","test",sdm_data_obj_1@models[["occurrence"]][["gam"]][["2"]]@evaluation[["test.dep"]]@threshold_based)
  
  
  GLM_stats1 <- as.data.frame(cbind("Prev.",sdm_data_obj_1@models[["occurrence"]][["glm"]][["1"]]@evaluation[["training"]]@statistics[["Prevalence"]]))
  GLM_stats2 <- as.data.frame(cbind("AUC",sdm_data_obj_1@models[["occurrence"]][["glm"]][["1"]]@evaluation[["training"]]@statistics[["AUC"]]))
  GLM_stats3 <- as.data.frame(cbind("COR",sdm_data_obj_1@models[["occurrence"]][["glm"]][["1"]]@evaluation[["training"]]@statistics[["COR"]][["cor"]]))
  GLM_stats4 <- as.data.frame(cbind("p.val.",sdm_data_obj_1@models[["occurrence"]][["glm"]][["1"]]@evaluation[["training"]]@statistics[["COR"]][["p.value"]]))
  GLM_stats5 <- as.data.frame(cbind("dev",sdm_data_obj_1@models[["occurrence"]][["glm"]][["1"]]@evaluation[["training"]]@statistics[["Deviance"]]))
  GLM_stats <- rbind("GLM_1", "training", GLM_stats1, GLM_stats2 ,GLM_stats3, GLM_stats4,GLM_stats5)
  
  GAM_stats1 <- cbind("Prev.","training",sdm_data_obj_1@models[["occurrence"]][["gam"]][["2"]]@evaluation[["training"]]@statistics[["Prevalence"]])
  GAM_stats2 <- cbind("AUC","training",sdm_data_obj_1@models[["occurrence"]][["gam"]][["2"]]@evaluation[["training"]]@statistics[["AUC"]])
  GAM_stats3 <- cbind("COR","training",sdm_data_obj_1@models[["occurrence"]][["gam"]][["2"]]@evaluation[["training"]]@statistics[["COR"]][["cor"]])
  GAM_stats4 <- cbind("p.val.","training",sdm_data_obj_1@models[["occurrence"]][["gam"]][["2"]]@evaluation[["training"]]@statistics[["COR"]][["p.value"]])
  GAM_stats5 <- cbind("dev","training",sdm_data_obj_1@models[["occurrence"]][["gam"]][["2"]]@evaluation[["training"]]@statistics[["Deviance"]])
  GAM_stats <- rbind("GAM_1", "training", GAM_stats1, GAM_stats2 ,GAM_stats3, GAM_stats4,GAM_stats5)
  
}