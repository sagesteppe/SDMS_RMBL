
library(tidyverse)
setwd('~/Documents/numberSppRmbl')

harmonizer <- read.csv('./data/processed/taxonomy_lookup_table.csv')
spp_modelled <- read.csv('./data/processed/spp2model.csv') %>% 
  mutate(name  = str_replace_all(scientificName, ' ', "_")) %>% 
  left_join(., harmonizer, by = c('name' = 'Original')) %>% 
  mutate(Harmonized = if_else(is.na(Harmonized), name, Harmonized)) 

# first we need to determine which species are at rmbl again, note after running the ML
# I went back and search through the Consortium of Southern Rockies herbaria and had 
# to add more species to the vascular plant checklist AGAIN. This and the spelling mismatches
# mean that the logistic regression for distance is truthfully inaccurate (the results probably
# much better), but that cake is baked. 

rmbl_spp <- read.csv('./data/processed/rmbl_wfo_synonyms.csv') %>%
  mutate(Node.1 = replace_na(Node.1, 'Vascular')) %>%
  mutate('PrtyFlrs' = case_when(
    Group != 'angiosperms' ~ 0,
    Group == 'angiosperms' & Order == 'Poales' ~ 0,
    TRUE ~ 1
  )) %>%
  filter(PrtyFlrs == 1) %>%
  mutate(name = str_replace_all(Species, " ", "_")) %>% 
  left_join(., harmonizer, by = c('name' = 'Original'))  %>% 
  mutate(Harmonized = if_else(is.na(Harmonized), name, Harmonized)) %>% 
  pull(Harmonized)

rmbl_spp <- c(rmbl_spp,  'Viola_renifolia', 
  'Erigeron_ursinus', 'Packera_werneriifolia', 'Packera_neomexicana', 'Artemisia_absinthium' ,
  'Artemisia_michauxiana', 'Boechera_grahamii', 'Boechera_spatifolia', 'Cirsium_clavatum',
  'Cirsium_clavatum', 'Delphinium_occidentale', 'Heterotheca_pumila', 'Mertensia_brevistyla',
  'Mertensia_franciscana', 'Populus_angustifolia', 'Potamogeton_gramineus', 'Potentilla_ovina',
  'Pseudognaphalium_macounii')

IDerror <- c('Ligusticum_apiifolium', 'Grindelia_inuloides', 'Taraxacum_alaskanum',
             'Cirsium_eriophorum', 'Saxifraga_platysepala', 'Trifolium_rusbyi',
             'Pedicularis_sylvatica', 'Myrcia_micropetala', 'Senecio_muirii', 
             'Acacia_delicatula', 'Erigeron_elatus','Phrynium_fasciculatum', 
             'Pseudognaphalium_viscosum', 'Rorippa_curvisiliqua', 'Rorippa_islandica', 
             'Callitriche_antarctica', 'Thelesperma_megapotamicum')  
# these are unequivocally # not at RMBL
synonymized <- c('Boechera_selbyi', 'Viola_nephrophylla', 'Gentiana_aquatica', 
                 'Valeriana_capitata', 'Polygonum_arenastrum', 'Pedicularis_scopulorum',
                 'Lupinus_caespitosus', 'Arabis_hirsuta', 'Ranunculus_reptans',
                 'Matricaria_matricarioides') # synonymized into a taxon left in the list. 
skipped <- c('Artemisia_tridentata', 'Populus_tremuloides') #  skipped due to size of datasets. 
typo <- c('Rorippa_curvisiliqua') # common outside range - typo of curvipes?

IDerror <- c(IDerror, synonymized, skipped)
rmbl_spp <- rmbl_spp[! rmbl_spp %in% c(IDerror)]
rmbl_sppV <- unique(rmbl_spp)

rm(synonymized, skipped, IDerror, rmbl_spp)

# a number of species not retrieved via the spatial search were modelled as test data, we need to
# remove them from the datasets 

spp_modelled <- read.csv('./data/processed/spp2model.csv') %>% 
  mutate(name  = str_replace_all(scientificName, ' ', "_")) %>% 
  filter(!name %in% c('Gymnocarpium_dryopteris', 'Cystopteris_fragilis',
                      'Cystopteris_fragilis', 'Cystopteris_reevesiana')) %>% 
  left_join(., harmonizer, by = c('name' = 'Original')) %>% 
  mutate(Harmonized = if_else(is.na(Harmonized), name, Harmonized)) %>% 
  select(name = Harmonized) 

spp_present <- read.csv('./data/processed/taxa_predicted_rmbl_stats_only.csv') %>% 
  filter(!name %in% c('Gymnocarpium_dryopteris', 'Cystopteris_fragilis',
                      'Cystopteris_fragilis', 'Cystopteris_reevesiana'))  %>% 
  left_join(., harmonizer, by = c('name' = 'Original')) %>% 
  mutate(Harmonized = if_else(is.na(Harmonized), name, Harmonized)) %>% 
  select(name = Harmonized, value, approach) %>% 
  mutate(Predicted = 1,
    Observed = if_else(name %in% rmbl_sppV, 1, 0) ) %>% 
  filter(name %in% spp_modelled$name)

spp_absent <- read.csv('./data/processed/taxa_not_predicted_rmbl.csv') %>% 
  left_join(., harmonizer, by = c('name' = 'Original')) %>% 
  mutate(Harmonized = if_else(is.na(Harmonized), name, Harmonized)) %>% 
  filter(!name %in% c('Cystopteris_reevesiana')) %>% 
  select(name = Harmonized, value, approach) %>% 
  mutate(Predicted = 0,
         Observed = if_else(name %in% rmbl_sppV, 1, 0) ) %>% 
  filter(name %in% spp_modelled$name)

data <- bind_rows(spp_absent, spp_present)

rmbl_sppV <- data.frame('Taxon' = rmbl_sppV)
miss <-anti_join(rmbl_sppV, data, by = c('Taxon' = 'name'))

lm <- data %>% filter(approach == 'LM') %>% 
  group_by(name) %>% 
  sample_n(1) %>% 
  ungroup() 

ml <- data %>% filter(approach == 'ML') 

rm(spp_absent, spp_present)

ml1 <- ml %>% mutate(across(where(is.numeric), ~ factor(.x, levels = c(0:1))))
lm1 <- lm %>% mutate(across(where(is.numeric), ~ factor(.x, levels = c(0:1))))

ml1 %>% filter(Predicted == 1 & Observed == 1) %>% nrow() # true presence

caret::confusionMatrix(ml1$Predicted, reference = ml1$Observed)
caret::confusionMatrix(lm1$Predicted, reference = lm1$Observed)

# We focused on maximizing specificity, the percentage of true negatives (i.e. species classified not present at the field site, which are truly missing) at the expense of other metrics. Given that the number of prospective species for a genomic database has already been greatly reduced by the spatial search, we focus on removing only the most unlikely species from the study area. 

# COMPARE TO JANES SPECIES RICHNESS FROM PLOTS

floral <- read.csv(
  '/home/sagesteppe/Documents/floral_observations/data/Bombus_flower_ranks_2015.csv') %>% 
  distinct(plant.species) %>% 
  mutate(plant.species = str_replace_all(plant.species, "\\.", "_")) %>% 
  left_join(., harmonizer, by = c('plant.species' = 'Original')) %>% 
  mutate(Harmonized = if_else(is.na(Harmonized), plant.species, Harmonized)) %>% 
  filter(!plant.species %in% c('Carex_sp', 'Salix_sp', 'Epilobium_sp'))

focal_taxa <- data %>% 
  select(name:Predicted) %>% 
  mutate(Observed = if_else(name %in% floral$Harmonized, 1, 0) )  %>% 
  filter(name %in% floral$Harmonized) %>% 
  mutate(across(.cols = c(Predicted, Observed), ~ factor(.x, levels = c(0:1))))

anti_join(floral, focal_taxa, by = c('Harmonized' = 'name'))

ftML <- filter(focal_taxa, approach == 'ML', value > 0.5) 
aj <- anti_join(floral, ftML, by = c('Harmonized' = 'name'))
ftML <- bind_rows(ftML, data.frame(name = aj$Harmonized, value = 0, approach = 'ML', 
                      Predicted = factor(0), Observed = factor(1))
)

ftLM <- filter(focal_taxa, approach == 'LM', value > 0.5)
aj <- anti_join(floral, ftLM, by = c('Harmonized' = 'name'))
ftLM <- bind_rows(ftLM, data.frame(name = aj$Harmonized, value = 0, approach = 'LM', 
                           Predicted = factor(0), Observed = factor(1))
)

anti_join( ftML, ftLM, by = 'name')

rm(aj, ftLM, ftML, floral, focal_taxa)

# now we have the proportion of species which are correct, but we need to add 
# to this the species not predicted to be present at the field site, make
# an evaluation table of it and share results. 
# also pipe these species into the script where the BIEN records are evaluated. 

modelled_taxa <- read.csv('./data/processed/modelled_taxa.csv') %>% 
  separate(modelled_taxa, into = c('genus', 'species', 'method')) %>% 
  count(method)

# Why were these taxa missed? 

rmbl_spp <- read.csv('./data/processed/rmbl_wfo_synonyms.csv') %>%
  mutate(Node.1 = replace_na(Node.1, 'Vascular')) %>% 
  mutate(name = str_replace_all(Species, " ", "_")) 
  
miss <- anti_join(rmbl_sppV, data, by = c('Taxon' = 'name')) %>% 
  left_join(., rmbl_spp %>% select(name, Family), by = c('Taxon' = 'name'))
  

filter(miss, Family == 'Orchidaceae') # 13 records 
13/554 * 100 # 2.3% of records. 

missNativity <- read.csv('./data/processed/missing_species-Nativity.csv')
filter(missNativity, Native == 0) # 41 # % of records.  

41/554 * 100


# suspect that R. idaeus is listed as non-native in database
# two orchids, six non-natives, 1 taxonomically conflicted,  - Sisyrnchium is missing !!! present as genus but missing species


# Stages of SDM evaluation - focus on three items - althought maybe change area of focus for two to be the whole field station??? but then that adversely impacts 3 which is what we actually care about... Just focus on three being the ACTUAL goal of this whole activity. 


# The quality of the models themselves in a typical evaluation. 



# The quality of the models themselves in predicting a large component of the alpha richness of a region

# - Of the 554 vascular plants with biotic pollination syndromes, the 493 ML ensembles accurately predicted the presence of 362 (65.3%), incorrectly predicted the presence of 64 (11.6%), incorrectly predicted 34 true presences (6.1%) as being absent, and correctly predicted the true absence of 33 (6.0%). The balanced accuracy of the ensembled models is 0.627 (Sensitivity = 0.340, Specificity 0.914), a P VALUE IS NOT REPORTED AS THE VALUES WERE MANUALLY PARSED INTO CLASSES BASED ON SUITABILITY PER UNIT AREA.  Of the Of the 554 vascular plants with biotic pollination syndromes, the 475 ML ensembles accurately predicted the presence of 286 (51.6%), incorrectly predicted the presence of 41 (14.3%), incorrectly predicted 93 true presences (16.8%) as being absent, and correctly predicted the true absence of 55 (9.9%). The balanced accuracy of the ensembled models is 0.664 (Sensitivity = 0.573, Specificity 0.754), a P VALUE IS NOT REPORTED AS THE VALUES WERE MANUALLY PARSED INTO CLASSES BASED ON SUITABILITY PER UNIT AREA. Of the 554 vascular plants with biotic pollination syndromes in the flora 13 (2.3%) were in the Orchid family and 41 (7.4%) are non-natives, both of which are restricted from the database. 
#  Of the 

# The quality of the models in predicting ecologically dominant plant taxa in a pollination network, 

 # - Of the 117 plant species identified to the species level across the spatial extents of all plots and duration of queen bee activity, the ML ensembles predicted the presence of 105 (89.7%) of them, and LM ensembles 102 (87.2). Of the missing species two (1.7%) are Orchids, six (5.1%) are non-native, and one (0.85%) is of contested taxonomic standing, all of which are restricted from the initial query database,

