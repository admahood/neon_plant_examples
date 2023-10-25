library(tidyverse)
library(neonUtilities)
library(neonPlantEcology)
library(neonPlants)

# data get =====================================================================

allDiv <- neonPlantEcology::npe_download_plant_div()

# looking at it with neonPlants ================================================

##extract 1m2 data from list of lists and some processing
data_1m2 <- allDiv[["div_1m2Data"]]
###get 10_100
data_10_100m2 <- allDiv[["div_10m2Data100m2Data"]]

data_stacked <- divStack(
  data_1m2 = data_1m2,
  data_10_100m2 = data_10_100m2)

all_spp <- unique(data_stacked$scientificName) %>% sort()
total_spp <- all_spp %>% length()

data_stacked %>%
  filter(!is.na(scientificName)) %>%
  group_by(eventID) %>%
  summarise(length(unique(scientificName)))

data_stacked %>%
  filter(eventID == "SRER.1.2016") %>%
  pull(taxonID)

# looking at it with neonPlantEcology ==========================================

lf <- npe_longform(allDiv, timescale ="subannual")

all_spp_npe <- unique(lf$scientificName) %>% sort()
total_spp_npe <- length(all_spp_npe)

diy <- npe_diversity_info(allDiv, scale = "site", timescale = "subannual")
diy |> dplyr::select(eventID, nspp_total) |> arrange(eventID) |> print()

di <- npe_diversity_info(allDiv, scale = "site", timescale= "all")
di$nspp_total

# turn a divStack into a community matrix
cm <- npe_community_matrix(allDiv)

data_stacked %>%
  dplyr::select(plotID, subplotID, eventID, taxonID, targetTaxaPresent) %>%
  transmute(row = paste(plotID, subplotID, eventID, sep = "_"),
            taxonID = taxonID, 
            present = ifelse(targetTaxaPresent == "Y",1,0)) %>%
  pivot_wider(values_from = present, names_from = taxonID, values_fill = 0, 
              values_fn = function(x)sum(x)) %>%
  tibble::column_to_rownames("row")

