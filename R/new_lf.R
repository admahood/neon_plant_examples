npe_longform <- function(neon_div_object,
                         trace_cover=0.5,
                         scale = "plot",
                         divDataType = "plantSpecies",
                         timescale = c("subannual", "annual", "all"),
                         fix_unks = FALSE){
  .datatable.aware <- TRUE
  requireNamespace("data.table")
  requireNamespace("dplyr")
  requireNamespace("dtplyr")
  requireNamespace("tidyverse")
  requireNamespace("tidyr")
  requireNamespace("stringr")
  requireNamespace("magrittr")
  
  if(scale == "plot"){
    cover <- neon_div_object$div_1m2Data %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(endDate = as.Date(endDate)) %>%
      dplyr::filter(divDataType == divDataType) %>%
      tidyr::replace_na(list(percentCover=trace_cover)) %>%
      dplyr::filter(taxonID != "") %>%
      dplyr::group_by(plotID, taxonID, eventID) %>%
      dplyr::summarise(cover = sum(percentCover, na.rm=TRUE)/ifelse(as.numeric(stringr::str_sub(eventID,8,11))< 2019, 8,6),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) %>%
      dplyr::ungroup() %>%
      tibble::as_tibble()
    
    traces <- neon_div_object$div_10m2Data100m2Data %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(endDate = as.Date(endDate)) %>%
      dplyr::filter(targetTaxaPresent == "Y") %>%
      dplyr::group_by(plotID, subplotID, taxonID, eventID) %>%
      dplyr::summarise(cover = trace_cover,
                       scientificName = first(scientificName),
                       nativeStatusCode = first(nativeStatusCode),
                       family = first(family)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(taxonID != "") %>%
      dplyr::group_by(plotID, taxonID, eventID) %>%
      dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(stringr::str_sub(eventID,8,11))< 2019, 8,6),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) %>%
      dplyr::ungroup() %>%
      tibble::as_tibble()
    
    full_on_cover <- dplyr::bind_rows(cover, traces) %>%
      dplyr::group_by(plotID, taxonID, eventID, nativeStatusCode, scientificName, family) %>%
      dplyr::summarise(cover = sum(cover)) %>%
      dplyr::ungroup()%>%
      tidyr::replace_na(list(family = "Unknown")) %>%
      dplyr::mutate(site = stringr::str_sub(plotID, 1,4),
                    subplotID = "plot")
    
    if(fix_unks) full_on_cover <- full_on_cover %>%  unk_fixer()
    
    if(timescale == "all") {
      full_on_cover <- full_on_cover %>%
        mutate(year = as.numeric(stringr::str_sub(eventID,8,11)))
      year_range <- unique(full_on_cover$year)%>%
        as.numeric %>%
        range %>%
        paste(collapse = "-")
      n_years <- length(unique(full_on_cover$year))
      full_on_cover <- full_on_cover %>%
        dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                        family, site, subplotID) %>%
        dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(year = year_range)
    }
    if(timescale == "annual") {
      full_on_cover <- full_on_cover %>%
        mutate(year = as.numeric(stringr::str_sub(eventID,8,11)))
      full_on_cover <- full_on_cover %>%
        dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                        family, site, subplotID,year) %>%
        dplyr::summarise(cover = max(cover, na.rm=T)) %>%
        dplyr::ungroup() 
    }
    return(full_on_cover)
  }
  
  if(scale == "site"){
    cover <- neon_div_object$div_1m2Data %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(endDate = as.Date(endDate)) %>%
      dplyr::filter(divDataType == divDataType) %>%
      tidyr::replace_na(list(percentCover=trace_cover)) %>%
      dplyr::filter(taxonID != "") %>%
      dplyr::group_by(plotID, taxonID, eventID) %>%
      dplyr::summarise(cover = sum(percentCover, na.rm=TRUE)/ifelse(as.numeric(stringr::str_sub(eventID,8,11))< 2019, 8,6),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) %>%
      dplyr::ungroup() %>%
      tibble::as_tibble()
    
    traces <- neon_div_object$div_10m2Data100m2Data %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(endDate = as.Date(endDate)) %>%
      dplyr::filter(targetTaxaPresent == "Y") %>%
      dplyr::group_by(plotID, subplotID, taxonID, eventID) %>%
      dplyr::summarise(cover = trace_cover,
                       scientificName = first(scientificName),
                       nativeStatusCode = first(nativeStatusCode),
                       family = first(family)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(taxonID != "") %>%
      dplyr::group_by(plotID, taxonID, eventID) %>%
      dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(stringr::str_sub(eventID,8,11))< 2019, 8,6),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) %>%
      dplyr::ungroup() %>%
      tibble::as_tibble()
    
    n_plots <- length(unique(cover$plotID))
    
    full_on_cover <- dplyr::bind_rows(cover, traces) %>%
      dtplyr::lazy_dt() %>%
      dplyr::group_by(plotID, taxonID, eventID, nativeStatusCode, scientificName, family) %>%
      dplyr::summarise(cover = sum(cover)) %>%
      dplyr::ungroup()%>%
      dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) %>%
      dplyr::group_by(site, taxonID, eventID, nativeStatusCode, scientificName, family) %>%
      dplyr::summarise(cover = sum(cover)/n_plots) %>%
      dplyr::mutate(subplotID = "site",
                    plotID = "site") %>%
      dplyr::ungroup() %>%
      tidyr::replace_na(list(family = "Unknown")) %>%
      tibble::as_tibble()
    
    if(fix_unks) full_on_cover <- full_on_cover %>%  unk_fixer()
    if(timescale == "all") {
      full_on_cover <- full_on_cover %>%
        mutate(year = as.numeric(stringr::str_sub(eventID,8,11)))
      year_range <- unique(full_on_cover$year)%>%
        as.numeric %>%
        range %>%
        paste(collapse = "-")
      n_years <- length(unique(full_on_cover$year))
      full_on_cover <- full_on_cover %>%
        dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                        family, site, subplotID) %>%
        dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(year = year_range)
    }
    if(timescale == "annual") {
      full_on_cover <- full_on_cover %>%
        mutate(year = as.numeric(stringr::str_sub(eventID,8,11)))
      full_on_cover <- full_on_cover %>%
        dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                        family, site, subplotID,year) %>%
        dplyr::summarise(cover = max(cover, na.rm=T)) %>%
        dplyr::ungroup() 
    }
    return(full_on_cover)
  }
  
  # cover 8 ===========
  cover8 <- neon_div_object$div_1m2Data %>%
    dtplyr::lazy_dt() %>%
    dplyr::filter(divDataType == divDataType) %>%
    tidyr::replace_na(list(percentCover=trace_cover)) %>%
    dplyr::select(plotID, subplotID, taxonID, eventID, cover = percentCover,
                  nativeStatusCode, scientificName, family) %>%
    dplyr::filter(taxonID != "") %>%
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1, 4)) %>%
    tidyr::replace_na(list(family = "Unknown")) %>%
    tibble::as_tibble()
  
  
  # 10m2,100m2 are given 0.5 (we can change later)
  # unique(x$div_10m2Data100m2Data$subplotID) # there are 12 subplots
  
  # traces8 (10m2) ==============
  traces8 <- neon_div_object$div_10m2Data100m2Data %>% 
    dtplyr::lazy_dt() %>%
    # dplyr::filter(targetTaxaPresent == "Y") %>%
    dplyr::group_by(plotID, subplotID, taxonID, eventID, scientificName,
                    nativeStatusCode, family) %>%
    dplyr::summarise(cover = trace_cover) %>%
    dplyr::ungroup() %>%
    dplyr::filter(taxonID != "",
                  subplotID != "31", # these are the 100m2 subplots under which two 1m2 and 10m2 pairs are nested
                  subplotID != "32",
                  subplotID != "40",
                  subplotID != "41")  %>%
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1, 4)) %>%
    tidyr::replace_na(list(family = "Unknown")) %>%
    tibble::as_tibble()
  
  # traces100s ========
  traces100s <- neon_div_object$div_10m2Data100m2Data %>%
    dtplyr::lazy_dt() %>%
    dplyr::filter(targetTaxaPresent == "Y") %>%
    dplyr::group_by(plotID, subplotID, taxonID, eventID, scientificName,
                    nativeStatusCode, family) %>%
    dplyr::summarise(cover = trace_cover) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) %>%
    dplyr::filter(taxonID != "",
                  subplotID == "31"| # these are the 100m2 subplots under which two 1m2 and 10m2 pairs are nested
                    subplotID == "32"|
                    subplotID == "40"|
                    subplotID == "41") %>%
    tidyr::replace_na(list(family = "Unknown")) %>%
    tibble::as_tibble()
  
  # aggregating at different spatial scales ------------------------------------
  cover8_1m2 <- cover8 %>%
    dplyr::group_by(plotID, subplotID, taxonID, eventID, nativeStatusCode, scientificName, family) %>%
    dplyr::summarise(cover = sum(cover)) %>%
    dplyr::ungroup()%>%
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4))
  if(fix_unks) cover8_1m2 <- unk_fixer(cover8_1m2)
  
  cover8_1m2_10m2 <- dplyr::bind_rows(cover8, traces8) %>%
    dplyr::group_by(plotID,subplotID, taxonID, eventID, nativeStatusCode, scientificName, family) %>%
    dplyr::summarise(cover = sum(cover)) %>%
    dplyr::ungroup()%>%
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4))
  if(fix_unks) cover8_1m2_10m2<-cover8_1m2_10m2 %>%  unk_fixer()
  
  cover4 <- cover8_1m2_10m2 %>%
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1,2)) %>%
    dplyr::bind_rows(traces100s) %>% # adding in the 100m2 subplots
    dplyr::group_by(plotID, subplotID, eventID, taxonID) %>%
    dplyr::summarise(cover = sum(cover), # this is summing together repeats from the rbinding
                     scientificName = first(scientificName),
                     nativeStatusCode = first(nativeStatusCode),
                     family = first(family),
                     site = first(site)) %>%
    dplyr::ungroup()
  if(fix_unks) cover4 <- cover4 %>%  unk_fixer()
  
  
  if(scale == "1m") full_on_cover <- cover8_1m2
  if(scale == "10m") full_on_cover <- cover8_1m2_10m2
  if(scale == "100m") full_on_cover <- cover4
  
  if(timescale == "all") {
    full_on_cover <- full_on_cover %>%
      mutate(year = as.numeric(stringr::str_sub(eventID,8,11)))
    year_range <- unique(full_on_cover$year)%>%
      as.numeric %>%
      range %>%
      paste(collapse = "-")
    n_years <- length(unique(full_on_cover$year))
    full_on_cover <- full_on_cover %>%
      dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                      family, site, subplotID) %>%
      dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(year = year_range)
  }
  if(timescale == "annual") {
    full_on_cover <- full_on_cover %>%
      mutate(year = as.numeric(stringr::str_sub(eventID,8,11)))
    full_on_cover <- full_on_cover %>%
      dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                      family, site, subplotID,year) %>%
      dplyr::summarise(cover = max(cover, na.rm=T)) %>%
      dplyr::ungroup() 
  }
  
  return(full_on_cover)
}


npe_heights <- function(neon_div_object){
  neon_div_object$div_1m2Data %>%
    dtplyr::lazy_dt() %>%
    dplyr::mutate(endDate = as.Date(endDate)) %>%
    dplyr::filter(divDataType == "plantSpecies") %>%
    dplyr::filter(taxonID != "") %>%
    dplyr::select(plotID, subplotID, eventID, scientificName, taxonID, heightPlantSpecies, percentCover, heightPlantOver300cm) %>%
    tibble::as_tibble()
}


# old lonfgorm

# npe_longform <- function(neon_div_object,
#                                trace_cover=0.5,
#                                scale = "plot",
#                                divDataType = "plantSpecies",
#                                dissolve_years = FALSE,
#                                fix_unks = FALSE){
#   .datatable.aware <- TRUE
#   requireNamespace("data.table")
#   requireNamespace("dplyr")
#   requireNamespace("dtplyr")
#   requireNamespace("tidyverse")
#   requireNamespace("tidyr")
#   requireNamespace("stringr")
#   requireNamespace("magrittr")
# 
#   if(scale == "plot"){
#     cover <- neon_div_object$div_1m2Data %>%
#       dtplyr::lazy_dt() %>%
#       dplyr::mutate(endDate = as.Date(endDate)) %>%
#       dplyr::filter(divDataType == divDataType) %>%
#       dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4))) %>%
#       tidyr::replace_na(list(percentCover=trace_cover)) %>%
#       dplyr::group_by(plotID, subplotID, taxonID, year) %>%
#       # dealing with the multiple bout issue by first getting the max cover
#       # per sampling effort
#       dplyr::summarise(cover = max(percentCover),
#                 nativeStatusCode = first(nativeStatusCode),
#                 scientificName = first(scientificName),
#                 family = first(family)) %>%
#       dplyr::ungroup()  %>%
#       dplyr::filter(taxonID != "") %>%
#       dplyr::group_by(plotID, taxonID, year) %>%
#       dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(year)<2019, 8,6),
#                 nativeStatusCode = first(nativeStatusCode),
#                 scientificName = first(scientificName),
#                 family = first(family)) %>%
#       dplyr::ungroup() %>%
#       tibble::as_tibble()
# 
#     traces <- neon_div_object$div_10m2Data100m2Data %>%
#       dtplyr::lazy_dt() %>%
#       dplyr::mutate(endDate = as.Date(endDate)) %>%
#       dplyr::filter(targetTaxaPresent == "Y") %>%
#       dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
#       dplyr::group_by(plotID, subplotID, taxonID, year) %>%
#       dplyr::summarise(cover = trace_cover,
#                 scientificName = first(scientificName),
#                 nativeStatusCode = first(nativeStatusCode),
#                 family = first(family)) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(taxonID != "") %>%
#       dplyr::group_by(plotID, taxonID, year) %>%
#       dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(year)<2019, 12,10),
#                 nativeStatusCode = first(nativeStatusCode),
#                 scientificName = first(scientificName),
#                 family = first(family)) %>%
#       dplyr::ungroup() %>%
#       tibble::as_tibble()
# 
#     full_on_cover <- dplyr::bind_rows(cover, traces) %>%
#       dplyr::group_by(plotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
#       dplyr::summarise(cover = sum(cover)) %>%
#       dplyr::ungroup()%>%
#       dplyr::mutate(site = stringr::str_sub(plotID, 1,4),
#              subplotID = "plot")
#     if(fix_unks) full_on_cover <- full_on_cover %>%  unk_fixer()
# 
#     if(dissolve_years) {
#       year_range <- unique(full_on_cover$year)%>%
#         as.numeric %>%
#         range %>%
#         paste(collapse = "-")
#       n_years <- length(unique(full_on_cover$year))
#       full_on_cover <- full_on_cover %>%
#         dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName, family, site, subplotID) %>%
#         dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) %>%
#         dplyr::ungroup() %>%
#         dplyr::mutate(year = year_range)
#     }
# 
#     return(full_on_cover)
#   }
#   if(scale == "site"){
#     cover <- neon_div_object$div_1m2Data %>%
#       dtplyr::lazy_dt() %>%
#       dplyr::mutate(endDate = as.Date(endDate)) %>%
#       dplyr::filter(divDataType == divDataType) %>%
#       dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
#       tidyr::replace_na(list(percentCover=trace_cover)) %>%
#       dplyr::group_by(plotID, subplotID, taxonID, year) %>%
#       # dealing with the multiple bout issue by first getting the max cover
#       # per sampling effort
#       dplyr::summarise(cover = max(percentCover),
#                        nativeStatusCode = first(nativeStatusCode),
#                        scientificName = first(scientificName),
#                        family = first(family)) %>%
#       dplyr::ungroup()  %>%
#       dplyr::filter(taxonID != "") %>%
#       dplyr::group_by(plotID, taxonID, year) %>%
#       dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(year)<2019, 8,6),
#                        nativeStatusCode = first(nativeStatusCode),
#                        scientificName = first(scientificName),
#                        family = first(family)) %>%
#       dplyr::ungroup() %>%
#       tibble::as_tibble()
# 
#     traces <- neon_div_object$div_10m2Data100m2Data %>%
#       dtplyr::lazy_dt() %>%
#       dplyr::mutate(endDate = as.Date(endDate)) %>%
#       dplyr::filter(targetTaxaPresent == "Y") %>%
#       dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
#       dplyr::group_by(plotID, subplotID, taxonID, year) %>%
#       dplyr::summarise(cover = trace_cover,
#                        scientificName = first(scientificName),
#                        nativeStatusCode = first(nativeStatusCode),
#                        family = first(family)) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(taxonID != "") %>%
#       dplyr::group_by(plotID, taxonID, year) %>%
#       dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(year)<2019, 12,10),
#                        nativeStatusCode = first(nativeStatusCode),
#                        scientificName = first(scientificName),
#                        family = first(family)) %>%
#       dplyr::ungroup() %>%
#       tibble::as_tibble()
# 
#     n_plots <- length(unique(cover$plotID))
# 
#     full_on_cover <- dplyr::bind_rows(cover, traces) %>%
#       dplyr::group_by(plotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
#       dplyr::summarise(cover = sum(cover)) %>%
#       dplyr::ungroup()%>%
#       dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) %>%
#       dplyr::group_by(site, taxonID, year, nativeStatusCode, scientificName, family) %>%
#       dplyr::summarise(cover = sum(cover)/n_plots) %>%
#       dplyr::mutate(subplotID = "site",
#              plotID = "site") %>%
#       dplyr::ungroup()
#     if(fix_unks) full_on_cover <- full_on_cover %>%  unk_fixer()
#     if(dissolve_years) {
#       year_range <- unique(full_on_cover$year)%>%
#         as.numeric %>%
#         range %>%
#         paste(collapse = "-")
#       n_years <- length(unique(full_on_cover$year))
#       full_on_cover <- full_on_cover %>%
#         dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName, family, site, subplotID) %>%
#         dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) %>%
#         dplyr::ungroup() %>%
#         dplyr::mutate(year = year_range)
#     }
#     return(full_on_cover)
#   }
#   
#   # cover 8 ===========
#   cover8 <- neon_div_object$div_1m2Data %>%
#     dtplyr::lazy_dt() %>%
#     dplyr::mutate(endDate = as.Date(endDate)) %>%
#     dplyr::filter(divDataType == divDataType) %>%
#     dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
#     # entries in the df with no values but species was there
#     # i.e. someone put the sp. code and forgot to fill in the number
#     # putting as trace cover value
#     tidyr::replace_na(list(percentCover=trace_cover)) %>%
#     dplyr::mutate(endDate = as.Date(endDate)) %>%
#     dplyr::filter(divDataType == divDataType) %>%
#     dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
#     dplyr::group_by(plotID, subplotID, taxonID, year) %>%
#     # dealing with the multiple bout issue by first getting the mean cover
#     # per sampling effort, without aggregating, then later we'll aggregate.
#     # that way, a fall-bloomer that isn't visible in spring, for example,
#     # will be given its full cover value for fall, but then a species
#     # that is there for both seasons will be averaged, if that makes sense
#     dplyr::summarise(cover = max(percentCover),
#               nativeStatusCode = first(nativeStatusCode),
#               scientificName = first(scientificName),
#               family = first(family)) %>%
#     dplyr::ungroup()  %>%
#     dplyr::filter(taxonID != "") %>%
#     dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1, 4)) %>%
#     tibble::as_tibble()
# 
# 
#   # 10m2,100m2 are given 0.5 (we can change later)
#   # unique(x$div_10m2Data100m2Data$subplotID) # there are 12 subplots
#   
#   # traces8 ==============
#   traces8 <- neon_div_object$div_10m2Data100m2Data %>%
#     dtplyr::lazy_dt() %>%
#     dplyr::mutate(endDate = as.Date(endDate)) %>%
#     dplyr::filter(targetTaxaPresent == "Y") %>%
#     dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
#     dplyr::group_by(plotID, subplotID, taxonID, year) %>%
#     dplyr::summarise(cover = trace_cover,
#               scientificName = first(scientificName),
#               nativeStatusCode = first(nativeStatusCode),
#               family = first(family)) %>%
#     dplyr::ungroup() %>%
#     dplyr::filter(taxonID != "",
#            subplotID != "31", # these are the 100m2 subplots under which two 1m2 and 10m2 pairs are nested
#            subplotID != "32",
#            subplotID != "40",
#            subplotID != "41")  %>%
#     dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1, 4)) %>%
#     tibble::as_tibble()
#   
#   # traces100s ========
#   traces100s <- neon_div_object$div_10m2Data100m2Data %>%
#     dtplyr::lazy_dt() %>%
#     dplyr::mutate(endDate = as.Date(endDate)) %>%
#     dplyr::filter(targetTaxaPresent == "Y") %>%
#     dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
#     dplyr::group_by(plotID, subplotID, taxonID, year) %>%
#     dplyr::summarise(cover = trace_cover,
#               scientificName = first(scientificName),
#               nativeStatusCode = first(nativeStatusCode),
#               family = first(family)) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) %>%
#     dplyr::filter(taxonID != "",
#            subplotID == "31"| # these are the 100m2 subplots under which two 1m2 and 10m2 pairs are nested
#              subplotID == "32"|
#              subplotID == "40"|
#              subplotID == "41") %>%
#     tibble::as_tibble()
# 
#   # aggregating at different scales ----------------------------------------------
#   cover8_1m2 <- cover8 %>%
#     dplyr::group_by(plotID, subplotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
#     dplyr::summarise(cover = sum(cover)) %>%
#     dplyr::ungroup()%>%
#     dplyr::mutate(site = stringr::str_sub(plotID, 1,4))
#   if(fix_unks) cover8_1m2 <- unk_fixer(cover8_1m2)
# 
#   cover8_1m2_10m2 <- dplyr::bind_rows(cover8, traces8) %>%
#     dplyr::group_by(plotID,subplotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
#     dplyr::summarise(cover = sum(cover)) %>%
#     dplyr::ungroup()%>%
#     dplyr::mutate(site = stringr::str_sub(plotID, 1,4))
#   if(fix_unks) cover8_1m2_10m2<-cover8_1m2_10m2 %>%  unk_fixer()
# 
#   cover4 <- cover8_1m2_10m2 %>%
#     dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1,2)) %>%
#     dplyr::bind_rows(traces100s) %>% # adding in the 100m2 subplots
#     dplyr::group_by(plotID, subplotID, year, taxonID) %>%
#     dplyr::summarise(cover = sum(cover), # this is summing together repeats from the rbinding
#               scientificName = first(scientificName),
#               nativeStatusCode = first(nativeStatusCode),
#               family = first(family),
#               site = first(site)) %>%
#     dplyr::ungroup()
#   if(fix_unks) cover4 <- cover4 %>%  unk_fixer()
# 
# 
#   if(scale == "1m") full_on_cover <- cover8_1m2
#   if(scale == "10m") full_on_cover <- cover8_1m2_10m2
#   if(scale == "100m") full_on_cover <- cover4
# 
#   if(dissolve_years) {
#     year_range <- unique(full_on_cover$year)%>%
#       as.numeric %>%
#       range %>%
#       paste(collapse = "-")
#     n_years <- length(unique(full_on_cover$year))
#     full_on_cover <- full_on_cover %>%
#       dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName, family, site, subplotID) %>%
#       dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) %>%
#       dplyr::ungroup() %>%
#       dplyr::mutate(year = year_range)
#   }
# 
#   return(full_on_cover)
# }