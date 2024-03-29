#!/usr/bin/Rscript
## Script name: get_wosis_crete.R
##
## Aim: Download the latest data from the ISRIC and WoSIS soil spatial service
## Open Geospatial Consortium (OGC) Web Feature Service (WFS)
## Author: Savvas Paragkamian
##
## Date Created: 2024-03-29

library(sf)
library(terra)
library(tidyverse)
library(httr)


#Specify the web address of the “latest” version of WoSIS:

wosis_wfs <- "WFS:https://maps.isric.org/mapserv?map=/map/wosis_latest.map"
layers.info <- st_layers(wosis_wfs)

# The layer "ms:wosis_latest:wosis_latest_profiles" contains the site information
wosis_latest <- "ms:wosis_latest_profiles"

#wosis_all <- st_read(dsn = wosis_wfs,
#                     layer=wosis_latest)

wosis_gr <- st_read(dsn = wosis_wfs,
                     layer=wosis_latest,
                     query= "select * from country_name = 'Gr'")

q = list(request = "GetCapabilities")
res$url

url <- parse_url(wosis_wfs)
url$query <- list(service = "wfs",
                  #version = "2.0.0", # facultative
                  request = "GetCapabilities")
request <- build_url(url)
request


test_gdal <- gdal_utils(
  util = "info",
  wosis_wfs
)
