#!/usr/bin/env Rscript

###############################################################################
# script name: search_jgi_gold.R
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to search JGI GOLD databse for biosamples
# in Crete.
#
###############################################################################
# OUTPUT:
#
###############################################################################
# usage:./isearch_jgi_gold.R
###############################################################################
library(sf)
library(tidyverse)

### read data
gold_Biosample <- read_excel("data/goldData.xlsx",sheet= "Biosample")

crete_shp <- sf::st_read("~/Documents/programming_projects/isd-crete/spatial_data/crete/crete.shp")

### modify the column names to remove spaces

names(gold_Biosample) = gsub(pattern = " ", replacement = "_", x = names(gold_Biosample))
gold_Biosample$LONGITUDE <- as.double(gold_Biosample$BIOSAMPLE_LONGITUDE)
gold_Biosample$LATITUDE <- as.double(gold_Biosample$BIOSAMPLE_LATITUDE)

### read as a sf object
gold_sf <- gold_Biosample |>
    filter(!is.na(LONGITUDE), LONGITUDE<181) |>
    st_as_sf(coords=c("LONGITUDE", "LATITUDE"),remove=F, crs="WGS84")

### global GOLD biosamples locations
gold_map<- ggplot()+ geom_point(gold_sf,mapping=aes(x=LONGITUDE,y=LATITUDE))
ggsave(plot=gold_map,filename= "figures/gold_map.png")

### terrestial Crete GOLD biosamples
gold_sf_crete <- st_intersection(gold_sf, crete_shp)

gold_map_crete <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_point(gold_sf_crete,mapping=aes(x=LONGITUDE,y=LATITUDE))

ggsave(plot=gold_map_crete,filename= "figures/gold_map_crete.png")

### text search for Crete samples, returns many false positives
gold_Biosample_crete <- gold_Biosample[grep("Crete", gold_Biosample$BIOSAMPLE_GEOGRAPHIC_LOCATION),]

