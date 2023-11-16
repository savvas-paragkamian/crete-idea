#!/usr/bin/Rscript

## Script name: get_species_info.R
##
## Aim: Retrieve and analyse species of Crete habitats
##
## Author: Savvas Paragkamian
##
## Date Created: 2023-11-16

library(sf)
library(tidyverse)
#library(readxl)
library(rredlist)
#library(taxize)
library(units)
#library(vegan)

# gbif download
# use this POLYGON((23.44952 34.65265,26.43362 34.65265,26.43362 35.7184,23.44952 35.7184,23.44952 34.65265))
# for Crete

# Red List IUCN

redlist_species <- read_delim("data/redlist_species_data/simple_summary.csv", delim=",")

redlist_points <- read_delim("data/redlist_species_spatial/points_data.csv", delim=",")

redlist_sf_point <- redlist_points %>% st_as_sf(coords=c("longitude", "latitude"),
                                                remove=F,
                                                crs="WGS84")

crete_shp <- sf::st_read("data/crete/crete.shp") 

redlist_polygons <- sf::st_read("data/redlist_species_polygons/data_0.shp")

sf_use_s2(FALSE)

redlist_polygons_valid <- st_make_valid(redlist_polygons)

redlist_polygons_crete <- st_crop(redlist_polygons_valid, xmin = 23.44952, ymin = 34.65265, xmax = 26.43362, ymax = 35.7184)

sf_use_s2(TRUE)

redlist_polygons_crete_a <- redlist_polygons_crete %>%
    left_join(redlist_species, by=c("SCI_NAME"="scientificName")) %>%
    mutate(area=st_area(.))

redlist_sf_point_crete <- st_intersection(redlist_sf_point,crete_shp) %>%
    left_join(redlist_species,
              by=c("sci_name"="scientificName"))

crete_redlist <- ggplot()+
    geom_sf(crete_shp, mapping=aes())+
    geom_sf(redlist_polygons_crete_a,
               mapping=aes(color=phylumName))+
    theme_bw()+
    facet_wrap(vars(redlistCategory))

ggsave(plot=crete_redlist,
       "figures/map_crete_redlist.png",
       height = 20,
       width=40,
       units="cm")

