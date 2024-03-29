#!/usr/bin/Rscript
## Script name: get_wosis_crete.R
##
## Aim: Download the latest data from the ISRIC and WoSIS soil spatial service
## World Soil Information service
## Open Geospatial Consortium (OGC) Web Feature Service (WFS)
## Author: Savvas Paragkamian
##
## Date Created: 2024-03-29

library(sf)
library(terra)
library(tidyverse)
library(httr)
library(RColorBrewer)

crete_shp <- sf::st_read("data/crete/crete.shp") 

#Specify the web address of the “latest” version of WoSIS:

wosis_wfs <- "WFS:https://maps.isric.org/mapserv?map=/map/wosis_latest.map"
layers.info <- st_layers(wosis_wfs)

# The layer "ms:wosis_latest:wosis_latest_profiles" contains the site information
wosis_latest <- "ms:wosis_latest_profiles"

#wosis_all <- st_read(dsn = wosis_wfs,
#                     layer=wosis_latest)

wosis_gr <- st_read(dsn = wosis_wfs,
                     layer=wosis_latest,
                     query= "select * from wosis_latest_profiles where country_name = 'Greece'")

wosis_crete <- st_intersection(wosis_gr, crete_shp)

st_write(wosis_crete, "results/wosis_crete/wosis_crete.shp")

wosis_crete |> st_drop_geometry() |> write_delim("results/wosis_crete/wosis_crete_metadata.tsv", delim="\t")

# plot

g_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    coord_sf(crs="WGS84") +
    theme_bw()

g_wosi <- g_base +
    geom_point(wosis_crete,
               mapping=aes(x=longitude,
                           y=latitude,
                           color=dataset_code),
               size=2,
               alpha=0.7) +
    scale_color_manual(values = c("cyan4","burlywood4", "darkgoldenrod1")) +
    theme(legend.position = 'bottom',legend.title=element_blank())

ggsave(paste0("figures/map_wosi_crete_soil.png",sep=""),
       plot=g_wosi, 
       height = 20, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

