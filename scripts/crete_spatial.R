#!/usr/bin/Rscript



library(sf)
library(tidyverse)


crete_shp <- sf::st_read("data/crete/crete.shp") 

metadata_long <- read_delim("results/ena_samples_attributes-crete.tsv", delim="\t") %>%
    mutate(VALUE=gsub("\\r(?!\\n)","", VALUE, perl=T)) %>%
    distinct(.)
# metadata to wide format

metadata_wide <- metadata_long %>% 
    dplyr::select(-c(UNITS)) %>%
    mutate(TAG=gsub(" ","_", TAG, perl=T)) %>%
    mutate(TAG=gsub("-","_", TAG, perl=T)) %>%
    mutate(TAG=gsub("\\(|\\)","", TAG, perl=T)) %>%
    filter(!(TAG %in% c("bio_material", "BioSampleModel"))) %>%
    pivot_wider(id_cols="file", names_from=TAG, 
                values_from=VALUE) %>%
    filter(ENA_STUDY!="ERP024063")

metadata_wide$elevation <- as.numeric(metadata_wide$geographic_location_elevation)
metadata_wide$geographic_location_latitude <- as.numeric(metadata_wide$geographic_location_latitude)
metadata_wide$latitude <- as.numeric(metadata_wide$latitude)
metadata_wide$geographic_location_longitude <- as.numeric(metadata_wide$geographic_location_longitude)
metadata_wide$longitude <- as.numeric(metadata_wide$longitude)
metadata_wide$lat_lon1 <- gsub(" N ",",",metadata_wide$lat_lon, perl=T)
metadata_wide$lat_lon1 <- gsub(" |E","",metadata_wide$lat_lon1, perl=T)

metadata_wide <- metadata_wide %>%
    separate(lat_lon1, into=c("lat","lon"), sep=",") %>%
    mutate(lat=as.numeric(lat), lon=as.numeric(lon))

metadata_wide$latitude <- coalesce(metadata_wide$lat, metadata_wide$geographic_location_latitude,metadata_wide$latitude)
metadata_wide$longitude <- coalesce(metadata_wide$lon, metadata_wide$geographic_location_longitude, metadata_wide$longitude)

metadata_wide_sf <- metadata_wide %>%
    dplyr::select(file,ENA_RUN,ENA_EXPERIMENT, ENA_STUDY,project_name, latitude, longitude,organism,environment_biome) %>%
    filter(!is.na(latitude)) %>%
    st_as_sf(coords=c("longitude", "latitude"),
             remove=F,
             crs="WGS84")

crete_terrestrial_metagenome <- st_intersection(metadata_wide_sf, crete_shp) %>%
    filter(grepl("metagenome",organism), !grepl("marine|sponge|seawater|sediment",organism))

length(unique(crete_terrestrial_metagenome$ENA_STUDY))
length(unique(crete_terrestrial_metagenome$project_name))

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(crete_terrestrial_metagenome$organism))


g_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    coord_sf(crs="WGS84") +
    theme_bw()

g_ena <- g_base +
    geom_point(crete_terrestrial_metagenome,
               mapping=aes(x=longitude,
                           y=latitude,
                           color=organism),
               size=2,
               alpha=0.7) +
    scale_color_manual(values = getPalette(colourCount)) +
    theme(legend.position = 'bottom',legend.title=element_blank())

ggsave(paste0("figures/map_ena_samples_crete.png",sep=""),
       plot=g_ena, 
       height = 20, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

