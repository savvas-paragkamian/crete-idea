#!/usr/bin/Rscript
## Script name: crete_redlist.R
##
## Aim: handle IUCN redlist assessments of Crete
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-07-14

library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(sf)

# Load data
crete_shp <- sf::st_read("data/crete/crete.shp") 
polygons <- read_sf("data/redlist_species_polygons/data_0.shp")

habitats <- read_delim("data/redlist_species_data/habitats.csv", delim=",")

# keep only the non marine taxa
non_marine_habitats <- habitats |> 
    filter(!grepl("marine|substrate", name, ignore.case = T))

simple_summary <- read_delim("data/redlist_species_data/simple_summary.csv", delim=",") |>
    filter(scientificName %in% unique(non_marine_habitats$scientificName))

threats <- read_delim("data/redlist_species_data/threats.csv", delim=",") |>
    filter(scientificName %in% unique(non_marine_habitats$scientificName))


# summary taxa
phyla_taxa <- simple_summary |>
    group_by(phylumName) |>
    summarise(species=n())

phyla_taxa_redlist <- simple_summary |>
    group_by(phylumName, redlistCategory) |>
    summarise(n=n(), .groups="keep") |>
    pivot_wider(names_from=redlistCategory, values_from=n) |>
    left_join(phyla_taxa)

write_delim(phyla_taxa_redlist, "results/iucn_phyla_taxa_redlist.tsv", delim="\t")

# summary threats

threats_summary <- threats |> 
    group_by(name) |> 
    summarise(n=n()) |>
    arrange(desc(n))

threats_summary_phyla <- threats |> 
    left_join(simple_summary) |>
    group_by(name,phylumName) |> 
    summarise(n=n(), .groups="keep") |>
    arrange(desc(n))

# maps 

g_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    coord_sf(crs="WGS84") +
    theme_bw()

g_iucn <- g_base +
    geom_sf(polygons, mapping=aes(fill=SCI_NAME),alpha=0.5, size=0.1, show.legend=F)+
    theme_bw()

ggsave("figures/iucn_taxa_map.png",
       g_iucn,
       height = 15, 
       width = 20,
       dpi = 600, 
       unit="cm",
       device="png")
