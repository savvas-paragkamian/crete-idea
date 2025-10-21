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

bbox_sf <- st_as_sfc(st_bbox(crete_shp))

#polygons <- read_sf("data/redlist_species_polygons/data_0.shp")

points <- read_delim("data/redlist_species_data_e48f862b-1a85-4d8e-a0e1-09fa0747a5e1/points_data.csv", delim=",")


points_sf <- points |>
    st_as_sf(coords=c("dec_long", "dec_lat"),
             remove=F,
             crs="WGS84")

points_sf_c <- st_intersection(points_sf, bbox_sf)


## other data from search summary of IUCN redlist
habitats <- read_delim("data/redlist_species_data_2a1d7856-4956-47ca-8b8c-8c079b336929/habitats.csv", delim=",")

# keep only the non marine taxa
#non_marine_habitats <- habitats |> 
#    filter(!grepl("marine|substrate", name, ignore.case = T))

simple_summary <- read_delim("data/redlist_species_data_2a1d7856-4956-47ca-8b8c-8c079b336929/simple_summary.csv", delim=",") 


#|>
#    filter(scientificName %in% unique(non_marine_habitats$scientificName))

threats <- read_delim("data/redlist_species_data_2a1d7856-4956-47ca-8b8c-8c079b336929/threats.csv", delim=",") 

# merge points and simple summary

taxonomy <- simple_summary |>
    distinct(scientificName,phylumName)

points_sf_c_s <- points_sf_c |>
    left_join(taxonomy, by=c("sci_name"="scientificName"))

# summary
points_sf_c_s_s <- points_sf_c_s |>
    distinct(sci_name, dec_long,dec_lat)

# summary taxa
phyla_taxa <- simple_summary |>
    group_by(phylumName) |>
    summarise(species=n())

# summary redlist categories
phyla_taxa_redlist <- simple_summary |>
    group_by(phylumName, redlistCategory, scopes) |>
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

write_delim(threats_summary_phyla, "results/threats_summary_phyla.tsv", delim="\t")
# maps 

colors <- c(
  "TRACHEOPHYTA" = "#228B22",      # Forest green for plants
  "ARTHROPODA" = "#FF4500",     # Orange-red for animals
  "PORIFERA" = "#8A2BE2",    # Blue-violet for chromists
  "MOLLUSCA" = "#8B4513",        # Saddle brown for fungi
  "BASIDIOMYCOTA" = "#FFD700"     # Gold for protozoa
)

g_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    coord_sf(crs="WGS84") +
    theme_bw() +
    theme(axis.title=element_blank(),
          legend.position = 'bottom',
          legend.title=element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          line = element_blank(),
          axis.text=element_blank(),
    )


g_iucn <- g_base +
    geom_sf(points_sf_c_s,
            mapping=aes(color=phylumName),
            alpha=0.5,
            size=0.1,
            show.legend=T)+
    scale_color_manual(values = colors) +
    theme_bw()+
    theme(legend.position = 'bottom')

ggsave("figures/iucn_taxa_map.png",
       g_iucn,
       height = 15, 
       width = 20,
       dpi = 600, 
       unit="cm",
       device="png")
