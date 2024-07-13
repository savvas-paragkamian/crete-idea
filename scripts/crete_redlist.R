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

# Load data

habitats <- read_delim("data/redlist_species_data/habitats.csv", delim=",")

# keep only the non marine taxa
non_marine_habitats <- habitats |> 
    filter(!grepl("marine|substrate", name, ignore.case = T))

simple_summary <- read_delim("data/redlist_species_data/simple_summary.csv", delim=",") |>
    filter(scientificName %in% unique(non_marine_habitats$scientificName))

# summary
phyla_taxa <- simple_summary |>
    group_by(phylumName) |>
    summarise(species=n())

phyla_taxa_redlist <- simple_summary |>
    group_by(phylumName, redlistCategory) |>
    summarise(n=n(), .groups="keep") |>
    pivot_wider(names_from=redlistCategory, values_from=n) |>
    left_join(phyla_taxa)

write_delim(phyla_taxa_redlist, "results/iucn_phyla_taxa_redlist.tsv", delim="\t")
