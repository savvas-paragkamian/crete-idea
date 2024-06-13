#!/usr/bin/Rscript
## Script name: crete_spatial.R
##
## Aim: handle spatial data of Crete
##
## Author: Savvas Paragkamian
##
## Date Created: 2023-11-16

library(sf)
library(terra)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

###################################### functions ############################

## create dir
create_dir <- function(dir_name){

    if (!dir.exists(dir_name)){
    dir.create(dir_name)
    print(paste0(dir_name, " directory created", sep=""))
    } else{
        print(paste0(dir_name, " directory exists", sep=""))
    }
}
convert_nested_l_df <- function(list_of_l_of_df){

    # transormation of nested list of lists to list of dataframes

    # this step merges lists of lists of dataframes to lists of dataframes.
    # The names of l1 are the concataned names of the original with . .
    # The next step is to create a new list with the right part of the name.
    # This way each element of the list contains dataframes of the same name.
    # In the last step the dataframes with the same name and subsequently the
    # same rows are merged. 
    # The final list is contains merged dataframes.
    l1 <- do.call('c', list_of_l_of_df)
    l2 <- split(l1,sub('.*\\.', '', names(l1)))
    results_df <- lapply(l2, function(x) Reduce(merge, x))
    return(results_df)


}
area_overlap_combination <- function(master_sf, list_sf, clabels){

    # this function takes as input a master sf, a named list of sf objects
    # and a labels factor. It calculates the area of overlap of the master_sf 
    # with each of the other sf objects per label. The output is a list of lists
    # of dataframes. 
    # Complementary functions are spatial_area_summary and convert_nested_l_df
    results <- list()

    for (i in seq_along(list_sf)){
    
        intersected <- st_intersection(master_sf,list_sf[[i]] )
        name_sf <- names(list_sf)[i]
        results[[name_sf]] <- list()
        
        for (j in seq_along(clabels)){

            label <- clabels[j]
            df_results <- spatial_area_summary(intersected,label)
            colnames(df_results) <- c(label, name_sf)
            print(j,i)

            results[[name_sf]][[label]] <- df_results
        }
    }
    return(results)
}

spatial_area_summary <- function(spatial_area,attribute){

    #this function takes a sf object and a grouping variable
    #and returns a dataframe of the area per name of the variable
    attribute <- as.character(attribute)
    spatial_area$area <- units::set_units(st_area(st_make_valid(spatial_area)),km^2)

    spatial_area_s <- spatial_area %>%
        st_drop_geometry() %>%
        group_by(.data[[attribute]]) %>%
        summarise(total_area=sum(area))

    return(spatial_area_s)
}

# crete borders
crete_shp <- sf::st_read("data/crete/crete.shp") 

############################### samplings #####################################
ena_crete_post <- read_delim("data/ena_post_crete.tsv", delim="\t")

arms_samples <- ena_crete_post |>
    filter(grepl("ARMS", sample_title)) |>
    st_as_sf(coords=c("lon", "lat"),
             remove=F,
             crs="WGS84")


metadata_long <- read_delim("results/ena_crete_long.tsv", delim="\t") %>%
    mutate(VALUE=gsub("\\r(?!\\n)","", VALUE, perl=T)) %>%
    distinct(.)
# metadata to wide format and column name tidy

metadata_wide <- metadata_long %>% 
    dplyr::select(-c(UNITS)) %>%
    mutate(TAG=gsub(" ","_", TAG, perl=T)) %>%
    mutate(TAG=gsub("-","_", TAG, perl=T)) %>%
    mutate(TAG=gsub("\\(|\\)","", TAG, perl=T)) %>%
    filter(!(TAG %in% c("ENA_FIRST_PUBLIC",
                        "ENA_LAST_UPDATE",
                        "bio_material",
                        "BioSampleModel"))) %>% # these contain duplicates
    pivot_wider(id_cols="file", names_from=TAG, 
                values_from=VALUE) %>%
    filter(ENA_STUDY!="ERP024063") # this is the ISD Crete project


metadata_wide$elevation <- as.numeric(metadata_wide$geographic_location_elevation)
metadata_wide$geographic_location_latitude <- as.numeric(metadata_wide$geographic_location_latitude)
metadata_wide$latitude <- as.numeric(metadata_wide$latitude)
metadata_wide$geographic_location_longitude <- as.numeric(metadata_wide$geographic_location_longitude)
metadata_wide$longitude <- as.numeric(metadata_wide$longitude)

# If locations are not in geographic_location_* they might be in lat_lon.
# Here we tranform lat lon to readable format
metadata_wide$lat_lon1 <- gsub(" N ",",",metadata_wide$lat_lon, perl=T)
metadata_wide$lat_lon1 <- gsub(" |E","",metadata_wide$lat_lon1, perl=T)

metadata_wide <- metadata_wide %>%
    separate(lat_lon1, into=c("lat","lon"), sep=",") %>%
    mutate(lat=as.numeric(lat), lon=as.numeric(lon))

# combine all the complementary columns
metadata_wide$latitude <- coalesce(metadata_wide$lat, metadata_wide$geographic_location_latitude,metadata_wide$latitude)
metadata_wide$longitude <- coalesce(metadata_wide$lon, metadata_wide$geographic_location_longitude, metadata_wide$longitude)

metadata_wide_sf <- metadata_wide %>%
    dplyr::select(file,ENA_RUN,ENA_EXPERIMENT, ENA_STUDY,environment_material, project_name, latitude, longitude,organism,environment_biome) %>%
    filter(!is.na(latitude)) %>%
    st_as_sf(coords=c("longitude", "latitude"),
             remove=F,
             crs="WGS84")

ena_sample_summary <- metadata_wide |>
    distinct(ENA_RUN,ENA_STUDY, organism, environment_feature, environment_biome, GOLD_Ecosystem_Classification) |>
    group_by(ENA_STUDY, organism) |>
    summarise(n=n(),.groups="keep") |>
    group_by(organism) |>
    summarise(studies=n(), samples=sum(n)) |>
    arrange(organism)

ena_sample_summary$organism <- factor(ena_sample_summary$organism, levels = c(sort(unique(ena_sample_summary$organism), decreasing = TRUE)))

ena_samples_bar <- ggplot() + 
    geom_col(ena_sample_summary, mapping=aes(x=samples, y=organism)) +
    theme_bw()+
    scale_x_continuous(breaks=seq(0,200,25))+
    xlab("Number of samples")+
    theme(axis.text.x = element_text(face="bold",
                                     size = 15),
          axis.title.x=element_text(face="bold",size=18),
          axis.ticks.y=element_blank(),axis.title.y=element_blank(),
          axis.text.y=element_text(size=10, hjust=0, vjust=0)) +
    theme(legend.position="bottom",
            panel.border = element_blank(),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank())

ggsave("figures/ena_samples_bar.png",
       plot=ena_samples_bar,
       device="png",
       height = 75,
       width = 38,
       units="cm")

# filter the terrestrial metagenomes only
crete_terrestrial_metagenome <- st_intersection(metadata_wide_sf, crete_shp) %>%
    filter(grepl("metagenome|bacterium",organism),
           #grepl("soil",environment_material),
           !grepl("marine|sponge|seawater|sediment",organism))


write_delim(crete_terrestrial_metagenome,"results/crete_terrestrial_metagenome.tsv", delim="\t")

length(unique(crete_terrestrial_metagenome$ENA_STUDY))
length(unique(crete_terrestrial_metagenome$project_name))

# Maps
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(crete_terrestrial_metagenome$organism))
colourCount_all = length(unique(metadata_wide_sf$organism))

g_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    coord_sf(crs="WGS84") +
    theme_bw()

ggsave(paste0("figures/map_crete_base.png",sep=""),
       plot=g_base, 
       height = 20, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")


g_arms <- g_base +
    geom_point(arms_samples,
               mapping=aes(x=lon,
                           y=lat),
               size=2,
               alpha=0.7) +
    theme(legend.position = 'bottom',legend.title=element_blank())

ggsave(paste0("figures/map_ena_arms.png",sep=""),
       plot=g_arms, 
       height = 20, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

g_ena <- g_base +
    geom_point(crete_terrestrial_metagenome,
               mapping=aes(x=longitude,
                           y=latitude,
                           color=organism),
               size=2,
               alpha=0.7) +
    scale_color_manual(values = getPalette(colourCount)) +
    theme(legend.position = 'bottom',legend.title=element_blank())

ggsave(paste0("figures/map_ena_terrestrial_metagenome_samples_crete.png",sep=""),
       plot=g_ena, 
       height = 20, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

g_ena_all <- g_base +
    geom_point(metadata_wide_sf,
               mapping=aes(x=longitude,
                           y=latitude,
                           color=organism),
               size=2,
               alpha=0.7) +
    scale_color_manual(values = getPalette(colourCount_all)) +
    theme(legend.position = 'bottom',legend.title=element_blank())

ggsave(paste0("figures/map_ena_all_samples_crete.png",sep=""),
       plot=g_ena_all, 
       height = 20, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

# EDAPHOBASE
edaphobase <- read_delim("data/20240603_edaphobase.csv", delim=";")

edaphobase$latitude <- edaphobase$`Latitude (Calculated single value)`
edaphobase$longitude <- edaphobase$`Longitude (Calculated single value)`

edaphobase_sf <- edaphobase |>
    filter(latitude!="field restricted",!is.na(longitude), !is.na(latitude)) |>
    mutate(latitude=as.numeric(latitude)) |>
    mutate(longitude=as.numeric(longitude)) |>
    distinct() |>
    st_as_sf(coords=c("longitude", "latitude"),
             remove=F,
             crs="WGS84")

edaphobase_crete_sf <- st_intersection(edaphobase_sf, crete_shp)

write_delim(st_drop_geometry(edaphobase_crete_sf),"results/edaphobase_crete.tsv", delim="\t")

edaphobase_crete_sf_samplings <- edaphobase_crete_sf |>
    distinct(Data.source,
             longitude,
             latitude) |>
    st_drop_geometry() |>
    mutate(reference=strtrim(Data.source,30))

colourCount = length(unique(edaphobase_crete_sf_samplings$reference))

g_edaphobase <- g_base +
    geom_point(edaphobase_crete_sf_samplings,
               position=position_jitter(h=0.015,w=0.015),
               mapping=aes(x=longitude,
                           y=latitude,
                           color=reference),
               size=2,
               alpha=0.5) +
    scale_color_manual(values = getPalette(colourCount)) +
    theme(legend.position = 'bottom',
          legend.text=element_text(size=5),
          legend.title=element_blank()) +
    guides(colour = guide_legend(ncol = 4))

#g_edaphobase <- g_edaphobase + facet_wrap(~Data.source)
ggsave(paste0("figures/map_edophobase_crete.png",sep=""),
       plot=g_edaphobase, 
       height = 15, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

############################### maps #####################################

crete_shp <- sf::st_read("data/crete/crete.shp")

clc_crete_shp <- st_read("data/clc_crete_shp/clc_crete_shp.shp")

clc_crete_shp$area <- units::set_units(st_area(clc_crete_shp),km^2)

clc_crete_shp_area <- clc_crete_shp |>
    st_drop_geometry() |>
    group_by(LABEL2) |>
    summarise(area=sum(area))

write_delim(clc_crete_shp_area, "results/clc_crete_label2_area.tsv", delim="\t")

### geology

crete_geology <- st_read("data/crete_geology/crete_geology.shp") |> st_make_valid()

crete_geology$area <- units::set_units(st_area(crete_geology),km^2)

crete_geology_area <- crete_geology |>
    st_drop_geometry() |>
    group_by(geology_na) |>
    summarise(area=sum(area))


write_delim(crete_geology_area, "results/crete_geology_area.tsv", delim="\t")

### protected areas
natura_crete <- sf::st_read("data/natura2000/natura2000_crete.shp")
wdpa_crete <- sf::st_read("data/wdpa_crete/wdpa_crete.shp") |>
    mutate(DESIG_ENG = gsub("Wildlife Refugee", "Wildlife Refuge", DESIG_ENG))

natura_crete_land <- st_intersection(natura_crete, crete_shp)
natura_crete_land_sci <- natura_crete_land %>% filter(SITETYPE=="B")
# raster DEM hangling
dem_crete <- rast("data/dem_crete/dem_crete.tif")
dem_crete_df <- as.data.frame(dem_crete, cells=T) |>
    filter(dem_crete>0) |>
    mutate(elevation=cut(dem_crete,
                             breaks=seq.int(from=0, to=2600, by=200),
                             dig.lab = 5 ))

dem_crete_area <- dem_crete_df |>
    group_by(elevation) |>
    summarise(area=0.0006*n())

write_delim(dem_crete_area, "results/dem_crete_area.tsv", delim="\t")
# raster bioclim 1 Annual Mean Temperature
bioclim1_crete <- rast("data/world_clim_crete/crete_wc2.1_30s_bio_1.tif")
bioclim1_crete_df <- as.data.frame(bioclim1_crete) 

# raster bioclim 12 Annual Mean Precipitation
bioclim12_crete <- rast("data/world_clim_crete/crete_wc2.1_30s_bio_12.tif")
bioclim12_crete_df <- as.data.frame(bioclim12_crete) 

# raster global aridity index
aridity_crete <- rast("data/crete_aridity_index.tif")
aridity_crete[aridity_crete[] == 0 ] = NA
aridity_crete_df <- as.data.frame(aridity_crete, cells=T)
aridity_crete_df$aridity <- aridity_crete_df$awi_pm_sr_yr*0.0001
aridity_crete_df$aridity_class <- cut(aridity_crete_df$aridity,
                                      breaks=c(0,0.03,0.2,0.5, 0.65,0.9),
                                      labels=c("Hyper Arid", "Arid", "Semi-Arid", "Dry sub-humid", "Humid"))

aridity_crete_area <- aridity_crete_df |>
    group_by(aridity_class) |>
    summarise(area=0.68*n())

write_delim(aridity_crete_area, "results/aridity_crete_area.tsv", delim="\t")
# raster desertification and Environmental Sensitive Areas Index
desertification_crete <- rast("data/crete_desertification_risk/esa3rdp_crete.tif")
desertification_crete_cat <- read_delim("data/crete_desertification_risk/esa3rdp_crete.tsv", delim="\t")
desertification_crete_df <- as.data.frame(desertification_crete,xy=T, cells=T)

desertification_crete_area <- desertification_crete_df |>
    group_by(ESA_12CL) |>
    summarise(area=4.3*n())

write_delim(desertification_crete_area, "results/desertification_crete_area.tsv", delim="\t")
# harmonised world soil database v2

hwsd2 <- rast("data/hwsd2_crete/hwsd2_crete.tif")
hwsd2[hwsd2[] == 7001 ] = NA
hwsd2_df <- as.data.frame(hwsd2,xy=T, cells=T) 
# hswd metadata
# with trimws the leading spaces are removed for the values.
HWSD2_wrb4 <- read_delim("data/hwsd2_crete/HWSD2_D_WRB4.tsv", delim="\t") |>
    mutate(VALUE=trimws(VALUE)) |>
    distinct(VALUE, CODE) 

HWSD2_SMU <- read_delim("data/hwsd2_crete/HWSD2_SMU.tsv", delim="\t") |>
    distinct(HWSD2_SMU_ID, WRB4) |>
    left_join(HWSD2_wrb4, by=c("WRB4"="CODE"))

hwsd2_df <- hwsd2_df |> 
    left_join(HWSD2_SMU, by=c("HWSD2"="HWSD2_SMU_ID")) 

hwsd2_area <- hwsd2_df |>
    group_by(VALUE) |>
    summarise(area=0.68*n())

write_delim(hwsd2_area, "results/hwsd2_area.tsv", delim="\t")

### protected areas
wdpa_crete$area <- units::set_units(st_area(wdpa_crete),km^2)

wdpa_crete_all <- data.frame(name="total protected", 
                             area=sum(wdpa_crete$area))

wdpa_crete_combine <- st_union(wdpa_crete) %>%
    st_make_valid() %>%
    st_as_sf() %>%
    filter(st_geometry_type(.) %in% c("MULTIPOLYGON"))

wdpa_crete_combine_area <- data.frame(name="total protected (no overlap)",
                                      area=sum(units::set_units(st_area(wdpa_crete_combine), km^2)))
crete_area <- data.frame(name="crete",
                         area=sum(units::set_units(st_area(crete_shp), km^2)))

protected_area <- wdpa_crete |> 
    group_by(DESIG_ENG) |>
    summarise(area=sum(area)) |>
    st_drop_geometry() |>
    dplyr::rename("name"="DESIG_ENG") |>
    bind_rows(crete_area,wdpa_crete_all, wdpa_crete_combine_area) |>
    arrange(area) |>
    mutate(area=round(area,2))

write_delim(protected_area, "results/protected_areas.tsv", delim="\t")

g_wdpa <- g_base +
    geom_sf(wdpa_crete, mapping=aes(fill=DESIG_ENG),alpha=0.5, size=0.1)+
    theme_bw()+
    theme(legend.position="bottom", legend.margin=margin())+
    guides(fill=guide_legend(nrow=5,byrow=TRUE, title="")) +
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank(),
          legend.key.size = unit(3, "mm"), 
          legend.text=element_text(size=7))

ggsave("figures/wdpa_protected_aread.png",
       g_wdpa,
       height = 15, 
       width = 20,
       dpi = 600, 
       unit="cm",
       device="png")

################ Global land use change hildap_GLOB-v1.0 #####################
### land cover change HILDA+
## crop data to crete
#path_hilda <- "/Users/talos/Documents/spatial_data/hildap_vGLOB-1.0_geotiff_wgs84/hildap_GLOB-v1.0_lulc-states/"
#output_directory <- "/Users/talos/Documents/spatial_data/hildap_vGLOB-1.0_geotiff_wgs84/hildap_GLOB-v1.0_lulc-states_crete/"
#hilda_files <- list.files(path_hilda)
#
#crete_bbox_polygon <- st_as_sf(st_as_sfc(st_bbox(crete_shp)))
#
#for (f in hilda_files) {
#    
#    if (grepl("*.tif$", f)) {
#        
#        #read_raster
#        path_raster <- paste0(path_hilda,f,sep="")
#        raster_tmp <- rast(path_raster)
#        
#        crete_raster <- terra::crop(raster_tmp, crete_bbox_polygon)
#        output_raster <- paste0(output_directory, "crete_",f,sep="")
#        terra::writeRaster(crete_raster, output_raster,overwrite=TRUE)
#
#        rm(path_raster,raster_tmp,crete_raster,output_raster)
#
#    }else{
#        
#        print(f, " not a tif")
#        next
#    }
#}
## Crete land use change
hilda_cat <- data.frame(hilda_id = c("11","22","33","44","55","66","77"),
                        hilda_name=c("urban","cropland","pasture/rangeland",
                                     "forest", "unmanaged grass/shrubland","sparse/no vegetation", "water"),
                        hilda_hex=c("#000000","#AE6120","#98BA6A","#07A07D","#BE81A3","#999999", "#1370A1"))
hilda_cat_v <- c("urban"="#000000",
                 "cropland"="#AE6120",
                 "pasture/rangeland"="#98BA6A",
                 "forest"="#07A07D",
                 "unmanaged grass/shrubland"="#BE81A3",
                 "sparse/no vegetation"="#999999",
                 "water"="#1370A1")


#### Hilda analysis
hilda_path <- "data/hildap_GLOB-v1.0_lulc-states_crete/"
hilda_id_names <- read_delim(paste0(hilda_path, "hilda_transitions_names.tsv", sep=""), delim="\t")
hilda_all <- list.files(hilda_path)
hilda_files <- hilda_all[grepl("*.tif", hilda_all)] 

create_dir("figures/hilda_crete")

for (i in 1:length(hilda_files)){
    print(i)
    filename <- hilda_files[i]
    
    raster_file <- rast(paste0(hilda_path,hilda_files[i],sep=""))
    
    raster_name <- paste0("hilda_",gsub(".*([0-9]{4}).*", "\\1", filename),sep="")
    
    raster_df <- terra::as.data.frame(raster_file, xy=TRUE, cells=TRUE)

    raster_df <- raster_df |>
        mutate(hilda_id=as.character(raster_df[,4])) |>
        filter(raster_df[,4]>0) |>
        left_join(hilda_cat, by=c("hilda_id"="hilda_id"))
    
    raster_df$hilda_name <- factor(raster_df$hilda_name, levels=as.character(unique(sort(raster_df$hilda_name))))
    
    g_hilda_map <- g_base +
        geom_raster(raster_df,
                    mapping=aes(x=x, y=y, fill=hilda_name)) +
        scale_fill_manual(values=hilda_cat_v) +
        guides(fill = guide_legend(nrow = 1)) +
        ggtitle(raster_name)+
        theme(axis.title=element_blank(),
              legend.position="bottom",
              legend.key.size = unit(4, "mm"), 
              legend.text=element_text(size=8),
              legend.title=element_blank())
    
#    ggsave(paste0("figures/hilda_crete/crete_",raster_name,"_map.png",sep=""),
#           plot=g_hilda_map,
#           height = 10, 
#           width = 20,
#           dpi = 300, 
#           units="cm",
#           device="png")
    
    hilda_sum <- zonal(cellSize(raster_file), raster_file, "sum") |> 
        mutate(area_m2=units::set_units(area,m^2)) |>
        mutate(area=units::set_units(area/10^6, km^2)) 
    
    hilda_sum <- hilda_sum |>
        mutate(hilda_id=as.character(hilda_sum[,1])) |>
        filter(hilda_sum[,1]>0) |>
        left_join(hilda_cat) 
    
    hilda_sum_g <- ggplot()+
        geom_col(hilda_sum,
                 mapping= aes(y=as.numeric(area),
                              x="",
                              fill = hilda_name),
                 position = position_stack(),
                 width = 0.2) +
        scale_fill_manual(values=hilda_cat_v) +
        scale_x_discrete(expand = expansion(add=c(0,0)))+
        scale_y_continuous(breaks=seq(0,9000,1000),
                           limits=c(0,8900),
                           expand = c(0,0))+
        ylab("Area sq. km") +
        xlab("") +
        theme_bw()+
        theme(legend.position='none',
              axis.ticks.x=element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(), #remove major gridlines
              panel.grid.minor = element_blank()) #remove minor gridlines
    
#    ggsave(paste0("figures/hilda_crete/crete_",raster_name,"_bar.png",sep=""),
#           plot=hilda_sum_g,
#           height = 10, 
#           width = 10,
#           dpi = 300, 
#           units="cm",
#           device="png")

    fig_hilda <- ggarrange(g_hilda_map,hilda_sum_g,
              labels = c("A", "B"),
              ncol = 2,
              nrow = 1,
              widths = c(0.85,0.15),
              font.label=list(color="black",size=15),
              common.legend = TRUE,
              legend="bottom") + bgcolor("white")
    
    ggsave(paste0("figures/hilda_crete/crete_",raster_name,".png",sep=""), 
           plot=fig_hilda, 
           height = 10, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")
}

## HILDA difference of land use, 1999-2019
##
hilda_1978 <- rast("data/hildap_GLOB-v1.0_lulc-states_crete/crete_hilda_plus_1978_states_GLOB-v1-0_wgs84-nn.tif")


hilda_2018 <- rast("data/hildap_GLOB-v1.0_lulc-states_crete/crete_hilda_plus_2018_states_GLOB-v1-0_wgs84-nn.tif")
## what is the transitions?
## the raster objects contain numeric values. The smart thing about this dataset
## is that I can use the numbers that are to show the category and create new
## numbers of the difference to symbolise the transitions.
## The oldest raster is the origin so transform it to the closest decade,
## 11 = 10, 22 = 20 etc. This can be accomplished with
## Modulus operation 44 - 44 %% 10
## Transform the latest raster to the units. so 44 = 4,
## Modulus operation 44 %% 10

hilda_1978_o <- app(hilda_1978, fun=function(i) i-i %% 10)
hilda_2018_o <- app(hilda_2018, fun=function(i) i %% 10)

hilda_1978_2018 <- hilda_1978_o + hilda_2018_o

hilda_1978_2018_sf <- st_as_sf(as.polygons(hilda_1978_2018,aggregate=F,values=T)) 

st_write(hilda_1978_2018_sf,
         "data/hilda_1978_2018/hilda_1978_2018.shp", 
         append=F,
         delete_layer=T,
         delete_dsn = TRUE) 



### Hilda summary
hilda_1978_2018_df <- terra::as.data.frame(hilda_1978_2018, xy=TRUE, cells=TRUE) |>
    filter(lyr.1>0)
#    left_join(hilda_cat, by=c("hilda_id"="hilda_id"))

hilda_sum <- zonal(cellSize(hilda_1978_2018), hilda_1978_2018, "sum") |> 
    mutate(crete_m2=units::set_units(area,m^2)) |>
    mutate(crete=units::set_units(area/10^6, km^2)) 

hilda_sum <- hilda_sum |>
    mutate(hilda_id_transition=hilda_sum[,1]) |>
    filter(hilda_sum[,1]>0) |> 
    left_join(hilda_id_names, by=c("hilda_id_transition"="hilda_id"))

write_delim(hilda_sum, "data/hilda_1978_2018/hilda_1978_2018.tsv", delim="\t")

### Natura2000 diffence of land use, 1998-2018
#natura_crete_land_hilda <- terra::mask(hilda_1998_2018, natura_crete_land)

#natura_crete_land_hilda_sum <- zonal(cellSize(natura_crete_land_hilda), natura_crete_land_hilda, "sum") |> 
#    mutate(natura2000_m2=units::set_units(area,m^2)) |>
#    mutate(natura2000=units::set_units(area/10^6, km^2)) |>
#    rename(hilda_id_transition=lyr.1) |>
#    filter(hilda_id_transition>0) |>
#    left_join(hilda_id_names, by=c("hilda_id_transition"="hilda_id")) |>
#    dplyr::select(-c(hilda_id_transition, area))

#write_delim(natura_crete_land_hilda_sum, "../results/natura_crete_land_hilda_sum.tsv", delim="\t")

################################### geology ################################
### the original data have EPSF:2001 and encoding = ISO 8859-7. I changed them to 
### EPSG:4326 and encoding UTF-8


## Desertification risk
#esa3rdp <- rast("data/ESDAC_CATENA_Desertification2018_RasterFiles/esa3rdp/w001001.adf")
#esa3rdp_wgs <- terra::project(esa3rdp, crs(metadata_aridity))

#esa3rdp_crete <- crop(esa3rdp_wgs, crete_shp)
#esa3rdp_attr <- data.frame(ESA=unique(esa3rdp_crete),
#                           ESA_id=seq(1,nrow(unique(esa3rdp_crete)),1))
#write_delim(esa3rdp_attr, "data/crete_desertification_risk/esa3rdp_crete.tsv", delim="\t")
#terra::writeRaster(esa3rdp_crete, "data/crete_desertification_risk/esa3rdp_crete.tif",overwrite=TRUE)
#esa3rdp_crete_r <- raster(esa3rdp_crete)

## Harmonised world soil database v2
#hwsd2 <- rast("data/hwsd2_crete/hwsd2_crete.tif")
#hwsd_r <- raster(hwsd2)
# hswd metadata
# with trimws the leading spaces are removed for the values.
#HWSD2_wrb4 <- read_delim("data/hwsd2_crete/HWSD2_D_WRB4.tsv", delim="\t") |>
#    mutate(VALUE=trimws(VALUE)) |>
#    distinct(VALUE, CODE) 

#HWSD2_SMU <- read_delim("data/hwsd2_crete/HWSD2_SMU.tsv", delim="\t") |>
#    distinct(HWSD2_SMU_ID, WRB4) |>
#    left_join(HWSD2_wrb4, by=c("WRB4"="CODE"))


## global_aridity_index-v3

#global_aridity_index <- rast("data/global_aridity_index-v3/Global-AI_ET0_v3_annual/ai_v3_yr.tif")
#crete_aridity <- terra::crop(global_aridity_index, crete_bbox_polygon)
#terra::writeRaster(crete_aridity, "results/crete_aridity_index.tif",overwrite=TRUE)


### Summary of overlaps of CLC with areas
list_sf <- list(natura2000 = natura_crete_land_sci,
                wildlife = wdpa_crete_wildlife)

clabels <- c("LABEL1", "LABEL2", "LABEL3")

### Here we calculate the overlap of areas with CLC
clc_area_ovelaps <- area_overlap_combination(clc_crete_shp, list_sf, clabels)

### data transformation
clc_area_ovelaps_df <- convert_nested_l_df(clc_area_ovelaps)

### Here we calculate the area of each LABEL for Crete
clc_crete_summary <- lapply(clabels, function(x) spatial_area_summary(clc_crete_shp, x))

### Here we merge the total area with the overlaps and 
### print the output

for (i in seq_along(clc_crete_summary)){
    
    merged <- clc_crete_summary[[i]] |>
        left_join(clc_area_ovelaps_df[[i]])

    write_delim(merged,
                paste0("results/clc_crete_", 
                       names(clc_area_ovelaps_df)[i],
                       ".tsv",sep=""),
                delim="\t")
}

clc_crete_label2 <- read_delim("results/clc_crete_LABEL2.tsv", delim="\t") |> 
    mutate(across(where(is.numeric), ~ round(.x,digits=2)))

