###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): MW, AC, LK, WRJ
###goal(s): Wrangling raw CND data such that it is ready for analysis
###date(s): July 2024
###note(s): 

# Housekeeping ------------------------------------------------------------

### load necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, googledrive, vegan, readxl, e1071, dplyr, splitstackshape)

# ### set google drive paths
# exc_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/1VakpcnFVckAYNggv_zNyfDRfkcGTjZxX")) |> 
#       ### updated this file after hierarchical NA-filling for some sites, plus correcting PISCO biomass estimates (i.e., previously at transect, not m2 level)
#       ### renamed previous version as 'harmonized_consumer_excretion_CLEANV1.csv'
#       dplyr::filter(name %in% c("harmonized_consumer_excretion_CLEAN.csv"))
# 
# strata_ids <- googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/1/folders/1CEgNtAnk4DuPNpR3lJN9IqpjWq0cM8F4")) %>%
#   dplyr::filter(name %in% c("strata_class.xlsx"))
# 
# ### combine file IDs
# harmonized_ids <- rbind(exc_ids, strata_ids)
# 
# ### for each raw data file, download it into the consumer folder
# for(k in 1:nrow(harmonized_ids)){
#   
#   ### download file (but silence how chatty this function is)
#   googledrive::with_drive_quiet(
#     googledrive::drive_download(file = harmonized_ids[k, ]$id, overwrite = T,
#                                 path = file.path("tier2", harmonized_ids[k, ]$name)) )
#   
#   ### print success message
#   message("Downloaded file ", k, " of ", nrow(harmonized_ids))
# }
# 
# ### cleans environment
# rm(list = ls()) 

### read in clean excretion and strata data from google drive
dt <- read.csv(file.path("tier2", "harmonized_consumer_excretion_CLEAN.csv"),stringsAsFactors = F,na.strings =".") |> 
  janitor::clean_names()
glimpse(dt)

strata_list <- readxl::read_excel(path = file.path("tier2", "strata_class.xlsx"),na=".") |> 
  ### remove decimals from numbered sites
  mutate(site = str_remove(site, "\\.0$"),
         subsite_level1 = str_remove(subsite_level1, "\\.0$"),
         subsite_level2 = str_remove(subsite_level2, "\\.0$"),
         subsite_level3 = str_remove(subsite_level3, "\\.0$"))

# set up data for summary statistics --------------------------------------

### replace NAs in subsite_level2 and subsite_level3 columns with "Not Available"
### to allow group_by function to go as far in sequence as it can for 
### each project without throwing NAs

dt1 <- dt |> 
  mutate(subsite_level1 = replace_na(subsite_level1, "Not Available"),
         subsite_level2 = replace_na(subsite_level2, "Not Available"),
         subsite_level3 = replace_na(subsite_level3, "Not Available"))

### check to see NA fixes incorporated
na_count_per_column <- sapply(dt1, function(x) sum(is.na(x)))
print(na_count_per_column)

# look into outliers for project-species combinations -----------
### removing 'biomass buster' sharks based on requests from PISCO, SBC, and MCR 

dt_og <- dt1 |> 
  group_by(project, habitat) |> 
  ### filtering out sharks and rays that are considered "biomass busters"
  mutate(mean_dmperind = mean(dmperind_g_ind, na.rm = TRUE),  
         sd_dmperind = sd(dmperind_g_ind, na.rm = TRUE),  
         lower_bound = mean_dmperind - 5 * sd_dmperind,  
         upper_bound = mean_dmperind + 5 * sd_dmperind,
         ### +/- 5 SD cutoff... rest of sharks and rays included
         outlier = dmperind_g_ind < lower_bound | dmperind_g_ind > upper_bound,
         sharkray = grepl("\\bshark\\b|\\bray\\b", common_name, ignore.case = TRUE),
         elasmo = class %in% c("Chondrichthyes", "Elasmobranchii")) |> 
  ungroup() |> 
  filter(!(outlier & (sharkray | elasmo))) |> #lose 268 biomass bustin' sharks and rays
  dplyr::select(-lower_bound, -upper_bound, -outlier, -sharkray, -elasmo)
glimpse(dt_og)

### check on the california moray to see if max size makes sense
test <- dt_og |> group_by(project,scientific_name) |> summarize(max_size = max(dmperind_g_ind, na.rm = TRUE), mean_size = mean(dmperind_g_ind, na.rm = TRUE))
### max reported weight for california moray is 80 lbs, so approximate dry weight of 20 lbs (or 9071 g)
dt_og1 <- dt_og |> 
      mutate(dmperind_g_ind = case_when(dmperind_g_ind > 9071 & scientific_name == "Gymnothorax mordax" ~ 9071,
                                        TRUE ~ dmperind_g_ind))

###########################################################################
# add vertebrate and invertebrate column ~ phylum -------------------------
###########################################################################

dt_mutate <- dt_og1 |> 
  ### classify each individual as either being a vertebrate or invertebrate
  mutate(vert_1 = if_else(phylum == "Chordata", "vertebrate", "invertebrate")) |> 
  mutate(vert2 = if_else(is.na(vert_1) & project == "CoastalCA", "vertebrate", vert_1)) |> 
  mutate(vert = if_else(is.na(vert2), "invertebrate", vert2)) |> 
  mutate(vertebrate_n = if_else(vert == "vertebrate" & dmperind_g_ind != 0, 1, 0),
         invertebrate_n = if_else(vert == "invertebrate" & dmperind_g_ind != 0, 1, 0)) |> 
  ### filtering out invertebrates - first manuscript focused on vertebrates (ie fishes)
  filter(vert == "vertebrate",
         ### removing invertebrate dominant projects, plus PIE given conversations with DB, NL, AS
         !project %in% c("NGA", "CCE", "PIE")) |> 
  ### removing some of the unnecessary columns      
  select(-vert, -vert_1, -vert2, -vertebrate_n, -invertebrate_n, -raw_filename, -row_num,-mean_dmperind,-sd_dmperind)

rm(dt, test, dt_og, dt_og1, dt1)
###########################################################################
# set up such that each row is a single individual ------------------------
###########################################################################
### want to save those zeros
fce_area<-read_csv("local_data/fce-area.csv") |> 
      mutate(subsite_level1 = as.character(subsite_level1),
             subsite_level2 = as.character(subsite_level2))

m_zero_save <- dt_mutate |> filter(project == "FCE" & density_num_m == 0)
m2_zero_save <- dt_mutate |> filter(project != "FCE"& density_num_m2 == 0)

fce <- dt_mutate |> filter(project=="FCE")
fce_no_zeros <- dt_mutate |> filter(density_num_m != 0)

# fce_plus_zero_save <- dt_mutate |> filter(project == "FCE"|density_num_m2 == 0)
# fce1 <- fce_plus_zero_save |> filter(project == "FCE")
fce_distance <- fce |> 
      left_join(fce_area, by = c("year", "month", "site", 
                                 "subsite_level1", "subsite_level2")) |> 
      mutate(area = if_else(
        is.na(area),
        mean(area, na.rm = TRUE),
        area
      )) |> 
      select(-calendar_yr)
# map(fce_distance, ~sum(is.na(.)))

fce_no_zeros <- fce_distance |> filter(density_num_m != 0)

dt_mutate_fce <- fce_no_zeros |> 
      # mutate(count = ceiling(density_num_m*area)) |>
      mutate(count = round(density_num_m*area)) |>
      expandRows(count = "count", drop = FALSE) |> 
      mutate(density_num_m = 1/area) |> 
      ungroup() |> 
      dplyr::select(-count,-area)

### filter out sites and instances where animals observed here to expand rows based on counts below 
### in this way, each row will become an observation of a single individual
no_fce_or_zeros <- dt_mutate |> filter(project != "FCE" & density_num_m2 != 0)
      
dt_mutate_1 <- no_fce_or_zeros |> 
      mutate(area = case_when(
            project == "CoastalCA" ~ 60,
            project == "SBC" ~ 40,
            project == "VCR" ~ 25,
            project == "MCR" & subsite_level3 == "1" ~ 50,
            project == "MCR" & subsite_level3 == "5" ~ 250
      )) |> 
      group_by(project) |> 
      # mutate(count = ceiling(density_num_m2*area)) |>
      mutate(count = round(density_num_m2*area),
             count = if_else(
                   count == 0,
                   1,
                   count
             )) |> 
      filter(count <10000,
             nind_ug_hr != 0) |> 
      ungroup() |> 
      group_by(project) |> 
      expandRows(count = "count", drop = FALSE) |> 
      ungroup() |> 
      group_by(project) |> 
      mutate(density_num_m2 = 1/area) |> 
      ungroup() |> 
      dplyr::select(-count,-area)

dt_mutate_2 <- rbind(dt_mutate_1,dt_mutate_fce,m_zero_save,m2_zero_save)
# rm(dt_mutate,dt_mutate_1,no_fce_or_zeros)
##########################################################################
### coding with AC to get the max size of each species in the population
### unique to step3 - we are not doing this for population or trophic levels
dt_mutate_3 <- dt_mutate_2 |> 
      group_by(project, habitat, year, month, site, subsite_level1,
               subsite_level2, subsite_level3, scientific_name) |>
      mutate(max_size = case_when(dmperind_g_ind != 0 ~ max(dmperind_g_ind),
      T ~ NA)) |> 
      ungroup()
##########################################################################
### here is where steps four-five begin to differ from steps one-three ---
##########################################################################
### summarize data at "finest" scale for each individual program (e.g.,
### transect or bout) - LK updates at June meeting allow appropriate
### resolution for PISCO datasets

### check for NAs
na_count_per_column <- sapply(dt_mutate_3, function(x) sum(is.na(x)))
print(na_count_per_column) #yay

dt_total <- dt_mutate_3 |> 
  group_by(project, habitat, year, month, 
           site, subsite_level1, subsite_level2, subsite_level3) |> 
  summarize(
    ### calculate total nitrogen supply at each sampling unit and then sum to get column with all totals
    total_nitrogen_m = sum(nind_ug_hr * density_num_m, na.rm = TRUE),
    total_nitrogen_m2 = sum(nind_ug_hr * density_num_m2, na.rm = TRUE),
    # total_nitrogen_m3 = sum(nind_ug_hr * density_num_m3, na.rm = TRUE),
    ### create column with total_nitrogen contribution for each program, regardless of units
    total_nitrogen = sum(total_nitrogen_m + total_nitrogen_m2, na.rm = TRUE),
    ### calculate total phosphorus supply at each sampling unit and then sum to get column with all totals
    total_phosphorus_m = sum(pind_ug_hr * density_num_m, na.rm = TRUE),
    total_phosphorus_m2 = sum(pind_ug_hr * density_num_m2, na.rm = TRUE),
    # total_phosphorus_m3 = sum(pind_ug_hr * density_num_m3, na.rm = TRUE),
    ### create column with total_phosphorus contribution for each program, regardless of units
    total_phosphorus = sum(total_phosphorus_m + total_phosphorus_m2, na.rm = TRUE),
    ### calculate total biomass at each sampling unit and then sum to get column with all totals
    total_bm_m = sum(dmperind_g_ind*density_num_m, na.rm = TRUE),
    total_bm_m2 = sum(dmperind_g_ind*density_num_m2, na.rm = TRUE),
    # total_bm_m3 = sum(dmperind_g_ind*density_num_m3, na.rm = TRUE),
    ### create column with total_biomass for each program, regardless of units
    total_biomass = sum(total_bm_m + total_bm_m2, na.rm = TRUE),
    # ### calculate species richness
    n_spp = n_distinct(scientific_name[dmperind_g_ind != 0]),
    ### calculate average community size metrics
    max_size = mean(max_size, na.rm = TRUE),
    # mean_size = mean(mean_size, na.rm = TRUE),
    # min_size = mean(min_size, na.rm = TRUE),
    ### calculate diversity indices of interest
    Species_Richness = length(unique(scientific_name[dmperind_g_ind != 0])),
    # Species_Shannon_Diversity_Index = diversity(x = table(scientific_name), index = "shannon"),
    Species_Inverse_Simpson_Diversity_Index = diversity(x = table(scientific_name[dmperind_g_ind != 0]), index = "invsimpson"),
    Trophic_Richness = length(unique(diet_cat[dmperind_g_ind != 0])),
    # Trophic_Shannon_Diversity_Index = diversity(x = table(diet_cat), index = "shannon"),
    Trophic_Inverse_Simpson_Diversity_Index = diversity(x = table(diet_cat[dmperind_g_ind != 0]), index = "invsimpson"),
    size_skewness = skewness(dmperind_g_ind[dmperind_g_ind != 0])) |> 
  ungroup() |>
  dplyr::select(-total_nitrogen_m, -total_nitrogen_m2,
                -total_phosphorus_m, -total_phosphorus_m2,
                -total_bm_m, -total_bm_m2, -n_spp) |>
  arrange(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3)

na_count_per_column <- sapply(dt_total, function(x) sum(is.na(x)))
print(na_count_per_column) #yay
# test <- dt_total |> filter(is.na(size_skewness))
### cleaning up infinite diversity metric values for instances where no fish were caught
dt_total_clean <- dt_total |> 
      mutate(Species_Inverse_Simpson_Diversity_Index = ifelse(Species_Richness == 0, 0, 
                                                              Species_Inverse_Simpson_Diversity_Index),
             Trophic_Inverse_Simpson_Diversity_Index = ifelse(Trophic_Richness == 0, 0, 
                                                              Trophic_Inverse_Simpson_Diversity_Index),
             max_size = if_else(Species_Richness == 0, 0, 
                                max_size),
             size_skewness = if_else(Species_Richness == 0, 0,
                                     size_skewness),
             size_skewness = case_when(is.na(size_skewness)~0,
                                       TRUE~size_skewness)) |> 
      rename(species_richness = Species_Richness,
             species_diversity = Species_Inverse_Simpson_Diversity_Index, 
             trophic_diversity = Trophic_Inverse_Simpson_Diversity_Index,
             trophic_richness = Trophic_Richness)
glimpse(dt_total_clean)

### check for NAs
na_count_per_column <- sapply(dt_total_clean, function(x) sum(is.na(x)))
print(na_count_per_column) #yay

### checking revised dataset
dt_total_clean |> ggplot(aes(x = species_richness, y = species_diversity, color = project)) + geom_abline() + geom_point() + facet_wrap(~project, scales = 'free')

###########################################################################
# add strata of interest to each project ----------------------------------
###########################################################################

### set up strata_list such that the "Not Available" aligns with dt_total
strata_list1 <- strata_list %>%
  mutate(subsite_level1 = replace_na(subsite_level1, "Not Available"),
         subsite_level2 = replace_na(subsite_level2, "Not Available"),
         subsite_level3 = replace_na(subsite_level3, "Not Available")) |> 
  ### LK suggested to remove subsite_level 3 given updates to PISCO 
  select(-subsite_level3) |> 
  distinct()

### join together the datasets of nutrient supply and biomass with strata
dt_total_strata <- left_join(dt_total_clean, 
                             strata_list1, 
                             by = c("project", "habitat", "site",
                                    "subsite_level1", "subsite_level2")) |> 
  unite("projecthabitat", project, habitat, sep = "-", remove = FALSE) |> 
  rename(strata = ecoregion_habitat) |> 
  select(project, habitat, projecthabitat, strata, year, month, site, subsite_level1, 
          subsite_level2, subsite_level3, everything())
      
### Check NAs
na_count_per_column <- sapply(dt_total_strata, function(x) sum(is.na(x)))
print(na_count_per_column) #yayay
glimpse(dt_total_strata)

###########################################################################
# set up individual projects/habitats for analyses and plotting -----------
###########################################################################

### Below I have separated each unique projecthabitat out to mutate new columns based on either
# the strata they want their data colored by (i.e., color = XXXX)and the level to which they want
# their data summarized (i.e., for FCE-estuary, I want summarized at subsite_level1..
# whereas SBC wants their data summarized at the site level. This approach sets up
# an easy way to map plots across all unique projecthabitats, instead of doing them
# individually

### CoastalCA-ocean
pisco_central <- dt_total_strata |> 
  filter(projecthabitat == "CoastalCA-ocean",
         site == "CENTRAL") |> #split pisco into central and southern
  mutate(group = subsite_level2,
         color = strata,
         units = 'm2',
         projecthabitat = "CoastalCA-ocean-CENTRAL") |> 
  ### added new resolution group wants considered for examination -> functionally the "site" for each project
  unite(color2, c(subsite_level2, color), sep = "-", remove = FALSE)

pisco_south <- dt_total_strata |> 
  filter(projecthabitat == "CoastalCA-ocean",
         site == "SOUTH") |> #split pisco into central and southern
  mutate(group = subsite_level2,
         color = strata,
         units = 'm2',
         projecthabitat = "CoastalCA-ocean-SOUTH") |> 
  ### added new resolution group wants considered for examination -> functionally the "site" for each project
  unite(color2, c(subsite_level2, color), sep = "-", remove = FALSE)

### FCE-estuary
fce <- dt_total_strata |> 
  filter(projecthabitat == "FCE-estuary") |>
  mutate(group = subsite_level1,
         color = strata,
         units = 'm') |> #grouped at subsite_level1
  ### added new resolution group wants considered for examination -> functionally the "site" for each project
  unite(color2, c(site, subsite_level1), sep = "-", remove = FALSE)

### MCR-ocean
mcr <- dt_total_strata |> 
  filter(projecthabitat == "MCR-ocean") |> 
  ### join site and subsite_level1 according DB request for grouping variable
  unite("group", site, subsite_level1, sep = "-", remove = FALSE) |>
  mutate(group = group,
         color = subsite_level1,
         units = 'm2') |> 
  ### added new resolution group wants considered for examination -> functionally the "site" for each project
  unite(color2, c(subsite_level1, site), sep = "-", remove = FALSE)

### SBC-ocean
sbc <- dt_total_strata |> 
  filter(projecthabitat == "SBC-ocean") |> 
  mutate(group = site,
         color = strata,
         units = 'm2') |> 
  ### added new resolution group wants considered for examination -> functionally the "site" for each project
  unite(color2, c(site, color), sep = "-", remove = FALSE)

### VCR-estuary
vcr <- dt_total_strata |> 
  filter(projecthabitat == "VCR-estuary") |> 
  mutate(group = subsite_level1,
         color = strata,
         units = 'm2') |> 
  ### added new resolution group wants considered for examination -> functionally the "site" for each project
  unite(color2, c(site, subsite_level1, color), sep = "-", remove = FALSE)

### binding everything back together, removing index row generated when saving out of R
## and arranging the data by date
dat_ready <- bind_rows(fce, mcr, pisco_central, pisco_south, sbc, vcr)

na_count_per_column <- sapply(dat_ready, function(x) sum(is.na(x)))
print(na_count_per_column) #yay

### tidy up working environment
rm(fce, mcr, pisco_central, pisco_south, sbc, vcr)

###########################################################################
# clean up dataset names for plotting and analysis ------------------------
###########################################################################

unique(dat_ready$projecthabitat)
label_mapping <- data.frame(
  projecthabitat = unique(dat_ready$projecthabitat),
  Project = c("FCE", "MCR", "PCCC", "PCCS",
              "SBC", "VCR")) 
print(label_mapping) #looks good

unique(dat_ready$color)
habitat_mapping <- data.frame(
  color = unique(dat_ready$color),
  Habitat = c(
              "Riverine", "Bay", #FCE
              "Back Reef", "Fore Reef", "Fringing Reef", #MCR
              "Marine Protected Area", "Reference", #PISCO-Central, PISCO-South, SBC-Beach, & SBC-Ocean
              "Seagrass", "Sand")) #VCR
print(habitat_mapping) #yayayay

dat_ready_2 <- dat_ready |> 
  left_join(label_mapping, by = "projecthabitat") |>
  left_join(habitat_mapping, by = "color") |>
  ### remove columns needed for joins up to this point
  select(-projecthabitat, -habitat, -project, -color, -site) |> 
  ### rename columns to be more representative/clean
  rename(site = color2,
         program = Project, 
         habitat = Habitat) |> 
  dplyr::select(program, habitat, site, year, month, everything())

glimpse(dat_ready_2)
unique(dat_ready_2$habitat)

dat_ready_3 <- dat_ready_2 |> 
  filter(!site %in% c("TB-5", "RB-17", "RB-19") ) |> 
  select(-strata, -subsite_level1, -subsite_level2, -subsite_level3, 
         -group, -units)

glimpse(dat_ready_3)
unique(dat_ready_3$site)
summary(dat_ready_3)

### summarize all sites measured within the dataset annualy, then across period of record
model_dt <- dat_ready_3 |> 
  group_by(program, habitat, site, year) |> 
      summarize(total_nitrogen_ann = mean(total_nitrogen),
                total_phosphorus_ann = mean(total_phosphorus),
                total_biomass_ann = mean(total_biomass),
                max_size_ann = mean(max_size),
                skew_size_ann = mean(size_skewness),
                species_richness_ann = mean(species_richness),
                species_diversity_ann = mean(species_diversity),
                trophic_richness_ann = mean(trophic_richness),
                trophic_diversity_ann = mean(trophic_diversity)) |> 
  ungroup() |> 
  group_by(program, habitat, site) |> 
  summarize(comm_mean_n = mean(total_nitrogen_ann),
            comm_sd_n = sd(total_nitrogen_ann),
            comm_cv_n = (sd(total_nitrogen_ann, na.rm = TRUE) / mean(total_nitrogen_ann, na.rm = TRUE)),
            ### nitrogen supply stability 
            comm_n_stability = 1/comm_cv_n,
            comm_mean_p = mean(total_phosphorus_ann),
            comm_sd_p = sd(total_phosphorus_ann),
            comm_cv_p = (sd(total_phosphorus_ann, na.rm = TRUE) / mean(total_phosphorus_ann, na.rm = TRUE)),
            ### phosphorus supply stability
            comm_p_stability = 1/comm_cv_p,
            comm_mean_bm = mean(total_biomass_ann),
            comm_sd_bm = sd(total_biomass_ann),
            comm_cv_bm = (sd(total_biomass_ann, na.rm = TRUE) / mean(total_biomass_ann, na.rm = TRUE)),
            ### biomass stability
            comm_bm_stability = 1/comm_cv_bm,
            comm_mean_max_ss = mean(max_size_ann),
            comm_cv_max_ss = (sd(max_size_ann, na.rm = TRUE) / mean(max_size_ann, na.rm = TRUE)),
            ### max size stability
            comm_max_size_stability = 1/comm_cv_max_ss,
            comm_mean_skew_ss = mean(skew_size_ann),
            comm_cv_skew_ss = (sd(skew_size_ann, na.rm = TRUE) / mean(skew_size_ann, na.rm = TRUE)),
            ### size structure stability
            comm_skew_size_stability = 1/comm_cv_skew_ss,
            mean_species_richness = mean(species_richness_ann),
            cv_species_richness = (sd(species_richness_ann, na.rm = TRUE) / mean(species_richness_ann, na.rm = TRUE)),
            ### species richness stability
            species_richness_stability = 1/cv_species_richness,
            mean_species_diversity = mean(species_diversity_ann),
            cv_species_diversity = (sd(species_diversity_ann, na.rm = TRUE) / mean(species_diversity_ann, na.rm = TRUE)),
            ### species diversity stability
            species_diversity_stability = 1/cv_species_diversity,
            mean_trophic_richness = mean(trophic_richness_ann),
            cv_trophic_richness = (sd(trophic_richness_ann, na.rm = TRUE) / mean(trophic_richness_ann, na.rm = TRUE)),
            ### trophic richness stability
            trophic_richness_stability = 1/cv_trophic_richness,
            mean_trophic_diversity = mean(trophic_diversity_ann),
            cv_trophic_diversity = (sd(trophic_diversity_ann, na.rm = TRUE) / mean(trophic_diversity_ann, na.rm = TRUE)),
            ### trophic diversity stability
            trophic_diversity_stability = 1/cv_trophic_diversity,
            raw_n_bm_differences = abs(comm_n_stability - comm_bm_stability),
            raw_n_p_differences = abs(comm_p_stability - comm_bm_stability))|> 
      ungroup()

model_dt_1 <- model_dt
### keeping everything for now, instead of filtering as previously
# model_dt_1 <- model_dt |> 
#       select(program, habitat, site, 
#              comm_mean_bm, comm_sd_bm, comm_bm_stability, 
#              comm_mean_n, comm_sd_n, comm_n_stability, 
#              comm_mean_p, comm_sd_p, comm_p_stability, 
#              comm_mean_max_ss, comm_mean_skew_ss,
#              mean_species_richness, 
#              mean_species_diversity, mean_trophic_richness, mean_trophic_diversity)

glimpse(model_dt_1)

# write_csv(model_dt_1, "2025-data/community-level-nutrient-stability.csv")

model_dt_1 |>
      ggplot(aes(comm_n_stability, comm_p_stability))+
      geom_point()+
      geom_abline()

model_dt_1 |>
      ggplot(aes(comm_bm_stability, comm_n_stability))+
      geom_point()+
      geom_abline()

### saving data at annual timesteps for summary figures and statistics
annual_dt <- dat_ready_3 |> 
      group_by(program, habitat, site, year) |> 
      summarize(total_nitrogen_ann = mean(total_nitrogen),
                total_phosphorus_ann = mean(total_phosphorus),
                total_biomass_ann = mean(total_biomass),
                max_size_ann = mean(max_size),
                skew_size_ann = mean(size_skewness),
                species_richness_ann = mean(species_richness),
                species_diversity_ann = mean(species_diversity),
                trophic_richness_ann = mean(trophic_richness),
                trophic_diversity_ann = mean(trophic_diversity)) |> 
      ungroup()
glimpse(annual_dt)
# write_csv(annual_dt, "local_data/annual-dt-for-summary11012024.csv")
# write_csv(annual_dt, "local_data/annual-dt-for-summary11142024.csv")
# write_csv(annual_dt, "2025-data/annual-dt-for-summary.csv")
