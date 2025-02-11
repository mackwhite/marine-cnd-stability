###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): Mack White
###goal(s): check out information you have for individual species
###date(s): February 2025
###note(s): 

# Housekeeping ------------------------------------------------------------

### load necessary libraries
### install.packages("librarian")

librarian::shelf(tidyverse, readr, readxl, rfishbase, fishflux)

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
test <- dt_og |> group_by(project,scientific_name) |> summarize(max_size = max(dmperind_g_ind), mean_size = mean(dmperind_g_ind))
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


spp_list <- dt_mutate |> 
      ungroup() |> 
      dplyr::select(scientific_name, common_name, kingdom, phylum, class, order, family, genus) |> 
      distinct() |> 
      mutate(species = scientific_name)

ests <- rfishbase::estimate(spp_list$species) |> 
      janitor::clean_names() |> 
      rename(scientific_name = species,
             tl_max = max_length_tl,
             sl_max = max_length_sl) |> 
      select(scientific_name, feeding_path, troph, k, depth_max, depth_min,
             temp_pref_min, temp_pref_mean, temp_pref_max,
             tl_max, sl_max)

spp_fill <- spp_list |> 
      left_join(ests, by = "scientific_name") |> 
      group_by(genus) |> 
      mutate(across(where(is.numeric),
                    ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) |> 
      ungroup() |> 
      group_by(family) |> 
      mutate(across(where(is.numeric),
                    ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) |> 
      ungroup() |> 
      group_by(order) |> 
      mutate(across(where(is.numeric),
                    ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) |> 
      ungroup() |> 
      mutate(across(where(is.numeric),
                    ~ifelse(is.nan(.), NA, .)))

get_mode <- function(x) {
      ux <- na.omit(unique(x))
      ux[which.max(tabulate(match(x, ux)))]
}

filled_data <- spp_fill |> 
      group_by(genus) |> 
      mutate(feeding_path = ifelse(is.na(feeding_path), get_mode(feeding_path), feeding_path)) |> 
      ungroup()  |> 
      group_by(family) |> 
      mutate(feeding_path = ifelse(is.na(feeding_path), get_mode(feeding_path), feeding_path)) |> 
      ungroup() |> 
      group_by(order) |> 
      mutate(feeding_path = ifelse(is.na(feeding_path), get_mode(feeding_path), feeding_path)) |> 
      ungroup()

lw <- growth_params(ests$scientific_name)

lw_distinct <- lw |> 
      rename(scientific_name = species) |> 
      group_by(scientific_name) |> 
      summarize(across(where(is.numeric), ~mean(., na.rm = TRUE))) |> 
      mutate(across(where(is.numeric), 
                    ~ifelse(is.nan(.), NA, .)))

spp_fill_test <- spp_fill |> 
      mutate(temp = 28.0)

library(dplyr)
library(purrr)
library(fishflux)  # Assuming this is the correct package

results <- spp_fill_test |> 
      filter(!is.na(family)) |> 
      mutate(
            metabolism_result = pmap(
                  list(family = family, temp = temp, troph_m = troph, troph_sd = 1e-10),
                  fishflux::metabolism
            )
      ) |> 
      unnest_wider(metabolism_result)

# View the results
head(results)

