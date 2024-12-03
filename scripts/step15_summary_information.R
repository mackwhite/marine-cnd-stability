###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): MW
###goal(s): wrangling raw cnd data for summary data
###date(s): november 2024 @ nceas
###note(s): 

# Housekeeping ------------------------------------------------------------

### load necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, googledrive, vegan, readxl, e1071, dplyr, splitstackshape)

### read in data ----
df <- read_csv("local_data/annual-dt-for-summary.csv")

summary <- df |> 
      group_by(program) |> 
      summarize(
            sites = n_distinct(site),
            habitats = n_distinct(habitat),
            year_min = min(year),
            year_max = max(year),
            biomass_m = mean(total_biomass_ann, na.rm = TRUE),
            biomass_sd = sd(total_biomass_ann, na.rm = TRUE),
            nsupply_m = mean(total_nitrogen_ann/1000, na.rm = TRUE),
            nsupply_sd = sd(total_nitrogen_ann/1000, na.rm = TRUE),
            spprich_m = mean(species_richness_ann, na.rm = TRUE),
            spprich_sd = sd(species_richness_ann, na.rm = TRUE),
            sppdiv_m = mean(species_diversity_ann, na.rm = TRUE),
            sppdiv_sd = sd(species_diversity_ann, na.rm = TRUE), 
            trophrich_m = mean(trophic_richness_ann, na.rm = TRUE), 
            trophrich_sd = sd(trophic_richness_ann, na.rm = TRUE),
            trophdiv_m = mean(trophic_diversity_ann, na.rm = TRUE),
            trophdiv_sd = sd(trophic_diversity_ann, na.rm = TRUE)
      ) |> 
      mutate(across(biomass_m:trophdiv_sd,\(x) round(x, digits = 2)),
             years = paste(year_min, year_max, sep = "-")) |>
      ungroup() |> 
      select(-year_min, -year_max) |> 
      select(program, habitats, sites, years, everything())

stability <- read_csv("local_data/community-level-nutrient-stability.csv")

summary2 <- stability |> 
      group_by(program) |> 
      summarize(
            stability_m = mean(comm_n_stability),
            stability_sd = sd(comm_n_stability)
      )

summary_all <- left_join(summary, summary2, by = "program")

# write_csv(summary_all, "output/ms-second-round/tables/site-summary-table.csv")

