###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): Mack White
###goal(s): visualize scaled dsr relationships across various scales of organization
###date(s): January 2025
###note(s): 

###########################################################################
# Housekeeping ------------------------------------------------------------
###########################################################################

### load necessary libraries
### install.packages("librarian")
# remotes::install_github('m-clark/mixedup')
librarian::shelf(tidyverse, readxl, MuMIn, sjPlot, lme4, corrplot, glmmTMB,
                 performance, ggeffects, ggpubr, parameters, ggstats, brms, mixedup, lterpalettefinder)

### read in necessary data ---

dat <- read_csv('2025-data/dsr-eco-org-raw-all.csv') |> 
      rename(trophic_group = troph_group,
             species = scientific_name,
             species_turnover = beta_time,
             species_synchrony = synch,
             trophic_turnover = troph_beta_time,
             trophic_synchrony = troph_synch,
             max_size = comm_mean_max_ss,
             max_size_stability = comm_max_size_stability,
             size_structure = comm_mean_skew_ss,
             size_structure_stability = comm_skew_size_stability,
             n_mean = comm_mean_n,
             n_stability = comm_n_stability,
             p_mean = comm_mean_p,
             p_stability = comm_p_stability,
             bm_mean = comm_mean_bm,
             bm_stability = comm_bm_stability,
             trophic_richness = mean_trophic_richness,
             trophic_diversity = mean_trophic_diversity,
             aggregate_n_bm_stability_diff = raw_n_bm_differences,
             aggregate_p_n_stability_diff = raw_n_p_differences
             ) |> 
      select(program, habitat, site,
             n_mean, n_stability,
             p_mean, p_stability,
             bm_mean, bm_stability,
             max_size, max_size_stability,
             size_structure, size_structure_stability,
             trophic_richness, trophic_richness_stability,
             trophic_diversity, trophic_diversity_stability,
             trophic_synchrony, trophic_turnover,
             aggregate_n_bm_stability_diff,
             aggregate_p_n_stability_diff) |> 
      distinct()

glimpse(dat)

dat_scaled <- dat |> 
      select(program, habitat, site, n_stability, p_stability, bm_stability,
             max_size_stability, size_structure_stability, trophic_richness_stability,
             trophic_diversity_stability, aggregate_n_bm_stability_diff, 
             aggregate_p_n_stability_diff, everything()) |> 
      mutate(across(n_stability:aggregate_p_n_stability_diff,\(x) scale(x, center = TRUE))) |> 
      group_by(program) |> 
      ## this is a function syntax
      mutate(across(n_mean:trophic_turnover,\(x) scale(x, center = TRUE))) |>
      ungroup()

glimpse(dat_scaled)

mod_df <- dat_scaled
### clean env (optional) ---
rm(dat_scaled, dat)
glimpse(mod_df)

### set color schemes ---
habitat_palette <- c("Overall"="#000000",
                     "Back Reef"="#fde725",
                     "Bay"="#addc30",
                     "Fore Reef"="#ff967d",
                     "Fringing Reef"="#21918c",
                     "Marine Protected Area"="#2c728e",
                     "MPA"="#2c728e",
                     "Reference"="#5ec962",
                     "Riverine"="#8b6b93",
                     "Sand"="#D2B48C",
                     "Seagrass"="#64a988")

program_palette <- c("Overall"="#000000", 
                     "FCE"="#64a988", 
                     "MCR"="#ff967d", 
                     'PCCC'="#2A788EFF", 
                     "PCCS"="#8b6b93",
                     'SBC'='#ff3f4c', 
                     "VCR"="#9b9254")

# disturbance_palette <- c("cold snap"="#44AA99",
#                          "drought"="#D55E00",
#                          "cyclone"="#117733",
#                          "COT"="#8b6b93",
#                          "heat wave"="#CC6677",
#                          "disease"="#999933")
# 
# thermal_palette <- c('thermal'='#CC6677',
#                      'nonthermal'='black')

# look at correlations ----------------------------------------------------
corr <- mod_df |> 
      select(n_stability, bm_stability, max_size_stability, size_structure_stability,
             trophic_richness_stability, trophic_diversity_stability, 
             n_mean, bm_mean, max_size, size_structure, trophic_richness, trophic_diversity,
             trophic_synchrony, trophic_turnover)

corr_matrix <- cor(corr, use = 'complete.obs')
corrplot(corr_matrix, method = "number", type = "lower", tl.col = "black", tl.srt = 45)

### important correlations to consider
## trophic richness stability & size structure stability 
## trophic richness stability & trophic diversity stability 
## trophic diversity & trophic richness 

### first run ---
global <- glmmTMB(n_stability ~ bm_stability + max_size_stability + size_structure_stability + 
                  trophic_richness_stability + trophic_diversity_stability + 
                  n_mean + bm_mean + max_size + size_structure + trophic_richness + 
                  trophic_diversity + trophic_synchrony + trophic_turnover + (1|program),
                  mod_df,
                  na.action = "na.fail",
                  family = gaussian(),
                  REML = FALSE
)

keep <- c("global", "mod_df") 
rm(list = setdiff(ls(), keep))

model_set <- dredge(global,
                   subset = !(`cond(trophic_richness_stability)`&&`cond(size_structure_stability)`) &
                   !(`cond(trophic_richness_stability)`&&`cond(trophic_diversity_stability)`) &
                   !(`cond(trophic_diversity)`&&`cond(trophic_richness)`)
                   ) |> 
      filter(delta < 2)

### second run ---
global <- glmmTMB(n_stability ~ bm_stability + max_size_stability + size_structure_stability + 
                        trophic_richness_stability + trophic_diversity_stability + 
                        trophic_synchrony + trophic_turnover + (1|program),
                  mod_df,
                  na.action = "na.fail",
                  family = gaussian(),
                  REML = FALSE
)

model_set2 <- dredge(global,
                    subset = !(`cond(trophic_richness_stability)`&&`cond(size_structure_stability)`) &
                          !(`cond(trophic_richness_stability)`&&`cond(trophic_diversity_stability)`)
) |> 
      filter(delta < 2)

### third run ---
global <- glmmTMB(n_stability ~ bm_stability + max_size_stability + size_structure_stability + 
                        trophic_richness_stability + trophic_diversity_stability + (1|program),
                  mod_df,
                  na.action = "na.fail",
                  family = gaussian(),
                  REML = FALSE
)

model_set3 <- dredge(global,
                     subset = !(`cond(trophic_richness_stability)`&&`cond(size_structure_stability)`) &
                           !(`cond(trophic_richness_stability)`&&`cond(trophic_diversity_stability)`)
) |> 
      filter(delta < 2)



bfm <- glmmTMB(n_stability ~ bm_stability + max_size_stability + trophic_richness_stability + (1|program),
               family = gaussian(), data = mod_df)
performance::performance(bfm)
summary(bfm)
