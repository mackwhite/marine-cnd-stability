###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): Mack White
###goal(s): visualize scaled dsr relationships across various scales of organization
###date(s): October 2024
###note(s):

###########################################################################
# Housekeeping ------------------------------------------------------------
###########################################################################

### load necessary libraries
### install.packages("librarian")
# remotes::install_github('m-clark/mixedup')
librarian::shelf(
      tidyverse,
      readxl,
      MuMIn,
      sjPlot,
      lme4,
      corrplot,
      performance,
      ggeffects,
      ggpubr,
      parameters,
      ggstats,
      brms,
      mixedup,
      lterpalettefinder
)

### read in necessary data ---

dat <- read_csv('local_data/dsr-eco-org-raw-all.csv') |>
      rename(
            Program = program,
            Trophic_Group = troph_group,
            Species = scientific_name,
            Habitat = habitat,
            Site = site
      ) |>
      select(
            Program,
            Habitat,
            Site,
            comm_mean_max_ss,
            comm_mean_skew_ss,
            comm_mean_bm,
            comm_sd_bm,
            comm_bm_stability,
            comm_mean_n,
            comm_sd_n,
            comm_n_stability,
            comm_mean_p,
            comm_sd_p,
            comm_p_stability,
            mean_species_richness,
            mean_species_diversity,
            mean_trophic_richness,
            mean_trophic_diversity,
            beta_time,
            synch,
            troph_beta_time,
            troph_synch
      ) |>
      distinct()

dat_scaled <- dat |>
      select(Program, Habitat, Site, comm_n_stability, everything()) |>
      mutate(comm_n_stability = scale(comm_n_stability)) |>
      group_by(Program) |>
      ## this is a function syntax
      mutate(across(comm_mean_max_ss:troph_synch, \(x) scale(x, center = TRUE))) |>
      ungroup()

glimpse(dat_scaled)

dat_ready <- dat_scaled
### clean env (optional) ---
rm(dat_scaled, dat)

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

disturbance_palette <- c("cold snap"="#44AA99",
                         "drought"="#D55E00",
                         "cyclone"="#117733",
                         "COT"="#8b6b93",
                         "heat wave"="#CC6677",
                         "disease"="#999933")

thermal_palette <- c('thermal'='#CC6677',
                     'nonthermal'='black')

### set priors following Lemoine (2019, Ecology)
pr = prior(normal(0, 1), class = 'b')

###########################################################################
# full models -------------------------------------------------------------
###########################################################################
test_corr <- dat_ready |> select(mean_species_richness,mean_species_diversity,
                                 mean_trophic_richness,mean_trophic_diversity,
                                 comm_mean_max_ss,comm_mean_skew_ss,
                                 beta_time,synch,
                                 troph_beta_time,troph_synch)
matrix <- cor(test_corr, use = 'complete.obs')
corrplot(matrix, method = "number", type = "lower", tl.col = "black", tl.srt = 45)

### following correlations to be aware of
# all of the richness and diversity metrics
# comm_mean_skew_ss ~ both richness metrics
# troph_synch ~ synch
# regardless, once synchrony or turnover is selected at a level of organization
# will quit dropping from the model selection process
rm(matrix,test_corr)
### round one ---
m1 <- brm(
      comm_n_stability ~ mean_trophic_richness + (mean_trophic_richness |
                                                        Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

# m2 <- brm(
#       comm_n_stability ~ comm_mean_max_ss + (comm_mean_max_ss | Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )
# 
# m3 <- brm(
#       comm_n_stability ~ comm_mean_skew_ss + (comm_mean_skew_ss |
#                                                     Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

m4 <- brm(
      comm_n_stability ~ synch + (synch | Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

saveRDS(m4, file = 'local_data/rds-single-synchrony.rds')

m5 <- brm(
      comm_n_stability ~ beta_time + (beta_time | Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

m6 <- brm(
      comm_n_stability ~ troph_synch + (troph_synch | Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

m7 <- brm(
      comm_n_stability ~ troph_beta_time + (troph_beta_time | Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

m8 <- brm(
      comm_n_stability ~ mean_species_richness + (mean_species_richness | Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)
saveRDS(m8, file = 'local_data/rds-single-richness.rds')

m9 <- brm(
      comm_n_stability ~ mean_species_diversity + (mean_species_diversity | Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

m10 <- brm(
      comm_n_stability ~ mean_trophic_diversity + (mean_trophic_diversity | Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)


# model_table <- performance::compare_performance(m1,m2,m3,m4,m5,m6,m7)
model_table_all <- performance::compare_performance(m1,m4,m5,m6,m7,m8,m9,m10)
model_selection <- model_table_all |>
      mutate(dWAIC = WAIC - min(WAIC))
# write_csv(model_selection, "output/ms-second-round/tables/brm-fullmodel-selection-table-roundone.csv")
write_csv(model_selection, "output/ms-second-round/tables/december-brm-fullmodel-selection-table-roundone.csv")
rm(list = setdiff(ls(), c("dat_ready", "pr", "palette", 'm4')))

### round two ---
m41 <- brm(
      comm_n_stability ~ mean_trophic_richness + synch + (mean_trophic_richness + synch |
                                                                Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

# m42 <- brm(
#       comm_n_stability ~ comm_mean_max_ss + synch + (comm_mean_max_ss + synch |
#                                                            Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )
# 
# m43 <- brm(
#       comm_n_stability ~ comm_mean_skew_ss + synch + (comm_mean_skew_ss + synch |
#                                                             Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

m45 <- brm(
      comm_n_stability ~ beta_time + synch + (beta_time + synch |
                                                    Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

m47 <- brm(
      comm_n_stability ~ troph_beta_time + synch + (troph_beta_time + synch |
                                                          Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

m48 <- brm(
      comm_n_stability ~ mean_species_richness + synch + (mean_species_richness + synch |
                                                          Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

m49 <- brm(
      comm_n_stability ~ mean_trophic_diversity + synch + (mean_trophic_diversity + synch |
                                                           Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

m40 <- brm(
      comm_n_stability ~ mean_species_diversity + synch + (mean_species_diversity + synch |
                                                            Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

model_table_all <- performance::compare_performance(m41,m4,m45,m47,m48,m49,m40)
model_selection <- model_table_all |>
      mutate(dWAIC = WAIC - min(WAIC))
# write_csv(model_selection, "output/ms-second-round/tables/brm-fullmodel-selection-table-roundtwo.csv")
write_csv(model_selection, "output/ms-second-round/tables/december-brm-fullmodel-selection-table-roundtwo.csv")
rm(list = setdiff(ls(), c("dat_ready", "pr", "palette",'m47')))

### round three - 

m471 <- brm(
      comm_n_stability ~ mean_trophic_richness + troph_beta_time + synch + (troph_beta_time + mean_trophic_richness + synch |
                                                                Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

# m472 <- brm(
#       comm_n_stability ~ comm_mean_max_ss + synch + troph_beta_time + (troph_beta_time + comm_mean_max_ss + synch |
#                                                            Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )
# 
# m473 <- brm(
#       comm_n_stability ~ comm_mean_skew_ss + synch + troph_beta_time + (troph_beta_time + comm_mean_skew_ss + synch |
#                                                             Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

m475 <- brm(
      comm_n_stability ~ beta_time + synch + troph_beta_time + (troph_beta_time + beta_time + synch |
                                                    Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)


m478 <- brm(
      comm_n_stability ~ mean_species_richness + troph_beta_time + synch + (troph_beta_time + mean_species_richness + synch |
                                                                Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

m479 <- brm(
      comm_n_stability ~ mean_trophic_diversity + synch + troph_beta_time + (troph_beta_time + mean_trophic_diversity + synch |
                                                                 Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

m470 <- brm(
      comm_n_stability ~ mean_species_diversity + synch + troph_beta_time + (troph_beta_time + mean_species_diversity + synch |
                                                                 Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

model_table <- performance::compare_performance(m47,m471,m475,m478,m479,m470)
model_selection_full <- model_table |>
      mutate(dWAIC = WAIC - min(WAIC))
# write_csv(model_selection_full, "output/ms-second-round/tables/brm-fullmodel-selection-table-roundthree.csv")
write_csv(model_selection_full, "output/ms-second-round/tables/december-brm-fullmodel-selection-table-roundthree.csv")
# rm(list = setdiff(ls(), c("dat_ready", "pr", "program_palette",'m475')))
rm(list = setdiff(ls(), c("dat_ready", "pr", "program_palette",'m470')))

### round four ---
# m4701 <- brm(
#       comm_n_stability ~ mean_trophic_richness + troph_beta_time + synch + beta_time + (beta_time + troph_beta_time + mean_trophic_richness + synch |
#                                                                                   Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )
# 
# m4752 <- brm(
#       comm_n_stability ~ comm_mean_max_ss + synch + troph_beta_time + beta_time + (beta_time + troph_beta_time + comm_mean_max_ss + synch |
#                                                                              Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

# m4753 <- brm(
#       comm_n_stability ~ comm_mean_skew_ss + synch + troph_beta_time + beta_time + (beta_time + troph_beta_time + comm_mean_skew_ss + synch |
#                                                                               Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

m4705 <- brm(
      comm_n_stability ~ mean_species_diversity + beta_time + synch + troph_beta_time + (mean_species_diversity+ troph_beta_time + beta_time + synch |
                                                                      Program),
      data = dat_ready,
      prior = pr,
      warmup = 1000,
      iter = 10000,
      chains = 4
)

# m4708 <- brm(
#       comm_n_stability ~ mean_species_richness + troph_beta_time + synch + beta_time + (beta_time + troph_beta_time + mean_species_richness + synch |
#                                                                                   Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

# m4709 <- brm(
#       comm_n_stability ~ mean_trophic_diversity + synch + troph_beta_time + beta_time + (beta_time + troph_beta_time + mean_trophic_diversity + synch |
#                                                                                    Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

# m470 <- brm(
#       comm_n_stability ~ mean_species_diversity + synch + troph_beta_time + beta_time + (beta_time + troph_beta_time + mean_species_diversity + synch |
#                                                                                    Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

model_table <- performance::compare_performance(m470,m4705)
# model_table <- performance::compare_performance(m475,m4750)
model_selection_full <- model_table |>
      mutate(dWAIC = WAIC - min(WAIC))
# write_csv(model_selection_full, "output/ms-second-round/tables/brm-fullmodel-selection-table-roundfour.csv")
write_csv(model_selection_full, "output/ms-second-round/tables/december-brm-fullmodel-selection-table-roundfour.csv")
rm(list = setdiff(ls(), c("dat_ready", "pr", "program_palette",'m4750')))
full_model <- m4705
saveRDS(full_model, file = 'local_data/rds-full-model.rds')

### round five ---

# m47501 <- brm(
#       comm_n_stability ~ mean_species_diversity + mean_trophic_richness + troph_beta_time + synch + beta_time + (beta_time + troph_beta_time + mean_trophic_richness + synch |
#                                                                                               Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

# m47502 <- brm(
#       comm_n_stability ~ mean_species_diversity + comm_mean_max_ss + synch + troph_beta_time + beta_time + (mean_species_diversity + beta_time + troph_beta_time + comm_mean_max_ss + synch |
#                                                                                          Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

# m47503 <- brm(
#       comm_n_stability ~ mean_species_diversity + comm_mean_skew_ss + synch + troph_beta_time + beta_time + (mean_species_diversity + beta_time + troph_beta_time + comm_mean_skew_ss + synch |
#                                                                                           Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

# m47508 <- brm(
#       comm_n_stability ~ mean_species_diversity + mean_species_richness + troph_beta_time + synch + beta_time + (beta_time + troph_beta_time + mean_species_richness + synch |
#                                                                                               Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

# m47509 <- brm(
#       comm_n_stability ~ mean_species_diversity + mean_trophic_diversity + synch + troph_beta_time + beta_time + (beta_time + troph_beta_time + mean_trophic_diversity + synch |
#                                                                                                Program),
#       data = dat_ready,
#       prior = pr,
#       warmup = 1000,
#       iter = 10000,
#       chains = 4
# )

# model_table <- performance::compare_performance(m4750,m47502,m47503)
# model_selection_full <- model_table |>
#       mutate(dWAIC = WAIC - min(WAIC))
# write_csv(model_selection_full, "output/ms-second-round/tables/brm-fullmodel-selection-table-roundfive.csv")
# rm(list = setdiff(ls(), c("dat_ready", "pr", "program_palette",'m4750')))


# ### visualize and save model predictions ---
# full_model <- m4750
# pp_check(full_model)
# full_model_re_slope <- mixedup::extract_random_coefs(full_model)
# full_model_re_slope_exp <- full_model_re_slope |>
#       mutate(nozero = map2_lgl(lower_2.5, upper_97.5, \(x, y) between(0, x, y)))
# full_model_fe_slope <- mixedup::extract_fixed_effects(full_model)
# summary(full_model)
# # save(full_model, file = "output/ms-second-round/models/fullmodel.RData")
# 
# # performance::performance(full_model)
# 
# ### visualizations for synchrony ---
# synchrony_re1 <- ggpredict(full_model,
#                            type = "re",
#                            terms = c('synch[-2:3 by=0.01]', 'Program'))
# synchrony_re <- as.data.frame(synchrony_re1)
# 
# synchrony_fe1 <- ggpredict(full_model,
#                            type = "fe",
#                            terms = c('synch[-2:3 by=0.01]'))
# synchrony_fe <- as.data.frame(synchrony_fe1) |>
#       mutate(group = 'Overall')
# 
# synchrony_all <- rbind(synchrony_re, synchrony_fe)
# 
# synchrony_all |>
#       mutate(group = factor(
#             group,
#             levels = c("Overall", "FCE", "MCR", "PCCC", "PCCS", "SBC", "VCR")
#       )) |>
#       ggplot(aes(x = x, y = predicted, color = group)) +
#       geom_smooth(method = "lm", linewidth = 1.5) +
#       labs(x = 'Population Synchrony', y = 'Aggregate Nitrogen Supply Stability', color = 'Program') +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(
#             axis.text.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.text.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.title.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.x = element_blank(),
#             axis.title.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.y = element_blank(),
#             legend.position = "none",
#             legend.background = element_blank(),
#             legend.key = element_rect(fill = 'white'),
#             legend.text = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             legend.title = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             )
#       )
# 
# ggsave("output/ms-second-round/plots/brm-fullmodel-population-synchrony.tiff", units = "in", width = 5,
#        height = 5, dpi =  600, compression = "lzw")
# ggsave("output/ms-second-round/plots/brm-fullmodel-population-synchrony.svg", units = "in", width = 5,
#        height = 5, dpi =  600)
# 
# ### visualizations for turnover ---
# trophic_turnover_re1 <- ggpredict(
#       full_model,
#       type = "re",
#       terms = c('troph_beta_time[-2:3 by=0.01]', 'Program')
# )
# trophic_turnover_re <- as.data.frame(trophic_turnover_re1)
# 
# trophic_turnover_fe1 <- ggpredict(full_model,
#                                   type = "fe",
#                                   terms = c('troph_beta_time[-2:3 by=0.01]'))
# trophic_turnover_fe <- as.data.frame(trophic_turnover_fe1) |>
#       mutate(group = 'Overall')
# 
# trophic_turnover_all <- rbind(trophic_turnover_re, trophic_turnover_fe)
# 
# trophic_turnover_all |>
#       mutate(group = factor(
#             group,
#             levels = c("Overall", "FCE", "MCR", "PCCC", "PCCS", "SBC", "VCR")
#       )) |>
#       ggplot(aes(x = x, y = predicted, color = group)) +
#       geom_smooth(method = "lm", linewidth = 1.5) +
#       labs(x = 'Trophic Group Turnover', y = 'Aggregate Nitrogen Supply Stability', color = 'Program') +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(
#             axis.text.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.text.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.title.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.x = element_blank(),
#             axis.title.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.y = element_blank(),
#             legend.position = "none",
#             legend.background = element_blank(),
#             legend.key = element_rect(fill = 'white'),
#             legend.text = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             legend.title = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             )
#       )
# 
# ggsave("output/ms-second-round/plots/brm-fullmodel-trophic-turnover.tiff", units = "in", width = 5,
#        height = 5, dpi =  600, compression = "lzw")
# ggsave("output/ms-second-round/plots/brm-fullmodel-trophic-turnover.svg", units = "in", width = 5,
#        height = 5, dpi =  600)
# 
# mech_fe <- full_model_fe_slope |>
#       rename(effect = term) |>
#       filter(effect != 'Intercept') |>
#       mutate(group = "Overall") |>
#       select(group, effect, value, se, lower_2.5, upper_97.5)
# 
# mech_re <- full_model_re_slope |>
#       filter(effect != "Intercept") |>
#       select(group, effect, value, se, lower_2.5, upper_97.5)
# 
# mech_slopes <- rbind(mech_fe, mech_re)
# glimpse(mech_slopes)
# 
# a <- mech_slopes |>
#       mutate(group = factor(
#             group,
#             levels = c("PCCC", "SBC", "VCR", "PCCS", "FCE", "MCR", "Overall")
#       )) |>
#       filter(effect == "synch") |>
#       ggplot(aes(x = value, y = group, color = group)) +
#       geom_point(size = 3) +
#       geom_errorbarh(aes(xmin = lower_2.5, xmax = upper_97.5),
#                      size = 1,
#                      height = 0) +
#       geom_vline(xintercept = 0, size = 1) +
#       # scale_x_continuous(labels = function(x) sprintf("%.1f", x),
#       #                    limits = c(-1.1,0.3),
#       #                    breaks = seq(-0.9,0.3, by = 0.3)) +
#       labs(x = 'Beta', y = 'Program') +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(
#             axis.text.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.text.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.title.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.x = element_blank(),
#             axis.title.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.y = element_blank(),
#             legend.position = "none",
#             legend.background = element_blank(),
#             legend.key = element_rect(fill = 'white'),
#             legend.text = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             legend.title = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             )
#       )
# ggsave("output/ms-second-round/plots/brm-fullmodel-populationsynchrony-slopecoeff.tiff", units = "in", width = 5,
#        height = 5, dpi =  600, compression = "lzw")
# ggsave("output/ms-second-round/plots/brm-fullmodel-populationsynchrony-slopecoeff.svg", units = "in", width = 5,
#        height = 5, dpi =  600)
# 
# b <- mech_slopes |>
#       mutate(group = factor(
#             group,
#             levels = c("PCCC", "SBC", "VCR", "PCCS", "FCE", "MCR", "Overall")
#       )) |>
#       filter(effect == "troph_beta_time") |>
#       ggplot(aes(x = value, y = group, color = group)) +
#       geom_point(size = 3) +
#       geom_errorbarh(aes(xmin = lower_2.5, xmax = upper_97.5),
#                      size = 1,
#                      height = 0) +
#       geom_vline(xintercept = 0, size = 1) +
#       scale_x_continuous(labels = function(x) sprintf("%.1f", x),
#                          limits = c(-1.0,0.3),
#                          breaks = seq(-0.9,0.3, by = 0.3)) +
#       labs(x = 'Beta', y = 'Program') +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(
#             axis.text.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.text.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.title.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.x = element_blank(),
#             axis.title.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.y = element_blank(),
#             legend.position = "none",
#             legend.background = element_blank(),
#             legend.key = element_rect(fill = 'white'),
#             legend.text = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             legend.title = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             )
#       )
# ggsave("output/ms-second-round/plots/brm-fullmodel-trophicturnover-slopecoeff.tiff", units = "in", width = 5,
#        height = 5, dpi =  600, compression = "lzw")
# ggsave("output/ms-second-round/plots/brm-fullmodel-trophicturnover-slopecoeff.svg", units = "in", width = 5,
#        height = 5, dpi =  600)
# 
# ### for the legend ---
# 
# # mech_slopes |>
# #       mutate(group = factor(
# #             group,
# #             levels = c("SBC", "PCCC", "VCR", "FCE", "PCCS", "MCR", "Overall")
# #       )) |>
# #       filter(effect == "troph_beta_time") |>
# #       ggplot(aes(x = value, y = group, color = group)) +
# #       # geom_point(size = 3) +
# #       geom_errorbarh(aes(xmin = lower_2.5, xmax = upper_97.5),
# #                      size = 1,
# #                      height = 0) +
# #       geom_vline(xintercept = 0, size = 1) +
# #       scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
# #       labs(x = 'Beta', y = 'Program') +
# #       theme_classic() +
# #       scale_color_manual(values = program_palette) +
# #       theme(
# #             axis.text.x = element_text(
# #                   face = "bold",
# #                   color = "black",
# #                   size = 12
# #             ),
# #             axis.text.y = element_text(
# #                   face = "bold",
# #                   color = "black",
# #                   size = 12
# #             ),
# #             axis.title.x = element_text(
# #                   face = "bold",
# #                   color = "black",
# #                   size = 14
# #             ),
# #             # axis.title.x = element_blank(),
# #             axis.title.y = element_text(
# #                   face = "bold",
# #                   color = "black",
# #                   size = 14
# #             ),
# #             # axis.title.y = element_blank(),
# #             legend.position = "right",
# #             legend.background = element_blank(),
# #             legend.key = element_rect(fill = 'white'),
# #             legend.text = element_text(
# #                   face = "bold",
# #                   color = "black",
# #                   size = 20
# #             ),
# #             legend.title = element_text(
# #                   face = "bold",
# #                   color = "black",
# #                   size = 14
# #             )
# #       ) +
# #       guides(
# #             color = guide_legend(override.aes = list(size = 8))
# #       )
# # 
# # ggsave("output/ms-second-round/plots/brm-fullmodel-legend.tiff", units = "in", width = 5,
# #        height = 5, dpi =  600, compression = "lzw")
# # ggsave("output/ms-second-round/plots/brm-fullmodel-legend.svg", units = "in", width = 5,
# #        height = 5, dpi =  600)
# 
# c <- mech_slopes |>
#       mutate(group = factor(
#             group,
#             levels = c("PCCC", "SBC", "VCR", "PCCS", "FCE", "MCR", "Overall")
#       )) |>
#       filter(effect == "beta_time") |>
#       ggplot(aes(x = value, y = group, color = group)) +
#       geom_point(size = 3) +
#       geom_errorbarh(aes(xmin = lower_2.5, xmax = upper_97.5),
#                      size = 1,
#                      height = 0) +
#       geom_vline(xintercept = 0, size = 1) +
#       # scale_x_continuous(labels = function(x) sprintf("%.1f", x),
#                          # limits = c(-1.1,0.3),
#                          # breaks = seq(-0.9,0.3, by = 0.3)) +
#       labs(x = 'Beta', y = 'Program') +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(
#             axis.text.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.text.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.title.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.x = element_blank(),
#             axis.title.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.y = element_blank(),
#             legend.position = "none",
#             legend.background = element_blank(),
#             legend.key = element_rect(fill = 'white'),
#             legend.text = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             legend.title = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             )
#       )
# 
# d <- mech_slopes |>
#       mutate(group = factor(
#             group,
#             levels = c("PCCC", "SBC", "VCR", "PCCS", "FCE", "MCR", "Overall")
#       )) |>
#       filter(effect == "mean_species_diversity") |>
#       ggplot(aes(x = value, y = group, color = group)) +
#       geom_point(size = 3) +
#       geom_errorbarh(aes(xmin = lower_2.5, xmax = upper_97.5),
#                      size = 1,
#                      height = 0) +
#       geom_vline(xintercept = 0, size = 1) +
#       # scale_x_continuous(labels = function(x) sprintf("%.1f", x),
#       # limits = c(-1.1,0.3),
#       # breaks = seq(-0.9,0.3, by = 0.3)) +
#       labs(x = 'Beta', y = 'Program') +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(
#             axis.text.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.text.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             axis.title.x = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.x = element_blank(),
#             axis.title.y = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             ),
#             # axis.title.y = element_blank(),
#             legend.position = "none",
#             legend.background = element_blank(),
#             legend.key = element_rect(fill = 'white'),
#             legend.text = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 12
#             ),
#             legend.title = element_text(
#                   face = "bold",
#                   color = "black",
#                   size = 14
#             )
#       )
# 
# test <- ggarrange(a,b,c,d)
# 
# count <- dat_ready |> group_by(Program) |> count()

# save(full_model, file = "local_data/full_model_11182024.RData")
# save(full_model, file = "local_data/full_model_11192024.RData")
saveRDS(full_model, file = 'local_data/rds-full-model.rds')
