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
librarian::shelf(tidyverse, readxl, MuMIn, sjPlot, lme4, corrplot, 
                 performance, ggeffects, ggpubr, parameters, ggstats, brms, mixedup, lterpalettefinder)

### read in necessary data ---

dat <- read_csv('local_data/dsr-eco-org-raw-all.csv') |> 
       rename(Program = program,
              Trophic_Group = troph_group,
              Species = scientific_name,
              Habitat = habitat,
              Site = site) |> 
       select(Program, Habitat, Site, 
             comm_mean_max_ss, comm_mean_skew_ss,
             comm_mean_bm, comm_sd_bm, comm_bm_stability,
             comm_mean_n, comm_sd_n, comm_n_stability,
             comm_mean_p, comm_sd_p, comm_p_stability,
             mean_species_richness, mean_species_diversity, 
             mean_trophic_richness, mean_trophic_diversity,
             beta_time, synch,
             troph_beta_time, troph_synch) |> 
      distinct()

dat_scaled <- dat |> 
      select(Program, Habitat, Site, comm_n_stability, everything()) |> 
      mutate(comm_n_stability = scale(comm_n_stability)) |> 
      group_by(Program) |> 
      ## this is a function syntax
      mutate(across(comm_mean_max_ss:troph_synch,\(x) scale(x, center = TRUE))) |>
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
pr = prior(normal(0,1), class = 'b')

###########################################################################
# exploratory regressions -------------------------------------------------
###########################################################################

dat_ready |>
      ggplot(aes(x = mean_species_richness, y = comm_n_stability, color = Program)) +
      geom_point(size = 3) +  # Adds the scatter plot points
      geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
      geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "solid") +  # Adds overall slope line
      labs(x = "Scaled Species Richness",
           y = "Scaled Aggregate Nitrogen Supply Stability (1/CV)") +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_y_continuous(breaks = -2:3) +
      scale_x_continuous(breaks = -2:2) +
      theme(axis.text.x = element_text(face = "bold", color = "black"),
            axis.text.y = element_text(face = "bold", color = "black"),
            axis.title.x = element_text(face = "bold", color = "black"),
            axis.title.y = element_text(face = "bold", color = "black"),
            legend.position = "right",
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"))

ggsave("output/ms-second-round/plots/raw-dsr-spprich-comm-n-stability.tiff", units = "in", width = 5,
       height = 5, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/raw-dsr-spprich-comm-n-stability.svg", units = "in", width = 5,
       height = 5, dpi =  600)

# species_diversity <- dat_ready |> 
#       ggplot(aes(x = mean_species_diversity, y = comm_n_stability, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Species Diversity",
#            y = "Scaled Aggregate Nitrogen Supply Stability (1/CV)") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# trophic_richness <- dat_ready |> 
#       ggplot(aes(x = mean_trophic_richness, y = comm_n_stability, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Trophic Richness",
#            y = "Scaled Aggregate Nitrogen Supply Stability (1/CV)") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# trophic_diversity <- dat_ready |> 
#       ggplot(aes(x = mean_trophic_diversity, y = comm_n_stability, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Trophic Diversity",
#            y = "Scaled Aggregate Nitrogen Supply Stability (1/CV)") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# p1 <- ggarrange(species_richness, species_diversity, trophic_diversity, trophic_richness,
#                 labels = c('a)','b)','c)','d)'),
#                 legend = 'bottom', common.legend = TRUE,
#                 ncol = 2, nrow = 2, align = "h")
# 
# spprich_popsynch <- dat_ready |> 
#       ggplot(aes(x = mean_species_richness, y = synch, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "solid") +  # Adds overall slope line
#       labs(x = "Scaled Species Richness",
#            y = "Scaled Population Synchrony") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       scale_x_continuous(breaks = -2:2) +
#       scale_y_continuous(breaks = -2:3) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_text(face = "bold", color = "black"),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# sppdiv_popsynch <- dat_ready |> 
#       ggplot(aes(x = mean_species_diversity, y = synch, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Species Diversity",
#            y = "Scaled Population Temporal Synchrony") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# spprich_popturnover <- dat_ready |> 
#       ggplot(aes(x = mean_species_richness, y = beta_time, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Species Richness",
#            y = "Scaled Population Temporal Turnover") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# sppdiv_popturnover <- dat_ready |> 
#       ggplot(aes(x = mean_species_diversity, y = beta_time, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Species Diversity",
#            y = "Scaled Population Temporal Turnover") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# trophrich_trophsynch <- dat_ready |> 
#       ggplot(aes(x = mean_trophic_richness, y = troph_synch, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Trophic Richness",
#            y = "Scaled Trophic Temporal Synchrony") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# trophdiv_trophsynch <- dat_ready |> 
#       ggplot(aes(x = mean_trophic_diversity, y = troph_synch, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Trophic Diversity",
#            y = "Scaled Trophic Temporal Synchrony") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# trophrich_trophturnover <- dat_ready |> 
#       ggplot(aes(x = mean_trophic_richness, y = troph_beta_time, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Trophic Richness",
#            y = "Scaled Trophic Temporal Turnover") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# trophdiv_trophturnover <- dat_ready |> 
#       ggplot(aes(x = mean_trophic_diversity, y = troph_beta_time, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Trophic Diversity",
#            y = "Scaled Trophic Temporal Turnover") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))
# 
# rm(species_diversity,species_richness,trophic_diversity,trophic_richness,
#    sppdiv_popsynch,sppdiv_popturnover,spprich_popsynch,spprich_popturnover,
#    trophdiv_trophsynch,trophdiv_trophturnover,trophrich_trophsynch,trophrich_trophturnover,
#    p1)

# dat_ready |>
#       ggplot(aes(x = mean_species_richness, y = comm_mean_skew_ss, color = Program)) +
#       geom_point() +  # Adds the scatter plot points
#       geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
#       labs(x = "Scaled Trophic Diversity",
#            y = "Scaled Trophic Temporal Turnover") +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black"),
#             axis.text.y = element_text(face = "bold", color = "black"),
#             axis.title.x = element_text(face = "bold", color = "black"),
#             axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.text = element_text(face = "bold", color = "black"),
#             legend.title = element_text(face = "bold", color = "black"))

###########################################################################
# diversity models --------------------------------------------------------
###########################################################################
test_corr <- dat_ready |> select(mean_species_richness,mean_species_diversity,mean_trophic_richness,mean_trophic_diversity,comm_mean_max_ss,comm_mean_skew_ss)
matrix <- cor(test_corr, use = 'complete.obs')
corrplot(matrix, method = "number", type = "lower", tl.col = "black", tl.srt = 45)
### following correlations to be aware of
# all of the richness and diversity metrics
# comm_mean_skew_ss ~ both richness metrics

### diversity -- 
m1 <- brm(comm_n_stability ~ mean_species_richness + (mean_species_richness|Program),
          data = dat_ready,
          # prior = pr, warmup = 100, iter = 1000, chains = 4)
          prior = pr, warmup = 1000, iter = 10000, chains = 4)

m2 <- brm(comm_n_stability ~ mean_species_diversity + (mean_species_diversity|Program),
          data = dat_ready,
          # prior = pr, warmup = 100, iter = 1000, chains = 4)
          prior = pr, warmup = 1000, iter = 10000, chains = 4)

m3 <- brm(comm_n_stability ~ mean_trophic_richness + (mean_trophic_richness|Program),
          data = dat_ready,
          # prior = pr, warmup = 100, iter = 1000, chains = 4)
          prior = pr, warmup = 1000, iter = 10000, chains = 4)
          # prior = pr, warmup = 10000, iter = 100000, chains = 4)

m4 <- brm(comm_n_stability ~ mean_trophic_diversity + (mean_trophic_diversity|Program),
          data = dat_ready,
          # prior = pr, warmup = 100, iter = 1000, chains = 4)
          prior = pr, warmup = 1000, iter = 10000, chains = 4)

div_model_table <- performance::compare_performance(m1,m2,m3,m4)
div_model_selection <- div_model_table |> 
      mutate(dWAIC = WAIC - min(WAIC))
write_csv(div_model_selection, "output/ms-second-round/tables/brm-dsr-selection-table.csv")

# adding other static metrics ---------------------------------------------
### round one ---
m5 <- brm(comm_n_stability ~ comm_mean_max_ss + (comm_mean_max_ss|Program),
          data = dat_ready,
          # prior = pr, warmup = 100, iter = 1000, chains = 4)
          prior = pr, warmup = 1000, iter = 10000, chains = 4)

m6 <- brm(comm_n_stability ~ comm_mean_skew_ss + (comm_mean_skew_ss|Program),
          data = dat_ready,
          # prior = pr, warmup = 100, iter = 1000, chains = 4)
          prior = pr, warmup = 1000, iter = 10000, chains = 4)

static_model_table <- performance::compare_performance(m1,m2,m3,m4,m5,m6)
static_model_selection <- static_model_table |> 
      mutate(dWAIC = WAIC - min(WAIC))
write_csv(static_model_selection, "output/ms-second-round/tables/brm-static-selection-roundone.csv")

### round two ---
rm(m1,m2,m4,m5,m6)
m31 <- brm(comm_n_stability ~ mean_trophic_richness + comm_mean_max_ss + (mean_trophic_richness + comm_mean_max_ss|Program),
          data = dat_ready,
          # prior = pr, warmup = 100, iter = 1000, chains = 4)
          prior = pr, warmup = 1000, iter = 10000, chains = 4)

### mean_trophic_richness and comm_mean_skew_ss are highly correlated unfortunately

static_model_table <- performance::compare_performance(m3,m31)
static_model_selection <- static_model_table |> 
      mutate(dWAIC = WAIC - min(WAIC))
write_csv(static_model_selection, "output/ms-second-round/tables/brm-static-selection-roundtwo.csv")
rm(m31)
### visualize and save model predictions ---
trophic_richness_a <- m3
pp_check(trophic_richness_a)
trophic_richness_re_slope <- mixedup::extract_random_coefs(trophic_richness_a)
trophic_richness_fe_slope <- mixedup::extract_fixed_effects(trophic_richness_a)

# species_richness_a <- m1
# pp_check(species_richness_a)
# species_richness_re_slope <- mixedup::extract_random_coefs(species_richness_a)
# species_richness_fe_slope <- mixedup::extract_fixed_effects(species_richness_a)
# 
# trophic_diversity_a <- m4
# pp_check(trophic_diversity_a)
# trophic_diversity_re_slope <- mixedup::extract_random_coefs(trophic_diversity_a)
# trophic_diversity_fe_slope <- mixedup::extract_fixed_effects(trophic_diversity_a)
# 
# species_diversity_a <- m2
# pp_check(species_diversity_a)
# species_diversity_re_slope <- mixedup::extract_random_coefs(species_diversity_a)
# species_diversity_fe_slope <- mixedup::extract_fixed_effects(species_diversity_a)

### trophic richness ---
trophic_richness_re1 <- ggpredict(trophic_richness_a, 
                              type = "re",
                              terms = c('mean_trophic_richness[-3:2 by=0.01]', 'Program')
                              )
trophic_richness_re <- as.data.frame(trophic_richness_re1)

trophic_richness_fe1 <- ggpredict(trophic_richness_a, 
                               type = "fe",
                               terms = c('mean_trophic_richness[-3:2 by=0.01]')
                               )

trophic_richness_fe <- as.data.frame(trophic_richness_fe1) |> 
      mutate(group = 'Overall')
      
trophic_richness_all <- rbind(trophic_richness_re, trophic_richness_fe)

trophic_richness_all |> 
      mutate(group = factor(group, levels = c("Overall", "FCE", "MCR", "PCCC", "PCCS", "SBC", "VCR"))) |> 
      ggplot(aes(x = x, y = predicted, color = group)) + 
      geom_smooth(method = "lm", linewidth = 1.5) + 
      labs(x = 'Trophic Group Richness',
           y = 'Aggregate Nitrogen Supply Stability',
           color = 'Program') +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            # axis.title.y = element_blank(),
            legend.position = "none",
            legend.background = element_blank(),
            legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

ggsave("output/ms-second-round/plots/brm-dsrmodel-trophic-richness.tiff", units = "in", width = 5,
       height = 5, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/brm-dsrmodel-trophic-richness.svg", units = "in", width = 5,
       height = 5, dpi =  600)

# ### species richness ---
# species_richness_re1 <- ggpredict(species_richness_a, 
#                                   type = "re",
#                                   terms = c('mean_species_richness[-3:2 by=0.01]', 'Program')
# )
# species_richness_re <- as.data.frame(species_richness_re1)
# 
# species_richness_fe1 <- ggpredict(species_richness_a, 
#                                   type = "fe",
#                                   terms = c('mean_species_richness[-3:2 by=0.01]')
# )
# 
# species_richness_fe <- as.data.frame(species_richness_fe1) |> 
#       mutate(group = 'Overall')
# 
# species_richness_all <- rbind(species_richness_re, species_richness_fe)
# 
# b <- species_richness_all |> 
#       mutate(group = factor(group, levels = c("Overall", "FCE", "MCR", "PCCC", "PCCS", "SBC", "VCR"))) |> 
#       ggplot(aes(x = x, y = predicted, color = group)) + 
#       geom_smooth(method = "lm", linewidth = 1.5) + 
#       labs(x = 'Scaled Species Richness',
#            y = 'Scaled Aggregate Nitrogen Supply Stability',
#            color = 'Program') +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
#             axis.text.y = element_text(face = "bold", color = "black", size = 12),
#             axis.title.x = element_text(face = "bold", color = "black", size = 14),
#             # axis.title.x = element_blank(),
#             axis.title.y = element_text(face = "bold", color = "black", size = 14),
#             # axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.background = element_blank(),
#             legend.key = element_rect(fill = 'white'),
#             legend.text = element_text(face = "bold", color = "black", size = 12),
#             legend.title = element_text(face = "bold", color = "black", size = 14))
# 
# ### trophic diversity ---
# trophic_diversity_re1 <- ggpredict(trophic_diversity_a, 
#                                   type = "re",
#                                   terms = c('mean_trophic_diversity[-3:2 by=0.01]', 'Program')
# )
# trophic_diversity_re <- as.data.frame(trophic_diversity_re1)
# 
# trophic_diversity_fe1 <- ggpredict(trophic_diversity_a, 
#                                   type = "fe",
#                                   terms = c('mean_trophic_diversity[-3:2 by=0.01]')
# )
# 
# trophic_diversity_fe <- as.data.frame(trophic_diversity_fe1) |> 
#       mutate(group = 'Overall')
# 
# trophic_diversity_all <- rbind(trophic_diversity_re, trophic_diversity_fe)
# 
# c <- trophic_diversity_all |> 
#       mutate(group = factor(group, levels = c("Overall", "FCE", "MCR", "PCCC", "PCCS", "SBC", "VCR"))) |> 
#       ggplot(aes(x = x, y = predicted, color = group)) + 
#       geom_smooth(method = "lm", linewidth = 1.5) + 
#       labs(x = 'Scaled Trophic Diversity',
#            y = 'Scaled Aggregate Nitrogen Supply Stability',
#            color = 'Program') +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
#             axis.text.y = element_text(face = "bold", color = "black", size = 12),
#             axis.title.x = element_text(face = "bold", color = "black", size = 14),
#             # axis.title.x = element_blank(),
#             axis.title.y = element_text(face = "bold", color = "black", size = 14),
#             # axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.background = element_blank(),
#             legend.key = element_rect(fill = 'white'),
#             legend.text = element_text(face = "bold", color = "black", size = 12),
#             legend.title = element_text(face = "bold", color = "black", size = 14))
# 
# ### species diversity ---
# species_diversity_re1 <- ggpredict(species_diversity_a, 
#                                    type = "re",
#                                    terms = c('mean_species_diversity[-3:2 by=0.01]', 'Program')
# )
# species_diversity_re <- as.data.frame(species_diversity_re1)
# 
# species_diversity_fe1 <- ggpredict(species_diversity_a, 
#                                    type = "fe",
#                                    terms = c('mean_species_diversity[-3:2 by=0.01]')
# )
# 
# species_diversity_fe <- as.data.frame(species_diversity_fe1) |> 
#       mutate(group = 'Overall')
# 
# species_diversity_all <- rbind(species_diversity_re, species_diversity_fe)
# 
# d <- species_diversity_all |> 
#       mutate(group = factor(group, levels = c("Overall", "FCE", "MCR", "PCCC", "PCCS", "SBC", "VCR"))) |> 
#       ggplot(aes(x = x, y = predicted, color = group)) + 
#       geom_smooth(method = "lm", linewidth = 1.5) + 
#       labs(x = 'Scaled Species Diversity',
#            y = 'Scaled Aggregate Nitrogen Supply Stability',
#            color = 'Program') +
#       theme_classic() +
#       scale_color_manual(values = program_palette) +
#       theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
#             axis.text.y = element_text(face = "bold", color = "black", size = 12),
#             axis.title.x = element_text(face = "bold", color = "black", size = 14),
#             # axis.title.x = element_blank(),
#             axis.title.y = element_text(face = "bold", color = "black", size = 14),
#             # axis.title.y = element_blank(),
#             legend.position = "right",
#             legend.background = element_blank(),
#             legend.key = element_rect(fill = 'white'),
#             legend.text = element_text(face = "bold", color = "black", size = 12),
#             legend.title = element_text(face = "bold", color = "black", size = 14))

### looking at slope coefficients for trophic richness model

div_fe <- trophic_richness_fe_slope |>
      rename(effect = term) |>
      filter(effect != 'Intercept') |>
      mutate(group = "Overall") |>
      select(group, effect, value, se, lower_2.5, upper_97.5)

div_re <- trophic_richness_re_slope |>
      filter(effect != "Intercept") |>
      select(group, effect, value, se, lower_2.5, upper_97.5)

div_slopes <- rbind(div_fe, div_re)
glimpse(div_slopes)

div_slopes |>
      mutate(group = factor(
            group,
            levels = c("VCR", "SBC", "PCCC", "FCE", "PCCS", "MCR", "Overall")
      )) |>
      filter(effect == "mean_trophic_richness") |>
      ggplot(aes(x = value, y = group, color = group)) +
      geom_point(size = 3) +
      geom_errorbarh(aes(xmin = lower_2.5, xmax = upper_97.5),
                     size = 1,
                     height = 0) +
      geom_vline(xintercept = 0, size = 1) +
      labs(x = 'Beta', y = 'Program') +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      theme(
            axis.text.x = element_text(
                  face = "bold",
                  color = "black",
                  size = 12
            ),
            axis.text.y = element_text(
                  face = "bold",
                  color = "black",
                  size = 12
            ),
            axis.title.x = element_text(
                  face = "bold",
                  color = "black",
                  size = 14
            ),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(
                  face = "bold",
                  color = "black",
                  size = 14
            ),
            # axis.title.y = element_blank(),
            legend.position = "none",
            legend.background = element_blank(),
            legend.key = element_rect(fill = 'white'),
            legend.text = element_text(
                  face = "bold",
                  color = "black",
                  size = 12
            ),
            legend.title = element_text(
                  face = "bold",
                  color = "black",
                  size = 14
            )
      )

ggsave("output/ms-second-round/plots/brm-dsrmodel-trophic-richness-slopecoeff.tiff", units = "in", width = 5,
       height = 5, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/brm-dsrmodel-trophic-richness-slopecoeff.svg", units = "in", width = 5,
       height = 5, dpi =  600)

### for the legend ---

div_slopes |>
      mutate(group = factor(
            group,
            levels = c("Overall", "FCE", "MCR", "PCCC", "PCCS", "SBC", "VCR")
      )) |>
      filter(effect == "mean_trophic_richness") |>
      ggplot(aes(x = value, y = group, color = group)) +
      # geom_point(size = 3) +
      geom_errorbarh(aes(xmin = lower_2.5, xmax = upper_97.5),
                     size = 1,
                     height = 0) +
      geom_vline(xintercept = 0, size = 1) +
      labs(x = 'Beta', y = 'Program') +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      theme(
            axis.text.x = element_text(
                  face = "bold",
                  color = "black",
                  size = 12
            ),
            axis.text.y = element_text(
                  face = "bold",
                  color = "black",
                  size = 12
            ),
            axis.title.x = element_text(
                  face = "bold",
                  color = "black",
                  size = 14
            ),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(
                  face = "bold",
                  color = "black",
                  size = 14
            ),
            # axis.title.y = element_blank(),
            legend.position = "right",
            legend.background = element_blank(),
            legend.key = element_rect(fill = 'white'),
            legend.text = element_text(
                  face = "bold",
                  color = "black",
                  size = 12
            ),
            legend.title = element_text(
                  face = "bold",
                  color = "black",
                  size = 14
            )
      )

ggsave("output/ms-second-round/plots/brm-dsrmodel-legend.tiff", units = "in", width = 5,
       height = 5, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/brm-dsrmodel-legend.svg", units = "in", width = 5,
       height = 5, dpi =  600)
