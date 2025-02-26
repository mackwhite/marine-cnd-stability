###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): Mack White, W. Ryan James
###goal(s): visualize modeled relationships
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

### read in necessary data ----
full_model = readRDS('local_data/rds-full-model.rds')

#summary stats -----
post = posterior_samples(full_model)

mean(post$`r_Program[MCR,beta_time]` + post$b_beta_time < 0)
mean(post$`r_Program[MCR,mean_species_diversity]` + post$b_mean_species_diversity < 0)

mean(post$`r_Program[PCCS,troph_beta_time]`+ post$b_troph_beta_time < 0)
mean(post$`r_Program[PCCS,beta_time]` + post$b_beta_time > 0)
mean(post$`r_Program[PCCS,mean_species_diversity]` + post$b_mean_species_diversity > 0)

mean(post$`r_Program[FCE,beta_time]` + post$b_beta_time > 0)
mean(post$`r_Program[FCE,mean_species_diversity]` + post$b_mean_species_diversity > 0)

mean(post$`r_Program[PCCC,synch]`+ post$b_synch < 0)
mean(post$`r_Program[PCCS,troph_beta_time]`+ post$b_troph_beta_time < 0)
mean(post$`r_Program[PCCC,beta_time]` + post$b_beta_time < 0)
mean(post$`r_Program[PCCC,mean_species_diversity]` + post$b_mean_species_diversity > 0)

mean(post$`r_Program[SBC,synch]`+ post$b_synch < 0)
mean(post$`r_Program[SBC,troph_beta_time]`+ post$b_troph_beta_time < 0)
mean(post$`r_Program[SBC,beta_time]` + post$b_beta_time < 0)
mean(post$`r_Program[SBC,mean_species_diversity]` + post$b_mean_species_diversity > 0)


mean(post$`r_Program[VCR,synch]`+ post$b_synch < 0)
mean(post$`r_Program[VCR,troph_beta_time]`+ post$b_troph_beta_time < 0)
mean(post$`r_Program[VCR,beta_time]` + post$b_beta_time < 0)
mean(post$`r_Program[VCR,mean_species_diversity]` + post$b_mean_species_diversity > 0)


# 
# order of sites 
stab = read_csv('local_data/dsr-eco-org-raw-all.csv') |> 
      rename(Program = program) |> 
      group_by(Program) |> 
      summarize(mn = mean(comm_n_stability),
                med = median(comm_n_stability)) |> 
      arrange(mn)

prog = c('VCR', 'SBC', 'PCCC', 'FCE', 'PCCS', 'MCR', 'Overall')

### set color schemes ----
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

# 
# pp_check(full_model)
# random effects
re95 = mixedup::extract_random_coefs(full_model, ci_level = c(0.95))
re80 = mixedup::extract_random_coefs(full_model, ci_level = c(0.8))
re90 = mixedup::extract_random_coefs(full_model, ci_level = c(0.9))
re50 = mixedup::extract_random_coefs(full_model, ci_level = c(0.5))

re_beta = left_join(re95, re80) |> 
      left_join(re90) |> 
      left_join(re50)|> 
      rename(term = effect,
             Program = group)

# fixed effects
fe95 = mixedup::extract_fixed_effects(full_model, ci_level = c(0.95))
fe80 = mixedup::extract_fixed_effects(full_model, ci_level = c(0.8))
fe90 = mixedup::extract_fixed_effects(full_model, ci_level = c(0.9))
fe50 = mixedup::extract_fixed_effects(full_model, ci_level = c(0.5))

fe_beta = left_join(fe95, fe80) |> 
      left_join(fe90) |> 
      left_join(fe50) |> 
      mutate(Program = 'Overall') 

# make data frame of all betas
df_beta = bind_rows(re_beta, fe_beta) |> 
      filter(term != 'Intercept') |> 
      mutate(Program = factor(Program, levels = prog),
             term = factor(term, levels = c('synch',
                                            'troph_beta_time',
                                            'mean_species_diversity',
                                            'beta_time'),
                           labels = c('Species Synchrony',
                                      'Trophic Turnover',
                                      'Species Evenness',
                                      'Species Turnover')))

# make equation data set
df_eq = bind_rows(re_beta, fe_beta) |> 
      select(term, Program, value) |> 
      pivot_wider(names_from = term, values_from = value) |> 
      # pivot_longer(mean_species_diversity:synch, names_to = 'term', values_to = 'beta')
      pivot_longer(mean_species_diversity:synch, names_to = 'term', values_to = 'beta')


# load and get min and max values
### read in necessary data ---
ov = read_csv('local_data/dsr-eco-org-raw-all.csv') |>
      # rename(Program = program,
      #       Trophic_Group = troph_group,
      #       Species = scientific_name,
      #       Habitat = habitat,
      #       Site = site) |>
      rename(Program = program) |> 
      select(Program,
             mean_species_diversity,
             beta_time,
             synch,
             troph_beta_time) |>
      distinct() |> 
      pivot_longer(mean_species_diversity:troph_beta_time, 
                   names_to = 'term', values_to = 'value') |> 
      group_by(Program,term) |> 
      mutate(scaled = scale(value)) |> 
      group_by(term) |> 
      slice(c(which.min(value), which.max(value))) |> 
      mutate(Program = 'Overall')

dat = read_csv('local_data/dsr-eco-org-raw-all.csv') |>
      # rename(Program = program,
      #       Trophic_Group = troph_group,
      #       Species = scientific_name,
      #       Habitat = habitat,
      #       Site = site) |>
      rename(Program = program) |> 
      select(Program,
             mean_species_diversity,
             beta_time,
             synch,
             troph_beta_time) |>
      distinct() |> 
      pivot_longer(mean_species_diversity:troph_beta_time, 
                   names_to = 'term', values_to = 'value') |> 
      group_by(Program,term) |> 
      mutate(scaled = scale(value)) |>  
      slice(c(which.min(value), which.max(value))) |> 
      ungroup() |> 
      bind_rows(ov)

dat_scaled = dat |>
      group_by(Program, term) |> 
      mutate(i = row_number()) |> 
      select(Program, term, scaled)

df = dat_scaled |> 
      left_join(df_eq) |> 
      mutate(stab = beta*scaled + Intercept) |> 
      left_join(dat) |> 
      mutate(Program = factor(Program, levels = prog),
             term = factor(term, levels = c('synch',
                                            'troph_beta_time',
                                            'mean_species_diversity',
                                            'beta_time'),
                           labels = c('Species Synchrony',
                                      'Trophic Turnover',
                                      'Species Evenness',
                                      'Species Turnover')))
      
a = ggplot(df|> filter(Program != 'Overall'), aes(value, stab, color = Program))+
      geom_line(linewidth = 1.75)+
      scale_color_manual(values = program_palette)+
      # facet_wrap(~term, nrow = 1, strip.position = 'right', scales = 'free')+
      facet_wrap(~term, ncol = 1, strip.position = 'right', scales = 'free')+
      labs(y = 'CND Stability', x = NULL)+
      theme_classic()+
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_text(face = "bold", color = "black", size = 12),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 12),
            #strip.text = element_text(face = "bold", color = "black", size = 12),
            strip.text = element_blank(),
            strip.background = element_blank(),
            legend.position = "bottom",
            #legend.background = element_blank(),
            # legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

# betas 
b = ggplot(df_beta|> filter(Program != 'Overall'), aes(Program, value, color = Program))+
      geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 0.75) +
      #geom_pointrange(aes(ymin = lower_5, ymax = upper_95), linewidth = 1.25)+
      geom_pointrange(aes(ymin = lower_10, ymax = upper_90), linewidth = 2)+
      # geom_pointrange(aes(ymin = lower_25, ymax = upper_75), linewidth = 2)+
      geom_pointrange(aes(ymin = lower_2.5, ymax = upper_97.5), linewidth = 1, size = .9)+
      labs(y = 'Beta', x = NULL)+
      scale_color_manual(values = program_palette)+
      coord_flip()+
      # facet_wrap(~term, nrow = 1, strip.position = 'top')+
      facet_wrap(~term, ncol = 1, strip.position = 'right')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
            axis.text.y = element_text(face = "bold", color = "black", size = 10),
            axis.title.x = element_text(face = "bold", color = "black", size = 12),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 12),
            strip.text = element_text(face = "bold", color = "black", size = 10),
            # axis.title.y = element_blank(),
            legend.position = "none",
            legend.background = element_blank(),
            #legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

# ggpubr::ggarrange(b,a, nrow=2, ncol=1, align = 'hv', common.legend = T, legend = 'bottom')
#ggsave('output/figs/fig4.png', dpi = 600, units= 'in', height = 10, width = 6)

plot = ggpubr::ggarrange(a,b,labels = c('a', 'b'), align = 'h', legend = 'none', label.x = -0.01)

f4 = annotate_figure(plot, 
                     top = text_grob(
                           bquote({R^2}[cond] == 0.67 ~ "," ~ {R^2}[mar] == 0.18),
                           #expression("R"^2*[cond] *"= 0.62," ~ "R"^2 * [mar]* "= 0.19"),
                           size = 14, face = 'bold'))

ggsave('output/figs/fig4_rev.png', plot = f4, dpi = 600, units= 'in', height = 7.5, width = 5.75)
