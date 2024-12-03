###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): MW
###goal(s): wrangling raw cnd data for summary data
###date(s): november 2024 @ nceas
###note(s): 

# Housekeeping ------------------------------------------------------------

### load necessary libraries
# install.packages("librarian")
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

### read in model data ----
synch = readRDS("local_data/rds-single-synchrony.rds")
rich = readRDS("local_data/rds-single-richness.rds")

# stats ----
p_synch = posterior_samples(synch)

mean(p_synch$`r_Program[PCCC,synch]` + p_synch$b_synch < 0)
mean(p_synch$`r_Program[SBC,synch]` + p_synch$b_synch < 0)
mean(p_synch$`r_Program[VCR,synch]` + p_synch$b_synch < 0)

p_sr = posterior_samples(rich)

mean(p_sr$`r_Program[PCCC,mean_species_richness]` + p_sr$b_mean_species_richness > 0)

mean(p_sr$`r_Program[SBC,mean_species_richness]` + p_sr$b_mean_species_richness > 0)

### set color scheme ---
program_palette = c("Overall"="#000000",
                     "FCE"="#64a988",
                     "MCR"="#ff967d",
                     'PCCC'="#2A788EFF",
                     "PCCS"="#8b6b93",
                     'SBC'='#ff3f4c',
                     "VCR"="#9b9254")


prog = c('VCR', 'SBC', 'PCCC', 'FCE', 'PCCS', 'MCR', 'Overall')

# species richness ----
# random effects
re95 = mixedup::extract_random_coefs(rich, ci_level = c(0.95))
re80 = mixedup::extract_random_coefs(rich, ci_level = c(0.8))

re_beta = left_join(re95, re80) |> 
      rename(term = effect,
             Program = group)

# fixed effects
fe95 = mixedup::extract_fixed_effects(rich, ci_level = c(0.95))
fe80 = mixedup::extract_fixed_effects(rich, ci_level = c(0.8))


fe_beta = left_join(fe95, fe80) |> 
      mutate(Program = 'Overall') 

# make data frame of all betas
df_beta = bind_rows(re_beta, fe_beta) |> 
      filter(term != 'Intercept') |> 
      mutate(Program = factor(Program, levels = prog))


# make equation data set
df_eq = bind_rows(re_beta, fe_beta) |> 
      select(term, Program, value) |> 
      pivot_wider(names_from = term, values_from = value)  |> 
      rename(beta = mean_species_richness)

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
             value = mean_species_richness) |>
      distinct() |> 
      group_by(Program) |> 
      mutate(scaled = scale(value)) |> 
      ungroup() |> 
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
             value = mean_species_richness) |>
      distinct() |> 
      group_by(Program) |> 
      mutate(scaled = scale(value)) |>  
      slice(c(which.min(value), which.max(value))) |> 
      ungroup() |> 
      bind_rows(ov)

dat_scaled = dat |>
      group_by(Program, scaled) |> 
      mutate(i = row_number()) |> 
      select(Program, scaled)

df = dat_scaled |> 
      left_join(df_eq) |> 
      mutate(stab = beta*scaled + Intercept) |> 
      left_join(dat) |> 
      mutate(Program = factor(Program, levels = prog))

raw = read_csv('local_data/dsr-eco-org-raw-all.csv') |>
      # rename(Program = program,
      #       Trophic_Group = troph_group,
      #       Species = scientific_name,
      #       Habitat = habitat,
      #       Site = site) |>
      rename(Program = program) |> 
      select(Program,comm_n_stability, 
             value = mean_species_richness) |> 
      mutate(stab = scale(comm_n_stability))


a = ggplot(df, aes(value, stab, color = Program))+
      geom_point(data = raw, aes(value, stab, color = Program), size = 2)+
      geom_line(linewidth = 1.75)+
      scale_color_manual(values = program_palette)+
      #facet_wrap(~term, ncol = 1, strip.position = 'right', scales = 'free')+
      labs(y = 'Aggregate nitrogen \nsupply stability', x = 'Species richness')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            #strip.text = element_text(face = "bold", color = "black", size = 12),
            strip.text = element_blank(),
            strip.background = element_blank(),
            legend.position = "none",
            #legend.background = element_blank(),
            # legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

# betas 
b = ggplot(df_beta, aes(Program, value, color = Program))+
      geom_hline(aes(yintercept = 0), size = 1)+
      #geom_pointrange(aes(ymin = lower_5, ymax = upper_95), linewidth = 1.25)+
      geom_pointrange(aes(ymin = lower_10, ymax = upper_90), linewidth = 2)+
      # geom_pointrange(aes(ymin = lower_25, ymax = upper_75), linewidth = 2)+
      geom_pointrange(aes(ymin = lower_2.5, ymax = upper_97.5), linewidth = 1, size = .9)+
      labs(y = 'Beta', x = NULL)+
      scale_color_manual(values = program_palette)+
      scale_y_continuous(limits = c(-0.26, 1.03))+
      coord_flip()+
      #facet_wrap(~term, ncol = 1, strip.position = 'right')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            strip.text = element_text(face = "bold", color = "black", size = 12),
            # axis.title.y = element_blank(),
            legend.position = "none",
            legend.background = element_blank(),
            #legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

ggpubr::ggarrange(a,b, align = 'h')

sr = ggpubr::ggarrange(a,b, align = 'h', legend = 'none')

# f4 = annotate_figure(plot, 
#                      top = text_grob(
#                            bquote({R^2}[cond] == 0.62 ~ "," ~ {R^2}[mar] == 0.19),
#                            #expression("R"^2*[cond] *"= 0.62," ~ "R"^2 * [mar]* "= 0.19"),
#                            size = 16, face = 'bold'))

# synchrony 
# random effects
re95 = mixedup::extract_random_coefs(synch, ci_level = c(0.95))
re80 = mixedup::extract_random_coefs(synch, ci_level = c(0.8))

re_beta = left_join(re95, re80) |> 
      rename(term = effect,
             Program = group)

# fixed effects
fe95 = mixedup::extract_fixed_effects(synch, ci_level = c(0.95))
fe80 = mixedup::extract_fixed_effects(synch, ci_level = c(0.8))


fe_beta = left_join(fe95, fe80) |> 
      mutate(Program = 'Overall') 

# make data frame of all betas
df_beta = bind_rows(re_beta, fe_beta) |> 
      filter(term != 'Intercept') |> 
      mutate(Program = factor(Program, levels = prog))


# make equation data set
df_eq = bind_rows(re_beta, fe_beta) |> 
      select(term, Program, value) |> 
      pivot_wider(names_from = term, values_from = value)  |> 
      rename(beta = synch)

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
             value = synch) |>
      distinct() |> 
      group_by(Program) |> 
      mutate(scaled = scale(value)) |> 
      ungroup() |> 
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
             value = synch) |>
      distinct() |> 
      group_by(Program) |> 
      mutate(scaled = scale(value)) |>  
      slice(c(which.min(value), which.max(value))) |> 
      ungroup() |> 
      bind_rows(ov)

dat_scaled = dat |>
      group_by(Program, scaled) |> 
      mutate(i = row_number()) |> 
      select(Program, scaled)

df = dat_scaled |> 
      left_join(df_eq) |> 
      mutate(stab = beta*scaled + Intercept) |> 
      left_join(dat) |> 
      mutate(Program = factor(Program, levels = prog))

raw = read_csv('local_data/dsr-eco-org-raw-all.csv') |>
      # rename(Program = program,
      #       Trophic_Group = troph_group,
      #       Species = scientific_name,
      #       Habitat = habitat,
      #       Site = site) |>
      rename(Program = program) |> 
      select(Program,comm_n_stability, 
             value = synch) |> 
      mutate(stab = scale(comm_n_stability))


a = ggplot(df, aes(value, stab, color = Program))+
      geom_point(data = raw, aes(value, stab, color = Program), size = 2)+
      geom_line(linewidth = 1.75)+
      scale_color_manual(values = program_palette)+
      #facet_wrap(~term, ncol = 1, strip.position = 'right', scales = 'free')+
      labs(y = 'Aggregate nitrogen \nsupply stability', x = 'Species synchrony')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            #strip.text = element_text(face = "bold", color = "black", size = 12),
            strip.text = element_blank(),
            strip.background = element_blank(),
            legend.position = "none",
            #legend.background = element_blank(),
            # legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

# betas 
b = ggplot(df_beta, aes(Program, value, color = Program))+
      geom_hline(aes(yintercept = 0), size = 1)+
      #geom_pointrange(aes(ymin = lower_5, ymax = upper_95), linewidth = 1.25)+
      geom_pointrange(aes(ymin = lower_10, ymax = upper_90), linewidth = 2)+
      # geom_pointrange(aes(ymin = lower_25, ymax = upper_75), linewidth = 2)+
      geom_pointrange(aes(ymin = lower_2.5, ymax = upper_97.5), linewidth = 1, size = .9)+
      labs(y = 'Beta', x = NULL)+
      scale_color_manual(values = program_palette)+
      coord_flip()+
      #facet_wrap(~term, ncol = 1, strip.position = 'right')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            strip.text = element_text(face = "bold", color = "black", size = 12),
            # axis.title.y = element_blank(),
            legend.position = "none",
            legend.background = element_blank(),
            #legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

syn = ggpubr::ggarrange(a,b, align = 'h', legend = 'none')


ggarrange(sr, syn, labels = c('a)', 'b)'), align = 'v', nrow =2)

# ggsave('output/figs/fig4.png', dpi = 600, units= 'in', height = 6, width = 6)

# removing overall average ----
# species richness 
# random effects
re95 = mixedup::extract_random_coefs(rich, ci_level = c(0.95))
re80 = mixedup::extract_random_coefs(rich, ci_level = c(0.8))

re_beta = left_join(re95, re80) |> 
      rename(term = effect,
             Program = group)

# fixed effects
fe95 = mixedup::extract_fixed_effects(rich, ci_level = c(0.95))
fe80 = mixedup::extract_fixed_effects(rich, ci_level = c(0.8))


fe_beta = left_join(fe95, fe80) |> 
      mutate(Program = 'Overall') 

# make data frame of all betas
df_beta = bind_rows(re_beta, fe_beta) |> 
      filter(term != 'Intercept') |> 
      mutate(Program = factor(Program, levels = prog))


# make equation data set
df_eq = bind_rows(re_beta, fe_beta) |> 
      select(term, Program, value) |> 
      pivot_wider(names_from = term, values_from = value)  |> 
      rename(beta = mean_species_richness)

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
             value = mean_species_richness) |>
      distinct() |> 
      group_by(Program) |> 
      mutate(scaled = scale(value)) |> 
      ungroup() |> 
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
             value = mean_species_richness) |>
      distinct() |> 
      group_by(Program) |> 
      mutate(scaled = scale(value)) |>  
      slice(c(which.min(value), which.max(value))) |> 
      ungroup() |> 
      bind_rows(ov)

dat_scaled = dat |>
      group_by(Program, scaled) |> 
      mutate(i = row_number()) |> 
      select(Program, scaled)

raw = read_csv('local_data/dsr-eco-org-raw-all.csv') |>
      # rename(Program = program,
      #       Trophic_Group = troph_group,
      #       Species = scientific_name,
      #       Habitat = habitat,
      #       Site = site) |>
      rename(Program = program) |> 
      select(Program, 
             value = mean_species_richness,
             stab = comm_n_stability) 

df = dat_scaled |> 
      left_join(df_eq) |> 
      mutate(pred = beta*scaled + Intercept,
             stab = pred*sd(raw$stab) + mean(raw$stab)) |> 
      left_join(dat) |> 
      mutate(Program = factor(Program, levels = prog))


a = ggplot(df |> filter(Program != 'Overall'), aes(value, stab, color = Program))+
      geom_point(data = raw, aes(value, stab, color = Program), size = 2)+
      geom_line(linewidth = 1.75)+
      scale_color_manual(values = program_palette)+
      #facet_wrap(~term, ncol = 1, strip.position = 'right', scales = 'free')+
      labs(y = 'CND Stability', x = 'Species Richness')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            #strip.text = element_text(face = "bold", color = "black", size = 12),
            strip.text = element_blank(),
            strip.background = element_blank(),
            legend.position = "none",
            #legend.background = element_blank(),
            # legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

# betas 
b = ggplot(df_beta|> filter(Program != 'Overall'), aes(Program, value, color = Program))+
      geom_hline(aes(yintercept = 0), linewidth = 1)+
      #geom_pointrange(aes(ymin = lower_5, ymax = upper_95), linewidth = 1.25)+
      geom_pointrange(aes(ymin = lower_10, ymax = upper_90), linewidth = 2)+
      # geom_pointrange(aes(ymin = lower_25, ymax = upper_75), linewidth = 2)+
      geom_pointrange(aes(ymin = lower_2.5, ymax = upper_97.5), linewidth = 1, size = .9)+
      labs(y = 'Beta', x = NULL)+
      scale_color_manual(values = program_palette)+
      scale_y_continuous(limits = c(-0.26, 1.03))+
      annotate('text', x = 2, y = 0.8, size = 4.5,
               label = bquote(atop(
                     {R^2}[cond] ~ '= 0.54',
                     {R^2}[mar] ~ '= 0.04 ')))+
      coord_flip()+
      #facet_wrap(~term, ncol = 1, strip.position = 'right')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            strip.text = element_text(face = "bold", color = "black", size = 12),
            # axis.title.y = element_blank(),
            legend.position = "none",
            legend.background = element_blank(),
            #legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

ggpubr::ggarrange(a,b, align = 'h')

sr = ggpubr::ggarrange(a,b, align = 'h', legend = 'none')

# f4 = annotate_figure(plot, 
#                      top = text_grob(
#                            bquote({R^2}[cond] == 0.62 ~ "," ~ {R^2}[mar] == 0.19),
#                            #expression("R"^2*[cond] *"= 0.62," ~ "R"^2 * [mar]* "= 0.19"),
#                            size = 16, face = 'bold'))

# synchrony 
# random effects
re95 = mixedup::extract_random_coefs(synch, ci_level = c(0.95))
re80 = mixedup::extract_random_coefs(synch, ci_level = c(0.8))

re_beta = left_join(re95, re80) |> 
      rename(term = effect,
             Program = group)

# fixed effects
fe95 = mixedup::extract_fixed_effects(synch, ci_level = c(0.95))
fe80 = mixedup::extract_fixed_effects(synch, ci_level = c(0.8))


fe_beta = left_join(fe95, fe80) |> 
      mutate(Program = 'Overall') 

# make data frame of all betas
df_beta = bind_rows(re_beta, fe_beta) |> 
      filter(term != 'Intercept') |> 
      mutate(Program = factor(Program, levels = prog))


# make equation data set
df_eq = bind_rows(re_beta, fe_beta) |> 
      select(term, Program, value) |> 
      pivot_wider(names_from = term, values_from = value)  |> 
      rename(beta = synch)

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
             value = synch) |>
      distinct() |> 
      group_by(Program) |> 
      mutate(scaled = scale(value)) |> 
      ungroup() |> 
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
             value = synch) |>
      distinct() |> 
      group_by(Program) |> 
      mutate(scaled = scale(value)) |>  
      slice(c(which.min(value), which.max(value))) |> 
      ungroup() |> 
      bind_rows(ov)

dat_scaled = dat |>
      group_by(Program, scaled) |> 
      mutate(i = row_number()) |> 
      select(Program, scaled)

raw = read_csv('local_data/dsr-eco-org-raw-all.csv') |>
      # rename(Program = program,
      #       Trophic_Group = troph_group,
      #       Species = scientific_name,
      #       Habitat = habitat,
      #       Site = site) |>
      rename(Program = program) |> 
      select(Program, 
             value = synch,
             stab = comm_n_stability) 

df = dat_scaled |> 
      left_join(df_eq) |> 
      mutate(pred = beta*scaled + Intercept,
             stab = pred*sd(raw$stab) + mean(raw$stab)) |> 
      left_join(dat) |> 
      mutate(Program = factor(Program, levels = prog))


a = ggplot(df |> filter(Program != 'Overall'), aes(value, stab, color = Program))+
      geom_point(data = raw, aes(value, stab, color = Program), size = 2)+
      geom_line(linewidth = 1.75)+
      scale_color_manual(values = program_palette)+
      #facet_wrap(~term, ncol = 1, strip.position = 'right', scales = 'free')+
      labs(y = 'CND Stability', x = 'Species Synchrony')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            #strip.text = element_text(face = "bold", color = "black", size = 12),
            strip.text = element_blank(),
            strip.background = element_blank(),
            legend.position = "none",
            #legend.background = element_blank(),
            # legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

# betas 
b = ggplot(df_beta|> filter(Program != 'Overall'), aes(Program, value, color = Program))+
      geom_hline(aes(yintercept = 0), linewidth = 1)+
      #geom_pointrange(aes(ymin = lower_5, ymax = upper_95), linewidth = 1.25)+
      geom_pointrange(aes(ymin = lower_10, ymax = upper_90), linewidth = 2)+
      # geom_pointrange(aes(ymin = lower_25, ymax = upper_75), linewidth = 2)+
      geom_pointrange(aes(ymin = lower_2.5, ymax = upper_97.5), linewidth = 1, size = .9)+
      labs(y = 'Beta', x = NULL)+
      scale_color_manual(values = program_palette)+
      annotate('text', x = 2, y = -0.8, size = 4.5,
               label = bquote(atop(
                     {R^2}[cond] ~ '= 0.61',
                     {R^2}[mar] ~ '= 0.10 ')))+
      coord_flip()+
      #facet_wrap(~term, ncol = 1, strip.position = 'right')+
      theme_classic()+
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 12),
            axis.text.y = element_text(face = "bold", color = "black", size = 12),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            # axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            strip.text = element_text(face = "bold", color = "black", size = 12),
            # axis.title.y = element_blank(),
            legend.position = "none",
            legend.background = element_blank(),
            #legend.key = element_rect(fill = 'white'),
            legend.text = element_text(face = "bold", color = "black", size = 12),
            legend.title = element_text(face = "bold", color = "black", size = 14))

syn = ggpubr::ggarrange(a,b, align = 'h', legend = 'none')


ggarrange(sr, syn, labels = c('a)', 'b)'), align = 'v', nrow =2)

ggsave('output/figs/fig4_noavg.png', dpi = 600, units= 'in', height = 6, width = 6.5)
