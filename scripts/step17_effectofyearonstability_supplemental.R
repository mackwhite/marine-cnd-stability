ann_dt <- read_csv("local_data/annual-dt-for-summary.csv")

program_palette <- c("Overall"="#000000", 
                     "FCE"="#64a988", 
                     "MCR"="#ff967d", 
                     'PCCC'="#2A788EFF", 
                     "PCCS"="#8b6b93",
                     'SBC'='#ff3f4c', 
                     "VCR"="#9b9254")

summary <- ann_dt |> 
      group_by(program, site) |> 
      summarize(years = n_distinct(year)) |> 
      rename(Program = program,
             Site = site)

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

dat_summary <- dat |> 
      select(Site, comm_n_stability)

all <- left_join(dat_summary, summary, by = "Site")  

model <- lm(comm_n_stability ~ years, data = all)
summary(model)$r.squared 
r2 <- summary(model)$r.squared

summary(model)

a <- all |>
      # filter(Program != "VCR") |> 
      ggplot(aes(x = years, y = comm_n_stability)) +
      geom_point(aes(color = Program), size = 2) +  # Adds the scatter plot points
      geom_smooth(method = "lm", size = 2, color = "black", linetype = "solid", se = FALSE) +
      # geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
      labs(x = "Years",
           y = "CND Stability") +
      theme_classic() +
      annotate('text',
               x = 7, y = 5,
               label = bquote({R^2} == 0.04),
               size = 5) +
      scale_color_manual(values = program_palette) +
      theme(axis.text.x = element_text(face = "bold", color = "black"),
            axis.text.y = element_text(face = "bold", color = "black"),
            axis.title.x = element_text(face = "bold", color = "black"),
            axis.title.y = element_text(face = "bold", color = "black"),
            legend.position = "right",
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"))

all_short <- all |> filter(Program != "VCR")

model <- lm(comm_n_stability ~ years, data = all_short)
summary(model)$r.squared 
r2 <- summary(model)$r.squared

b <- all |>
      filter(Program != "VCR") |>
      ggplot(aes(x = years, y = comm_n_stability)) +
      geom_point(aes(color = Program), size = 2) +  # Adds the scatter plot points
      geom_smooth(method = "lm", size = 2, color = "black", linetype = "solid", se = FALSE) +
      # geom_smooth(method = "lm", se = FALSE) +  # Adds linear model lines for each program
      labs(x = "Years",
           y = "CND Stability") +
      theme_classic() +
      annotate('text',
               x = 12.7, y = 5,
               label = bquote({R^2} == 0.004),
               size = 5) +
      scale_color_manual(values = program_palette) +
      theme(axis.text.x = element_text(face = "bold", color = "black"),
            axis.text.y = element_text(face = "bold", color = "black"),
            axis.title.x = element_text(face = "bold", color = "black"),
            axis.title.y = element_text(face = "bold", color = "black"),
            legend.position = "right",
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"))

ggarrange(a, b,
          labels = c('a)','b)'),
          ncol = 2, vjust = 1, align = "h",
          common.legend = TRUE, legend = 'right')

ggsave("output/figs/supplemental-effect-of-year-on-stability.png", units = "in", width = 10,
       height = 5, dpi =  600)
