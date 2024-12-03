###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): Mack White
###goal(s): create time series of community nutrient supply, jitterbox plot, and map
###date(s): November 2024
###note(s): 

###########################################################################
# Housekeeping ------------------------------------------------------------
###########################################################################

### load necessary libraries
### install.packages("librarian")
# remotes::install_github('m-clark/mixedup')
librarian::shelf(tidyverse, readxl, MuMIn, sjPlot, lme4, corrplot, 
                 performance, ggeffects, ggpubr, ggridges, parameters, ggstats, brms, mixedup, lterpalettefinder)

### read in necessary data ---
ann_dt <- read_csv("local_data/annual-dt-for-summary.csv")
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


# map insets --------------------------------------------------------------

ann_dt |> 
      filter(program == "SBC") |> 
      ggplot(aes(x = year, y = total_nitrogen_ann, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-2~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_y_continuous(breaks = c(0,2000,4000,6000,8000,10000,12000,14000)) +
      scale_x_continuous(breaks = c(2000,2005,2010,2015,2020)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/figs/map_insets/sbc-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

ann_dt |> 
      filter(program == "MCR") |> 
      ggplot(aes(x = year, y = total_nitrogen_ann, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-2~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_y_continuous(breaks = c(0,1000,2000,3000,4000,5000))+
      scale_x_continuous(breaks = c(2006,2010,2014,2018,2022)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/figs/map_insets/mcr-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

ann_dt |> 
      filter(program == "VCR") |> 
      ggplot(aes(x = year, y = total_nitrogen_ann, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-2~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_y_continuous(breaks = c(0,250,500,750,1000,1250))+
      scale_x_continuous(breaks = c(2012,2014,2016,2018)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/figs/map_insets/vcr-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

ann_dt|> 
      filter(program == "FCE") |> 
      ggplot(aes(x = year, y = total_nitrogen_ann, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-1~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_y_continuous(breaks = c(0,1500,3000,4500,6000,7500,9000,10500,12000,13500))+
      scale_x_continuous(breaks = c(2005,2008,2011,2014,2017,2020,2023)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/figs/map_insets/fce-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

ann_dt |> 
      filter(program == "PCCC") |> 
      ggplot(aes(x = year, y = total_nitrogen_ann, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-2~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_y_continuous(breaks = c(0,2000,4000,6000,8000)) +
      scale_x_continuous(breaks = c(2000,2005,2010,2015,2020)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/figs/map_insets/pccc-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

ann_dt |> 
      filter(program == "PCCS"&total_nitrogen_ann <15000) |> 
      ggplot(aes(x = year, y = total_nitrogen_ann, group = site, color = program)) +
      geom_line(alpha = 0.8, linewidth = 1) +
      labs(x = 'Year',
           y = expression(bold('Aggregate Nitrogen Supply (μg '*~m^-2~""~hr^-1*')'))) +
      theme_classic() +
      scale_color_manual(values = program_palette) +
      scale_y_continuous(breaks = c(0,3000,6000,9000,12000))+
      scale_x_continuous(breaks = c(2000,2005,2010,2015,2020)) +
      theme(
            axis.text = element_text(face = "bold", size = 12, color = "black"),
            axis.title.y = element_text(face = "bold", size = 14, color = "black"),
            axis.title.x = element_blank(),
            axis.line = element_line("black"),
            legend.position = "none",
            legend.text = element_text(face = "bold", size = 14, color = "black"),
            legend.title = element_text(face = "bold", size = 14, color = "black"),
            panel.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 12, color = "black"))

ggsave("output/figs/map_insets/pccs-timeseries.png", units = "in", width = 5,
       height = 5, dpi =  600)

# density plots -----------------------------------------------------------

### nitrogen supply ---
ann_dt |> 
      filter(total_nitrogen_ann < 25000) |>
      filter(total_nitrogen_ann > 0) |>
      ggplot(aes(x = log1p(total_nitrogen_ann), color = program, fill = program, group = program)) +
      geom_density(aes(color = program), size = 1.5, alpha = 0.3) + 
      geom_density(aes(fill = program), alpha = 0.3) +
      scale_color_manual(values = program_palette) + 
      scale_fill_manual(values = program_palette) + 
      labs(title = "Aggregate Nitrogen Supply") +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line("black"),
            axis.text.x = element_text(color = "black", face = "bold", size = 18),
            axis.text.y = element_text(color = "black", face = "bold", size = 18),
            # axis.text.y = element_text(color = "black", face = "bold", size = 18),
            plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
            legend.position = "none")

ggsave("output/ms-second-round/plots/figure2/nsupply.tiff", units = "in", width = 12,
       height = 3, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/figure2/nsupply.svg", units = "in", width = 12,
       height = 3, dpi =  600)

### biomass ---
ann_dt |> 
      filter(total_nitrogen_ann < 25000) |>
      filter(total_nitrogen_ann > 0) |>
      ggplot(aes(x = log1p(total_biomass_ann), color = program, fill = program, group = program)) +
      geom_density(aes(color = program), size = 1.5, alpha = 0.3) + 
      geom_density(aes(fill = program), alpha = 0.3) +
      scale_color_manual(values = program_palette) + 
      scale_fill_manual(values = program_palette) + 
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
      labs(title = "Aggregate Biomass") +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line("black"),
            axis.text.x = element_text(color = "black", face = "bold", size = 18),
            axis.text.y = element_text(color = "black", face = "bold", size = 18),
            # axis.text.y = element_text(color = "black", face = "bold", size = 18),
            plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
            legend.position = "none")

ggsave("output/ms-second-round/plots/figure2/biomass.tiff", units = "in", width = 12,
       height = 3, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/figure2/biomass.svg", units = "in", width = 12,
       height = 3, dpi =  600)

### max size ---
ann_dt |> 
      filter(total_nitrogen_ann < 25000) |>
      filter(total_nitrogen_ann > 0) |>
      ggplot(aes(x = log1p(max_size_ann), color = program, fill = program, group = program)) +
      geom_density(aes(color = program), size = 1.5, alpha = 0.3) + 
      geom_density(aes(fill = program), alpha = 0.3) +
      scale_color_manual(values = program_palette) + 
      scale_fill_manual(values = program_palette) + 
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
      labs(title = "Max Size") +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line("black"),
            axis.text.x = element_text(color = "black", face = "bold", size = 18),
            axis.text.y = element_text(color = "black", face = "bold", size = 18),
            # axis.text.y = element_text(color = "black", face = "bold", size = 18),
            plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
            legend.position = "none")

ggsave("output/ms-second-round/plots/figure2/maxsize.tiff", units = "in", width = 12,
       height = 3, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/figure2/maxsize.svg", units = "in", width = 12,
       height = 3, dpi =  600)

### size structure ---
ann_dt |> 
      filter(total_nitrogen_ann < 25000) |>
      filter(total_nitrogen_ann > 0) |>
      ggplot(aes(x = log1p(skew_size_ann), color = program, fill = program, group = program)) +
      geom_density(aes(color = program), size = 1.5, alpha = 0.3) + 
      geom_density(aes(fill = program), alpha = 0.3) +
      scale_color_manual(values = program_palette) + 
      scale_fill_manual(values = program_palette) + 
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
      labs(title = "Size Structure") +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line("black"),
            axis.text.x = element_text(color = "black", face = "bold", size = 18),
            axis.text.y = element_text(color = "black", face = "bold", size = 18),
            # axis.text.y = element_text(color = "black", face = "bold", size = 18),
            plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
            legend.position = "none")

ggsave("output/ms-second-round/plots/figure2/sizestructure.tiff", units = "in", width = 12,
       height = 3, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/figure2/sizestructure.svg", units = "in", width = 12,
       height = 3, dpi =  600)

### species richness ---
ann_dt |> 
      filter(total_nitrogen_ann < 25000) |>
      filter(total_nitrogen_ann > 0) |>
      ggplot(aes(x = log1p(species_richness_ann), color = program, fill = program, group = program)) +
      geom_density(aes(color = program), size = 1.5, alpha = 0.3) + 
      geom_density(aes(fill = program), alpha = 0.3) +
      scale_color_manual(values = program_palette) + 
      scale_fill_manual(values = program_palette) + 
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
      labs(title = "Species Richness") +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line("black"),
            axis.text.x = element_text(color = "black", face = "bold", size = 18),
            axis.text.y = element_text(color = "black", face = "bold", size = 18),
            # axis.text.y = element_text(color = "black", face = "bold", size = 18),
            plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
            legend.position = "none")

ggsave("output/ms-second-round/plots/figure2/speciesrichness.tiff", units = "in", width = 12,
       height = 3, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/figure2/speciesrichness.svg", units = "in", width = 12,
       height = 3, dpi =  600)

### species diversity ---
ann_dt |> 
      filter(total_nitrogen_ann < 25000) |>
      filter(total_nitrogen_ann > 0) |>
      ggplot(aes(x = log1p(species_diversity_ann), color = program, fill = program, group = program)) +
      geom_density(aes(color = program), size = 1.5, alpha = 0.3) + 
      geom_density(aes(fill = program), alpha = 0.3) +
      scale_color_manual(values = program_palette) + 
      scale_fill_manual(values = program_palette) +
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
      labs(title = "Species Diversity") +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line("black"),
            axis.text.x = element_text(color = "black", face = "bold", size = 18),
            axis.text.y = element_text(color = "black", face = "bold", size = 18),
            # axis.text.y = element_text(color = "black", face = "bold", size = 18),
            plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
            legend.position = "none")

ggsave("output/ms-second-round/plots/figure2/speciesdiversity.tiff", units = "in", width = 12,
       height = 3, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/figure2/speciesdiversity.svg", units = "in", width = 12,
       height = 3, dpi =  600)

### trophic richness ---
ann_dt |> 
      filter(total_nitrogen_ann < 25000) |>
      filter(total_nitrogen_ann > 0) |>
      ggplot(aes(x = log1p(trophic_richness_ann), color = program, fill = program, group = program)) +
      geom_density(aes(color = program), size = 1.5, alpha = 0.3) + 
      geom_density(aes(fill = program), alpha = 0.3) +
      scale_color_manual(values = program_palette) + 
      scale_fill_manual(values = program_palette) + 
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
      labs(title = "Trophic Richness") +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line("black"),
            axis.text.x = element_text(color = "black", face = "bold", size = 18),
            axis.text.y = element_text(color = "black", face = "bold", size = 18),
            # axis.text.y = element_text(color = "black", face = "bold", size = 18),
            plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
            legend.position = "none")

ggsave("output/ms-second-round/plots/figure2/trophicrichness.tiff", units = "in", width = 12,
       height = 3, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/figure2/trophicrichness.svg", units = "in", width = 12,
       height = 3, dpi =  600)

### trophic diversity ---
ann_dt |> 
      filter(total_nitrogen_ann < 25000) |>
      filter(total_nitrogen_ann > 0) |>
      ggplot(aes(x = log1p(trophic_diversity_ann), color = program, fill = program, group = program)) +
      geom_density(aes(color = program), size = 1.5, alpha = 0.3) + 
      geom_density(aes(fill = program), alpha = 0.3) +
      scale_color_manual(values = program_palette) + 
      scale_fill_manual(values = program_palette) + 
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
      labs(title = "Trophic Diversity") +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line("black"),
            axis.text.x = element_text(color = "black", face = "bold", size = 18),
            axis.text.y = element_text(color = "black", face = "bold", size = 18),
            # axis.text.y = element_text(color = "black", face = "bold", size = 18),
            plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
            legend.position = "none")

ggsave("output/ms-second-round/plots/figure2/trophicdiversity.tiff", units = "in", width = 12,
       height = 3, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/figure2/trophicdiversity.svg", units = "in", width = 12,
       height = 3, dpi =  600)

### legend ---
ann_dt |> 
      filter(total_nitrogen_ann < 25000) |>
      filter(total_nitrogen_ann > 0) |>
      ggplot(aes(x = log1p(trophic_diversity_ann), color = program, fill = program, group = program)) +
      geom_density(aes(color = program), size = 1.5, alpha = 0.3) + 
      geom_density(aes(fill = program), alpha = 0.3) +
      scale_color_manual(values = program_palette) + 
      scale_fill_manual(values = program_palette) + 
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
      labs(title = "Trophic Diversity", color = 'Program', fill = 'Program') +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_line("black"),
            axis.text.x = element_text(color = "black", face = "bold", size = 18),
            axis.text.y = element_text(color = "black", face = "bold", size = 18),
            # axis.text.y = element_text(color = "black", face = "bold", size = 18),
            plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
            legend.position = "bottom",
            legend.text = element_text(color = "black", face = "bold", size = 18),
            legend.title = element_text(color = "black", face = "bold", size = 18))

ggsave("output/ms-second-round/plots/figure2/legend.tiff", units = "in", width = 12,
       height = 3, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/figure2/legend.svg", units = "in", width = 12,
       height = 3, dpi =  600)

### exploratory figures ---

ann_dt |> 
      filter(program == "PCCC") |> 
      ggplot(aes(x = log1p(total_nitrogen_ann), y = factor(year), fill = habitat)) +
      geom_density_ridges() +
      scale_fill_manual(values = c("orange", "blue")) + # Adjust colors for your programs
      theme_minimal() +
      labs(x = "Total Annual Nitrogen Supply",
           y = "year",
           fill = "Program")

ann_dt |> 
      filter(program == "PCCS") |> 
      ggplot(aes(x = log1p(total_nitrogen_ann), y = factor(year), fill = habitat)) +
      geom_density_ridges() +
      scale_fill_manual(values = c("orange", "blue")) + # Adjust colors for your programs
      theme_minimal() +
      labs(x = "Total Annual Nitrogen Supply",
           y = "year",
           fill = "Habitat")

ann_dt |> 
      filter(program == "FCE") |> 
      ggplot(aes(x = log1p(total_nitrogen_ann), y = factor(year), fill = habitat)) +
      geom_density_ridges() +
      scale_fill_manual(values = c("orange", "blue")) + # Adjust colors for your programs
      theme_minimal() +
      labs(x = "Total Annual Nitrogen Supply",
           y = "year",
           fill = "Habitat")

ann_dt |> 
      filter(program == "MCR") |> 
      ggplot(aes(x = log1p(total_nitrogen_ann), y = factor(year), fill = habitat)) +
      geom_density_ridges() +
      scale_fill_manual(values = c("orange", "blue", "purple")) + # Adjust colors for your programs
      theme_minimal() +
      labs(x = "Total Annual Nitrogen Supply",
           y = "year",
           fill = "Habitat")

ann_dt |> 
      filter(program == "SBC") |> 
      ggplot(aes(x = log1p(total_nitrogen_ann), y = factor(year), fill = habitat)) +
      geom_density_ridges() +
      ### not sufficient number of each points each year for MPAs
      scale_fill_manual(values = c("orange", "blue")) + 
      theme_minimal() +
      labs(x = "Total Annual Nitrogen Supply",
           y = "year",
           fill = "Habitat")

ann_dt |> 
      filter(program == "VCR") |> 
      ggplot(aes(x = log1p(total_nitrogen_ann), y = factor(year), fill = habitat)) +
      geom_density_ridges() +
      ### not sufficient number of each points each year for MPAs
      scale_fill_manual(values = c("orange", "blue")) + 
      theme_minimal() +
      labs(x = "Total Annual Nitrogen Supply",
           y = "year",
           fill = "Habitat")

