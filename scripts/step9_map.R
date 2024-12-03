###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): Mack White
###goal(s): create a map for figure one
###date(s): November 2024
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
map_dt <- read_csv("../cndwg_website/data/LTER_Site_Coords_w_information.csv") |> 
      filter(Site %in% c("FCE", "MCR", "SBC", "VCR")) |> 
      rename(program = Site,
             lat = Latitude,
             long = Longitude) |> 
      select(program, lat, long)
glimpse(dt)

pisco <- data.frame(
      program = c("PCCC", "PCCS"),
      lat = c(38.47409, 35.5000),
      long = c(-123.2485, -121.019020)
)

dt <- rbind(map_dt, pisco)

### set up color palette ---
palette <- c("Overall"="#000000", "FCE"="#64a988", "MCR"="#ff967d", 'PCCC'="#2A788EFF", "PCCS"="#8b6b93",
             'SBC'='#ff3f4c', "VCR"="#9b9254")

### install necessary packages ---
librarian::shelf(tidyverse, readr, janitor, zoo, 
                 lubridate, openintro, maps, ggmap,
                 ggthemes, shapefiles, broom, sf, ggspatial, 
                 GISTools)

### set theme ---
theme_set(theme_minimal())

### set up google map key ---
api_secret <- 'AIzaSyB_Q-Ow1klpX9jblm7M2614k5KCVYUXTZM'
register_google(key = api_secret)
has_google_key()

### read in the map
site_map <- get_map(
      ### determine bounding box: https://www.openstreetmap.org/#map=5/25.304/-69.412
      c(left = -188.7, bottom = -25.0, right = -24.9, top = 55.0),
      maptype = 'satellite',
      source = 'google',
      api_key = api_secret
)

ggmap(site_map) +
      geom_jitter(data = dt, aes(
                 x = long,
                 y = lat,
                 color = program),
                 size = 4,
                 position = 'jitter') +
      scale_color_manual(values = palette) +
      theme_map() + 
      theme(legend.background = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(face = 'bold',
                                       color = 'white',
                                       size = 16),
            legend.position = c(0.001, 0.350))
      # annotation_north_arrow(location = 'br',
      #                        style = north_arrow_fancy_orienteering(text_col = 'white',
      #                                                               fill = 'white',
      #                                                               line_col = 'white',
      #                                                               text_face = "bold",
      #                                                               text_size = 18)) +
      # annotation_scale(location = 'br', width_hint = 0.5, text_cex = 1.25,
      #                  text_face = "bold", text_col = 'white') +
      # coord_sf(crs = 4326) +
      # theme_map() +
      # theme(legend.background = element_blank(),
      #       legend.title = element_blank(),
      #       legend.text = element_text(face = 'bold',
      #                                  color = 'white',
      #                                  size = 16),
      #       legend.position = c(0.001, 0.275))

### save for publication
ggsave("output/ms-second-round/plots/figure1/test-map.tiff", units = "in", width = 10,
       height = 8, dpi =  600, compression = "lzw")
ggsave("output/ms-second-round/plots/figure1/test-map.svg", units = "in", width = 10,
       height = 8, dpi =  600)
# ggsave("output/ms-second-round/plots/figure1/map-scalebar-arrow.svg", units = "in", width = 10,
#        height = 8, dpi =  600)
