###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): Mack White, Adrian Stier, Nate Lemoine
###goal(s): join cnd, population, community, and turnover/synchrony datasets
###date(s): July 2024
###note(s): 

# Housekeeping ------------------------------------------------------------

### load necessary libraries
### install.packages("librarian")

librarian::shelf(tidyverse, readr, readxl)

###########################################################################
# read in csv files for joins ---------------------------------------------
###########################################################################

### nutrient and diversity calculations ---

cnd_pop <- read_csv('local_data/species-level-nutrient-stability.csv') |> 
      mutate(program = as.factor(program),
             habitat = as.factor(habitat),
             site = as.factor(site))

cnd_troph <- read_csv('local_data/trophic-level-nutrient-stability.csv') |> 
      mutate(program = as.factor(program),
             habitat = as.factor(habitat),
             site = as.factor(site))

cnd_comm <- read_csv("local_data/community-level-nutrient-stability.csv") |> 
      mutate(program = as.factor(program),
             habitat = as.factor(habitat),
             site = as.factor(site))

### turnover and synchrony ---

turnsynch_pop <- read_csv("local_data/population-turnover-synchrony.csv") |> 
      mutate(program = as.factor(program),
             habitat = as.factor(habitat),
             site = as.factor(site))

turnsynch_troph <- read_csv("local_data/trophic-turnover-synchrony.csv") |> 
      mutate(program = as.factor(program),
             habitat = as.factor(habitat),
             site = as.factor(site)) |> 
      rename(troph_beta_time = beta_time,
             troph_synch = synch)

comm_cnd_div_mech <- left_join(cnd_comm, turnsynch_pop, by = c("program", "habitat", "site"))

all <- left_join(comm_cnd_div_mech, cnd_troph, by = c("program", "habitat", "site"))

df_raw <- left_join(all, cnd_pop, by = c("program", "habitat", "site", "troph_group"))
# write_csv(df_raw, "local_data/dsr-eco-org-raw-10172024.csv")
# write_csv(df_raw, "local_data/dsr-eco-org-raw-10292024.csv")
write_csv(df_raw, "local_data/dsr-eco-org-raw.csv")

df_raw_wtroph_dynamics <- left_join(df_raw, turnsynch_troph, by = c("program", "habitat", "site"))
# write_csv(df_raw_wtroph_dynamics, "local_data/dsr-eco-org-raw-10292024-all.csv")
write_csv(df_raw_wtroph_dynamics, "local_data/dsr-eco-org-raw-all.csv")

df_raw_wtroph_dynamics |> 
      ggplot(aes(x=troph_synch, y = synch)) +
      geom_point() +
      geom_smooth(method = "lm") + 
      facet_wrap(~program, scales = "free")

df_raw_wtroph_dynamics |> 
      ggplot(aes(x=troph_beta_time, y = beta_time)) +
      geom_point() +
      geom_smooth(method = "lm") + 
      facet_wrap(~program, scales = "free")