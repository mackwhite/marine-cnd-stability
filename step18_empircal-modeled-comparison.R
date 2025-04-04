###project: LTER Marine Consumer Nutrient Dynamic Synthesis Working Group
###author(s): Mack White, Adrian Stier, Nate Lemoine
###goal(s): join cnd, population, community, and turnover/synchrony datasets
###date(s): July 2024
###note(s): 

# Housekeeping ------------------------------------------------------------

### load necessary libraries
### install.packages("librarian")

librarian::shelf(tidyverse, readr, lme4, ggpubr, performance, 
                 scales, ggstats, ggeffects, visreg, mgcv, MuMIn, glmmTMB, corrplot,
                 broom, purrr)

### read in necessary data ----
emp <- read_csv('local_data/empirical_excretion_kelp.csv') |> 
      rename(order = TAXON_ORDER,
             family = TAXON_FAMILY,
             genus = TAXON_GENUS,
             species = TAXON_SPECIES,
             wetmass_g = WM_g,
             diet_cat = `Functional group`,
             nind_umol_hr = `exc_rate_NH4_umol_hour-1`,
             region = Region) |> 
      mutate(scientific_name = paste(genus, species, sep = " "),
             nind_ug_hr = 18.042*nind_umol_hr) |> 
      select(scientific_name, wetmass_g, nind_ug_hr, nind_umol_hr, order, family, genus, species, diet_cat, region)
glimpse(emp)

species_list <- emp |> select(scientific_name, order, family, genus, diet_cat, region) |> distinct()
# write_csv(species_list, "../../../to be filed/kelp_empirical_species_list.csv")

# dt <- read.csv(file.path("tier2", "harmonized_consumer_excretion_CLEAN.csv"),stringsAsFactors = F,na.strings =".") |> 
#       janitor::clean_names()
# glimpse(dt)
# 
# dt_test1 <- dt |> 
#       select(project, scientific_name, order, family, genus, diet_cat) |>
#       filter(project %in% c('CoastalCA', 'SBC')) |> 
#       distinct() |> 
#       select(-project) |> 
#       distinct()
# write_csv(dt_test1, "../../../to be filed/kelp_modeled_species_list.csv")

### read in empirical species list with diet categorized for models ----
kelp_diet <- read_csv('../../../to be filed/kelp_empirical_species_list_updated.csv')
glimpse(kelp_diet)      
kelp <- emp |> select(-diet_cat)
glimpse(kelp)
kelp_all <- kelp |> left_join(kelp_diet) |> distinct()

### clean up environment
keep <- c("kelp_all")
rm(list = setdiff(ls(), keep))

### format model data ----
kelp1 <- kelp_all |> 
      rename(nind_ug_hr_emp = nind_ug_hr,
             nind_umol_hr_emp = nind_umol_hr) |> 
      mutate(drymass_g = wetmass_g*dm_conv,
             phylum = "Chordata") |> 
      mutate(temp = case_when(
            region == "Central California" ~ 13,
            region == "Southern California" ~ 17
      ))
glimpse(kelp1)

### generate coefficients from vanni and mcintyre model ----
kelp2 <- kelp1 |> 
      ### vertebrate coefficient classification
      mutate(vert_coef = if_else(phylum == "Chordata", 0.7804, 0),
             vert_coef_upper = if_else(phylum == "Chordata", 0.7804 + 0.0655, 0),
             vert_coef_lower = if_else(phylum == "Chordata", 0.7804 - 0.0655, 0)) |> 
      ### diet coefficient classification
      mutate(diet_coef = case_when(
            diet_cat == "algae_detritus" ~ -0.0389,
            diet_cat == "invert" ~ -0.2013,
            diet_cat == "fish" ~ -0.0537,
            diet_cat == "fish_invert" ~ -0.1732,
            diet_cat == "algae_invert" ~ 0,
            TRUE ~ NA)) |> 
      mutate(diet_coef_upper = case_when(
            diet_cat == "algae_detritus" ~ -0.0389 + 0.0765,
            diet_cat == "invert" ~ -0.2013 + 0.0771,
            diet_cat == "fish" ~ -0.0537 + 0.2786,
            diet_cat == "fish_invert" ~ -0.1732 + 0.1384,
            diet_cat == "algae_invert" ~ 0,
            TRUE ~ NA)) |> 
      mutate(diet_coef_lower = case_when(
            diet_cat == "algae_detritus" ~ -0.0389 - 0.0765,
            diet_cat == "invert" ~ -0.2013 - 0.0771,
            diet_cat == "fish" ~ -0.0537 - 0.2786,
            diet_cat == "fish_invert" ~ -0.1732 - 0.1384,
            diet_cat == "algae_invert" ~ 0,
            TRUE ~ NA)) |> 
      ### temperature coefficient classification
      mutate(temp_coef = 0.0246,
             temp_coef_upper = 0.0246 + 0.0014,
             temp_coef_lower = 0.0246 - 0.0014) |> 
      ### dry mass coefficient classification
      mutate(dm_coef = 0.6840,
             dm_coef_upper = 0.6840 + 0.0177,
             dm_coef_lower = 0.6840 - 0.0177) |> 
      ### intercept coefficient classification
      mutate(int_coef = 1.4610,
             int_coef_upper = 1.4610 + 0.0897,
             int_coef_lower = 1.4610 - 0.0897)
glimpse(kelp2)

### calculate excretion rates ----
kelp3 <- kelp2 |> 
      mutate(n10 = int_coef + dm_coef*(log10(drymass_g)) + temp_coef*temp - diet_coef + vert_coef,
             n10_lower = int_coef_lower + dm_coef_lower*(log10(drymass_g)) + temp_coef_lower*temp - diet_coef_lower + vert_coef_lower,
             n10_upper = int_coef_upper + dm_coef_upper*(log10(drymass_g)) + temp_coef_upper*temp - diet_coef_upper + vert_coef_upper) |> 
      mutate(nind_ug_hr = 10^n10,
             nind_ug_hr_lower = 10^n10_lower,
             nind_ug_hr_upper = 10^n10_upper)

### format model data for visualization ----
kelp4 <- kelp3 |> 
      select(scientific_name, wetmass_g, drymass_g, diet_cat,
             nind_ug_hr_emp, nind_ug_hr, nind_ug_hr_lower, nind_ug_hr_upper,
             phylum, order, family, genus, species) |> 
      mutate(
            n10_emp = log10(nind_ug_hr_emp),
            n10_mod = log10(nind_ug_hr),
            n10_mod_low = log10(nind_ug_hr_lower),
            n10_mod_upp = log10(nind_ug_hr_lower),
            n10_dm = log10(drymass_g)
      ) |> 
      group_by(scientific_name) |> mutate(n = n()) |> ungroup() |> 
      filter(n > 2) |> 
      select(-n) |> 
      rename(diet = diet_cat) |> 
      group_by(scientific_name) |> mutate(species_n = n()) |> ungroup() |> 
      group_by(genus) |> mutate(genus_n = n()) |> ungroup() |> 
      group_by(family) |> mutate(family_n = n()) |> ungroup()

kelp4 |> 
      ggplot(aes(x = n10_dm)) +
      geom_point(aes(y = n10_emp), color = "blue") +
      geom_point(aes(y = n10_mod), color = "black") +
      facet_wrap(~scientific_name, scale = "free")

kelp4 |> 
      ggplot(aes(x = n10_dm)) +
      geom_point(aes(y = n10_emp), color = "blue") +
      geom_point(aes(y = n10_mod), color = "black") +
      facet_wrap(~scientific_name, scale = "free")

test <- lm(n10_emp ~ n10_mod, data = kelp4)
summary(test)
kelp4 |> 
      ggplot(aes(x = n10_mod, y = n10_emp)) +
      geom_point()
      geom_point(aes(y = n10_emp), color = "blue") +
      geom_point(aes(y = n10_mod), color = "black") +
      facet_wrap(~scientific_name, scale = "free")

### generate r2 values for all species ----

### species-level evaluation ---
speciesR2 <- kelp4 |> 
      group_by(scientific_name, species_n) |> 
      nest() |> 
      mutate(
            model = map(data, ~ lm(n10_emp ~ n10_mod, data = .x)),
            model_summary = map(model, glance)
      ) |> 
      unnest(model_summary) |> 
      select(scientific_name, r.squared, p.value, species_n) |> 
      ungroup()

species_summary <- speciesR2 |> 
      summarize(r2 = mean(r.squared),
                p = mean(p.value))

### genus-level evaluation ---
genusR2 <- kelp4 |> 
      group_by(genus, genus_n) |> 
      nest() |> 
      mutate(
            model = map(data, ~ lm(n10_emp ~ n10_mod, data = .x)),
            model_summary = map(model, glance)
      ) |> 
      unnest(model_summary) |> 
      select(genus, r.squared, p.value, genus_n) |> 
      ungroup()

genus_summary <- genusR2 |> 
      summarize(r2 = mean(r.squared),
                p = mean(p.value))

### family-level evaluation ---
familyR2 <- kelp4 |> 
      group_by(family, family_n) |> 
      nest() |> 
      mutate(
            model = map(data, ~ lm(n10_emp ~ n10_mod, data = .x)),
            model_summary = map(model, glance)
      ) |> 
      unnest(model_summary) |> 
      select(family, r.squared, p.value, family_n) |> 
      ungroup()

family_summary <- familyR2 |> 
      summarize(r2 = mean(r.squared),
                p = mean(p.value))

keep <- c("speciesR2", "genusR2", "familyR2", "kelp4",
          "species_summary", "genus_summary", "family_summary")
rm(list = setdiff(ls(), keep))

### review output
print(species_summary)
print(genus_summary)
print(family_summary)

### evaluate correlation ----
mod <- glmmTMB(
      n10_emp ~ n10_mod, family = gaussian(), data = kelp4
)
summary(mod)
performance::performance(mod)
model_validation <- ggemmeans(mod, 'n10_mod')
plot(model_validation)

model_validation |> 
      ggplot(aes(x = x, y = predicted)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                  fill = '#003153', alpha = 0.4) +
      annotate('text', 
               x = 1.76, y = 5.4,
               label = bquote({R^2} == 0.76),
               size = 6) +
      annotate('text', 
               x = 1.89, y = 5.12,
               label = bquote(italic(p) < 2 %*% 10^-16),
               size = 6) +
      geom_line(size = 1.5, color = "#003153") +
      theme_classic() +
      scale_y_continuous(breaks = c(1,2,3,4,5), limits = c(0.68,5.41)) +
      scale_x_continuous(breaks = c(2,3,4,5), limits = c(1.5,5.5)) +
      ylab(expression(bold("Empirical Log" [10] * " N Excretion (" * mu * "g" %.% ind^-1 %.% hr^-1 * ")"))) +
      xlab(expression(bold("Modeled Log" [10] * " N Excretion (" * mu * "g" %.% ind^-1 %.% hr^-1 * ")"))) +
      theme(axis.text.x = element_text(face = "bold", color = "black", size = 14),
            axis.text.y = element_text(face = "bold", color = "black", size = 14),
            axis.title.x = element_text(face = "bold", color = "black", size = 14),
            axis.title.y = element_text(face = "bold", color = "black", size = 14),
            legend.position = "none",
            legend.text = element_text(face = "bold", color = "black"),
            legend.title = element_text(face = "bold", color = "black"))

ggsave("output/figs/supplemental-model-validation.png", units = "in", width = 6,
       height = 5, dpi =  600)

