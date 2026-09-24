###### Code for TPC and Physiology comparisons ####### 
### Created by: Maya Powell
#### Last updated on: Oct 2 2025

### Install Packages #####
## if these packages are not yet installed, install them 
## great for updates or new users 
#if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented')

####Read in required libraries#####
##### Include Versions of libraries
#library(lubridate)
library(tidyverse)
library(here)
library(PNWColors)
library(viridis)
library(car)
library(dplyr)
library(ggplot2)
library(powerjoin)
library(forcats)
library(car)
library(emmeans)
library(see)
library(performance)
library(purrr)
library(rlang)
library(tibble)
library(here)
library(psych)
library(corrplot)
library(ggpubr)
library(vegan)
#remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(multcomp)

###### Initial Data Read In ########
#TPC data with only the seven species that we have physio data for
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4.csv"))
topt_df <- topt_df %>% filter(sample_ID != "B08_TPC") %>% select(-species,-full_species,-SA_cm2)
#topt_df <- topt_df %>% dplyr::select(rmax, topt, e, PR, frag_ID)
phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta_all.csv"))
topt_df <- topt_df %>% left_join(phys_meta, by = "frag_ID")

#additional setup
sp_cols <- c(
  "Acropora hyacinthus" = '#d8aedd',
  "Echinopora lamellosa" = '#ba7999',
  "Favites complanata"   = '#dd4124',
  "Montipora aequituberculata" = '#ed8b00',
  "Montipora vietnamensis" = '#efbc82',
  "Pachyseris rugosa" = '#edd746',
  "Pocillopora eydouxi" = '#d0e2af',
  "Porites cylindrica" = '#45681e',
  "Porites rus" = '#7bbcd5',
  "Turbinaria frondens" = '#00496f'
)

se_fun <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

###Correlation coefficient plot####
#generated below - only need once
# avg_AFDW <- read.csv(here("Data", "Physiology", "Average_Ash_Free_Dry_Weight.csv")) #afdw_mg_cm2
# avg_AFDW <- avg_AFDW %>% mutate(afdw_log = log(afdw_mg_cm2)) %>%
#   mutate(dw_log = log(dw_mg_cm2)) %>%
#   dplyr::select(frag_ID, dw_mg_cm2, dw_log, afdw_mg_cm2,afdw_log)
# chla_avg <- read.csv(here("Data", "Physiology", "Chla_avg.csv")) #chla_ug_cm2_mean
# chla_avg <- chla_avg %>% mutate(chla_log = log(chla_ug_cm2_mean)) %>%
#   mutate(chla_sym_log = log(chla_pg_sym)) %>%
#   dplyr::select(frag_ID, chla_ug_cm2_mean,chla_log, chla_pg_sym, chla_sym_log)
# avg_sym <- read.csv(here("Data", "Physiology", "Average_Sym_Density.csv"))
# avg_sym <- avg_sym %>% filter(frag_ID != "C07") %>% filter(frag_ID != "D10") #remove crazy outliers!
# avg_sym <- avg_sym %>% mutate(sym_log = log(sym_cm2)) %>% dplyr::select(frag_ID, sym_cm2,sym_log)
# prot <- read.csv(here("Data", "Physiology", "protein_all_summary.csv")) #prot_ug_cm2
# prot <- prot %>% mutate(prot_log = log(prot_ug_cm2)) %>% dplyr::select(frag_ID, prot_ug_cm2,prot_log)
# physio_list <- list(avg_AFDW,chla_avg,avg_sym,prot)
# all_physio <- physio_list %>% reduce(left_join)
# write_csv(all_physio, here("Data", "Physiology", "all_physio_data.csv"))
all_physio <- read_csv(here("Data", "Physiology", "all_physio_data.csv"))

#generate full dataframe
# result <- Reduce(function(x, y) merge(x, y, all = TRUE), list(topt_df,avg_AFDW,chla_avg,avg_sym,prot))
# #drop NA data from respo data because we donʻt have all reps
# all_data <- result %>% drop_na(rmax)
# write_csv(result, here("Data", "Physiology", "all_data_concatenated.csv"))

###Correlation plots between thermal performance and physio data #####
####Significant topt and rmax relationship plots based on corr plots####
#NP
#temp_at_NPR1 chla_ug_cm2
#temp_at_NPR1 chla_pg_sym
#rmax chla_ug_cm2
#rmax chla_pg_sym
#topt chla_ug_cm2
#topt chla_pg_sym
#topt prot_ug_cm2

#GP
#rmax chla_ug_cm2
#rmax chla_pg_sym
#breadth prot_ug_cm2

#R
#topt prot_ug_cm2
#breadth chla_ug_cm2
#bredth afdw_mg_cm2

#np rmax and chla pg sym
np_rmax_chla_sym_scatter <- ggplot(np_data) +
  geom_point(aes(y = rmax, x = chla_pg_sym, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~full_species, scales = "free")+
  geom_smooth(aes(y = rmax, x = chla_pg_sym, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Net Photosynthesis",
       x = expression("Chlorophyll a" ~ (pg ~ symbiont^{-1})), 
       y = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
       color = "Species")
np_rmax_chla_sym_scatter
#ggsave(here("Output", "Physiology", "np_rmax_chla_sym_scatter.pdf"), np_rmax_chla_sym_scatter, h = 5, w = 10)

#np rmax and chla ug cm
np_rmax_chla_scatter <- ggplot(np_data) +
  geom_point(aes(y = rmax, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Net Photosynthesis",
       x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
       color = "Species")
np_rmax_chla_scatter
#ggsave(here("Output", "Physiology", "np_rmax_chla_scatter.pdf"), np_rmax_chla_scatter, h = 5, w = 10)

#np topt and chla pg sym
np_topt_chla_sym_scatter <- ggplot(np_data) +
  geom_point(aes(y = topt, x = chla_pg_sym, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~full_species, scales = "free")+
  geom_smooth(aes(y = topt, x = chla_pg_sym, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Net Photosynthesis",
       x = expression("Chlorophyll a" ~ (pg ~ symbiont^{-1})), 
       y = "Thermal optimum (°C)",
       color = "Species")
np_topt_chla_sym_scatter
#ggsave(here("Output", "Physiology", "np_topt_chla_sym_scatter.pdf"), np_topt_chla_sym_scatter, h = 5, w = 10)

#np topt and chla ug cm
np_topt_chla_scatter <- ggplot(np_data) +
  geom_point(aes(y = topt, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Net Photosynthesis",
       x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = "Thermal optimum (°C)",
       color = "Species")
np_topt_chla_scatter
#ggsave(here("Output", "Physiology", "np_topt_chla_scatter.pdf"), np_topt_chla_scatter, h = 5, w = 10)

#np topt and prot_ug_cm2
np_topt_prot_scatter <- ggplot(np_data) +
  geom_point(aes(y = topt, x = prot_ug_cm2, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = topt, x = prot_ug_cm2, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Net Photosynthesis",
       x = expression("Protein Content" ~ (mu*g ~ cm^{-2})),
       y = "Thermal optimum (°C)",
       color = "Species")
np_topt_prot_scatter

#np p:r and chla ug cm sym
np_pr_chla_sym_scatter <- ggplot(np_data) +
  geom_point(aes(y = temp_at_NPR1, x = chla_pg_sym, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = temp_at_NPR1, x = chla_pg_sym, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Net Photosynthesis",
       x = expression("Chlorophyll a" ~ (pg ~ symbiont^{-1})),
       y = "Temp at NP:R = 1 (°C)",
       color = "Species")
np_pr_chla_sym_scatter

#np p:r and chla ug cm
np_pr_chla_scatter <- ggplot(np_data) +
  geom_point(aes(y = temp_at_NPR1, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = temp_at_NPR1, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Net Photosynthesis",
       x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = "Temp at NP:R = 1 (°C)",
       color = "Species")
np_pr_chla_scatter

#gp rmax and chla pg sym
gp_rmax_chla_sym_scatter <- ggplot(gp_data) +
  geom_point(aes(y = rmax, x = chla_pg_sym, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~full_species, scales = "free")+
  geom_smooth(aes(y = rmax, x = chla_pg_sym, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Gross Photosynthesis",
       x = expression("Chlorophyll a" ~ (pg ~ symbiont^{-1})), 
       y = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
       color = "Species")
gp_rmax_chla_sym_scatter
#ggsave(here("Output", "Physiology", "gp_rmax_chla_sym_scatter.pdf"), gp_rmax_chla_sym_scatter, h = 5, w = 10)

#gp rmax and chla ug cm
gp_rmax_chla_scatter <- ggplot(gp_data) +
  geom_point(aes(y = rmax, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Gross Photosynthesis",
       x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
       color = "Species")
gp_rmax_chla_scatter

#r topt and prot ug cm
r_topt_prot_scatter <- ggplot(r_data) +
  geom_point(aes(y = topt, x = prot_ug_cm2, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~full_species, scales = "free")+
  geom_smooth(aes(y = topt, x = prot_ug_cm2, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Respiration",
       x = expression("Protein" ~ (mu*g ~ cm^{-2})),
       y = "Thermal optimum (°C)",
       color = "Species")
r_topt_prot_scatter
#ggsave(here("Output", "Physiology", "r_topt_prot_scatter.pdf"), r_topt_prot_scatter, h = 5, w = 10)

#r breadth chla
r_breadth_chla_scatter <- ggplot(r_data) +
  geom_point(aes(y = breadth, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = breadth, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Respiration",
       x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = "Breadth (°C)",
       color = "Species")
r_breadth_chla_scatter

#r breadth biomass
r_breadth_afdw_scatter <- ggplot(r_data) +
  geom_point(aes(y = breadth, x = afdw_mg_cm2, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = breadth, x = afdw_mg_cm2, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(title = "Respiration",
       x = expression("Biomass" ~ (mg ~ cm^{-2})),
       y = "Breadth (°C)",
       color = "Species")
r_breadth_afdw_scatter

#temp_at_NPR1 chla_ug_cm2
#temp_at_NPR1 chla_pg_sym
#rmax chla_ug_cm2
#rmax chla_pg_sym
#topt chla_ug_cm2
#topt chla_pg_sym
#topt prot_ug_cm2

all_corr_topt_physio <- ggarrange(np_topt_chla_scatter, np_topt_chla_sym_scatter, np_topt_prot_scatter,
                                  np_rmax_chla_scatter, np_rmax_chla_sym_scatter, np_pr_chla_scatter,
                                  np_pr_chla_sym_scatter, gp_rmax_chla_scatter, gp_rmax_chla_sym_scatter, 
                                  r_topt_prot_scatter,r_breadth_chla_scatter,r_breadth_afdw_scatter,
                          common.legend = T, legend = "right",
                          ncol = 3, nrow=4, labels = c("A","B","C","D","E","F","G","H","I","J","K","L"), font.label = list(size = 30, color = "black"))
ggsave(here("Output", "Physiology", "all_corr_tpc_physio.pdf"), all_corr_topt_physio, h = 20, w = 20)

##### Intraspecific variation ####
#### Effect size plots #####
all_data <- read_csv(here("Data", "Physiology", "all_data_concatenated.csv"))

gp_data <- all_data %>% filter(PR == "GrossPhoto")
np_data <- all_data %>% filter(PR == "NetPhoto")
r_data <- all_data %>% filter(PR == "Respiration")
models<- np_data %>%
  nest(.by = species) %>% # nest all the data by species
  mutate(fit = map(data, ~lm(topt~chla_ug_cm2_mean, data = .)))
#column of dataframes is called data 
# "." is like "i" with the for loop
#
models

#effect size function:
get_effect_sizes <- function(data, response, predictor, group_var) {
  form <- as.formula(paste0("scale(", response, ") ~ scale(", predictor, ")"))
  
  data %>%
    filter(!is.na(.data[[response]]), !is.na(.data[[predictor]])) %>%
    nest(.by = all_of(group_var)) %>%
    mutate(
      n     = map_int(data, nrow),
      fit   = map(data, ~ lm(form, data = .x)),
      coeffs = map(fit, tidy, conf.int = TRUE)
    ) %>%
    dplyr::select(-data, -fit) %>%
    unnest(coeffs) %>%
    filter(str_detect(term, "^scale\\(")) %>%   # keep only the slope term
    mutate(response = response, predictor = predictor)
}

predictors_to_test <- c("chla_ug_cm2_mean", "chla_pg_sym","prot_ug_cm2", "sym_cm2", "afdw_mg_cm2")

np_species_effects_topt <- map_dfr(
  predictors_to_test,
  ~ get_effect_sizes(np_data, response = "topt", predictor = .x,
                     group_var = "full_species"))
np_species_effects_rmax <- map_dfr(
  predictors_to_test,
  ~ get_effect_sizes(np_data, response = "topt", predictor = .x,
                     group_var = "full_species"))
np_species_effects_breadth <- map_dfr(
  predictors_to_test,
  ~ get_effect_sizes(np_data, response = "breadth", predictor = .x,
                     group_var = "full_species"))

#significant intraspecific correlations:
#topt
#Porites rus scale(chla_ug_cm2_mean) 0.009783428
#Montipora vietnamensis scale(chla_ug_cm2_mean) 0.015100501

#rmax
# Porites cylindrica	5	scale(afdw_mg_cm2)	0.002121706
# 7	Echinopora lamellosa	5	scale(chla_ug_cm2_mean)	0.006317979	
# 2	Montipora vietnamensis	4	scale(chla_ug_cm2_mean) 0.009980522
# 17	Echinopora lamellosa	5	scale(chla_pg_sym) 0.033929838
# 13	Favites complanata	5	scale(chla_pg_sym)	0.033995514	-1.68313340	-0.1292694	rmax	chla_pg_sym

#none for e

#breadth
# Acropora hyacinthus	5	scale(sym_cm2) 0.007918908
# Porites rus	5	scale(chla_pg_sym) 0.025449297	

gp_species_effects <- map_dfr(
  predictors_to_test,
  ~ get_effect_sizes(gp_data, response = "breadth", predictor = .x,
                     group_var = "full_species"))

#topt
#Acropora hyacinthus scale(sym_cm2) 0.03264724

#rmax
# Porites cylindrica	5	scale(afdw_mg_cm2) 0.001454162	-1.26498890	-0.7122578	rmax	afdw_mg_cm2
# Favites complanata	5	scale(chla_pg_sym) 0.010793493	-1.49217072	-0.4209702	rmax	chla_pg_sym
# Pocillopora eydouxi	5	scale(chla_pg_sym) 0.018157119	0.30383361	1.5730813	rmax	chla_pg_sym
# Echinopora lamellosa	5	scale(chla_ug_cm2_mean)	0.022140472	0.25292542	1.6064630	rmax	chla_ug_cm2_mean

#ct max
# Acropora hyacinthus	5	scale(prot_ug_cm2) 0.003581437	-1.35174628	-0.60671842	ctmax	prot_ug_cm2
# Turbinaria frondens	5	scale(prot_ug_cm2) 0.034759498	0.12231175	1.68726312	ctmax	prot_ug_cm2
# Pocillopora eydouxi	5	scale(sym_cm2) 0.043648256	0.04763586	1.73035899	ctmax	sym_cm2

#e
# Acropora hyacinthus	5	scale(sym_cm2)0.01079117	0.42101418	1.492139267	e	sym_cm2
# Porites rus	5	scale(sym_cm2)0.02078394	0.26955915	1.595673314	e	sym_cm2
# Porites rus	5	scale(chla_ug_cm2_mean)0.02390969	0.23217158	1.619765975	e	chla_ug_cm2_mean
# Porites rus	5	scale(prot_ug_cm2)0.04206668	0.06016947	1.723282436	e	prot_ug_cm2
# Pachyseris rugosa	5	scale(chla_ug_cm2_mean)0.04890456	-1.75236011	-0.007930897	e	chla_ug_cm2_mean

#breadth
#Porites cylindrica	5	scale(prot_ug_cm2) 0.03805736	-1.7042516	-0.09333619	breadth	prot_ug_cm2
#Acropora hyacinthus	5	scale(afdw_mg_cm2)	0.04276406	-1.7264309	-0.05460691	breadth	afdw_mg_cm2

# r_species_effects <- map_dfr(
#   predictors_to_test,
#   ~ get_effect_sizes(r_data, response = "topt", predictor = .x,
#                      group_var = "full_species"))
#none with large enough sample sizes

##COULD ALSO LOOK AT THIS WITH MORPHOLOGY

#Plots
#pull significant intraspecific correlations for NP:
# -	Topt vs chla: Mvie, Prus
# -	Rmax vs chla: Elam, Mvie
# -	Rmax vs chla_sym: Elam, Fcom
# -	Rmax vs afdw: Pcyl
# -	Breadth vs sym: Ahya
# -	Breadth vs chla_sym: Prus

np_topt_chla_sp_eff <- get_effect_sizes(np_data, response  = "topt", predictor = "chla_ug_cm2_mean", group_var = "full_species")
np_rmax_chla_sp_eff <- get_effect_sizes(np_data, response  = "rmax", predictor = "chla_ug_cm2_mean", group_var = "full_species")
np_rmax_chla_sym_sp_eff <- get_effect_sizes(np_data, response  = "rmax", predictor = "chla_pg_sym", group_var = "full_species")
np_rmax_afdw_sp_eff <- get_effect_sizes(np_data, response  = "rmax", predictor = "afdw_mg_cm2", group_var = "full_species")
np_breadth_chla_sym_sp_eff <- get_effect_sizes(np_data, response  = "breadth", predictor = "chla_pg_sym", group_var = "full_species")
np_breadth_sym_sp_eff <- get_effect_sizes(np_data, response  = "breadth", predictor = "sym_cm2", group_var = "full_species")
all_np_sig_eff_intra <- rbind(np_topt_chla_sp_eff,np_rmax_chla_sp_eff,np_rmax_chla_sym_sp_eff,
                              np_rmax_afdw_sp_eff,np_breadth_chla_sym_sp_eff,np_breadth_sym_sp_eff)
all_np_sig_eff_intra <- all_np_sig_eff_intra %>%
  mutate(sig_color = if_else(p.value < 0.05, as.character(full_species), "ns")) %>%
  mutate(facet = as.factor(paste(response,"vs",predictor))) %>%
  mutate(facet_nice = case_when(
    facet == "topt vs chla_ug_cm2_mean" ~ "Topt vs Chl a",
    facet == "rmax vs chla_ug_cm2_mean" ~ "Rmax vs Chl a",
    facet == "rmax vs chla_pg_sym" ~ "Rmax vs Chl a per sym",
    facet == "rmax vs afdw_mg_cm2" ~ "Rmax vs Biomass",
    facet == "breadth vs chla_pg_sym" ~ "Breadth vs Chl a per sym",
    facet == "breadth vs sym_cm2" ~ "Breadth vs Sym density"
    )) |>
  mutate(facet_nice=fct_relevel(facet_nice,c("Rmax vs Chl a","Rmax vs Chl a per sym","Rmax vs Biomass",
                                             "Topt vs Chl a","Breadth vs Chl a per sym","Breadth vs Sym density")))
  
  

color_values <- c(sp_cols, "ns" = "grey85")

np_topt_chla_effect_plot <- all_np_sig_eff_intra %>%
  ggplot(aes(x = estimate, y = full_species, color = sig_color)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), width = 0.2) +
  geom_point(size = 3) +
  scale_color_manual(values = color_values, breaks = names(sp_cols), name = "Species",
    labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')"))) +
  labs(x = "Standardized effect size") +
  theme_classic(base_size = 20) +
  facet_wrap(facet_nice~., nrow=2, axis.labels = "margins")+
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_text(face = "italic"))
np_topt_chla_effect_plot

ggsave(here("Output", "Physiology", "effect_size_np_tpc_physio_intraspecific.pdf"), 
       np_topt_chla_effect_plot, h = 8, w = 12)

topt_intra_effect_np <- np_species_effects %>%
  #ggplot(aes(x = estimate, y = morphology, color = morphology)) +
  ggplot(aes(x = estimate, y = full_species, color = full_species)) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.3, size = 2) +
  scale_color_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.7) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), width = 0.2) +
  labs(x = "Standardized effect size on thermal optimum", y = "") +
  facet_wrap(~predictor, scales = "free_x", nrow = 1) +
  theme_bw(base_size = 22) +
  theme(legend.position = "right", axis.title.y = element_blank(), axis.text.y = element_blank())
  #theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_text(face = "italic"))
topt_intra_effect_np

ggsave(here("Output","Physiology","physio_np_topt_effectsize_morphology_meanse.pdf"), multi_effect_plot, height = 5, width = 15)


##np topt and chla
np_topt_chla_mod <- lm(topt~chla_ug_cm2_mean, data = np_data)
Anova(np_topt_chla_mod)
summary(np_topt_chla_mod) 
coeffs<-tidy(np_topt_chla_mod)
#Residual standard error: 0.8391 on 47 degrees of freedom
#Multiple R-squared:  0.1288,	Adjusted R-squared:  0.1103 
#F-statistic: 6.949 on 1 and 47 DF,  p-value: 0.01133 *
check_model(np_topt_chla_mod)

#np rmax and chla
np_rmax_chla_mod <- lm(rmax~chla_ug_cm2_mean, data = np_data)
Anova(np_rmax_chla_mod)
summary(np_rmax_chla_mod) 
# Residual standard error: 0.2411 on 47 degrees of freedom
# Multiple R-squared:  0.1463,	Adjusted R-squared:  0.1281 
# F-statistic: 8.055 on 1 and 47 DF,  p-value: 0.00668 **
check_model(np_rmax_chla_mod)

#gp rmax and chla
gp_rmax_chla_mod <- lm(rmax~chla_pg_sym, data = gp_data)
Anova(gp_rmax_chla_mod)
summary(gp_rmax_chla_mod) 
# Residual standard error: 0.3211 on 47 degrees of freedom
# Multiple R-squared:  0.1187,	Adjusted R-squared:  0.09992 
# F-statistic: 6.329 on 1 and 47 DF,  p-value: 0.01535
check_model(gp_rmax_chla_mod)

#r topt and prot
r_topt_prot_mod <- lm(topt~prot_ug_cm2, data = r_data)
Anova(r_topt_prot_mod)
summary(r_topt_prot_mod) 
# Residual standard error: 0.8712 on 14 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.3958,	Adjusted R-squared:  0.3526 
# F-statistic: 9.171 on 1 and 14 DF,  p-value: 0.009029 **
check_model(r_topt_prot_mod)

####P:R ratios

#read in data
respo_constant_temps <- read_csv(here("Data", "RespoFiles","TPC", "respo_constant_temps.csv"))
phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta_all.csv"))
respo_constant_temps <- respo_constant_temps %>% left_join(phys_meta, by = "frag_ID") %>%
  filter(frag_ID != "B03") #remove extreme outlier
respo_pared_temps <- respo_constant_temps %>% 
  filter(temp_c_value == "24.5" | temp_c_value == "28" | temp_c_value == "31" | temp_c_value == "34")

#calculate NP:R inflection point
#create models
npr_mods <- respo_constant_temps|>
  nest(.by = frag_ID) %>%
  mutate(n   = map_int(data, nrow),
         fit = map(data, ~ lm(NPR ~ temp_c_value, data = .x)))

npr_mods_sp <- respo_constant_temps|>
  nest(.by = full_species) %>%
  mutate(n   = map_int(data, nrow),
         fit = map(data, ~ lm(NPR ~ temp_c_value, data = .x)))
#pull out coefficients
npr_coeffs <- npr_mods |>
  mutate(coeffs = map(fit, tidy)) |>
  dplyr::select(frag_ID, coeffs) |>
  unnest(coeffs) |>
  dplyr::select(frag_ID, term, estimate) |>
  pivot_wider(names_from = term, values_from = estimate) |>
  rename(intercept = `(Intercept)`, slope = temp_c_value)

npr_coeffs_sp <- npr_mods_sp |>
  mutate(coeffs = map(fit, tidy)) |>
  dplyr::select(full_species, coeffs) |>
  unnest(coeffs) |>
  dplyr::select(full_species, term, estimate) |>
  pivot_wider(names_from = term, values_from = estimate) |>
  rename(intercept = `(Intercept)`, slope = temp_c_value)

#create lines for plot
npr_preds_sp <- npr_mods_sp |>
  mutate(temp_seq = map(data, ~ tibble(temp_c_value = seq(min(20), max(40),length.out = 100))), #changed temp seq higher
    preds = map2(fit, temp_seq, ~ .y %>% mutate(NPR = predict(.x, newdata = .y)))) |>
  dplyr::select(full_species, preds) |>
  unnest(preds)

id_species <- respo_constant_temps |> dplyr::select(frag_ID, full_species)
npr_preds_frag <- npr_mods |>
  mutate(temp_seq = map(data, ~ tibble(temp_c_value = seq(min(20), max(40),length.out = 100))), #changed temp seq higher
         preds = map2(fit, temp_seq, ~ .y %>% mutate(NPR = predict(.x, newdata = .y)))) |>
  dplyr::select(frag_ID, preds) |>
  unnest(preds) |>
  left_join(id_species, relationship = "many-to-many")

#get intercepts
npr_intercepts_sp <- npr_coeffs_sp |>
  mutate(temp_at_NPR1 = (1 - intercept) / slope) #|>
  #left_join(respo_constant_temps) |>
  #group_by(full_species) |>
  #mutate(full_species = fct_reorder(full_species, temp_at_NPR1)) |>
  #filter(temp_c_value == "28") #just take one temp data since otherwise it duplicates
npr_intercepts <- npr_coeffs |>
  mutate(temp_at_NPR1 = (1 - intercept) / slope) |>
  left_join(respo_constant_temps) |>
  group_by(full_species) |>
  mutate(full_species = fct_reorder(full_species, temp_at_NPR1)) |>
  filter(temp_c_value == "28") |> #just take one temp data since otherwise it duplicates
  dplyr::select(frag_ID, full_species, temp_at_NPR1) 

#calculate NP:R percent less than 1
npr_pct_less_1 <- respo_constant_temps |>
  group_by(frag_ID) |>
  summarise(n_total = sum(!is.na(NPR)),
            n_below = sum(NPR <1, na.rm = T),
            prop_below_1 = n_below/n_total,
            pct_below_1 = 100 * prop_below_1,
            .groups = "drop") |>
  left_join(respo_constant_temps) |>
  group_by(full_species) |>
  mutate(full_species = fct_reorder(full_species, prop_below_1)) |>
  filter(temp_c_value == "28") |>  #just take one temp data since otherwise it duplicates
  dplyr::select(frag_ID, full_species, n_total,n_below,prop_below_1,pct_below_1) 

#save data to include later:
npr_metrics <- left_join(npr_pct_less_1, npr_intercepts)
npr_metrics[13,7] <- NA #remove extreme outlier of NP:R1 = 60
write.csv(npr_metrics, here("Data/Physiology/npr_metrics.csv"), row.names=FALSE)
npr_metrics <- read.csv(here("Data/Physiology/npr_metrics.csv"))

##### P:R plots #####
PR_plot <- ggplot() +
  geom_jitter(data = respo_constant_temps, aes(x = temp_c_value, y = NPR, color = full_species), width = 0.15, alpha = 0.8) +
  theme_classic(base_size = 22) +
  #facet_wrap(~full_species, ncol = 5, scales = "free") +
  facet_grid(~factor(full_species, levels = c("Favites complanata","Porites cylindrica","Pachyseris rugosa",
                                                     "Echinopora lamellosa","Acropora hyacinthus","Montipora aequituberculata",
                                                     "Montipora vietnamensis","Pocillopora eydouxi","Porites rus","Turbinaria frondens")))+
  facet_wrap(~full_species, ncol = 5, scales = "free") +
  geom_line(data = npr_preds_sp, aes(x = temp_c_value, y = NPR, color = full_species, group = full_species), inherit.aes = FALSE)+
  theme(legend.position = "none", strip.text = element_text(face = "italic"))+
  #geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(data = npr_intercepts, aes(xintercept = temp_at_NPR1), linetype = "dashed", color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_label(data = npr_intercepts, aes(x = temp_at_NPR1, y = 1.5, label = round(temp_at_NPR1,2)), 
             color = "black", size = 6, inherit.aes = FALSE) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #xlim(24, 36) +
  #ylim(0,4.2)+
  scale_y_continuous(breaks = scales::breaks_width(0.5)) +
  labs(y = "NP:R",
       x = "Temperature (°C)")
PR_plot
ggsave(here("Output", "Physiology", "PR_temp_species.pdf"), PR_plot, h = 8, w = 20)

npr_preds_sp %>% filter(full_species == "Acropora hyacinthus") %>%
  ggplot(aes(temp_c_value, NPR)) + geom_line()

#look at average NP:R where it becomes <1
PR_inflect <- ggplot(data = npr_intercepts) +
  geom_jitter(aes(x = temp_at_NPR1, y = full_species, color = full_species)) +
  theme_bw(base_size = 22) +
  #theme(axis.text.x = element_blank())+
  theme(legend.position = "none")+
  #facet_wrap(~species, ncol = 5, scales = "free") +
  #geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  labs(color = "Species",
       y = "",
       x = "Temperature (°C) \nat NP:R = 1")
PR_inflect

#Proportion of NP:R points less than 1
PR_less_1 <- ggplot(data = npr_pct_less_1) +
  geom_point(aes(x = prop_below_1, y = full_species, color = full_species)) +
  theme_bw(base_size = 22) +
  theme(legend.position = "none")+
  #theme(axis.text.x = element_blank())+
  #facet_wrap(~species, ncol = 5, scales = "free") +
  #geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  labs(color = "Species",
       y = "",
       x = "Proportion of NP:R < 1")
PR_less_1

#filled area plot for prop < 1 >
# Counts of points above/below 1 per species
pr_counts <- respo_constant_temps %>%
  group_by(full_species) %>%
  summarise(below_1 = sum(NPR < 1, na.rm = TRUE),
            above_1 = sum(NPR >= 1, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_longer(cols = c(below_1, above_1), names_to = "category", values_to = "n") %>%
  mutate(fill_group = if_else(category == "above_1", as.character(full_species), "below_1"))

# Order species by % below 1, for a cleaner visual gradient across the plot
species_order <- pr_counts %>%
  filter(category == "below_1") %>%
  arrange(n) %>%
  pull(full_species)

pr_counts <- pr_counts %>%
  mutate(fill_group = factor(fill_group, levels = c("below_1", names(sp_cols))))

fill_values <- c(sp_cols, "below_1" = "grey70")

filled_area_plot <- ggplot(pr_counts, aes(x = full_species, y = n, fill = fill_group)) +
  geom_col(position = "fill", width = 0.9) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(
    values = fill_values,
    breaks = names(sp_cols)) +
  scale_x_discrete(labels = function(x) paste0(x)) +
  labs(x = NULL, y = "NP:R > 1", fill = NULL) +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),legend.position = "none")

filled_area_plot

PR_area_1 <- ggarrange(PR_plot,filled_area_plot,nrow = 2, ncol = 1, 
                       labels = c("A","B"), heights = c(2,1),
                       font.label = list(size = 30, color = "black"))
ggsave(here("Output", "Physiology", "NPR_1_Pct_Temp.pdf"), PR_area_1, h = 15, w = 20)

np_allparams_jitter_pca <- ggarrange(np_rmax_plot,np_topt_plot,np_e_plot,np_breadth_plot,
                                     np_pca_spp_arrows,tpc_schematic,
                                     common.legend = T, legend = "right", ncol = 3, nrow=2,
                                     labels = c("A","B","C", "D","E","F"), 
                                     font.label = list(size = 30, color = "black"))

#####Ordination plots#####
#load metadata
phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta_all.csv"))
#load topt data
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4.csv"))
topt_df <- topt_df %>% filter(sample_ID != "B08_TPC") %>% 
  dplyr::select(-ctmin,-ctmax,-eh,-q10,-thermal_tolerance,-skewness,-thermal_safety_margin) %>%
  drop_na(e) #cleanup dataframe so things will run, all parameters taken out have too many NAs or infinity values
topt_matrix <- topt_df %>% dplyr::select(rmax:frag_ID) #generate data for matrix
#separate out data, make sure to put frag_ID as rownames so you can re-join with metadata later
topt_r_data <- topt_matrix %>% filter(PR == "Respiration") %>% dplyr::select(-PR) %>% column_to_rownames("frag_ID")
topt_np_data <- topt_matrix %>% filter(PR == "NetPhoto") %>% dplyr::select(-PR) %>% column_to_rownames("frag_ID")
topt_gp_data <- topt_matrix %>% filter(PR == "GrossPhoto") %>% dplyr::select(-PR) %>% column_to_rownames("frag_ID")
#make them matrices
topt_r <- as.matrix(topt_r_data)
topt_np <- as.matrix(topt_np_data)
topt_gp <- as.matrix(topt_gp_data)
#quick glance at data to look for strong patterns
pairs(x = topt_gp, gap = 0, cex.labels = 0.5) #look similar even between metrics
#scale across variable types
topt_r <- scale(topt_r)
topt_np <- scale(topt_np)
topt_g <- scale(topt_gp)
#generate pca data
pca_topt_r <- prcomp(topt_r)
pca_topt_np <- prcomp(topt_np)
pca_topt_gp <- prcomp(topt_gp)
#collapse PCA data
pc_axes_r <- as.data.frame(pca_topt_r$x)
pc_axes_np <- as.data.frame(pca_topt_np$x)
pc_axes_gp <- as.data.frame(pca_topt_gp$x)
#add frag_ID back
pc_axes_r$frag_ID <- rownames(pc_axes_r) 
pc_axes_np$frag_ID <- rownames(pc_axes_np) 
pc_axes_gp$frag_ID <- rownames(pc_axes_gp) 
#put it back with metadata
pca_r <- pc_axes_r %>% left_join(phys_meta, by = "frag_ID")
pca_np <- pc_axes_np %>% left_join(phys_meta, by = "frag_ID")
pca_gp <- pc_axes_gp %>% left_join(phys_meta, by = "frag_ID")

topt_fit_r <- envfit(pca_topt_r$x[, c("PC1", "PC2")], topt_r_data, permutations = 999, na.rm = TRUE)
topt_fit_np <- envfit(pca_topt_np$x[, c("PC1", "PC2")], topt_np_data, permutations = 999, na.rm = TRUE)
topt_fit_gp <- envfit(pca_topt_gp$x[, c("PC1", "PC2")], topt_gp_data, permutations = 999, na.rm = TRUE)

#topt_fit #can look at p-values of data to see what is significant in driving differences

topt_scores_r <- as.data.frame(scores(topt_fit_r, display = "vectors")) %>% mutate(variable = rownames(.))
topt_scores_np <- as.data.frame(scores(topt_fit_np, display = "vectors")) %>% mutate(variable = rownames(.))
topt_scores_gp <- as.data.frame(scores(topt_fit_gp, display = "vectors")) %>% mutate(variable = rownames(.))

#species
gp_pca_spp_arrows <- ggplot(pca_gp, aes(x = PC1, y = PC2, color = full_species, fill = full_species)) +
  geom_point(size = 3) +
  #scale_color_manual(values = morph_colors, name = "Morphology", labels = c("Branching/Tabular", "Massive/Submassive", "Encrusting/Plating")) +
  #scale_fill_manual(values = morph_colors, name = "Morphology", labels = c("Branching/Tabular", "Massive/Submassive", "Encrusting/Plating")) +
  scale_color_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  scale_fill_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  stat_ellipse(geom = "polygon", alpha = 0.1)+
  #stat_ellipse(aes(group = perf_imperf), level = 0.95, alpha = 0.5, color = "black", linewidth = 0.8) +
  geom_segment(data = topt_scores_r,
               aes(x = 0, y = 0, xend = PC1*6, yend = PC2*6),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  labs(title = "Gross Photosynthesis") +
  # geom_text(data = topt_scores_r,
  #          aes(x = PC1 * 3.5, y = PC2 * 3.5, label = variable),
  #          color = "black", size = 5, inherit.aes = FALSE) +
  theme_classic(base_size = 22)

gp_pca_spp_arrows

np_pca_spp_arrows <- ggplot(pca_np, aes(x = PC1, y = PC2, color = full_species, fill = full_species)) +
  geom_point(size = 3, alpha = 0.1) +
  #scale_color_manual(values = morph_colors, name = "Morphology", labels = c("Branching/Tabular", "Massive/Submassive", "Encrusting/Plating")) +
  #scale_fill_manual(values = morph_colors, name = "Morphology", labels = c("Branching/Tabular", "Massive/Submassive", "Encrusting/Plating")) +
  scale_color_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  scale_fill_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  stat_ellipse(geom = "polygon", alpha = 0.1, linewidth = 0.2)+
  stat_ellipse(level = 0.0001, geom = "point", shape = 21) +
  #stat_ellipse(aes(group = perf_imperf), level = 0.95, alpha = 0.5, color = "black", linewidth = 0.8) +
  geom_segment(data = topt_scores_r,
               aes(x = 0, y = 0, xend = PC1*4, yend = PC2*4),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  #labs(title = "Net Photosynthesis") +
  # geom_text(data = topt_scores_r,
  #          aes(x = PC1 * 3.5, y = PC2 * 3.5, label = variable),
  #          color = "black", size = 5, inherit.aes = FALSE) +
  theme_classic(base_size = 22)

np_pca_spp_arrows

ggsave(here("Output", "Physiology", "gp_topt_params_sp_arrows_nolabels.pdf"), gp_pca_spp_arrows, h = 8, w = 12)

#stats for beta dispersion and diversity
####TOPT
#generate distance matrix
np_dist <- vegdist(pca_topt_np$x, method = "euclidean")

set.seed(8)
#beta dispersion comparisons
bet.np <- betadisper(np_dist,pca_np$species)
anova(bet.np) 
#species: p= 0.5685
#plot(bet.phys)
#permutest(bet.phys, pairwise = TRUE, permutations = 999)

np_perm <- adonis2(np_dist ~ species, data = pca_np, permutations = 999)
np_perm
#species: R2 = 0.37422, F = 2.4585, p = 0.001 ***
pw_comps_np <- pairwise.adonis2(np_dist ~ species, data=pca_np, permutations = 999)

#table of pairwise comparison output
pw_list_np <- pw_comps_np[names(pw_comps_np) != "parent_call"]
pw_table_np <- do.call(rbind, lapply(names(pw_list_np), function(nm) {
  df <- as.data.frame(pw_list_np[[nm]])
  df$Term <- rownames(df)
  df$comparison <- nm
  df[, c("comparison", "Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
}))

write.csv(pw_table_np, here("Output/Physiology/PERMANOVA_ord_pairwise_comparisons_np_tpc_params.csv"))

#gp
#generate distance matrix
gp_dist <- vegdist(pca_topt_gp$x, method = "euclidean")

#beta dispersion comparisons
bet.gp <- betadisper(gp_dist,pca_gp$species)
anova(bet.gp) 
#species: p= 0.3234
#morphology: p = 0.7414
#plot(bet.phys)
#permutest(bet.phys, pairwise = TRUE, permutations = 999)

gp_perm <- adonis2(gp_dist ~ morphology, data = pca_gp, permutations = 999)
gp_perm
#species: R2 = 0.40042, F = 2.8939, p = 0.005 **
#morphology: R2 =  0.07046,   F = 1.7434 p = 0.173

pairwise.adonis2(gp_dist ~ morphology, data=pca_gp, permutations = 999)

#######All TPC plots####
#use NP dataset for now
all_data <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4.csv"))
gp_data <- all_data %>% filter(PR == "GrossPhoto")
np_data <- all_data %>% filter(PR == "NetPhoto")
r_data <- all_data %>% filter(PR == "Respiration")

#write for loop - for column, test species
columns <- c("rmax","topt","e","breadth")

# Optional: nicer y-axis labels per variable
y_labels <- list(
  rmax              = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
  topt              = "Thermal optimum (°C)",
  breadth           = "Breadth (°C)",
  e                 = expression("Activation energy" ~ (cal ~ mol^{-1}))
  # afdw_mg_cm2       = expression("Tissue Biomass" ~ (mg ~ cm^{-2})),
  # chla_ug_cm2_mean  = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})),
  # chla_pg_sym       = expression("Chlorophyll a" ~ (pg ~ cm^{-2})),
  # sym_cm2           = expression("Symbiont density" ~ (cells ~ cm^{-2})),
  # prot_ug_cm2       = expression("Protein" ~ (mu*g ~ cm^{-2}))
)

analyze_var_sp <- function(var, data) {
  var_sym <- sym(var)
  
  message("Processing variable: ", var)
  
  ## 1. Summarise by species for this variable
  summary_tbl <- data %>%
    group_by(full_species) %>%
    summarise(
      n    = sum(!is.na(!!var_sym)),
      mean = mean(!!var_sym, na.rm = TRUE),
      se   = se_fun(!!var_sym),
      .groups = "drop"
    )
  
  ## 2. Fit model: var ~ full_species
  form <- new_formula(lhs = expr(!!var_sym), rhs = expr(full_species))
  mod  <- lm(form, data = data)
  
  anova_tbl <- car::Anova(mod) %>%
    as.data.frame() %>%
    rownames_to_column("term")
  
  p_full <- anova_tbl %>%
    filter(term == "full_species") %>%
    pull(`Pr(>F)`) %>%
    first()
  
  performance::check_model(mod)
  
  ## 3. emmeans if significant
  emm_obj   <- NULL
  emm_pairs <- NULL
  
  if (!is.na(p_full) && p_full < 0.05) {
    emm_obj   <- emmeans::emmeans(mod, ~ full_species)
    emm_pairs <- pairs(emm_obj)
  }
  
  ## 3b. NEW — generate compact letter display and attach to summary_tbl
  letters_tbl <- if (!is.null(emm_obj)) {
    multcomp::cld(emm_obj, Letters = letters, adjust = "tukey") %>%
      as.data.frame() %>%
      transmute(full_species, group = trimws(.group))
  } else {
    # not significant overall -> no meaningful letters; leave blank
    summary_tbl %>% transmute(full_species, group = "")
  }
  
  summary_tbl <- summary_tbl %>% left_join(letters_tbl, by = "full_species")
  
  ## 4. Plot jitter + means + SE for this variable
  df_ordered <- data %>%
    left_join(summary_tbl, by = "full_species") %>%
    mutate(full_species = fct_reorder(full_species, mean))
  
  y_lab <- y_labels[[var]]
  if (is.null(y_lab)) y_lab <- var
  
  # reorder summary_tbl to match df_ordered's factor levels, so the text
  # layer's x positions line up with the jitter/point layers
  summary_tbl <- summary_tbl %>%
    mutate(full_species = factor(full_species, levels = levels(df_ordered$full_species)))
  
  p <- ggplot() +
    geom_jitter(
      data = df_ordered,
      aes(x = full_species, y = !!var_sym, color = full_species),
      width = 0.15, alpha = 0.8
    ) +
    geom_errorbar(
      data = summary_tbl,
      aes(x = full_species, ymin = mean - se, ymax = mean + se),
      width = 0.2, linewidth = 0.6
    ) +
    geom_point(
      data = summary_tbl,
      aes(x = full_species, y = mean),
      size = 2
    ) +
    stat_summary(
      data = df_ordered,
      aes(x = full_species, y = !!var_sym, label = group),
      fun = max, geom = "text", vjust = -0.8
    ) +
    theme_bw(base_size = 22) +
    theme(
      legend.position = "right",
      axis.text.x     = element_blank(),
      axis.title.x    = element_blank()
    ) +
    scale_color_manual(
      values = sp_cols,
      labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')"))
    ) +
    labs(
      x     = "Species",
      color = "Species",
      y     = y_lab
    )
  
  # ggsave(
  #   filename = here("Output", "Physiology", paste0(var, "_species_jitter.pdf")),
  #   plot     = p,
  #   h        = 5,
  #   w        = 10
  # )
  
  list(
    variable   = var,
    summary    = summary_tbl,
    model      = mod,
    anova      = anova_tbl,
    emm_obj    = emm_obj,
    emm_pairs  = emm_pairs,
    plot       = p
  )
}

# Run the pipeline for each variable using purrr
np_results <- columns %>%
  set_names() %>%
  map(~ analyze_var_sp(.x, np_data))

gp_results <- columns %>%
  set_names() %>%
  map(~ analyze_var_sp(.x, gp_data))

r_results <- columns %>%
  set_names() %>%
  map(~ analyze_var_sp(.x, r_data))

# Examples of accessing outputs:
# results[["rmax"]]$summary   # summary table
# results[["rmax"]]$anova     # ANOVA table
# results[["rmax"]]$emm_pairs # emmeans pairs (if significant)
# results[["rmax"]]$plot      # ggplot object

##Look through plots and annotate and check models

###NP PLOTS
#np rmax
np_rmax_plot <- np_results[["rmax"]]$plot +
  ylim(0.4,1.7)
np_rmax_plot
np_results[["rmax"]]$anova #0.0172721
check_model(np_results[["rmax"]]$model) #good
leveneTest(np_results[["rmax"]]$model)
np_results[["rmax"]]$emm_pairs
# Favites complanata - Montipora aequituberculata       0.5440 0.143 39   3.808  0.0155

#np topt
np_topt_plot <- np_results[["topt"]]$plot 
np_results[["topt"]]$anova #0.1096963
check_model(np_results[["topt"]]$model) #looks a bit weird but ns and transforms didn't help
leveneTest(np_results[["topt"]]$model)

#np e
np_e_plot <- np_results[["e"]]$plot + ylim(-0.4,1.2)
np_e_plot
np_results[["e"]]$anova #0.0006690892
check_model(np_results[["e"]]$model) #good
leveneTest(np_results[["e"]]$model)
np_results[["e"]]$emm_pairs

#np breadth
np_breadth_plot <- np_results[["breadth"]]$plot 
np_results[["breadth"]]$anova #0.5425902
check_model(np_results[["breadth"]]$model) #good
leveneTest(np_results[["breadth"]]$model)

###GP PLOTS
#gp rmax
gp_rmax_plot <- gp_results[["rmax"]]$plot +
  ylim(0.6,2.8)
gp_rmax_plot
gp_results[["rmax"]]$anova #0.008926322
check_model(gp_results[["rmax"]]$model) #good
leveneTest(gp_results[["rmax"]]$model)
gp_results[["rmax"]]$emm_pairs

#gp topt
gp_topt_plot <- gp_results[["topt"]]$plot 
check_model(gp_results[["topt"]]$model) #same as np topt
leveneTest(gp_results[["topt"]]$model)
gp_results[["topt"]]$anova #0.1013523

#gp e
gp_e_plot <- gp_results[["e"]]$plot + 
  ylim(0,1)
gp_e_plot
gp_results[["e"]]$anova #0.01694122
check_model(gp_results[["e"]]$model) #good
leveneTest(gp_results[["e"]]$model)
gp_results[["e"]]$emm_pairs

#gp breadth
gp_breadth_plot <- gp_results[["breadth"]]$plot 
check_model(gp_results[["breadth"]]$model)
leveneTest(gp_results[["breadth"]]$model)
gp_results[["breadth"]]$anova #0.09421194

##R PLOTS
#r rmax
r_rmax_plot <- r_results[["rmax"]]$plot
r_rmax_plot
r_results[["rmax"]]$anova #0.786303

#r topt
r_topt_plot <- r_results[["topt"]]$plot 
r_results[["topt"]]$anova #0.4577987

#r ctmax
r_ctmax_plot <- r_results[["ctmax"]]$plot
r_ctmax_plot
r_results[["ctmax"]]$anova #0.7630376

#r e
r_e_plot <- r_results[["e"]]$plot
r_e_plot
r_results[["e"]]$anova #0.2652668

#r breadth
r_breadth_plot <- r_results[["breadth"]]$plot 
r_results[["breadth"]]$anova #0.5881213

##NP plots with ordination
#add tpc schematic
tpc_schematic <- readRDS(here("Output/Okinawa_Map/tpc_schematic.rds"))
#put all plots together
np_allparams_jitter_pca <- ggarrange(np_rmax_plot,np_topt_plot,np_e_plot,np_breadth_plot,
                                 np_pca_spp_arrows,tpc_schematic,
                                 common.legend = T, legend = "right", ncol = 3, nrow=2,
                                 labels = c("A","B","C", "D","E","F"), 
                                 font.label = list(size = 30, color = "black"))
np_allparams_jitter_pca

ggsave(here("Output", "Physiology", "np_tpc_params_jitter_pca.pdf"), 
       np_allparams_jitter_pca, h = 10, w = 20)


# Examples of accessing outputs:
# results[["rmax"]]$summary   # summary table
# results[["rmax"]]$anova     # ANOVA table
# results[["rmax"]]$emm_pairs # emmeans pairs (if significant)
# results[["rmax"]]$plot      # ggplot object

#### Tile plot of correlations between physio and tpc data####
#read in data with both tpc and physio together
all_data <- read_csv(here("Data", "Physiology", "all_data_concatenated.csv"))
npr_metrics <- read.csv(here("Data/Physiology/npr_metrics.csv"))
all_data <- all_data |> 
  left_join(npr_metrics)
gp_data <- all_data %>% filter(PR == "GrossPhoto")
np_data <- all_data %>% filter(PR == "NetPhoto") #just use NP data for this analysis
r_data <- all_data %>% filter(PR == "Respiration")

physio_vars <- c("chla_ug_cm2_mean", "chla_pg_sym", "sym_cm2", "afdw_mg_cm2", "prot_ug_cm2")
tpc_vars    <- c("rmax", "topt", "e", "breadth", "temp_at_NPR1")

#np correlations
np_cor_pvals <- expand_grid(var1 = physio_vars, var2 = tpc_vars) |>
  mutate(fit      = map2(var1, var2, ~ lm(reformulate(.x, .y), data = np_data)),
         tidy_fit = map(fit, tidy)) |>
  mutate(p.value = map_dbl(tidy_fit, ~ .x |> filter(term != "(Intercept)") |> pull(p.value))) |>
  dplyr::select(var1, var2, p.value)

#create correlation matrix between variables
np_physio_cor <- np_data |>
  dplyr::select(
    rmax, topt, e, breadth, temp_at_NPR1,
    afdw_mg_cm2,
    chla_ug_cm2_mean,
    chla_log,
    chla_pg_sym,
    sym_cm2,
    prot_ug_cm2) |>
  na.omit() |>
  cor()

#pivot
np_cor_data <- np_physio_cor |>
  as.data.frame() |>
  rownames_to_column("var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "correlation") |> 
  left_join(np_cor_pvals, by = c("var1", "var2")) |>
  mutate(star = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "")) |>
  mutate(corr_star = paste(round(correlation,2), star, sep = "\n"))

np_cor_data <- np_cor_data |>
  filter(var1 %in% c("chla_ug_cm2_mean","chla_pg_sym","sym_cm2","afdw_mg_cm2","prot_ug_cm2")) |>
  filter(var2 %in% c("rmax", "topt","e","breadth","temp_at_NPR1")) |>
  mutate(physio = case_when(
    var1 == "chla_ug_cm2_mean" ~ "Chl a",
    var1 == "chla_pg_sym" ~ "Chl a per sym",
    var1 == "sym_cm2" ~ "Sym density",
    var1 == "afdw_mg_cm2" ~ "Biomass",
    var1 == "prot_ug_cm2" ~ "Protein")) |>
  mutate(tpc = case_when(
    var2 == "rmax" ~ "Rmax",
    var2 == "topt" ~ "Topt",
    var2 == "e" ~ "e",
    var2 == "breadth" ~ "breadth",
    #var2 == "pct_below_1" ~ "% NP:R < 1",
    var2 == "temp_at_NPR1" ~ "°C at NP:R = 1")) |>
  mutate(physio=fct_relevel(physio,c("Chl a","Chl a per sym","Sym density","Biomass","Protein"))) |>
  mutate(tpc=fct_relevel(tpc,c("breadth","e","Topt","Rmax","°C at NP:R = 1")))

np_physio_cor_plot <- ggplot(np_cor_data, aes(x = physio, y = tpc)) +
  geom_tile(aes(fill = correlation), color = "white",linewidth = 1) +
  #geom_text(aes(label = round(correlation, 2))) +
  geom_text(aes(label = corr_star)) +
  scale_fill_gradient2(low = "orangered3", mid = "white", high = "dodgerblue3", midpoint = 0) +
  labs(title = "Net Photosynthesis",
       x = "", y = "", fill = "Correlation") +
  theme_classic(base_size = 22) +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
np_physio_cor_plot

#gp correlations
gp_cor_pvals <- expand_grid(var1 = physio_vars, var2 = tpc_vars) |>
  mutate(fit      = map2(var1, var2, ~ lm(reformulate(.x, .y), data = gp_data)),
         tidy_fit = map(fit, tidy)) |>
  mutate(p.value = map_dbl(tidy_fit, ~ .x |> filter(term != "(Intercept)") |> pull(p.value))) |>
  dplyr::select(var1, var2, p.value)

#create correlation matrix between variables
gp_physio_cor <- gp_data |>
  dplyr::select(
    rmax, topt, e, breadth,
    afdw_mg_cm2,
    chla_ug_cm2_mean,
    chla_log,
    chla_pg_sym,
    sym_cm2,
    prot_ug_cm2) |>
  na.omit() |>
  cor()

#pivot
gp_cor_data <- gp_physio_cor |>
  as.data.frame() |>
  rownames_to_column("var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "correlation") |> 
  left_join(gp_cor_pvals, by = c("var1", "var2")) |>
  mutate(star = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "")) |>
  mutate(corr_star = paste(round(correlation,2), star, sep = "\n"))

gp_cor_data <- gp_cor_data |>
  filter(var1 %in% c("chla_ug_cm2_mean","chla_pg_sym","sym_cm2","afdw_mg_cm2","prot_ug_cm2")) |>
  filter(var2 %in% c("rmax", "topt","e","breadth")) |>
  mutate(physio = case_when(
    var1 == "chla_ug_cm2_mean" ~ "Chl a",
    var1 == "chla_pg_sym" ~ "Chl a per sym",
    var1 == "sym_cm2" ~ "Sym density",
    var1 == "afdw_mg_cm2" ~ "Biomass",
    var1 == "prot_ug_cm2" ~ "Protein")) |>
  mutate(tpc = case_when(
    var2 == "rmax" ~ "Rmax",
    var2 == "topt" ~ "Topt",
    var2 == "e" ~ "e",
    var2 == "breadth" ~ "breadth")) |>
    #var2 == "pct_below_1" ~ "% gp:R < 1",
    #var2 == "temp_at_NPR1" ~ "°C at gp:R = 1")) |>
  mutate(physio=fct_relevel(physio,c("Chl a","Chl a per sym","Sym density","Biomass","Protein"))) |>
  mutate(tpc=fct_relevel(tpc,c("breadth","e","Topt","Rmax")))

gp_physio_cor_plot <- ggplot(gp_cor_data, aes(x = physio, y = tpc)) +
  geom_tile(aes(fill = correlation), color = "white",linewidth = 1) +
  #geom_text(aes(label = round(correlation, 2))) +
  geom_text(aes(label = corr_star)) +
  scale_fill_gradient2(low = "orangered3", mid = "white", high = "dodgerblue3", midpoint = 0) +
  labs(title = "Gross Photosynthesis",
       x = "", y = "", fill = "Correlation") +
  theme_classic(base_size = 22) +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gp_physio_cor_plot

#R correlations
r_cor_pvals <- expand_grid(var1 = physio_vars, var2 = tpc_vars) |>
  mutate(fit      = map2(var1, var2, ~ lm(reformulate(.x, .y), data = r_data)),
         tidy_fit = map(fit, tidy)) |>
  mutate(p.value = map_dbl(tidy_fit, ~ .x |> filter(term != "(Intercept)") |> pull(p.value))) |>
  dplyr::select(var1, var2, p.value)

#create correlation matrix between variables
r_physio_cor <- r_data |>
  dplyr::select(
    rmax, topt, e, breadth,
    afdw_mg_cm2,
    chla_ug_cm2_mean,
    chla_log,
    chla_pg_sym,
    sym_cm2,
    prot_ug_cm2) |>
  na.omit() |>
  cor()

#pivot
r_cor_data <- r_physio_cor |>
  as.data.frame() |>
  rownames_to_column("var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "correlation") |> 
  left_join(r_cor_pvals, by = c("var1", "var2")) |>
  mutate(star = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "")) |>
  mutate(corr_star = paste(round(correlation,2), star, sep = "\n"))

r_cor_data <- r_cor_data |>
  filter(var1 %in% c("chla_ug_cm2_mean","chla_pg_sym","sym_cm2","afdw_mg_cm2","prot_ug_cm2")) |>
  filter(var2 %in% c("rmax", "topt","e","breadth","temp_at_rR1")) |>
  mutate(physio = case_when(
    var1 == "chla_ug_cm2_mean" ~ "Chl a",
    var1 == "chla_pg_sym" ~ "Chl a per sym",
    var1 == "sym_cm2" ~ "Sym density",
    var1 == "afdw_mg_cm2" ~ "Biomass",
    var1 == "prot_ug_cm2" ~ "Protein")) |>
  mutate(tpc = case_when(
    var2 == "rmax" ~ "Rmax",
    var2 == "topt" ~ "Topt",
    var2 == "e" ~ "e",
    var2 == "breadth" ~ "breadth")) |> 
    #var2 == "pct_below_1" ~ "% r:R < 1",
    #var2 == "temp_at_rR1" ~ "°C at r:R = 1")) |>
  mutate(physio=fct_relevel(physio,c("Chl a","Chl a per sym","Sym density","Biomass","Protein"))) |>
  mutate(tpc=fct_relevel(tpc,c("breadth","e","Topt","Rmax")))

r_physio_cor_plot <- ggplot(r_cor_data, aes(x = physio, y = tpc)) +
  geom_tile(aes(fill = correlation), color = "white",linewidth = 1) +
  #geom_text(aes(label = round(correlation, 2))) +
  geom_text(aes(label = corr_star)) +
  scale_fill_gradient2(low = "orangered3", mid = "white", high = "dodgerblue3", midpoint = 0) +
  labs(title = "Respiration",
       x = "", y = "", fill = "Correlation") +
  theme_classic(base_size = 22) +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
r_physio_cor_plot

#for looking at models if u want
# mod <- lm(temp_at_NPR1 ~ afdw_mg_cm2, data = np_data)
# Anova(mod)
# check_model(mod)

r_gp_physio_cor <- ggarrange(gp_physio_cor_plot, r_physio_cor_plot,
                                  common.legend = F,
                                  ncol = 2, nrow=1, labels = c("A","B"), font.label = list(size = 30, color = "black"))
ggsave(here("Output", "Physiology", "gp_r_tpc_physio_correlations.pdf"), r_gp_physio_cor, h = 8, w = 20)

ggsave(here("Output", "Physiology", "np_tpc_physio_correlations.pdf"), np_physio_cor_plot, h = 8, w = 10)

