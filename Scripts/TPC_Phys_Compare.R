###### Code for TPC and Physiology comparisons ####### 
### Created by: Maya Powell
#### Last updated on: 2025-07-26

############## Introduction to code/script ####################
## this script will help us process the raw data gathered during respirometry runs. 
## need to change for specific project/experimental variables 

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

###### Initial Data Read In ########
#TPC data with only the seven species that we have physio data for
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4_7sp.csv"))

#additional setup
sp_cols <- c(
  "Acropora hyacinthus" = '#ba7999',
  "Echinopora lamellosa" = '#dd4124',
  "Favites complanata"   = '#ed8b00',
  "Montipora aequituberculata" = '#edd746',
  "Pocillopora eydouxi" = '#45681e',
  "Porites rus" = '#7bbcd5',
  "Turbinaria frondens" = '#00496f'
)


##### dry weight #####
avg_DW <- read.csv(here("Data", "Physiology", "Average_Dry_Weight.csv"))
avg_DW <- avg_DW %>% rename_at('species_long', ~'full_species')

library(powerjoin)
power_full_join(df1, df2, by = "Species", conflict = coalesce_xy)

DW_topt <- power_full_join(avg_DW, topt_df, by = "frag_ID", conflict = coalesce_xy) %>% 
  drop_na() %>%
  filter(PR != "Respiration")

#scatterplot of dw and topt params

#topt
dw_topt_scatter <- ggplot(DW_topt) +
  geom_point(aes(x = topt, y = dw_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free_y")+
  ylim(25,250)+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), x = "Thermal optimum (°C)", color = "Species")
dw_topt_scatter

dw_topt_scatter <- dw_topt_scatter +
  #species-specific regression lines
  geom_smooth(aes(x = topt, y = dw_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(x = topt, y = dw_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_topt_scatter_reglines.pdf"), dw_topt_scatter, h = 8, w = 10)

#rmax
dw_rmax_scatter <- ggplot(DW_topt) +
  geom_point(aes(x = rmax, y = dw_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  ylim(25,250)+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       x = expression("Rate Maximum" ~ (mu*mol ~ cm^{-2} ~ h^{-1})), 
       color = "Species")
dw_rmax_scatter

dw_rmax_scatter <- dw_rmax_scatter +
  #species-specific regression lines
  geom_smooth(aes(x = rmax, y = dw_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(x = rmax, y = dw_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_rmax_scatter_reglines.pdf"), dw_rmax_scatter, h = 8, w = 10)

#e
dw_e_scatter <- ggplot(DW_topt) +
  geom_point(aes(x = e, y = dw_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  ylim(25,250)+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       x = ("Activation energy (eV)"), 
       color = "Species")
dw_e_scatter

dw_e_scatter <- dw_e_scatter +
  #species-specific regression lines
  geom_smooth(aes(x = e, y = dw_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(x = e, y = dw_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_e_scatter_reglines.pdf"), dw_e_scatter, h = 8, w = 10)

#breadth
dw_breadth_scatter <- ggplot(DW_topt) +
  geom_point(aes(x = breadth, y = dw_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  ylim(25,250)+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       x = expression("Thermal performance breadth (°C)"), 
       color = "Species")
dw_breadth_scatter

dw_breadth_scatter <- dw_breadth_scatter +
  #species-specific regression lines
  geom_smooth(aes(x = breadth, y = dw_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(x = breadth, y = dw_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_breadth_scatter_reglines.pdf"), dw_breadth_scatter, h = 8, w = 10)


##### chlorophyll a #####
chla_avg <- read.csv(here("Data", "Physiology", "Chla_avg.csv"))

library(powerjoin)

chla_topt <- power_full_join(chla_avg, topt_df, by = "frag_ID", conflict = coalesce_xy) %>% 
  drop_na() %>%
  filter(PR != "Respiration")

#scatterplot of chla and topt params

#topt
chla_topt_scatter <- ggplot(chla_topt) +
  geom_point(aes(x = topt, y = chla_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free_y")+
  ylim(25,250)+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), x = "Thermal optimum (°C)", color = "Species")
chla_topt_scatter

chla_topt_scatter <- chla_topt_scatter +
  #species-specific regression lines
  geom_smooth(aes(x = topt, y = chla_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(x = topt, y = chla_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_topt_scatter_reglines.pdf"), chla_topt_scatter, h = 8, w = 10)

#rmax
chla_rmax_scatter <- ggplot(chla_topt) +
  geom_point(aes(x = rmax, y = chla_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  ylim(25,250)+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       x = expression("Rate Maximum" ~ (mu*mol ~ cm^{-2} ~ h^{-1})), 
       color = "Species")
chla_rmax_scatter

chla_rmax_scatter <- chla_rmax_scatter +
  #species-specific regression lines
  geom_smooth(aes(x = rmax, y = chla_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(x = rmax, y = chla_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_rmax_scatter_reglines.pdf"), chla_rmax_scatter, h = 8, w = 10)

#e
chla_e_scatter <- ggplot(chla_topt) +
  geom_point(aes(x = e, y = chla_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  ylim(25,250)+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       x = ("Activation energy (eV)"), 
       color = "Species")
chla_e_scatter

chla_e_scatter <- chla_e_scatter +
  #species-specific regression lines
  geom_smooth(aes(x = e, y = chla_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(x = e, y = chla_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_e_scatter_reglines.pdf"), chla_e_scatter, h = 8, w = 10)

#breadth
chla_breadth_scatter <- ggplot(chla_topt) +
  geom_point(aes(x = breadth, y = chla_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  ylim(25,250)+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       x = expression("Thermal performance breadth (°C)"), 
       color = "Species")
chla_breadth_scatter

chla_breadth_scatter <- chla_breadth_scatter +
  #species-specific regression lines
  geom_smooth(aes(x = breadth, y = chla_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(x = breadth, y = chla_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_breadth_scatter_reglines.pdf"), chla_breadth_scatter, h = 8, w = 10)

