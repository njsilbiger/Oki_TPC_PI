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

###### Initial Data Read In ########
#TPC data with only the seven species that we have physio data for
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4_7sp.csv"))

#additional setup
sp_cols <- c(
  "Acropora hyacinthus" = '#ba7999',
  "Echinopora lamellosa" = '#dd4124',
  "Favites complanata"   = '#ed8b00',
  "Montipora aequituberculata" = '#edd746',
  "Montipora vietnamensis" = '#89689d',
  "Pachyseris rugosa" = '#d0e2af',
  "Pocillopora eydouxi" = '#45681e',
  "Porites cylindrica" = '#f2af4a',
  "Porites rus" = '#7bbcd5',
  "Turbinaria frondens" = '#00496f'
)


##### dry weight #####
avg_DW <- read.csv(here("Data", "Physiology", "Average_Dry_Weight.csv"))
avg_DW <- avg_DW %>% rename_at('species_long', ~'full_species')

DW_topt <- power_full_join(avg_DW, topt_df, by = "frag_ID", conflict = coalesce_xy) %>% 
  drop_na() %>%
  filter(PR != "Respiration")

#scatterplot of dw and topt params

#topt
dw_topt_scatter <- ggplot(DW_topt) +
  geom_point(aes(y = topt, x = dw_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free_y")+
  labs(x = expression("Dry weight" ~ (mg ~ cm^{-2})), y = "Thermal optimum (°C)", color = "Species")
dw_topt_scatter

dw_topt_scatter <- dw_topt_scatter +
  #species-specific regression lines
  geom_smooth(aes(y = topt, x = dw_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = topt, x = dw_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

by(DW_topt, DW_topt$PR, function(d) summary(lm(topt ~ dw_mg_cm2, data = d)))
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(topt ~ dw_mg_cm2, data = d)))

#topt
#DW_topt$PR: NetPhoto
#dw_mg_cm2    0.006035   0.002387   2.528   0.0172 *  significant
#none individual sig dif

ggsave(here("Output", "Physiology", "dryweight_topt_scatter_reglines.pdf"), dw_topt_scatter, h = 4, w = 12)

#rmax
dw_rmax_scatter <- ggplot(DW_topt) +
  geom_point(aes(y = rmax, x = dw_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       x = expression("Rate Maximum" ~ (mu*mol ~ cm^{-2} ~ h^{-1})), 
       color = "Species")
dw_rmax_scatter

dw_rmax_scatter <- dw_rmax_scatter +
  #species-specific regression lines
  # geom_smooth(aes(y = rmax, x = dw_mg_cm2, color = full_species),
  #             method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = rmax, x = dw_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_rmax_scatter_reg.pdf"), dw_rmax_scatter, h = 4, w = 12)

by(DW_topt, DW_topt$PR, function(d) summary(lm(rmax ~ dw_mg_cm2, data = d))) #ns
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(rmax ~ dw_mg_cm2, data = d)))
#none sig

#e
dw_e_scatter <- ggplot(DW_topt) +
  geom_point(aes(y = e, x = dw_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       x = ("Activation energy (eV)"), 
       color = "Species")
dw_e_scatter

dw_e_scatter <- dw_e_scatter +
  #species-specific regression lines
  geom_smooth(aes(y = e, x = dw_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = e, x = dw_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_e_scatter_reglines.pdf"), dw_e_scatter, h = 4, w = 12)

by(DW_topt, DW_topt$PR, function(d) summary(lm(e ~ dw_mg_cm2, data = d)))
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(e ~ dw_mg_cm2, data = d)))

#e
#DW_topt$PR: GrossPhoto
#dw_mg_cm2   -0.0019061  0.0007072  -2.695   0.0136 *  
#none significant


#breadth
dw_breadth_scatter <- ggplot(DW_topt) +
  geom_point(aes(y = breadth, x = dw_mg_cm2, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  labs(y = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       x = expression("Thermal performance breadth (°C)"), 
       color = "Species")
dw_breadth_scatter

dw_breadth_scatter <- dw_breadth_scatter +
  #species-specific regression lines
  geom_smooth(aes(y = breadth, x = dw_mg_cm2, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = breadth, x = dw_mg_cm2, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_breadth_scatter_reglines.pdf"), dw_breadth_scatter, h = 8, w = 10)

by(DW_topt, DW_topt$PR, function(d) summary(lm(breadth ~ dw_mg_cm2, data = d))) #ns
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(breadth ~ dw_mg_cm2, data = d)))

#breadth
#DW_topt$PR: NetPhoto
#dw_mg_cm2   0.005614   0.002442   2.299   0.0289 * 

#breadth
#NetPhoto
#Turbinaria frondens
#dw_mg_cm2   0.035475   0.008652    4.10   0.0262 *

##### chlorophyll a #####
chla_avg <- read.csv(here("Data", "Physiology", "Chla_avg.csv"))

chla_topt <- power_full_join(chla_avg, topt_df, by = "frag_ID", conflict = coalesce_xy) %>% 
  drop_na() %>%
  filter(PR != "Respiration")

#scatterplot of chla and topt params

#topt
chla_topt_scatter <- ggplot(chla_topt) +
  geom_point(aes(y = topt, x = chla_ug_cm2_mean, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free_y")+
  #ylim(25,250)+
  labs(x = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})), 
       y = "Thermal optimum (°C)", color = "Species")
chla_topt_scatter

chla_topt_scatter <- chla_topt_scatter +
  #species-specific regression lines
  #geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, color = full_species),
  #            method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "chla_topt_scatter_reg.pdf"), chla_topt_scatter, h = 4, w = 12)

by(chla_topt, chla_topt$PR, function(d) summary(lm(topt ~ chla_ug_cm2_mean, data = d))) #ns
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(topt ~ chla_ug_cm2_mean, data = d))) #ns

#rmax
chla_rmax_scatter <- ggplot(chla_topt) +
  geom_point(aes(y = rmax, x = chla_ug_cm2_mean, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  labs(x = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})), 
       y = expression("Rate Maximum" ~ (mu*mol ~ cm^{-2} ~ h^{-1})), 
       color = "Species")
chla_rmax_scatter

chla_rmax_scatter <- chla_rmax_scatter +
  #species-specific regression lines
  geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "chla_rmax_scatter_reglines.pdf"), chla_rmax_scatter, h = 4, w = 12)

by(chla_topt, chla_topt$PR, function(d) summary(lm(rmax ~ chla_ug_cm2_mean, data = d))) #ns
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(rmax ~ chla_ug_cm2_mean, data = d)))

#rmax chla
#NetPhoto Echinopora lamellosa
#chla_ug_cm2_mean  0.13595    0.01915   7.098  0.00575 **

#e
chla_e_scatter <- ggplot(chla_topt) +
  geom_point(aes(y = e, x = chla_ug_cm2_mean, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  labs(x = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})), 
       y = ("Activation energy (eV)"), 
       color = "Species")
chla_e_scatter

chla_e_scatter <- chla_e_scatter +
  #species-specific regression lines
  geom_smooth(aes(y = e, x = chla_ug_cm2_mean, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = e, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "chla_e_scatter_reglines.pdf"), chla_e_scatter, h = 4, w = 12)

by(chla_topt, chla_topt$PR, function(d) summary(lm(e ~ chla_ug_cm2_mean, data = d))) #ns
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(e ~ chla_ug_cm2_mean, data = d))) #ns


#breadth
chla_breadth_scatter <- ggplot(chla_topt) +
  geom_point(aes(y = breadth, x = chla_ug_cm2_mean, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  labs(x = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})), 
       y = expression("Thermal performance breadth (°C)"), 
       color = "Species")
chla_breadth_scatter


chla_breadth_scatter <- chla_breadth_scatter +
  #species-specific regression lines
  geom_smooth(aes(y = breadth, x = chla_ug_cm2_mean, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = breadth, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "chla_breadth_scatter_reglines.pdf"), chla_breadth_scatter, h = 4, w = 12)

by(chla_topt, chla_topt$PR, function(d) summary(lm(breadth ~ chla_ug_cm2_mean, data = d))) #ns
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(breadth ~ chla_ug_cm2_mean, data = d))) #ns
