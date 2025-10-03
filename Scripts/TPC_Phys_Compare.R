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


#DW chla scatter
avg_DW <- read.csv(here("Data", "Physiology", "Average_Dry_Weight.csv"))
avg_DW <- avg_DW %>% rename_at('species_long', ~'full_species')
avg_DW <- avg_DW %>% mutate(dw_log = log(dw_mg_cm2))

chla_avg <- read.csv(here("Data", "Physiology", "Chla_avg.csv"))
chla_avg <- chla_avg %>% mutate(chla_log = log(chla_ug_cm2_mean))

chla_dw <- power_full_join(chla_avg, avg_DW, by = "frag_ID", conflict = coalesce_xy) %>% 
  drop_na()

dw_chla_scat <- ggplot(chla_dw) +
  geom_point(aes(x = dw_mg_cm2, y = chla_ug_cm2_mean, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  labs(x = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       y = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})), 
       color = "Species") +
  geom_smooth(aes(x = dw_mg_cm2, y = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black", linewidth = 1.1)
dw_chla_scat

ggsave(here("Output", "Physiology", "dryweight_chla_scatter.pdf"), dw_chla_scat, h = 6, w = 10)

##### dry weight #####
avg_DW <- read.csv(here("Data", "Physiology", "Average_Dry_Weight.csv"))
avg_DW <- avg_DW %>% rename_at('species_long', ~'full_species')

DW_topt <- power_full_join(avg_DW, topt_df, by = "frag_ID", conflict = coalesce_xy) %>% 
  drop_na() %>%
  filter(PR != "Respiration")

#scatterplot of dw and topt params

#topt
dw_topt_scatter <- ggplot(filter(DW_topt, PR == "NetPhoto")) +
  geom_point(aes(x = dw_log, y = topt, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #facet_wrap(~PR, scales = "free_y")+
  labs(x = expression("log(Dry weight" ~ (mg ~ cm^{-2})~")"), y = "Thermal optimum (°C)", color = "Species")
dw_topt_scatter

dw_topt_scatter <- dw_topt_scatter +
  #species-specific regression lines
  #geom_smooth(aes(y = topt, x = dw_log, color = full_species),
  #            method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = topt, x = dw_log, group = 1),
              method = "lm", se = FALSE, color = "black", linewidth = 1.1)

by(DW_topt, DW_topt$PR, function(d) summary(lm(topt ~ dw_log, data = d))) #*
library(performance)
check_model(lm(topt ~ dw_log, data = filter(DW_topt, PR == "NetPhoto")))
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(topt ~ dw_mg_cm2, data = d)))

#topt net photo
#dw_log        0.5526     0.2563   2.156   0.0395 * 

#topt
#DW_topt$PR: NetPhoto
#dw_mg_cm2    0.006035   0.002387   2.528   0.0172 *  significant
#none individual sig dif

ggsave(here("Output", "Physiology", "dryweight_topt_log_np.pdf"), dw_topt_scatter, h = 6, w = 10)

#rmax
dw_rmax_scatter <- ggplot(filter(DW_topt, PR == "NetPhoto")) +
  geom_point(aes(y = rmax, x = dw_log, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #facet_wrap(~PR, scales = "free")+
  labs(x = expression("log(Dry weight" ~ (mg ~ cm^{-2})~")"), 
       y = expression("Rmax" ~ (mu*mol ~ cm^{-2} ~ h^{-1})), 
       color = "Species")
dw_rmax_scatter

dw_rmax_scatter <- dw_rmax_scatter +
  #species-specific regression lines
  # geom_smooth(aes(y = rmax, x = dw_mg_cm2, color = full_species),
  #             method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = rmax, x = dw_log, group = 1),
              method = "lm", se = FALSE, color = "black", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_rmax_log.pdf"), dw_rmax_scatter, h = 4, w = 12)

by(DW_topt, DW_topt$PR, function(d) summary(lm(rmax ~ dw_log, data = d))) #ns
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(rmax ~ dw_mg_cm2, data = d)))
#none sig

#e
dw_e_scatter <- ggplot(filter(DW_topt, PR == "NetPhoto")) +
  geom_point(aes(y = e, x = dw_log, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #facet_wrap(~PR, scales = "free")+
  labs(x = expression("log(Dry weight" ~ (mg ~ cm^{-2})~")"), 
       y = ("E (eV)"), 
       color = "Species")
dw_e_scatter

dw_e_scatter <- dw_e_scatter +
  #species-specific regression lines
  #geom_smooth(aes(y = e, x = dw_mg_cm2, color = full_species),
  #            method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = e, x = dw_log, group = 1),
              method = "lm", se = FALSE, color = "black", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_e_log.pdf"), dw_e_scatter, h = 4, w = 12)

by(DW_topt, DW_topt$PR, function(d) summary(lm(e ~ dw_log, data = d)))
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(e ~ dw_mg_cm2, data = d)))

#e DW_topt$PR: GrossPhoto
#dw_log      -0.16421    0.07198  -2.281  0.03307 * 

#e
#DW_topt$PR: GrossPhoto
#dw_mg_cm2   -0.0019061  0.0007072  -2.695   0.0136 *  
#none significant


#breadth
dw_breadth_scatter <- ggplot(filter(DW_topt, PR == "NetPhoto")) +
  geom_point(aes(y = breadth, x = dw_mg_cm2, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  coord_transform(x = "log")+
  #facet_wrap(~PR, scales = "free")+
  geom_smooth(aes(y = breadth, x = dw_mg_cm2, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(x = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       y = expression("Breadth (°C)"), 
       color = "Species")
dw_breadth_scatter

dw_breadth_scatter <- dw_breadth_scatter +
  #species-specific regression lines
  #geom_smooth(aes(y = breadth, x = dw_mg_cm2, color = full_species),
  #            method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = breadth, x = dw_log, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_breadth_log.pdf"), dw_breadth_scatter, h = 4, w = 12)

by(DW_topt, DW_topt$PR, function(d) summary(lm(breadth ~ dw_log, data = d))) #ns
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(breadth ~ dw_mg_cm2, data = d)))

#breadth DW_topt$PR: NetPhoto
#dw_log        0.6504     0.2504   2.598   0.0146 *  

#breadth
#DW_topt$PR: NetPhoto
#dw_mg_cm2   0.005614   0.002442   2.299   0.0289 * 

#breadth
#NetPhoto
#Turbinaria frondens
#dw_mg_cm2   0.035475   0.008652    4.10   0.0262 *


library(ggpubr)
dw_topt_log <- ggarrange(dw_topt_scatter, dw_rmax_scatter, dw_e_scatter, dw_breadth_scatter, 
                           nrow = 2, ncol = 2, legend = "right", common.legend = TRUE)
dw_topt_log

ggsave(here("Output","TPC","Graphs","np_topt_plots.pdf"), dw_topt_log, h = 8, w = 12)


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
