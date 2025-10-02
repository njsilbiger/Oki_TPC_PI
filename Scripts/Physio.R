######################################
##### Coral physiology Okinawa 2025###
######################################
#Maya Powell
#Created August 2025
#Last edited October 1 2025

#load libraries
library(tidyverse)
library(here)
library(dplyr)
library(car)
library(PNWColors)

#load data for all
meta <- read_csv(here("Data","Physiology","Physio_meta.csv"))
#assign colors to species
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

# pnw_palettes <- list(
#   Starfish = rbind(c('#24492e', '#015b58', '#2c6184', '#59629b', '#89689d', '#ba7999', '#e69b99'),c(7,4,5,3,1,6,2)),
#   Shuksan = rbind(c('#33271e', '#74677e', '#ac8eab', '#d7b1c5', '#ebbdc8', '#f2cec7', '#f8e3d1', '#fefbe9'),c(7,2,4,1,6,3,5,8)),
#   Bay = rbind(c('#00496f', '#0f85a0', '#edd746', '#ed8b00', '#dd4124'),c(4,1,3,5,2)),
#   Winter = rbind(c('#2d2926', '#33454e', '#537380', '#81a9ad', '#ececec'),c(1,4,5,2,3)),
#   Lake = rbind(c('#362904', '#54450f', '#45681e', '#4a9152', '#64a8a8', '#85b6ce', '#cde5f9', '#eef3ff'),c(4,8,7,2,6,1,3,5)),
#   Sunset = rbind(c('#41476b', '#675478', '#9e6374', '#c67b6f', '#de9b71', '#efbc82', '#fbdfa2'),c(3,5,1,7,4,6,2)),
#   Shuksan2 = rbind(c('#5d74a5', '#b0cbe7', '#fef7c7', '#eba07e', '#a8554e'),c(2,4,1,5,3)),
#   Cascades = rbind(c("#2d4030","#516823","#dec000","#e2e260","#677e8e","#88a2b9"),c(4,1,5,2,6,3)),
#   Sailboat = rbind(c('#6e7cb9', '#7bbcd5', '#d0e2af', '#f5db99', '#e89c81', '#d2848d'),c(1,4,6,2,5,3)),
#   Moth = rbind(c('#4a3a3b', '#984136', '#c26a7a', '#ecc0a1', '#f0f0e4'),c(4,1,2,3,5)),
#   Spring = rbind(c('#d8aedd', '#bf9bdd', '#cb74ad', '#e69e9c', '#ffc3a3', '#fbe4c6'),c(1,5,2,4,3,6)),
#   Mushroom = rbind(c('#4f412b', '#865a3c', '#ba783e', '#e69c4c', '#fbcc74', '#fffbda'),c(6,1,4,2,3,5)),
#   Sunset2 = rbind(c('#1d457f', '#61599d', '#c36377', '#eb7f54', '#f2af4a'),c(5,1,2,4,3)),
#   Anemone = rbind(c("#009474" ,"#11c2b5" ,"#72e1e1", "#f1f4ee" ,"#efddcf", "#dcbe9b" ,"#b0986c"),c(3, 5, 1 ,7, 2, 6, 4))
# )
##### Dry weight#####

#load data
avg_DW <- read.csv(here("Data", "Physiology", "Average_Dry_Weight.csv"))

#code to calculate average dry weight per fragment
#needed initially but not after writing average dry weight datasheet
# #read in metadata
# phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta.csv"))
# 
# DW <- read.csv(here("Data", "Physiology", "Dry_Weight_Oki2025.csv"))
# 
# DW <- DW %>% left_join(phys_meta, by = "frag_ID")
# 
# #calculate dry weight
# DW <- DW %>%
#   mutate(dw_g_mL = dry_pan - pan_weight,
#          dw_total_g = dw_g_mL*blastate_vol_mL,
#          dw_g_cm2 = dw_total_g/SA_cm2,
#          dw_mg_cm2 = dw_g_cm2*1000)
# 
# #average across the reps
# avg_DW <- DW %>%
#   group_by(frag_ID) %>%
#   summarise(dw_g_mL = mean(dw_g_mL),
#             dw_total_g = mean(dw_total_g),
#             dw_g_cm2 = mean(dw_g_cm2),
#             dw_mg_cm2 = mean(dw_mg_cm2)) %>%
#   left_join(phys_meta, by = "frag_ID")
# 
# write.csv(avg_DW, here("Data", "Physiology", "Average_Dry_Weight.csv"), row.names = FALSE)

#now generate mean and se dataframe of dw
se_fun <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

dw_summary <- avg_DW %>%
  group_by(species_long) %>%
  summarise(
    n    = sum(!is.na(dw_mg_cm2)),
    mean = mean(dw_mg_cm2, na.rm = TRUE),
    se   = se_fun(dw_mg_cm2),
    .groups = "drop")


#reorder based on mean values
means <- avg_DW %>% group_by(species) %>% summarize(m = mean(dw_mg_cm2, na.rm = TRUE), .groups = "drop")

avg_DW_ordered <- avg_DW %>% left_join(dw_summary, by = "species_long") %>% mutate(species_long = fct_reorder(species_long, mean))  # ascending by mean

dw_plot <- ggplot() +
  geom_jitter(data = avg_DW_ordered, aes(x = species_long, y = dw_mg_cm2, color = species_long), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = dw_summary, aes(x = species_long, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = dw_summary, aes(x = species_long, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1, face = "italic")) +
  scale_color_manual(values = sp_cols) + 
  ylim(25,325)+
  labs(x = "Species", 
       y = expression("Dry weight" ~ (mg ~ cm^{-2})))
dw_plot

ggsave(here("Output", "Physiology", "dryweight_species_jitter.pdf"), dw_plot, h = 8, w = 6)

#summarize mean, sd, se
avg_DW_spp <- avg_DW %>% 
  group_by(species) %>% 
  summarise(avg_dw_mg_cm2 = mean(dw_mg_cm2),
            sd_dw_mg_cm2 = sd(dw_mg_cm2),
            se_dw_mg_cm2 = sd(dw_mg_cm2)/sqrt(n())
  )

#quickly see if there is a dif between species
DW.mod.spp <- lm(dw_mg_cm2~species_long, data = avg_DW)
Anova(DW.mod.spp)
summary(DW.mod.spp)

emm_pairs <- pairs(emm_obj)
emm_obj <- emmeans::emmeans(DW.mod.spp, ~ species_long)
emm_tbl <- as_tibble(summary(emm_obj, infer = TRUE)) %>% mutate(species_long = fct_reorder(species_long, emmean))

emm_plot <- ggplot(emm_tbl, aes(x = emmean, y = species_long, color = species_long)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.2) +
  theme_bw(base_size = 22) +
  theme(legend.position = "none", axis.text.y = element_text(face = "italic"))+
  scale_color_manual(values = sp_cols, ) + 
  labs(x = "Emmean", y = "Species")
emm_plot

ggsave(here("Output", "TPC", "Graphs","dw_emmeans.pdf"),
       device = "pdf", height = 8, width = 8, emm_plot)

#yes there is - Prus, Tfro, and Fcom all higher

#and sensitivity
DW.mod.stress <- lm(dw_mg_cm2~stress, data = avg_DW)
Anova(DW.mod.stress)
summary(DW.mod.stress)
#based on graph, looks like Peyd and Tfro should be switched
#we'll see what happens with the rest of the TPC data!


##### Chlorophyll a #####

#read in data
#data is corrected for blank rates and surface area and blanks are removed
chla_plate1 <- read.csv(here("Data", "ChlaSpec", "oki_chla_plate1_summary.csv"))
chla_plate2 <- read.csv(here("Data", "ChlaSpec", "oki_chla_plate2_summary.csv"))
chla_plate3 <- read.csv(here("Data", "ChlaSpec", "oki_chla_plate3_summary.csv"))
chla_plate4 <- read.csv(here("Data", "ChlaSpec", "oki_chla_plate4_summary.csv"))

# concatenate all plates
chla_all_corr <- bind_rows(chla_plate1,chla_plate2,chla_plate3,chla_plate4)
chla_all_corr <- chla_all_corr %>% rename(frag_ID = Sample.Name)

#add metadata
SA_frag_final <- read.csv(here("Data","SurfaceArea","SA_fragments_final.csv"))
chla_meta <- SA_frag_final %>%
  filter(Frag_type == "physio") %>%
  select(frag_ID, full_species, species, calc_SA)

chla_data <- left_join(chla_all_corr, chla_meta, by = "frag_ID") %>% 
  rename(SA_cm2 = calc_SA)

#save data
write_csv(chla_data, here("Data", "Physiology", "Chla_avg.csv"))

#read in updated chla data
chla_avg <- read.csv(here("Data", "Physiology", "Chla_avg.csv"))
sp_keep <- c("Acropora hyacinthus","Echinopora lamellosa", "Favites complanata","Montipora aequituberculata","Pocillopora eydouxi","Porites rus","Turbinaria frondens")

chla_avg_7sp <- chla_avg %>% filter(full_species %in% sp_keep)


#now generate mean and se dataframe of chla
se_fun <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

#summarize by species
chla_summary <- chla_avg_7sp %>%
  group_by(full_species) %>%
  summarise(
    n    = sum(!is.na(chla_ug_cm2_mean)),
    mean = mean(chla_ug_cm2_mean, na.rm = TRUE),
    se   = se_fun(chla_ug_cm2_mean),
    .groups = "drop")

#reorder based on mean values
means <- chla_avg_7sp %>% group_by(species) %>% summarize(m = mean(chla_ug_cm2_mean, na.rm = TRUE), .groups = "drop")

chla_avg_7sp_ordered <- chla_avg_7sp %>% 
  left_join(chla_summary, by = "full_species") %>% 
  mutate(full_species = fct_reorder(full_species, mean))  # ascending by mean

chla_plot <- ggplot() +
  geom_jitter(data = chla_avg_7sp_ordered, aes(x = full_species, y = chla_ug_cm2_mean, color = full_species), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = chla_summary, aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6)+
  geom_point(data = chla_summary, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22)+
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1, face = "italic")) +
  scale_color_manual(values = sp_cols) +
  labs(x = "Species", 
       y = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})))
chla_plot

ggsave(here("Output", "Physiology", "chla_species_jitter.pdf"), chla_plot, h = 8, w = 6)

#summarize mean, sd, se
chla_avg_7sp_spp <- chla_avg_7sp %>% 
  group_by(species) %>% 
  summarise(chla_avg_7sp_mg_cm2 = mean(chla_ug_cm2_mean),
            sd_chla_ug_cm2_mean = sd(chla_ug_cm2_mean),
            se_chla_ug_cm2_mean = sd(chla_ug_cm2_mean)/sqrt(n())
  )

#quickly see if there is a dif between species
chla.mod.spp <- lm(chla_ug_cm2_mean~full_species, data = chla_avg_7sp)
Anova(chla.mod.spp)
summary(chla.mod.spp)

emm_obj <- emmeans::emmeans(chla.mod.spp, ~ full_species)
emm_pairs <- pairs(emm_obj)
emm_tbl <- as_tibble(summary(emm_obj, infer = TRUE)) %>% mutate(full_species = fct_reorder(full_species, emmean))

emm_plot <- ggplot(emm_tbl, aes(x = emmean, y = full_species, color = full_species)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.2) +
  theme_bw(base_size = 22) +
  theme(legend.position = "none", axis.text.y = element_text(face = "italic"))+
  scale_color_manual(values = sp_cols, ) + 
  labs(x = "Emmean", y = "Species")
emm_plot

ggsave(here("Output", "TPC", "Graphs","chla_emmeans.pdf"),
       device = "pdf", height = 8, width = 8, emm_plot)

