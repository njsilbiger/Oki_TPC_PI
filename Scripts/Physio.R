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

#assign colors
# pnw_palette("Bay",7,type="continuous")
# '#24492e', '#015b58', '#2c6184', '#59629b', '#89689d', '#ba7999', '#e69b99'
# '#00496f', '#0f85a0', '#edd746', '#ed8b00', '#dd4124'
#pal <- c('#00496f','#7bbcd5','#45681e','#edd746', '#ed8b00','#dd4124','#ba7999')
sp_cols <- c(
  "Acropora hyacinthus" = '#ba7999',
  "Echinopora lamellosa" = '#dd4124',
  "Favites complanata"   = '#ed8b00',
  "Montipora aequituberculata" = '#edd746',
  "Pocillopora eydouxi" = '#45681e',
  "Porites rus" = '#7bbcd5',
  "Turbinaria frondens" = '#00496f'
)

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
