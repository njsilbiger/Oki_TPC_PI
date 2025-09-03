####################################
### Coral physiology Okinawa 2025
####################################
#Maya Powell
#Last edited August 29th 2025

#load libraries
library(tidyverse)
library(here)
library(dplyr)
library(car)

##################
## Dry weight
###################

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

##boxplot by species
#and add facet by morphology or stress (sensitive/tolerant)
dw_box <- ggplot(avg_DW, aes(x = species, y=dw_mg_cm2, fill = species))+
  geom_boxplot()+
  geom_jitter(alpha=0.8, width=0.2)+
  theme_classic(base_size = 22)+
  labs(x = "Coral Species", y = "Dry Weight (mg/cm2)")+
  theme(legend.position = "none") +
  #facet_wrap(.~stress, scales = "free")+
  ylim(25,325)
  #scale_fill_manual(values = c("tan","brown"))
dw_box

ggsave(here("Output", "Physiology", "dryweight_box.pdf"), dw_box, h = 5, w = 8)

#summarize mean, sd, se
avg_DW_spp <- avg_DW %>% 
  group_by(species) %>% 
  summarise(avg_dw_mg_cm2 = mean(dw_mg_cm2),
            sd_dw_mg_cm2 = sd(dw_mg_cm2),
            se_dw_mg_cm2 = sd(dw_mg_cm2)/sqrt(n())
  )

#quickly see if there is a dif between species
DW.mod.spp <- lm(dw_mg_cm2~species, data = avg_DW)
Anova(DW.mod.spp)
summary(DW.mod.spp)

#yes there is - Prus, Tfro, and Fcom all higher

#and sensitivity
DW.mod.stress <- lm(dw_mg_cm2~stress, data = avg_DW)
Anova(DW.mod.stress)
summary(DW.mod.stress)
#based on graph, looks like Peyd and Tfro should be switched
#we'll see what happens with the rest of the TPC data!
