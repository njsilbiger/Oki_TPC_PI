#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(here)
library(dplyr)

#############################
### READ IN DATA
#############################
SA_standard_CMO <- read_csv(here("Data","SurfaceArea", "SA_standards_CMO.csv"))
#add DMB standards here
SA_frag <- read_csv(here("Data", "SurfaceArea","SA_fragments.csv"))
#once DMB standards in - do two separate calibrations and then concatenate them down below
#SA_frag_CMO <- SA_frag %>% subset(dipper == "CMO")
#SA_frag_DMB <- SA_frag %>% subset(dipper == "DMB")

###Calculate standards wax weight
SA_standard_CMO$Wax_weight <- SA_standard_CMO$weight_wWax - SA_standard_CMO$weight

model<-lm(data = SA_standard_CMO, Wax_weight ~ SA_cm2)
eq <- function(model){
  coefs <- coef(model)
  r2 <- summary(model)$r.squared
  paste0("y = ", round(coefs[1], 2), " + ", round(coefs[2], 2), "x\n",
         "R² = ", round(r2, 3))
}
wax.plot <- SA_standard_CMO%>% 
  ggplot(aes(x = SA_cm2, y = Wax_weight)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = expression("Dowel surface area (cm"^2*")"),
       y = "Wax weight (g)") +
  annotate("text",x = 5, y = max(SA_standard_CMO$Wax_weight), label = eq(model), hjust = 0)
wax.plot

###Calculate wax weight
SA_frag$Wax_weight <- SA_frag$weight_wWax - SA_frag$weight

### Extract coefficients from reg line formula for SA
a <- coef(model)[1]  # intercept
b <- coef(model)[2]  # slope

SA_frag$calc_SA <- (SA_frag$Wax_weight - a) / b
write.csv(SA_frag, here("Data","SurfaceArea","SA_fragments_final.csv"), row.names = FALSE)

#MP 9/26/2025 
#need to update with DMB curve for TPC #2 data
#all of the SA is set for PI and physio data
#shout out to the wax dipper extraordinaires!!

##############################################
### PUT SA DATA INTO OTHER METADATA SHEETS
##############################################

#import SA final dataframe
SA_frag_final <- read.csv(here("Data","SurfaceArea","SA_fragments_final.csv"))

#####TPC data####
#change to regular dataframe not est for TPC only later
SA_frag_final_est <- read.csv(here("Data","SurfaceArea","SA_fragments_final_est.csv"))

frag_data_TPC <- read_csv(here("Data", "RespoFiles", "TPC", "Fragment_Measurements_TPC.csv"))
SA_frag_final_TPC <- SA_frag_final_est %>% filter(Frag_type == "TPC") %>% select(frag_ID, calc_SA)

#join with SA and remove old SA dummy variable rename for it to work easily in scripts
frag_data_TPC <- frag_data_TPC %>% 
  left_join(SA_frag_final_TPC, by = "frag_ID") %>% 
  select(-SA_cm2) %>% 
  rename(SA_cm2 = calc_SA)

write.csv(frag_data_TPC, here("Data", "RespoFiles", "TPC", "Fragment_Measurements_TPC.csv"), row.names = FALSE)

#####PI data#####
frag_data_PI <- read_csv(here("Data", "RespoFiles", "PI", "Fragment_Measurements_PI.csv"))
SA_frag_final_PI <- SA_frag_final %>% filter(Frag_type == "PI") %>% select(frag_ID, calc_SA)

#join with SA and remove old SA dummy variable and rename for it to work easily in scripts
frag_data_PI <- frag_data_PI %>% 
  left_join(SA_frag_final_PI, by = "frag_ID") %>% 
  select(-SA_cm2) %>% 
  rename(SA_cm2 = calc_SA)

write.csv(frag_data_PI, here("Data", "RespoFiles", "PI", "Fragment_Measurements_PI.csv"), row.names = FALSE)

###Physiology data####
frag_data_phys <- read_csv(here("Data", "Physiology", "Physio_meta.csv"))
SA_frag_final_phys <- SA_frag_final %>% filter(Frag_type == "physio") %>% select(frag_ID, calc_SA)

#join with SA and remove old SA dummy variable and rename for it to work easily in scripts
frag_data_phys <- frag_data_phys %>% 
  left_join(SA_frag_final_phys, by = "frag_ID") %>% 
  #select(-SA_cm2) %>% 
  rename(SA_cm2 = calc_SA)

write.csv(frag_data_phys, here("Data", "Physiology", "Physio_meta.csv"), row.names = FALSE)
