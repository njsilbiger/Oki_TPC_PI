#############################
### LOAD LIBRARIES
#############################
library(tidyverse)
library(here)


#############################
### READ IN DATA
#############################
SA_standard <- read_csv(here("Data", "SA_standards.csv"))
SA_frag <- read_csv(here("Data", "SA_fragments.csv"))

###Calculate standards wax weight
SA_standard$Wax_weight <- SA_standard$weight_wWax - SA_standard$weight

model<-lm(data = SA_standard, waxweight_g ~ SA_cm2)
eq <- function(model){
  coefs <- coef(model)
  r2 <- summary(model)$r.squared
  paste0("y = ", round(coefs[1], 2), " + ", round(coefs[2], 2), "x\n",
         "R² = ", round(r2, 3))
}
wax.plot <- SA_standard%>% 
  ggplot(aes(x = SA_cm2, y = waxweight_g)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = expression("Dowel surface area (cm"^2*")"),
       y = "Wax weight (g)") +
  annotate("text",x = 5, y = max(waxSA$waxweight_g), label = eq(model), hjust = 0)
wax.plot

###Calculate wax weight
SA_frag$Wax_weight <- SA_frag$weight_wWax - SA_frag$weight

### Extract coefficients from reg line formula for SA
a <- coef(model)[1]  # intercept
b <- coef(model)[2]  # slope

SA_frag$calc_SA <- (SA_frag$Wax_weight - a) / b
write.csv(SA_frag, here("Data","SA_fragments_final.csv"))



