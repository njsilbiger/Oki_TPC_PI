######Making a map of Okinawa ####
##Also then added temperature data
#Maya Powell
#October 2nd, 2025

#####Load packages####
library(here)
library(tidyverse)
library(ggplot2)
library(sf)
library(lubridate)
library(dplyr)

#map using shape file from this link:
#https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v2_3.html

#data saved in Data/Okinawa_Map folder for hi-res data

hires <- st_read(here("Data", "Okinawa_Map", "N03-19_47_190101.shp"), quiet = TRUE) #hires map of okinawa prefectures

pref_outline <- hires %>% 
  st_union() %>% #just take biggest shape for each island so cuts out prefecture smaller internal outlines
  st_as_sf() #make it a shapefile again

oki_outline <- ggplot() + #generate map of just outline
  geom_sf(data = pref_outline, fill = "honeydew4", color = "black", linewidth = 0.2) +
  coord_sf(xlim = c(127.5, 128.5), ylim = c(26, 27), expand = FALSE) +
  #labs(title = "Okinawa Island, Japan") +
  theme_minimal()
oki_outline

ggsave(here("Output", "Okinawa_map", "oki_outline.pdf"), oki_outline, h = 8, w = 8)

crs_target <- st_crs(pref_outline) #assign coordinate system to be the same as the outline
oist <- st_sf(name = "OIST", geometry = st_sfc(st_point(c(127.83015620094221, 26.465355941024466)), crs = crs_target))
oist_mss <- st_sf(name = "OIST MSS", geometry = st_sfc(st_point(c(127.87022582794472, 26.510131894538446)), crs = crs_target))
afuso <- st_sf(name = "Afuso Reef", geometry = st_sfc(st_point(c(127.88984, 26.51454)), crs = crs_target))

oki_outline_labels <- ggplot() +
  geom_sf(data = pref_outline, fill = "honeydew4", color = "black", linewidth = 0.2) +
  #labs(title = "Okinawa Island, Japan") +
  theme_bw() +
  geom_sf(data = oist, shape = 21, fill = "firebrick3", size = 4, stroke = 0.5) + #add labels on map
  geom_sf_text(data = oist, aes(label = name), nudge_x = -0.05, fontface = "bold", size = 5) +
  geom_sf(data = oist_mss, shape = 21, fill = "firebrick3", size = 4, stroke = 0.5) +
  geom_sf_text(data = oist_mss, aes(label = name), nudge_x = -0.09, fontface = "bold", size = 5) +
  geom_sf(data = afuso, shape = 21, fill = "cornflowerblue", size = 4, stroke = 0.5) +
  geom_sf_text(data = afuso, aes(label = name), nudge_y = 0.03, fontface = "bold", size = 5) +
  coord_sf(xlim = c(127.5, 128.5), ylim = c(26, 27), expand = FALSE) + #make sure to set boundary for map after adding labels because coord system will make map big if not
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) 
oki_outline_labels

ggsave(here("Output", "Okinawa_map", "oki_outline_labels.pdf"), oki_outline_labels, h = 8, w = 8)

##### Temperature data

#read in temp data from Tilt 2 (at coral collection site)
temp <- read_csv(here("Data", "TiltMeterData", "Tilt2_Temperature.csv"))

#add column to denote temperature points above Topt
#Fcom highest mean Topt for NP = 29.87
#Maeq lowest mean Topt for NP = 28.294
#temp <- temp %>% mutate(hot = Temperature > 28.294) #works if doing geom_point but I think geom_line looks better
#note which values are above Topt for np
temp_np <- temp %>% arrange(DateTime) %>% 
  mutate(hot = Temperature >= 28.294,
         temp_hot  = ifelse(hot,  Temperature, NA_real_),
         temp_cool = ifelse(!hot, Temperature, NA_real_))

#plot np topt temperature graph
temp_topt_np <- ggplot(temp_np, aes(x = DateTime)) +
  geom_line(aes(y = temp_cool), linewidth = 0.6, color = "black") +
  geom_line(aes(y = temp_hot),  linewidth = 0.8, color = "firebrick3") +
  geom_hline(yintercept = 28.294, linetype = "dashed") +
  geom_hline(yintercept = 29.87, linetype = "dashed") +
  annotate("text", x = ymd("2025-08-08"), y = 28.294,
           label = "Minimum~mean~T['opt']", parse = TRUE, vjust = -0.4, size = 7, color = "grey20") +
  annotate("text", x = ymd("2025-08-08"), y = 29.87,
           label = "Maximum~mean~T['opt']", parse = TRUE, vjust = -0.4, size = 7, color = "grey20") +
  #geom_vline(xintercept = )
  labs(x = "Date", y = "Temperature (°C)") +
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylim(27,30.2)+
  scale_x_datetime(date_breaks = "1 day", 
                   minor_breaks = NULL,
                   date_labels = "%b %d")

temp_topt_np

ggsave(here("Output", "Okinawa_map", "temp_topt_np.pdf"), temp_topt_np, h = 8, w = 10)

#generate same dataframe but for topt values with gp
#Tfro highest mean Topt for GP = 30.572
#Mvie lowest mean Topt for GP = 29.552

temp_gp <- temp %>% arrange(DateTime) %>% 
  mutate(hot = Temperature >= 29.552,
         temp_hot  = ifelse(hot,  Temperature, NA_real_),
         temp_cool = ifelse(!hot, Temperature, NA_real_))

#plot gp topt temperature graph
temp_topt_gp <- ggplot(temp_gp, aes(x = DateTime)) +
  geom_line(aes(y = temp_cool), linewidth = 0.6, color = "black") +
  geom_line(aes(y = temp_hot),  linewidth = 0.6, color = "firebrick3") +
  geom_hline(yintercept = 29.552, linetype = "dashed") +
  geom_hline(yintercept = 30.572, linetype = "dashed") +
  annotate("text", x = ymd("2025-08-04"), y = 29.552,
           label = "Minimum~mean~T['opt']", parse = TRUE, vjust = -0.4, size = 7, color = "grey20") +
  annotate("text", x = ymd("2025-08-04"), y = 30.572,
           label = "Maximum~mean~T['opt']", parse = TRUE, vjust = -0.4, size = 7, color = "grey20") +
  #geom_vline(xintercept = )
  labs(x = "Date", y = "Temperature (°C)") +
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylim(27,30.7)+
  scale_x_datetime(date_breaks = "1 day", 
                   minor_breaks = NULL,
                   date_labels = "%b %d")

temp_topt_gp

ggsave(here("Output", "Okinawa_map", "temp_topt_gp.pdf"), temp_topt_gp, h = 8, w = 10)


#topt fake curve schematic

# Thermal performance curve schematic (fake data + annotations)
library(grid)  # for arrow()

set.seed(1)

# Thermal performance curve (bounded Beta shape) + full annotations
library(ggplot2)
library(dplyr)
library(grid)   # for unit() in arrows

# ---- Parameters ----
CTmin <- 20
CTmax <- 35
Topt  <- 29
Rmax  <- 1.5

temp <- seq(CTmin, CTmax, by = 0.1)
s    <- (temp - CTmin) / (CTmax - CTmin)             # map temperature to [0, 1]
s    <- pmin(pmax(s, 0), 1)
s0   <- (Topt - CTmin) / (CTmax - CTmin)             # desired mode location in [0,1]

# Choose skew (left-skew = longer cold tail) by setting shape2 > shape1
shape2 <- 3
shape1 <- (s0*shape2 - 2*s0 + 1) / (1 - s0)          # ensure mode = s0; requires shape1, shape2 > 1
while (shape1 <= 1) {
  shape2 <- shape2 + 1
  shape1 <- (s0*shape2 - 2*s0 + 1) / (1 - s0)
}

# Beta-shaped curve; scale to peak at Rmax
curve_beta <- dbeta(s, shape1 = shape1, shape2 = shape2)
rate_mean  <- Rmax * curve_beta / max(curve_beta)

df <- tibble(temp, rate_mean)

# ---- Peak (use discrete max for exact plotting coords) ----
i_peak <- which.max(rate_mean)
Topt_eff <- temp[i_peak]          # should be ~ Topt
Rmax_eff <- rate_mean[i_peak]     # should be ~ Rmax

#breadth
half_height <- 0.5 * Rmax_eff

left_idx  <- which(temp <= Topt_eff)
right_idx <- which(temp >= Topt_eff)

T_low  <- temp[left_idx][  which.min(abs(rate_mean[left_idx]  - half_height)) ]
T_high <- temp[right_idx][ which.min(abs(rate_mean[right_idx] - half_height)) ]
breadth <- T_high - T_low


# ---- E (activation energy) schematic: a rising-limb segment ----
E_x1 <- max(CTmin, Topt - 7)
E_x2 <- max(CTmin + 0.1, Topt - 3)
E_y1 <- approx(temp, rate_mean, xout = E_x1)$y
E_y2 <- approx(temp, rate_mean, xout = E_x2)$y

# ---- Plot ----
p <- ggplot(df, aes(x = temp, y = rate_mean)) +
  # curve + points
  geom_line(aes(y = rate_mean), linewidth = 1.2) +
  
  # Rmax: mark peak, arrow, and label
  annotate("point", x = Topt_eff, y = Rmax_eff, size = 3) +
  annotate("segment",
           x = Topt_eff, xend = Topt_eff + 2,
           y = Rmax_eff, yend = Rmax_eff,
           arrow = arrow(length = unit(6, "pt"))) +
  annotate("text",
           x = Topt_eff + 2.3, y = Rmax_eff, hjust = 0, vjust = 0,
           label = "R['max']",
           parse = TRUE, size = 8) +
  
  # Topt: vertical dashed line + parsed label near x-axis
  geom_vline(xintercept = Topt, linetype = 2) +
  annotate("text",
           x = Topt, y = 0, vjust = 0.5, hjust = 2,
           label = "T['opt']",
           parse = TRUE, size = 8) +
  
  # Breadth (FWHM): double-headed arrow at half height + ticks + label
  annotate("segment",
           x = T_low, xend = T_high,
           y = half_height, yend = half_height,
           arrow = arrow(ends = "both", length = unit(6, "pt"))) +
  annotate("text",
           x = (T_low + T_high)/2, y = half_height + 0.08*Rmax_eff,
           label = "breadth",
           parse = TRUE, hjust = 1, vjust = 1, size = 8) +
  
  # E (activation energy) schematic on rising limb + label
  annotate("segment",
           x = E_x1, xend = E_x2,
           y = E_y1, yend = E_y2,
           linewidth = 1.1, arrow = arrow(length = unit(6, "pt"))) +
  annotate("text",
           x = E_x1 - 0.2, y = (E_y1 + E_y2)/2,
           label = "italic('E')",
           parse = TRUE, hjust = -0.5, vjust = 0.5, size = 8) +
  annotate("text", x = CTmin, y = 0.2, vjust = 1.6, hjust = -0.05,
           label = "CT['min']", parse = TRUE, size = 8, color = "slategray4") +
  annotate("text", x = CTmax, y = 0.2, vjust = 1.6, hjust = -0.01,
           label = "CT['max']", parse = TRUE, size = 8, color = "slategray4") +
  labs(x = expression("Temperature ("*degree*C*")"),
       y = expression("Physiological Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1}))) +
  coord_cartesian(ylim = c(0, Rmax*1.15)) +
  theme_classic(base_size = 22)+
  xlim(20,36)

p

ggsave(here("Output","Okinawa_Map","tpc_schematic.pdf"), p, h = 8, w = 10)
