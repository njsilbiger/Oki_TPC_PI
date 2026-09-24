######Making a map of Okinawa ####
##Also then added temperature data
#Maya Powell
#September 2026

#####Load packages####
library(here)
library(tidyverse)
library(ggplot2)
library(sf)
library(lubridate)
library(rmapshaper)
library(patchwork)
library(png)
library(grid)
library(ggpubr)
library(ggspatial)

#map using shape file from this link:
#https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v2_3.html

#data saved in Data/Okinawa_Map folder for hi-res data
#only need to make shape files of outlines once - takes a long time
#THESE FILES ARE NOT ON THIS GITHUB BUT YOU CAN FIND THEM AT THE LINK ABOVE IF YOU NEED THE FULL SHAPE FILES
# hires <- st_read(here("Data", "Okinawa_Map", "N03-19_47_190101.shp"), quiet = TRUE) #hires map of okinawa prefectures
# jp <- st_read(here("Data", "Okinawa_Map", "N03-19_190101.shp"), quiet = TRUE) #hires map of japan
# 
# pref_outline <- hires %>%
#   st_make_valid() %>%
#   st_union() %>% #just take biggest shape for each island so cuts out prefecture smaller internal outlines
#   st_as_sf() %>% #make it a shapefile again
#   ms_simplify(keep = 0.05, keep_shapes = TRUE)
# 
# jp_outline <- jp %>%
#   st_make_valid() %>%
#   st_union() %>% #just take biggest shape for each island so cuts out prefecture smaller internal outlines
#   st_as_sf() %>% #make it a shapefile again
#   ms_simplify(keep = 0.01, keep_shapes = TRUE)
# 
# st_write(pref_outline, here("Data/Okinawa_Map/Okinawa_outline.shp"), delete_layer = TRUE)
# st_write(jp_outline, here("Data/Okinawa_Map/Japan_outline.shp"), delete_layer = TRUE)

#read shape files in of just outlines
pref_outline <- st_read(here("Data", "Okinawa_Map", "Okinawa_outline.shp"), quiet = TRUE) #hires map of okinawa prefectures
jp <- st_read(here("Data", "Okinawa_Map", "Japan_outline.shp"), quiet = TRUE) #hires map of japan

oki_outline <- ggplot() + #generate map of just outline
  geom_sf(data = pref_outline, fill = "honeydew4", color = "black", linewidth = 0.2) +
  coord_sf(xlim = c(127.5, 128.5), ylim = c(26, 27), expand = FALSE) +
  ggspatial::annotation_scale(location = 'tl') + 
  #labs(title = "Okinawa Island, Japan") +
  theme_classic(base_size = 22)
oki_outline

#ggsave(here("Output", "Okinawa_map", "oki_outline.pdf"), oki_outline, h = 8, w = 8)

crs_target <- st_crs(pref_outline) #assign coordinate system to be the same as the outline
crs_target <- st_crs(jp) #assign coordinate system to be the same as the outline
oist <- st_sf(name = "OIST", geometry = st_sfc(st_point(c(127.83015620094221, 26.465355941024466)), crs = crs_target))
oist_mss <- st_sf(name = "OIST MSS", geometry = st_sfc(st_point(c(127.87022582794472, 26.510131894538446)), crs = crs_target))
afuso <- st_sf(name = "Afuso Reef", geometry = st_sfc(st_point(c(127.88984, 26.51454)), crs = crs_target))
japan <- st_sf(name = "Japan", geometry = st_sfc(st_point(c(127,42)), crs = crs_target))
oki <- st_sf(name = "Okinawa", geometry = st_sfc(st_point(c(127.7,26.9)), crs = crs_target))

oki_outline_labels <- ggplot() +
  geom_sf(data = pref_outline, fill = "honeydew4", color = "black", linewidth = 0.2) +
  #labs(title = "Okinawa Island, Japan") +
  theme_classic(base_size = 15) +
  #geom_sf(data = oist, shape = 21, fill = "firebrick3", size = 4, stroke = 0.5) + #add labels on map
  #geom_sf_text(data = oist, aes(label = name), nudge_x = -0.05, fontface = "bold", size = 5) +
  #geom_sf(data = oist_mss, shape = 21, fill = "firebrick3", size = 8, stroke = 0.5) +
  #geom_sf_text(data = oist_mss, aes(label = name), nudge_x = -0.15, fontface = "bold", size = 8) +
  geom_sf(data = afuso, shape = 21, fill = "cornflowerblue", size = 4, stroke = 0.5) +
  geom_sf_text(data = afuso, aes(label = name), nudge_y = 0.03, nudge_x = -0.07, fontface = "bold", size = 5) +
  geom_sf_text(data = oki, aes(label = name), fontface = "bold", size = 10) +
  ggspatial::annotation_scale(location = 'tl') + 
  coord_sf(xlim = c(127.5, 128.5), ylim = c(26, 27), expand = FALSE) + #make sure to set boundary for map after adding labels because coord system will make map big if not
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) 
oki_outline_labels

#ggsave(here("Output", "Okinawa_map", "oki_outline_afuso.pdf"), oki_outline_labels, h = 8, w = 8)

jp_outline <- ggplot() + #generate map of just outline
  geom_sf(data = jp, fill = "honeydew4", color = "black", linewidth = 0.2) +
  #labs(title = "Okinawa Island, Japan") +
  theme_classic(base_size = 10) +
  geom_sf(data = afuso, shape = 0, size = 9, stroke = 1) +
  geom_sf_text(data = japan, aes(label = name), fontface = "bold", size = 5) +
  coord_sf(xlim = c(122, 150), ylim = c(22, 48), expand = F) +
  ggspatial::annotation_scale(location = 'tl') + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) 
jp_outline

#ggsave(here("Output", "Okinawa_map", "jp_outline.pdf"), jp_outline, h = 8, w = 8)

ok_jp_inset <- oki_outline_labels + inset_element(jp_outline, 0.56, 0.01, 1, 0.46)
ok_jp_inset

ggsave(here("Output", "Okinawa_map", "oki_jp_map_inset.pdf"), ok_jp_inset, h = 8, w = 8)

##### Temperature data

#read in temp data from Tilt 2 (at coral collection site)
temp <- read_csv(here("Data", "TiltMeterData", "Tilt2_Temperature.csv"))
aqua_temp <- read_csv(here("Data/TiltMeterData/2025_NOAA_Satellite_Temp.csv"))
aqua_temp <- aqua_temp %>% 
  filter(timestamp > as_date("2025-01-01"),
         timestamp < as_date("2025-12-31"))

#temp plot
temp_plot <- ggplot(temp, aes(y = Temperature, x = DateTime)) +
  geom_line(color = "cornflowerblue")+
  #geom_hline(yintercept = 28.294, linetype = "dashed") +
  #geom_hline(yintercept = 29.87, linetype = "dashed") +
  labs(x = "Date", y = "Temperature (°C)") +
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylim(27,30.2)+
  scale_x_datetime(date_breaks = "2 days", 
                   minor_breaks = NULL,
                   date_labels = "%b %d")

temp_plot
ggsave(here("Output", "Okinawa_map", "temp_plot.pdf"), temp_plot, h = 8, w = 10)

max(temp$Temperature) #29.99
mean(temp$Temperature) #29.16
min(temp$Temperature) #27.06

#yearly satellite temp
noaa_temp_plot <- ggplot(aqua_temp) +
  geom_line(aes(y = satellite_temperature_noaa, x = timestamp))+
  #geom_line(aes(y = dhw_noaa, x = timestamp))+
  #geom_hline(yintercept = 28.294, linetype = "dashed") +
  #geom_hline(yintercept = 29.87, linetype = "dashed") +
  labs(x = "Date", y = "Temperature (°C)") +
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #ylim(27,30.2)+
  scale_x_datetime(date_breaks = "1 month", 
                   minor_breaks = NULL,
                   date_labels = "%b")
noaa_temp_plot
ggsave(here("Output", "Okinawa_map", "noaa_2025_temp_plot.pdf"), noaa_temp_plot, h = 4, w = 12)

max(aqua_temp$satellite_temperature_noaa) #30.28
mean(aqua_temp$satellite_temperature_noaa) #25.62
min(aqua_temp$satellite_temperature_noaa) #20.20


##Create figure 1 for paper

map_temp <- ggarrange(oki_outline_labels,temp_plot,ncol = 2, nrow = 1, labels = c("A","B"), font.label = list(size = 25, color = "black"))

map_temp_temp <- ggarrange(map_temp, noaa_temp_plot, ncol = 1, nrow = 2, labels = c("A", "C"), font.label = list(size = 25, color = "black"))
map_temp_temp

ggsave(here("Output", "Okinawa_map", "map_temp_noaatemp.pdf"), map_temp_temp, h = 8, w = 10)

##add species photos to first plot
all_sp_pics <- readPNG(here("Output", "Physiology","Okinawa2025_AllSpecies.png"))
all_sp_pics <- rasterGrob(all_sp_pics, interpolate = TRUE)

map_sp <- ok_jp_inset + all_sp_pics + plot_layout(ncol = 2, widths = c(2,1))
map_sp_temp <- ggarrange(map_sp, noaa_temp_plot, ncol = 1, nrow = 2, labels = c("A", "C"), font.label = list(size = 25, color = "black"))
map_sp_temp

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
           parse = TRUE, size = 8, color = "red") +
  
  # Topt: vertical dashed line + parsed label near x-axis
  geom_vline(xintercept = Topt, linetype = 2) +
  annotate("text",
           x = Topt, y = 0, vjust = 0.5, hjust = 1.3,
           label = "T['opt']",
           parse = TRUE, size = 8, color = "red") +
  
  # Breadth (FWHM): double-headed arrow at half height + ticks + label
  annotate("segment",
           x = T_low, xend = T_high,
           y = half_height, yend = half_height,
           arrow = arrow(ends = "both", length = unit(6, "pt"))) +
  annotate("text",
           x = 29, y = 0.85,
           label = "breadth",
           parse = TRUE, size = 8, color = "red") +
  
  # E (activation energy) schematic on rising limb + label
  annotate("segment",
           x = 23.8, xend = 27.3,
           y = 0.5, yend = 1.4,
           linewidth = 1.1, arrow = arrow(length = unit(6, "pt"))) +
  annotate("text",
           x = 24, y = 1,
           label = "e",
           parse = TRUE, hjust = -0.5, vjust = 0.5, size = 8, color = "red") +
  # annotate("text", x = CTmin, y = 0.2, vjust = 1.6, hjust = -0.05,
  #          label = "CT['min']", parse = TRUE, size = 8, color = "slategray4") +
  # annotate("text", x = CTmax, y = 0.2, vjust = 1.6, hjust = -0.01,
  #          label = "CT['max']", parse = TRUE, size = 8, color = "slategray4") +
  labs(x = expression("Temperature ("*degree*C*")"),
       y = expression("Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1}))) +
  coord_cartesian(ylim = c(0, Rmax*1.15)) +
  theme_classic(base_size = 22)+
  xlim(20,36)

p

ggsave(here("Output","Okinawa_Map","tpc_schematic.pdf"), p, h = 8, w = 10)
saveRDS(p, here("Output/Okinawa_Map/tpc_schematic.rds"))

#simple plot
schematic <- ggplot(df, aes(x = temp, y = rate_mean)) +
  # curve + points
  geom_line(aes(y = rate_mean), linewidth = 1.2) +
  # Topt: vertical dashed line + parsed label near x-axis
  geom_vline(xintercept = Topt, linetype = 2) +
  labs(x = expression("Temperature ("*degree*C*")"),
       y = expression("Physiological Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1}))) +
  theme_classic(base_size = 22)

ggsave(here("Output","Okinawa_Map","tpc_simple_schematic.pdf"), p, h = 8, w = 10)
