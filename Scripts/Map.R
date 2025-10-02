######Making a map of Okinawa ####
#Maya Powell
#October 2nd, 2025

#####Load packages####
library(here)
library(tidyverse)
library(ggplot2)
library(sf)

#map using shape file from this link:
#https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v2_3.html

#data saved in Data/Okinawa_Map folder for hi-res data

hires <- st_read(here("Data", "Okinawa_Map", "N03-19_47_190101.shp"), quiet = TRUE) #hires map of okinawa prefectures

pref_outline <- hires %>% 
  st_union() %>% #just take biggest shape for each island so cuts out prefecture smaller outlines
  st_as_sf() #make it a shapefile again

oki_outline <- ggplot() + #generate map of just outline
  geom_sf(data = pref_outline, fill = "honeydew4", color = "black", linewidth = 0.2) +
  coord_sf(xlim = c(127.5, 128.5), ylim = c(26, 27), expand = FALSE) +
  labs(title = "Okinawa Island, Japan") +
  theme_minimal()
oki_outline

ggsave(here("Output", "Okinawa_map", "oki_outline.pdf"), oki_outline, h = 8, w = 8)

crs_target <- st_crs(pref_outline) #assign coordinate system to be the same as the outline
oist <- st_sf(name = "OIST", geometry = st_sfc(st_point(c(127.83015620094221, 26.465355941024466)), crs = crs_target))
oist_mss <- st_sf(name = "OIST MSS", geometry = st_sfc(st_point(c(127.87022582794472, 26.510131894538446)), crs = crs_target))
afuso <- st_sf(name = "Afuso Reef", geometry = st_sfc(st_point(c(127.88984, 26.51454)), crs = crs_target))

oki_outline_labels <- ggplot() +
  geom_sf(data = pref_outline, fill = "honeydew4", color = "black", linewidth = 0.2) +
  labs(title = "Okinawa Island, Japan") +
  theme_bw() +
  geom_sf(data = oist, shape = 21, fill = "darkred", size = 4, stroke = 0.5) + #add labels on map
  geom_sf_text(data = oist, aes(label = name), nudge_x = -0.05, fontface = "bold", size = 5) +
  geom_sf(data = oist_mss, shape = 21, fill = "darkred", size = 4, stroke = 0.5) +
  geom_sf_text(data = oist_mss, aes(label = name), nudge_x = -0.1, fontface = "bold", size = 5) +
  geom_sf(data = afuso, shape = 21, fill = "cornflowerblue", size = 4, stroke = 0.5) +
  geom_sf_text(data = afuso, aes(label = name), nudge_y = 0.03, fontface = "bold", size = 5) +
  coord_sf(xlim = c(127.5, 128.5), ylim = c(26, 27), expand = FALSE) + #make sure to set boundary for map after adding labels because coord system will make map big if not
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) 
oki_outline_labels

ggsave(here("Output", "Okinawa_map", "oki_outline_labels.pdf"), oki_outline_labels, h = 8, w = 8)


