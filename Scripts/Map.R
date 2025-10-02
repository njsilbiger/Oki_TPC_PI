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

hires <- st_read(here("Data", "Okinawa_Map", "N03-19_47_190101.shp"), quiet = TRUE)  
pref_outline <- hires %>%
  st_union() %>%
  st_as_sf()

oki_outline <- ggplot() +
  geom_sf(data = pref_outline, fill = "honeydew4", color = "black", linewidth = 0.2) +
  coord_sf(xlim = c(127.5, 128.5), ylim = c(26, 27), expand = FALSE) +
  labs(title = "Okinawa Island, Japan") +
  theme_minimal()
oki_outline

ggsave(here("Output", "Okinawa_map", "oki_outline.pdf"), oki_outline, h = 8, w = 8)

