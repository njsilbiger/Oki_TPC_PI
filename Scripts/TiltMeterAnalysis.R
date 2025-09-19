## Load Libraries ###

library(tidyverse)
library(here)
library(climaemet)

#### Load Data ###

Tilt1Data <- read_csv(here("TiltMeterData", "Tilt1_CurrentData.csv"))
view(Tilt1Data)

### Using Lubridate to set DateTime Column"

tilt_fun<-function(Tilt1Data = Tilt1Data, plotname){

  Tilt1Data$DateTime <- ymd_hms(Tilt1Data$DateTime)
  Tilt1Data$Speed_m <- Tilt1Data$Speed * 0.01
  Tilt1Data$VelocityN_m <- Tilt1Data$Velocity_N * 0.01
  Tilt1Data$VelocityE_m <- Tilt1Data$Velocity_E * 0.01
  
  Tilt1Hourly <- Tilt1Data %>% #Created a dataframe that has the mean for everyhour
    group_by(DateTime = floor_date(DateTime, "min")) %>% #New datetime that is rounded to the start of each hour
    summarise( #Averaging the remaining columns by hour
      Speed = mean(Speed_m, na.rm = TRUE), 
      Heading = mean(Heading, na.rm = TRUE), 
      Velocity_N = mean(VelocityN_m, na.rm = TRUE), 
      Velocity_E = mean(VelocityE_m, na.rm = TRUE), 
      .groups = "drop") #Give you a clean data frame and drops the groups
  
  
  ### Making a windrose plot of the current data
  
  ggwindrose(
    speed = Tilt1Data$Speed_m, 
    direction = Tilt1Data$Heading,
    legend_title =  "Current speed (m/s)", 
    speed_cuts = seq(0, 0.15,.025)) +
    
    
    scale_fill_manual (values = c("blue", "cyan", "green", "yellow", "orange", "red", "darkred"), 
                       drop =FALSE)
  
  ggsave(here("Output","TiltMeterData",plotname))
  
}


tilt_fun(Tilt1Data, plotname = "tilt1.png")

Tilt2Data<-read_csv(here("TiltMeterData","Tilt2_CurrentData.csv"))
tilt_fun(Tilt2Data, plotname = "tilt2.png")

Tilt3Data <- read_csv(hear("TiltMeterData", "Tilt3_CurrentData.csv"))
tilt_fun(Tilt3Data, plotname = "tilt3.png")

Tilt4Data <- read_csv(here("TiltMeterData", "Tilt4_CurrentData.csv"))
tilt_fun(Tilt4Data, plotname = "tilt4.png")

Tilt5Data <- read_csv(here("TiltMeterData", "Tilt5_CurrentData.csv"))
tilt_fun(Tilt5Data, plotname = "tilt5.png")

Tilt6Data <- read_csv(here("TiltMeterData", "Tilt6_CurrentData.csv"))
tilt_fun(Tilt6Data, plotname = "tilt6.png")

Tilt8Data <- read_csv(here("TiltMeterData", "Tilt8_CurrentData.csv"))
tilt_fun(Tilt8Data, plotname = "tilt8.png")


table(cut(Tilt1Data$Speed_m, seq(0, 0.15, 0.025), include.lowest = TRUE)) ##Showing me the different bins and how many data points I have in each one

  ### Trying the code from Shaun to R ###
# Convert to m/s
Tilt1Data <- Tilt1Data %>%
  mutate(
    Speed = Speed / 100,   # cm/s → m/s
    Heading = Heading %% 360  # wrap into [0,360)
  )

## --- Define Sectors & Speed Bins ---
numSectors <- 12
edgesHeading <- seq(0, 360, length.out = numSectors + 1)

speedBins_cm <- seq(0, 15, by = 2.5)
speedBins <- speedBins_cm / 100  # m/s

# Custom colors (same order as bins)
speedColors <- c("#0000FF",   # blue
                 "#0080FF",   # lighter blue
                 "#00FFFF",   # cyan
                 "#80FF00",   # green
                 "#FF8000",   # orange
                 "#FF0000")   # red

## --- Bin the data ---
Tilt1Binned <- Tilt1Data %>%
  mutate(
    dir_bin = cut(Heading,
                  breaks = edgesHeading,
                  include.lowest = TRUE,
                  labels = FALSE),
    spd_bin = cut(Speed,
                  breaks = speedBins,
                  include.lowest = TRUE,
                  labels = FALSE)
  ) %>%
  filter(!is.na(dir_bin), !is.na(spd_bin)) %>%
  count(dir_bin, spd_bin)

## --- Convert to percent ---
Tilt1Binned <- Tilt1Binned %>%
  mutate(percent = n / sum(n) * 100)

## --- Prepare for plotting ---
# Midpoint of each direction bin
Tilt1Binned <- Tilt1Binned %>%
  mutate(
    dir_center = (edgesHeading[dir_bin] + edgesHeading[dir_bin + 1]) / 2,
    spd_bin = as.integer(spd_bin)
  )

## --- Plot as polar bar chart ---
ggplot(Tilt1Binned, aes(x = dir_center, y = percent, fill = factor(spd_bin))) +
  geom_col(width = 360/numSectors, color = "black", linewidth = 0.3) +
  scale_x_continuous(
    limits = c(0, 360),
    breaks = seq(0, 315, 45),
    labels = c("E","NE","N","NW","W","SW","S","SE") # adjust orientation if needed
  ) +
  scale_fill_manual(
    values = speedColors,
    labels = paste0(head(speedBins, -1), "–", tail(speedBins, -1), " m/s"),
    name = "Speed (m/s)"
  ) +
  coord_polar(start = pi/2, direction = -1) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(title ="Tilt1 Current Flow Rose (1-min data)")
