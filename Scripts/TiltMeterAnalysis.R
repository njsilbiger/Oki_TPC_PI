## Load Libraries ###

library(tidyverse)
library(here)
library(climaemet)

#### Load Data ###

Tilt1Data <- read_csv(here("TiltMeterData", "Tilt1_CurrentData.csv"))
view(Tilt1Data)



### Tilt Meter 7 Data Condensing ###
##Tilt 7 was set to take a reading at every second we will use the lubridate function to condense it to one minute intervals

Tilt7_CurrentData <- read_csv(here("TiltMeterData", "Tilt7_CurrentData.csv")) %>%

  group_by(DateTime = floor_date(DateTime, "min")) %>% #New datetime that is rounded to the start of each hour
  summarise( #Averaging the remaining columns by hour
    Speed = mean(Speed, na.rm = TRUE), 
    Heading = mean(Heading, na.rm = TRUE), 
    Velocity_N = mean(Velocity_N, na.rm = TRUE), 
    Velocity_E = mean(Velocity_E, na.rm = TRUE), 
    .groups = "drop") #Give you a clean data frame and drops the groups


write.csv(Tilt7_CurrentData,  file ="TiltMeterData/Tilt7_CurrentData.csv") #Overwriting the CSV file with the new condensed one

## Tilt 8 was taking data every 20 seconds 
Tilt8_CurrentData <- read_csv(here("TiltMeterData", "Tilt8_CurrentData.csv")) %>%
  
  group_by(DateTime = floor_date(DateTime, "min")) %>% #New datetime that is rounded to the start of each hour
  summarise( #Averaging the remaining columns by hour
    Speed = mean(Speed, na.rm = TRUE), 
    Heading = mean(Heading, na.rm = TRUE), 
    Velocity_N = mean(Velocity_N, na.rm = TRUE), 
    Velocity_E = mean(Velocity_E, na.rm = TRUE), 
    .groups = "drop") #Give you a clean data frame and drops the groups


write.csv(Tilt8_CurrentData,  file ="TiltMeterData/Tilt8_CurrentData.csv")



### Using Lubridate to set DateTime Column"

tilt_fun<-function(Tilt4Data = Tilt4Data, plotname){ #Creating a function to run all tilt meter data easier

  declination <- -5.85 # West declination in degrees to adjust magnetic north (which the sensors take) to true north
  
  Tilt4Data$DateTime <- ymd_hms(Tilt4Data$DateTime) #Setting the DateTime Column
  Tilt4Data$Speed_m <- Tilt4Data$Speed * 0.01 #Converting cm/s to m/s
  Tilt4Data$VelocityN_m <- Tilt4Data$Velocity_N * 0.01
  Tilt4Data$VelocityE_m <- Tilt4Data$Velocity_E * 0.01
  Tilt4Data$Heading_True <- (Tilt4Data$Heading + declination) %% 360 # Create a new column for True Heading
  

  
  ### Making a windrose plot of the current data


  
  ggwindrose(
    speed = Tilt4Data$Speed_m, 
    direction = Tilt4Data$Heading_True,
    n_directions = 8,
    speed_cuts = seq(0, 0.50,.05), Inf) + 
    #stack_revers = TRUE) + #If we want the higher speeds on the outside of the plot
    
    

    scale_fill_manual(
      values = c( "skyblue","blue", "cyan" ,"green","yellow", "tan", "darkorange1","hotpink","magenta", "red", "red4"),
      drop = FALSE
    ) +
    
    labs(fill = "Current Speed (m/s)") +
      theme_minimal() +
    theme(axis.title = element_blank(), 
          text = element_text(face = "bold"))
          #panel.background = element_rect(fill = "transparent", color = NA),  # Plot panel
         # plot.background = element_rect(fill = "transparent", color = NA),   # Entire plot
         # legend.background = element_rect(fill = "transparent", color = NA), # Legend
          #legend.box.background = element_rect(fill = "transparent", color = NA))
      
  
  
  
  ggsave(here("Output","TiltMeterData",plotname))
  
}

Tilt1Data<-read_csv(here("TiltMeterData","Tilt1_CurrentData.csv"))
tilt_fun(Tilt1Data, plotname = "tilt1.png")


Tilt2Data<-read_csv(here("TiltMeterData","Tilt2_CurrentData.csv"))
tilt_fun(Tilt2Data, plotname = "tilt2.png")

Tilt3Data <- read_csv(here("TiltMeterData", "Tilt3_CurrentData.csv"))


tilt_fun(Tilt3Data, plotname = "tilt3.png")

Tilt4Data <- read_csv(here("TiltMeterData", "Tilt4_CurrentData.csv"))
tilt_fun(Tilt4Data, plotname = "tilt4.png")

Tilt5Data <- read_csv(here("TiltMeterData", "Tilt5_CurrentData.csv"))

tilt_fun(Tilt5Data, plotname = "tilt5.png")

Tilt6Data <- read_csv(here("TiltMeterData", "Tilt6_CurrentData.csv"))
tilt_fun(Tilt6Data, plotname = "tilt6.png")

Tilt7Data <- read_csv(here("TiltMeterData", "Tilt7_CurrentData.csv"))
tilt_fun(Tilt7Data, plotname = "tilt7.png")


Tilt8Data <- read_csv(here("TiltMeterData", "Tilt8_CurrentData.csv"))
tilt_fun(Tilt8Data, plotname = "tilt8.png")

view(Tilt4Data)
table(cut(Tilt4Data$Speed *0.01 , seq(0, .5,.05), include.lowest = TRUE)) ##Showing me the different bins and how many data points I have in each one

mean(Tilt1Data$Speed) * 0.01
mean(Tilt2Data$Speed) * 0.01
mean(Tilt3Data$Speed) * 0.01
mean(Tilt4Data$Speed) * 0.01
mean(Tilt5Data$Speed) * 0.01
mean(Tilt6Data$Speed) * 0.01
mean(Tilt7Data$Speed) * 0.01
mean(Tilt8Data$Speed) * 0.01




#### The rest of the script is testing different plots ####

# Pre-bin the data by direction + speed bins
rose_data <- 
  Tilt1Data %>%
  mutate(speed_bin = cut(Speed_m, breaks = seq(0, 0.15, .025), include.lowest = TRUE)) %>%
  group_by(Heading, speed_bin) %>%
  tally() %>%
  ungroup() %>%
  mutate(percent = 100 * n / sum(n))   # normalize counts to %

# Plot with ggplot (manual, but very flexible)
ggplot(rose_data, aes(x = Heading, y = percent, fill = speed_bin)) +
  geom_bar(stat = "identity", width = 10, color = "black") +  # width controls bin size
  coord_polar(start = 0) +
  scale_y_continuous(
    limits = c(0, 60),                    # fix max radius
    breaks = seq(0, 60, 10),              # breaks at 10%
    labels = paste0(seq(0, 60, 10), "%")  # percent labels
  ) +
  scale_fill_manual(values = c("blue", "cyan", "green", "yellow", "orange", "red", "darkred")) +
  theme_minimal()


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
