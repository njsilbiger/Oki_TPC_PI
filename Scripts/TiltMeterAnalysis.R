## Load Libraries ###

library(tidyverse)
library(here)
library(climaemet)
library(rstatix)
library(cowplot)
library(ggpubr)

####Condensing Data IGNORE and continue from Load Data#####
##Only overwrote the orginal CSV because it was too big to upload to GitHub 

Tilt7_CurrentData <- read_csv(here("TiltMeterData", "Tilt7_CurrentData.csv")) %>% #Was taking recordings every second
  
  group_by(DateTime = floor_date(DateTime, "min")) %>% #New datetime that is rounded to the start of each hour
  summarise( #Averaging the remaining columns by hour
    Speed = mean(Speed, na.rm = TRUE), 
    Heading = mean(Heading, na.rm = TRUE), 
    Velocity_N = mean(Velocity_N, na.rm = TRUE), 
    Velocity_E = mean(Velocity_E, na.rm = TRUE), 
    .groups = "drop") #Give you a clean data frame and drops the groups


write_csv(Tilt7_CurrentData,  file ="TiltMeterData/Tilt7_CurrentData.csv") #Overwriting the CSV file with the new condensed one


### Tilt 4 was taking data every 20 seconds ###
Tilt4_CurrentData <- read_csv(here("TiltMeterData", "Tilt4_CurrentData.csv")) %>%
  
  group_by(DateTime = floor_date(DateTime, "min")) %>% #New datetime that is rounded to the start of each hour
  summarise( #Averaging the remaining columns by hour
    Speed = mean(Speed, na.rm = TRUE), 
    Heading = mean(Heading, na.rm = TRUE), 
    Velocity_N = mean(Velocity_N, na.rm = TRUE), 
    Velocity_E = mean(Velocity_E, na.rm = TRUE), 
    .groups = "drop") #Give you a clean data frame and drops the groups

write_csv(Tilt4_CurrentData,  file ="TiltMeterData/Tilt4_CurrentData_Condensed.csv") ##Writing a new CSV file

##Importing Tide Data from .txt file and save filtered data for appropriate dates
#https://www.gsi.go.jp/kanshi/tide_data_21_e.html

TideData_txtFile <- "Data/TiltMeterData/21_2021-2025_hour.txt" #Reading in the txt file

TideData1 <- read.table(TideData_txtFile, header = FALSE)
view(TideData)

colnames(TideData1) <- c("Number", "Date", "Hour", "Sealevel_mm", "DatumConstant_mm", "FixedPointHeight_mm") #Changing the column names
head(TideData1)
TideData1 <- TideData1[,-1] #Removing the first row which was #No from the data?
head(TideData1)

TideDataFiltered <- TideData1 %>% #Filtering out the dates that data was not recorded during 
  filter(between(Date, "2025/08/01", "2025/08/15"))
View(TideData1Filtered)

write_csv(TideData1Filtered, here("Data/TiltMeterData/TideData_2025.08.csv"))

### Tilt 8 was taking data every 20 seconds 
Tilt8_CurrentData <- read_csv(here("TiltMeterData", "Tilt8_CurrentData.csv")) %>%
  
  group_by(DateTime = floor_date(DateTime, "min")) %>% #New datetime that is rounded to the start of each hour
  summarise( #Averaging the remaining columns by hour
    Speed = mean(Speed, na.rm = TRUE), 
    Heading = mean(Heading, na.rm = TRUE), 
    Velocity_N = mean(Velocity_N, na.rm = TRUE), 
    Velocity_E = mean(Velocity_E, na.rm = TRUE), 
    .groups = "drop") #Give you a clean data frame and drops the groups

write_csv(Tilt8_CurrentData,  file ="TiltMeterData/Tilt8_CurrentData.csv")

#### Load Data ####
 #Tilt Meter Data
Tilt1Data <- read_csv(here("Data", "TiltMeterData", "Tilt1_CurrentData.csv"))
Tilt2Data <- read_csv(here("Data", "TiltMeterData", "Tilt2_CurrentData.csv"))
Tilt3Data <- read_csv(here("Data","TiltMeterData", "Tilt3_CurrentData.csv"))
Tilt4Data <- read_csv(here("Data","TiltMeterData", "Tilt4_CurrentData_Condensed.csv"))
Tilt5Data <- read_csv(here("Data","TiltMeterData", "Tilt5_CurrentData.csv"))
Tilt6Data <- read_csv(here("Data","TiltMeterData", "Tilt6_CurrentData.csv"))
Tilt7Data <- read_csv(here("Data","TiltMeterData", "Tilt7_CurrentData.csv"))
Tilt7Data <- Tilt7Data[,-1] # Had to remove a column that was created giving each entry a point
Tilt8Data <- read_csv(here("Data","TiltMeterData", "Tilt8_CurrentData.csv"))
Tilt8Data <- Tilt8Data[,-1]

#Load Temperature Data

Tilt1Temp <- read_csv(here("Data/TiltMeterData/Tilt1_TemperatureData.csv"))
Tilt2Temp <- read_csv(here("Data/TiltMeterData/Tilt2_TemperatureData.csv"))
Tilt3Temp <- read_csv(here("Data/TiltMeterData/Tilt3_TemperatureData.csv"))
Tilt4Temp <- read_csv(here("Data/TiltMeterData/Tilt4_TemperatureData.csv"))
Tilt5Temp <- read_csv(here("Data/TiltMeterData/Tilt5_TemperatureData.csv"))
Tilt6Temp <- read_csv(here("Data/TiltMeterData/Tilt6_TemperatureData.csv"))
Tilt7Temp <- read_csv(here("Data/TiltMeterData/Tilt7_TemperatureData.csv"))
Tilt8Temp <- read_csv(here("Data/TiltMeterData/Tilt8_TemperatureData.csv"))

#Load Tide Data, from Okinawa Tide Station 

TideData <- read_csv(here("Data/TiltMeterData/TideData_2025.08.csv"))

#### Organizing the Tilt Speed Data ####

declination <- -5.85 # West declination in degrees to adjust magnetic north (which the sensors take) to true north
#declination is the same for all sensors 

Tilt_Data_fun <- function(df, tilt_name){  #Function to add a tilt name column and separate data and time column 
df %>%
    
separate(col = DateTime, #Choose the column that you want separated
             into = c("Date", "Time"), # Separate it into two columns "Data" and "Time"
             sep = " ", #What are the two separated by 
             remove = FALSE) %>% #Keep the orginial DateTime Column    
    
mutate(DateTime = lubridate::ymd_hms(gsub(" UTC", "", DateTime), tz = "UTC"), #Set datetime,a dding in the UTC timezone so wer are removing taht
       Date = lubridate::ymd(Date), #Set Date
       Time = lubridate::hms(Time), #Set Time
       Speed_m = Speed * 0.01, #Change cm/s -> m/s
       VelocityN_m = Velocity_N * 0.01, #Change cm/s -> m/s
       VelocityE_m = Velocity_E* 0.01, #Change cm/s -> m/s
       Heading_True = (Heading + declination) %% 360, #Adjusting magnetic north Heading to true north using the declination of the area
       TiltMeter = tilt_name) %>%
    
    select(TiltMeter, everything()) #Moving TiltMeter# Column to be the first column
    
}

Tilt1DataClean <- Tilt_Data_fun(Tilt1Data, "Tilt1")
Tilt2DataClean <- Tilt_Data_fun(Tilt2Data, "Tilt2")
Tilt3DataClean <- Tilt_Data_fun(Tilt3Data, "Tilt3")
Tilt4DataClean <- Tilt_Data_fun(Tilt4Data, "Tilt4")
Tilt5DataClean <- Tilt_Data_fun(Tilt5Data, "Tilt5")
Tilt6DataClean <- Tilt_Data_fun(Tilt6Data, "Tilt6")
Tilt7DataClean <- Tilt_Data_fun(Tilt7Data, "Tilt7")
Tilt8DataClean <- Tilt_Data_fun(Tilt8Data, "Tilt8")

####Summary Stats of all TiltMeter Speed Data ####

JoinedData <- rbind(Tilt1DataClean, Tilt2DataClean, Tilt3DataClean, Tilt4DataClean, Tilt5DataClean,  Tilt6DataClean,  Tilt7DataClean, Tilt8DataClean)
view(JoinedData) #Joined the data so that we can look at min and max values for all tilt meters


## Visualizing the summary stats to see where to filter the data off at ##

summary(JoinedData$Speed_m)
is_outlier(JoinedData$Speed_m) %>%
  table()
lower_bound <- quantile(JoinedData$Speed_m, 0.001)
lower_bound #0.1% of the data lower quantile number
Upper_bound <- quantile(JoinedData$Speed_m, 0.999)
Upper_bound #99.9% of the data upper bound limit (0.2501416 )


boxplot.stats(JoinedData$Speed_m)

sum(JoinedData$Speed_m >= .5) #Counting the number of points above .5
# = 5

sum(JoinedData$Speed_m > 0.2501416) #Number of points above the upper bound, 0.25, 145 points

table(cut(JoinedData$Speed_m , seq(0, .25,.05), include.lowest = TRUE)) ##Showing me the different bins and how many data points I have in each one

Means <- JoinedData %>%
  group_by(TiltMeter) %>%
  summarise(mean(Speed_m))
Means #Giving the average speed m/s of each tilt meter

#Visualize the data 

ggplot(JoinedData, aes(x = DateTime, y = Speed_m)) +
  geom_point() +
   facet_wrap(~TiltMeter)


### Windrose plot of Tilt Meter Speed Data ####

tilt_fun<-function(df, plotname){ #Creating a function to run all tilt meter data easier

    ggwindrose(
    speed = df$Speed_m, 
    direction = df$Heading_True,
    n_directions = 8,
    speed_cuts = seq(0, 0.25,.05), Inf) + 
    #stack_revers = TRUE) + #If we want the higher speeds on the outside of the plot
  

    scale_fill_manual(
      values = c( "skyblue","orange" ,"green2", "blue","red",  "darkred"),
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

tilt_fun(df= Tilt1DataClean, plotname = "tilt1b.jpeg") #Saving as jpeg to upload into googleEarth with a white background
tilt_fun(Tilt2DataClean, plotname = "tilt2b.jpeg")
tilt_fun(Tilt3DataClean, plotname = "tilt3b.jpeg")
tilt_fun(Tilt4DataClean, plotname = "tilt4b.jpeg")
tilt_fun(Tilt5DataClean, plotname = "tilt5b.jpeg")
tilt_fun(Tilt6DataClean, plotname = "tilt6b.jpeg")
tilt_fun(Tilt7DataClean, plotname = "tilt7b.jpeg")
tilt_fun(Tilt8DataClean, plotname = "tilt8b.jpeg")




#### Stick Plot of Tilt Meter Speed Data ####
# --- velocity components VelocityE_m, VelocityN_m set teh arrow tile, arrows points in the true direction fo the current
# --- length of the arrow reflects the speed
# --- Tip: If you want the y-axis labeled in real m/s, keep VelocityN_m as is (no scale factor on yend), and only scale the xend.

# We want to plot the data in hours instead of every minute so that we can compare with the tide data #

#### Plotting Tide Data ####


Velocity_TidePlot_fun <- function(df, TideData, plotname, combinedplot){
  
Newdf <- df %>%
group_by(DateTime = floor_date(DateTime, "hour")) %>% #New datetime that is rounded to the start of each hour
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop")

speed_limits <- c(0, 0.25) #Setting the min/max speeds so that each graph has the same scale (m/s)
time_limits <- as.POSIXct(range(Newdf$DateTime, na.rm = TRUE))

TideData$DateTime <- as.POSIXct(TideData$DateTime, tz = "UTC")

plot1 <- ggplot(Newdf) +
 
  geom_segment(aes( x = DateTime, y = 0, 
                    xend = DateTime + VelocityE_m, #Makes the arrows lean left/right
                    yend = VelocityN_m, #Makes the arrows point up/down
                    color = Speed_m), 
               arrow = arrow(length = unit(0.2, "cm")),
               linewidth = 0.5) +
 
  scale_color_viridis_c(limits = speed_limits, 
                        breaks = seq(0, 0.25,.05), 
                        oob = scales::squish, 
                        direction = -1) +
  scale_x_datetime(limits =  time_limits) +
  labs(
     
    x = "Time", 
    y = "Velocity", 
    color = "Speed (m/s)") +
  theme_minimal()  
 
 ggsave(here("Output","TiltMeterData",plotname))

 TidePlot <- ggplot(TideData, aes(DateTime, Sealevel_m)) +
   geom_line(color = "steelblue") +
   scale_x_datetime(limits = time_limits) 
   
 
 
 combined <- plot_grid(plot1, TidePlot, ncol = 1, align = "v", rel_heights =  c(2,1))


 ggsave(here::here("Output", "TiltMeterData", combinedplot))
 
 return(combined) #View the plot

}


Velocity_TidePlot_fun(Tilt1DataClean, TideData, plotname = "Tilt1Velocity.png", combinedplot = "Tilt1TideCombo.png")
Velocity_TidePlot_fun(Tilt2DataClean, TideData, plotname = "Tilt2Velocity.png", combinedplot = "Tilt2TideCombo.png")
Velocity_TidePlot_fun(Tilt3DataClean, TideData, plotname = "Tilt3Velocity.png", combinedplot = "Tilt3TideCombo.png")
Velocity_TidePlot_fun(Tilt4DataClean, TideData, plotname = "Tilt4Velocity.png", combinedplot = "Tilt4TideCombo.png")
Velocity_TidePlot_fun(Tilt5DataClean, TideData, plotname = "Tilt5Velocity.png", combinedplot = "Tilt5TideCombo.png")
Velocity_TidePlot_fun(Tilt6DataClean, TideData, plotname = "Tilt6Velocity.png", combinedplot = "Tilt6TideCombo.png")
Velocity_TidePlot_fun(Tilt7DataClean, TideData, plotname = "Tilt7Velocity.png", combinedplot = "Tilt7TideCombo.png")
Velocity_TidePlot_fun(Tilt8DataClean, TideData, plotname = "Tilt8Velocity.png", combinedplot = "Tilt8TideCombo.png")
 

## -- Adding on the Tide Data ## 
 TideData <- TideData %>%
   mutate(DateTime = ymd_hms(paste(Date, sprintf("%02d:00:00", Hour))), #Combining Data and time into one column 
          Sealevel_m = Sealevel_mm / 1000) #Convert mm to m 

 
 
 TideData <- TideData %>%
   filter(between(DateTime, ymd_hms("2025-08-01 11:00:00"), ymd_hms("2025-08-15 15:00:00"))) #Filtering out the dates that data was not recorded during 
 
 
 TidePlot <- ggplot(TideData, aes(DateTime, Sealevel_m)) +
   geom_line(color = "steelblue") +
   scale_x_datetime(limits = time_limits) +
   labs(title = "Tide Plot")
 
 
 plot_grid(p1, TidePlot, ncol = 1, align = "v", rel_heights =  c(2,1))
 
 

 #ggarrange(p1, TidePlot, nrow = 2) #Another way to add the two graphs together
 

#### Organizing Temperature Data ####
Tilt_Temp_fun <- function(df){ 
df %>%
separate(col = DateTime, #Choose the column that you want separated
         into = c("Date", "Time"), # Separate it into two columns "Data" and "Time"
         sep = " ", #What are the two separated by 
         remove = FALSE) %>% #Keep the orginial DateTime Column    
  mutate(DateTime = lubridate::ymd_hms(DateTime), #Set datetime
         Date = lubridate::ymd(Date), #Set Date
         Time = lubridate::hms(Time))
}

Tilt1TempData <- Tilt_Temp_fun(Tilt1Temp)
Tilt2TempData <- Tilt_Temp_fun(Tilt2Temp)
Tilt3TempData <- Tilt_Temp_fun(Tilt3Temp)
Tilt4TempData <- Tilt_Temp_fun(Tilt4Temp)
Tilt5TempData <- Tilt_Temp_fun(Tilt5Temp)
Tilt6TempData <- Tilt_Temp_fun(Tilt6Temp)
Tilt7TempData <- Tilt_Temp_fun(Tilt7Temp)
Tilt8TempData <- Tilt_Temp_fun(Tilt8Temp)

####Plotting Temperature Data ####
Tilt_TempPlot_fun <- function(df, plotname){ 

  ggplot(df, aes(DateTime, Temperature_C)) +
  geom_point()
  

ggsave(here("Output", "TiltMeterData", plotname))
}

Tilt_TempPlot_fun(Tilt1TempData, plotname = "Tilt1TempPlot.png")
Tilt_TempPlot_fun(Tilt2TempData, plotname = "Tilt2TempPlot.png")
Tilt_TempPlot_fun(Tilt3TempData, plotname = "Tilt3TempPlot.png")
Tilt_TempPlot_fun(Tilt4TempData, plotname = "Tilt4TempPlot.png")
Tilt_TempPlot_fun(Tilt5TempData, plotname = "Tilt5TempPlot.png")
Tilt_TempPlot_fun(Tilt6TempData, plotname = "Tilt6TempPlot.png")
Tilt_TempPlot_fun(Tilt7TempData, plotname = "Tilt7TempPlot.png")
Tilt_TempPlot_fun(Tilt8TempData, plotname = "Tilt8TempPlot.png")




#CorrectedSealevel_mm = TideData$Sealevel_mm + TideData$DatumConstant_mm - TideData$FixedPointHeight_mm
#(CorrectedSealevel_mm/1000)

### The rest of the script is testing different plots ####

# Pre-bin the data by direction + speed bins
rose_data <- 
  Tilt1DataClean %>%
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
Tilt1DataClean <- Tilt1DataClean %>%
  mutate(

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

