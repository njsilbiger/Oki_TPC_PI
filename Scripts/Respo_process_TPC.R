###### Respo Code for TPC Light and Dark Runs ####### 
### Created by: Nyssa Silbiger
#### Last updated on: 2025-07-26

############## Introduction to code/script ####################
## this script will help us process the raw data gathered during respirometry runs. 
## need to change for specific project/experimental variables 

### Install Packages #####
## if these packages are not yet installed, install them 
## great for updates or new users 
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented')
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix')
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra')
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') devtools::install_github('colin-olito/LoLinR')
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron')
if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse')
if ("here" %in% rownames(installed.packages()) == 'FALSE') install.packages('here')
if ("patchwork" %in% rownames(installed.packages()) == 'FALSE') install.packages('patchwork')
if ("PNWColors" %in% rownames(installed.packages()) == 'FALSE') install.packages('PNWColors')

#Read in required libraries
##### Include Versions of libraries
library(segmented)
library(plotrix)
library(gridExtra)
library(LoLinR)
library(lubridate)
library(chron)
library(patchwork)
library(tidyverse)
library(here)
library(PNWColors)
library(ggrepel)
library(reshape2)
library(viridis)
library(car)
library(future)
library(furrr)

############# now it's time to code ############
################################################
# get the file path

#set the path to all of the raw oxygen datasheets
## these are saved onto the computer in whatever file path/naming scheme you saved things to 
path.p<-here("Data","RespoFiles","TPC","RawO2") #the location of all your respirometry files
#you can change to individual run folders if needed

# bring in all of the individual files
filenames_final<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders

#basename above removes the subdirectory name from the file, re-name as file.names.full
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) 

#empty chamber volume
ch.vol <- 475 #mL #of small chambers 

######### Load and tidy files ###############
############################################
#Load your respiration data file, with all the times, water volumes(mL), #not doing dry weight just SA
#RespoMeta <- read_csv(here("Data","RespoFiles","Respo_Metadata_SGDDilutions_Cabral_Varari.csv"))
BioData <- read_csv(here("Data","RespoFiles","TPC","Fragment_Measurements_TPC.csv"))

RespoMeta <- read_csv(here("Data","RespoFiles","TPC","TPC_meta.csv"))
#View(BioData)
## try first with prelim fake data to make sure script runs
## then switch to real calculated data after getting volumes and weight and surface area


# join the data together
Sample_Info <- left_join(RespoMeta, BioData)
#View(Sample.Info)

##### Make sure times are consistent ####
# make start and stop times real times, so that we can join the respo output and sample_info data frames
Sample_Info <- Sample_Info %>% 
  #drop_na(sample_ID) %>% 
  unite(date,start_time,col="start_time",remove=F, sep=" ") %>% 
  unite(date,stop_time,col="stop_time",remove=F, sep=" ") %>%
  mutate(start_time = mdy_hms(start_time)) %>% 
  mutate(stop_time = mdy_hms(stop_time)) %>% 
  mutate(date = mdy(date))

#view(Sample_Info)

#generate a 4 column dataframe with specific column names
# data is in umol.L.sec

n_temp_levels<-9 # number of unique light levels

RespoR <- tibble(.rows =length(filenames_final)*n_temp_levels*2, # temp * 2 for light and dark runs
                 sample_ID = NA,
                 Intercept = NA,
                 umol.L.sec = NA,
                 Temp.C = NA,
                 temp_c_level = NA,
                 temp_c_value = NA,
                 run_block = NA,
                 light_dark = NA)

######### Create a for loop! ###############
############################################

###forloop##### 
for(i in 1:length(filenames_final)) {
  FRow <- as.numeric(which(Sample_Info$FileID_csv==filenames_final[i])) # stringsplit this renames our file
  
  Respo.Data1 <- read_csv(skip=1,file.path(path.p, paste0(file.names.full[i]))) %>% # reads in each file in list
    dplyr::select(Date, Time, Value, Temp) %>% # keep only what we need: Time stamp per 1sec, Raw O2 value per 1sec, in situ temp per 1sec
    unite(Date,Time,col="Time",remove=T, sep = " ") %>%
    drop_na() %>% 
    mutate(Time = mdy_hms(Time)) #%>% # convert time
  #mutate(help = i) ##if stuck in forloop with error from filter, can check RespoR and see at what row the forloop stopped working  
  
  ## cut the data by start and stop times from metadata
  #Use start time of each light step from the metadata to separate data by light stop
  
  oxy_subsets <- Sample_Info[FRow,] %>%
    pmap(function(temp_c_value, light_dark, start_time, stop_time, ...) {
      data <- Respo.Data1  %>%
        filter(Time >= start_time & Time <= stop_time) %>%
        mutate(sec = row_number()) %>%# add an id for each row to help remove the first few mins
        mutate(light_dark = light_dark,
               temp_c_value = temp_c_value,
               sec = sec) %>%
        filter(sec > 60)  %>%# delete the first 2 mins of data assuming freq of 2 Hz
        mutate(row_number = row_number()) %>%
        filter(row_number %% 10 == 0) %>%  # keep every 10th row only to thin the data
        select(-row_number) %>%
        mutate(sec2 = row_number())  #update the row numbers
      #return(subset)
    }) 
  
  
  # Combine into one long dataframe with ID labels
  combined_oxy <- bind_rows(oxy_subsets)
  
  # Get the filename without the .csv
  rename<- sub("_O2.csv","", filenames_final[i])
  
  ### plot and export the thinned data ####
  p1<- ggplot(combined_oxy, aes(x = sec, y = Value)) +
    geom_point(aes(color = light_dark)) +
    labs(
      x = 'Time (seconds)',
      y = expression(paste(' O'[2],' (',mu,'mol/L)')),
      title = "original"
    )+
    facet_wrap(~temp_c_value, scales = "free_y")
  
  
  ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
  #option to add multiple outputs method= c("z", "e "pc")
  
  # Define function for fitting LoLinR regressions to be applied to all intervals for all samples
  fit_reg <- function(data) {
    rankLocReg(xall = data$sec2, yall = data$Value, 
               alpha = 0.2, method = "pc", verbose = FALSE)
  }
  
  # Setup for parallel processing
  future::plan(multisession)
  
  # Map LoLinR function onto all intervals of each sample's thinned dataset
  df <- combined_oxy %>%
    select(sec2, Value, temp_c_value, Temp, light_dark)%>%
    mutate(sec2 = as.numeric(sec2))%>%
    nest_by(temp_c_value, light_dark) %>%
    ungroup()%>%
    mutate(regs = furrr::future_map(data, fit_reg), # run the LOLinR fit in parallel
           Temp.C = map_dbl(map(data, "Temp"), mean),# get the mean temperature
           RegStats =map(regs, function(x){ # extract the intercept and slope for the parameters
             x$allRegs %>%
               slice(1) %>%
               select(Intercept = b0,
                      umol.L.sec = b1)
           }) )
  
  
  #  Plot regression diagnostics
  df <- df %>% 
    mutate(temp_light = paste(temp_c_value,light_dark))
  
  for(j in 1:length(df$temp_light)){
    pdf(paste0(here("Output","TPC"),"/",rename,"_",j,".pdf" ))
    plot(df$regs[[j]])
    dev.off() 
  }
  
  
  df<-df %>%
    select(temp_c_value,Temp.C,light_dark,RegStats ) %>%
    unnest(RegStats) %>%
    mutate(sample_ID = rename) %>%
    left_join(Sample_Info[FRow,] %>%
                select(sample_ID, temp_c_level, temp_c_value, run_block, light_dark))  # make sure the light value (or whatever other metadata you want) is in the final DF
  
  ################################
  # fill in all the O2 consumption and rate data
  
  RespoR[FRow,"Temp.C"]<-df$Temp.C
  RespoR[FRow, "light_dark"]<-df$light_dark
  RespoR[FRow,"sample_ID"]<-df$sample_ID
  RespoR[FRow,"Intercept"]<-df$Intercept
  RespoR[FRow,"umol.L.sec"]<-df$umol.L.sec
  RespoR[FRow,"temp_c_level"]<-df$temp_c_level
  RespoR[FRow,"temp_c_value"]<-df$temp_c_value
  RespoR[FRow,"run_block"]<-df$run_block
  
}  


######### end of for loop - celebrate victory of getting through that ###############
############################################

#export raw data and read back in as a failsafe 
#this allows me to not have to run the for loop again !!!!!
write_csv(RespoR, here("Data","RespoFiles","TPC","Respo_R.csv"))  

##### 

RespoR <- read_csv(here("Data","RespoFiles","TPC","Respo_R.csv"))

######### Calculate Respiration rate ###############
############################################

RespoR2 <- RespoR %>%
  #drop_na(FileID_csv) %>% # drop NAs
  left_join(Sample_Info) %>% # Join the raw respo calculations with the metadata
  mutate(Ch.Volume.mL = ch.vol) %>% # measured volume of chambers with coral + stand + stirbar displacement
  mutate(Ch.Volume.L = Ch.Volume.mL * 0.001) %>% # mL to L conversion
  mutate(umol.sec = umol.L.sec*Ch.Volume.L) %>% #Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
  mutate_if(sapply(., is.character), as.factor)  #convert character columns to factors

#Account for blank rate by sample run Block (if we do at least one blank per block)

#View(RespoR)

####### normalize the respo rates to the blanks ##### 

RespoR_Normalized <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(temp_c_value, run_block,light_dark) %>%
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE)) %>% # if you have multiple blanks per run take the average
  ungroup() %>% 
  #dplyr::select(Light_level, run_block, blank.rate, date) %>% # this is what I will use to join the blanks back with the raw data
  right_join(RespoR2) %>% # join blanks with the respo data
  mutate(umol.sec.corr = umol.sec - blank.rate, # subtract the blank rates from the raw rates   
         mmol.cm2.hr = 0.001*(umol.sec.corr*3600)/SA_cm2, # convert to mmol cm-2 hr-1
         mmol.cm2.hr_uncorr = 0.001*(umol.sec*3600)/SA_cm2) %>% 
  filter(blank == 0) %>% # remove the Blank data
  #ungroup() %>%
  dplyr::select(date, species, sample_ID, frag_ID, light_dark, run_block, SA_cm2, run_block, mmol.cm2.hr, chamber_channel, 
                Temp.C, mmol.cm2.hr_uncorr, temp_c_level, temp_c_value) #keep only what we need

RespoR_PR <- RespoR_Normalized %>%
  select(-mmol.cm2.hr_uncorr, -Temp.C) %>% # remove to pivot
  pivot_wider(names_from = light_dark, values_from = mmol.cm2.hr) %>% 
  rename(Respiration = DARK , NetPhoto = LIGHT) %>% # rename the columns
  mutate(Respiration = -1 * Respiration) %>%  # Make respiration positive
  mutate(GrossPhoto = Respiration + NetPhoto) %>% 
  pivot_longer(cols = Respiration:GrossPhoto, names_to = "PR", values_to = "Values")


write_csv(RespoR_PR,here("Data","RespoFiles","TPC","PnR_rates.csv")) # export all the uptake rates

#dev.off() # may need if plot doesn't run?
PR_plot <- RespoR_PR %>% 
  ggplot(aes(x = temp_c_value, y = Values, group=frag_ID, color = species)) +
  geom_point() +
  geom_line() +
  facet_wrap(~PR) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

# go back and delete from raw data the weird noisy data drops from O2 probes 1:4

ggsave(here("Output", "TPC", "PR_boxplots.pdf"),
       device = "pdf", height = 8, width = 8, PR_plot)


#######################
### making a df for just blank data for future use in plots ### 
Blank_only <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(temp_c_level, temp_c_value, run_block, light_dark) %>%
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE))

#############################
write_csv(RespoR_Normalized , here("Data","RespoFiles","TPC","Respo_RNormalized_AllTPCRates.csv"))  


## Plot the blanks across treatments to make sure nothing is funky
Blank_only %>%
  ggplot(aes(x = temp_c_value, blank.rate)) +
  geom_point()

#  basic plot of rates versus temp before you make the real TPC curve

RespoR_Normalized %>%
  ggplot(aes(x = temp_c_value, y = mmol.cm2.hr, color = species))+
  geom_point() +
  facet_wrap(~light_dark)

### run an nls model for PI curve and extract Ik for each species ###