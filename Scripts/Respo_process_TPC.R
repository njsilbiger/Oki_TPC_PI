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
if ("nls.multstart" %in% rownames(installed.packages()) == 'FALSE') install.packages('nls.multstart')
if ("rTPC" %in% rownames(installed.packages()) == 'FALSE') install.packages('rTPC')

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
library(dplyr)

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
#View(RespoMeta)
## try first with prelim fake data to make sure script runs
## then switch to real calculated data after getting volumes and weight and surface area

# species names
sp_names <- read_csv(here("Data", "species_names.csv"))

# join the data together
BioData <- BioData %>% 
  full_join(sp_names) %>%
  dplyr::select(-full_species, -species_ID)

RespoMeta <- RespoMeta %>% 
  dplyr::select(-notes)

Sample_Info <- left_join(RespoMeta, BioData)
#View(Sample_Info)

##### Make sure times are consistent ####
# make start and stop times real times, so that we can join the respo output and sample_info data frames
Sample_Info <- Sample_Info %>% 
  #drop_na(sample_ID) %>% 
  unite(date,start_time,col="start_time",remove=F, sep=" ") %>% 
  unite(date,stop_time,col="stop_time",remove=F, sep=" ") %>%
  mutate(start_time = mdy_hms(start_time)) %>% 
  mutate(stop_time = mdy_hms(stop_time)) %>% 
  mutate(date = mdy(date))

#View(Sample_Info)
write_csv(Sample_Info, here("Data","RespoFiles","TPC","Sample_Info.csv"))
Sample_Info <- read_csv(here("Data","RespoFiles","TPC","Sample_Info.csv"))

#generate a 4 column dataframe with specific column names
# data is in umol.L.sec

n_temp_levels<-9 # number of unique temperature levels

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
        arrange(Time) %>%
        mutate(t_sec = as.numeric(difftime(Time, first(Time), units = "secs"))) %>% #keep everything in seconds
        mutate(light_dark = light_dark,
               temp_c_value = temp_c_value) %>%
        filter(t_sec > 120) %>%                          # drop first 2 min (120 s)
        filter(row_number() %% 10 == 0)                  # keep every 10th row
      # now t_sec increments by ~20 s for your kept rows
      #   filter(Time >= start_time & Time <= stop_time) %>%
      #   mutate(sec = row_number()) %>%# add an id for each row to help remove the first few mins
      #   mutate(light_dark = light_dark,
      #          temp_c_value = temp_c_value,
      #          sec = sec) %>%
      #   filter(sec > 60)  %>%# delete the first 2 mins of data assuming freq of 2 Hz
      #   mutate(row_number = row_number()) %>%
      #   filter(row_number %% 10 == 0) %>%  # keep every 10th row only to thin the data
      #   dplyr::select(-row_number) %>%
      #   mutate(sec2 = row_number())  #update the row numbers
      # #return(subset)
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
    rankLocReg(xall = data$t_sec, yall = data$Value, #changed to secs in time instead of sec2 variable
               alpha = 0.2, method = "pc", verbose = FALSE)
  }
  
  # Setup for parallel processing
  future::plan(multisession)
  
  # Map LoLinR function onto all intervals of each sample's thinned dataset
  df <- combined_oxy %>%
    dplyr::select(t_sec, Value, temp_c_value, Temp, light_dark)%>%
    #mutate(sec2 = as.numeric(sec2))%>%
    nest_by(temp_c_value, light_dark) %>%
    ungroup()%>%
    mutate(regs = furrr::future_map(data, fit_reg), # run the LOLinR fit in parallel
           Temp.C = map_dbl(map(data, "Temp"), mean),# get the mean temperature
           RegStats =map(regs, function(x){ # extract the intercept and slope for the parameters
             x$allRegs %>%
               slice(1) %>%
               dplyr::select(Intercept = b0,
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
    dplyr::select(temp_c_value,Temp.C,light_dark,RegStats ) %>%
    unnest(RegStats) %>%
    mutate(sample_ID = rename) %>%
    left_join(Sample_Info[FRow,] %>%
                dplyr::select(sample_ID, temp_c_level, temp_c_value, run_block, light_dark))  # make sure the light value (or whatever other metadata you want) is in the final DF
  
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

#View(RespoR)
######### end of for loop - celebrate victory of getting through that ###############
############################################

#export raw data and read back in as a failsafe 
#this allows me to not have to run the for loop again !!!!!
write_csv(RespoR, here("Data","RespoFiles","TPC","Respo_R.csv"))  
#MP ran for TPC 1-6 (all) on Sept 29th 2025
##### 

######### Calculate Respiration rate ###############
############################################

RespoR <- read_csv(here("Data","RespoFiles","TPC","Respo_R.csv"))
Sample_Info <- read_csv(here("Data","RespoFiles","TPC","Sample_Info.csv"))
ch.vol <- 475 #mL #of small chambers 

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
         umol.cm2.hr = (umol.sec.corr*3600)/SA_cm2,
         umol.cm2.hr_uncorr = (umol.sec*3600)/SA_cm2) %>%
         #mmol.cm2.hr = 0.001*(umol.sec.corr*3600)/SA_cm2, # convert to mmol cm-2 hr-1
         #mmol.cm2.hr_uncorr = 0.001*(umol.sec*3600)/SA_cm2) %>% 
  filter(blank == 0) %>% # remove the Blank data
  #ungroup() %>%
  dplyr::select(date, species, sample_ID, frag_ID, light_dark, run_block, SA_cm2, run_block, umol.cm2.hr, chamber_channel, 
                Temp.C, umol.cm2.hr_uncorr, temp_c_level, temp_c_value) #keep only what we need

RespoR_PR <- RespoR_Normalized %>%
  dplyr::select(-umol.cm2.hr_uncorr, -Temp.C) %>% # remove to pivot
  pivot_wider(names_from = light_dark, values_from = umol.cm2.hr) %>% 
  rename(Respiration = DARK , NetPhoto = LIGHT) %>% # rename the columns
  mutate(Respiration = -1 * Respiration) %>%  # Make respiration positive
  mutate(GrossPhoto = Respiration + NetPhoto) %>% 
  pivot_longer(cols = Respiration:GrossPhoto, names_to = "PR", values_to = "Values") #values still in umol.cm2.hr


write_csv(RespoR_PR,here("Data","RespoFiles","TPC","PnR_rates.csv")) # export all the uptake rates
RespoR_PR <- read_csv(here("Data","RespoFiles","TPC","PnR_rates.csv"))

#dev.off() # may need if plot doesn't run?
PR_plot <- RespoR_PR %>% 
  ggplot(aes(x = temp_c_value, y = Values, group=frag_ID, color = species)) +
  geom_point() +
  geom_line() +
  facet_wrap(~PR*species, scales = "free") +
  theme_bw() +
  labs(x = "Temperature °C", y = "umol.cm2.hr")+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

ggsave(here("Output", "TPC", "Graphs","PR_boxplots.pdf"),
       device = "pdf", height = 8, width = 8, PR_plot)

GP_plot <- RespoR_PR %>% filter(PR == "GrossPhoto") %>%
  ggplot(aes(x = temp_c_value, y = Values, group=frag_ID, color = species, shape = run_block)) +
  geom_point() +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

NP_plot <- RespoR_PR %>% filter(PR == "NetPhoto") %>%
  ggplot(aes(x = temp_c_value, y = Values, group=frag_ID, color = species, shape = run_block)) +
  geom_point() +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

Resp_plot <- RespoR_PR %>% filter(PR == "Respiration") %>%
  ggplot(aes(x = temp_c_value, y = Values, group=frag_ID, color = species, shape = run_block)) +
  geom_point() +
  geom_line() +
  facet_wrap(~species, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

#######################
### making a df for just blank data for future use in plots ### 
Blank_only <- RespoR2 %>% 
  filter(blank == 1) %>% # grab the blanks
  group_by(temp_c_level, temp_c_value, run_block, light_dark) %>%
  #dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  summarise(blank.rate = mean(umol.sec, na.rm = TRUE))

#############################
write_csv(RespoR_Normalized , here("Data","RespoFiles","TPC","Respo_RNormalized_AllTPCRates.csv"))  
RespoR_Normalized <- read_csv(here("Data","RespoFiles","TPC","Respo_RNormalized_AllTPCRates.csv"))

## Plot the blanks across treatments to make sure nothing is funky
Blank_only %>%
  ggplot(aes(x = temp_c_value, blank.rate)) +
  geom_point()

#  basic plot of rates versus temp before you make the real TPC curve

RespoR_Normalized %>%
  ggplot(aes(x = temp_c_value, y = umol.cm2.hr, color = species))+
  geom_point() +
  facet_wrap(~light_dark)

##########################################################
### run an nls model for TPC curves for each species ###

# sharpeschoolhigh_1981
# https://padpadpadpad.github.io/rTPC/articles/rTPC.html

# first run Topt by species, then if not working, do it by FragID
# is I get an error, starting values are off, so play around with starting values that would make sense

# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(here)

# read in data
df <- read_csv(here("Data","RespoFiles","TPC","PnR_rates.csv"))
df <- df %>% mutate(full_ID = paste(sample_ID,run_block,temp_c_value,PR)) %>% mutate(id_block = paste(sample_ID,run_block,PR))
#log_df <- df %>% mutate(logValues = log(Values))


#create outlier IQR function
is_outlier_iqr <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  (x < (q1 - 1.5 * iqr)) | (x > (q3 + 1.5 * iqr))
}

#add column to flag outliers
#just grouped by species and PR, could also do more finescale grouping
df_iqr <- df %>%
  group_by(PR, species) %>%
  mutate(outlier_iqr = if (n() >= 4) is_outlier_iqr(Values) else FALSE) %>%
  ungroup()

#look at graph of IQR values to remove
ggplot(df_iqr, aes(temp_c_value, Values, color = outlier_iqr)) +
  geom_point() +
  facet_wrap(PR ~ species, scales = "free_y") +
  theme_bw()

#now remove those from runs where <4 data points were removed by IQR removal
df_ir_check <- df_iqr %>% filter(outlier_iqr)
bad_ids <- c("C06_TPC RUN4 Respiration", "D04_TPC RUN6 NetPhoto", "I08_TPC RUN6 NetPhoto")

df_iqr <- df %>%
  mutate(outlier_manual = id_block %in% bad_ids) %>% #manual flags
  group_by(PR, species) %>%
  mutate(outlier_iqr = if (n() >= 4) is_outlier_iqr(Values) else FALSE) %>%
  ungroup() %>%
  mutate(outlier_any = outlier_iqr | outlier_manual)   

#check updated graph of IQR values to remove
ggplot(df_iqr, aes(temp_c_value, Values, color = outlier_any)) +
  geom_point() +
  facet_wrap(PR ~ species, scales = "free_y") +
  theme_bw()

#drop IQR outliers
df_clean_iqr <- df_iqr %>% filter(!outlier_any)
write_csv(df_clean_iqr, here("Data","RespoFiles","TPC","PnR_clean_iqr.csv"))
#removes 66 data points

#choose model for rTPC: 'sharpschoolhigh_1981'

#generate empty dataframes to fill in data for topt
topt_df <- tibble(rmax = as.numeric(),
                  topt = as.numeric(),
                  ctmin = as.numeric(),
                  ctmax = as.numeric(),
                  e = as.numeric(),
                  eh = as.numeric(),
                  q10 = as.numeric(),
                  thermal_safety_margin = as.numeric(),
                  thermal_tolerance = as.numeric(),
                  breadth = as.numeric(),
                  skewness = as.numeric(),
                  frag_ID = as.character(),
                  PR = as.character())

#and predicted data
preds_all <- tibble(
  PR = as.character(),
  frag_ID = as.character(),
  temp_c_value = as.numeric(),
  .fitted = numeric()
)

new_data <- tibble(temp_c_value = c(24.5, 26, 27, 28, 29, 30, 31, 32, 34))

for (j in unique(df_clean_iqr$PR)){
  pr = j
  PR_df <- df_clean_iqr %>%
    filter(PR == pr)
  
  for(i in unique(PR_df$frag_ID)){
    id = i
    my_df <- PR_df %>%
      filter(frag_ID == id)
    
    # get start vals
    start_vals <- get_start_vals(my_df$temp_c_value, my_df$Values, model_name = 'sharpeschoolhigh_1981')
    
    #issues with setting start values with some parameters
    #skip these
    if (anyNA(start_vals) || any(!is.finite(start_vals))) {
      message(sprintf("Skipping %s %s: could not set start values", id, pr))
      next
    }
    
    # get limits
    low_lims <- get_lower_lims(my_df$temp_c_value, my_df$Values, model_name = 'sharpeschoolhigh_1981')
    upper_lims <- get_upper_lims(my_df$temp_c_value, my_df$Values, model_name = 'sharpeschoolhigh_1981')
    
    # fit model
    fit <- nls_multstart(Values~sharpeschoolhigh_1981(temp = temp_c_value, r_tref,e,eh,th, tref = 15),
                         data = my_df,
                         iter = 500,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = low_lims,
                         upper = upper_lims,
                         supp_errors = 'Y')
    
    fit
    
    # calculate additional traits
    topt_params <- calc_params(fit) %>%
      # round for easy viewing
      mutate_all(round, 2) %>% 
      mutate(PR = pr,
             frag_ID = id)
    
    topt_df <- topt_df %>%
      rbind(topt_params)
    
    #generate predictions
    preds_id <- augment(fit, newdata = new_data) %>%
      transmute(PR = pr,
                frag_ID = id,
                temp_c_value = new_data$temp_c_value,
                .fitted = .fitted)
    
    preds_all <- bind_rows(preds_all, preds_id)
  }
}

#Values that this code skipped because couldn't set start values:
# Skipping C01 Respiration: could not set start values
# Skipping F01 Respiration: could not set start values
# Skipping E01 Respiration: could not set start values
# Skipping J01 Respiration: could not set start values
# Skipping I01 Respiration: could not set start values
# Skipping B01 Respiration: could not set start values
# Skipping D01 Respiration: could not set start values
# Skipping E02 Respiration: could not set start values
# Skipping F04 Respiration: could not set start values
# Skipping D02 Respiration: could not set start values
# Skipping C02 Respiration: could not set start values
# Skipping J02 Respiration: could not set start values
# Skipping B02 Respiration: could not set start values
# Skipping E03 Respiration: could not set start values
# Skipping I03 Respiration: could not set start values
# Skipping A03 Respiration: could not set start values
# Skipping H04 Respiration: could not set start values
# Skipping F03 Respiration: could not set start values
# Skipping J03 Respiration: could not set start values
# Skipping G05 Respiration: could not set start values
# Skipping I06 Respiration: could not set start values
# Skipping E06 Respiration: could not set start values
# Skipping F06 Respiration: could not set start values
# Skipping B07 Respiration: could not set start values
# Skipping A07 Respiration: could not set start values
# Skipping E07 Respiration: could not set start values
# Skipping J07 Respiration: could not set start values
# Skipping C07 Respiration: could not set start values
# Skipping I07 Respiration: could not set start values
# Skipping H07 Respiration: could not set start values
# Skipping B08 Respiration: could not set start values
# Skipping E08 Respiration: could not set start values
# Skipping D04 Respiration: could not set start values
# Skipping J08 Respiration: could not set start values
# Skipping I08 Respiration: could not set start values
# Skipping G08 Respiration: could not set start values
# Skipping A06 NetPhoto: could not set start values

#add species names to predictions and topt dfs
BioData <- read_csv(here("Data","RespoFiles","TPC","Fragment_Measurements_TPC.csv"))

preds_all <- preds_all %>% left_join(BioData, by = "frag_ID")
topt_df <- topt_df %>% left_join(BioData, by = "frag_ID")

#save data files
#clean = outliers removed
write_csv(topt_df, here("Data","RespoFiles","TPC","Topt_data_clean.csv"))
write_csv(preds_all, here("Data","RespoFiles","TPC","Preds_data_clean.csv"))

#read in dataframes and generate prediction dfs for each metric to graph
PnR_clean_iqr <- read_csv(here("Data","RespoFiles","TPC","PnR_clean_iqr.csv"))
preds_all <- read_csv(here("Data","RespoFiles","TPC","Preds_data_clean.csv"))
preds_gp <- preds_all %>% filter(PR == "GrossPhoto")
preds_np <- preds_all %>% filter(PR == "NetPhoto")
preds_resp <- preds_all %>% filter(PR == "Respiration")

######plots of predicted TPCs with data#####

#gross photo
gp_pred_plot <- PnR_clean_iqr %>% filter(PR == "GrossPhoto") %>% ggplot(aes(temp_c_value, Values)) +
  geom_point(alpha = 0.7) +
  geom_line(data = preds_gp,
            aes(temp_c_value, .fitted, group = frag_ID),
            linewidth = 0.6, color = "blue") +
  facet_wrap(~ species, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(x = "Temperature (ºC)",
       y = "Gross Photosynthesis (umol.cm2.hr)",
       title = "Thermal performance: gross photosynthesis by species")

gp_pred_plot

ggsave(here("Output", "TPC", "Graphs", "gp_predicted_plot_clean.pdf"),
       device = "pdf", height = 8, width = 6, gp_pred_plot)

#net photo
np_pred_plot <- PnR_clean_iqr %>% filter(PR == "NetPhoto") %>% ggplot(aes(temp_c_value, Values)) +
  geom_point(alpha = 0.7) +
  geom_line(data = preds_np,
            aes(temp_c_value, .fitted, group = frag_ID),
            linewidth = 0.6, color = "blue") +
  facet_wrap(~ species, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(x = "Temperature (ºC)",
       y = "Net Photosynthesis",
       title = "Thermal performance: net photosynthesis by species")

np_pred_plot

ggsave(here("Output", "TPC", "Graphs","np_predicted_plot_clean.pdf"),
       device = "pdf", height = 8, width = 6, np_pred_plot)

#respiration 
resp_pred_plot <- PnR_clean_iqr %>% filter(PR == "Respiration") %>% ggplot(aes(temp_c_value, Values)) +
  geom_point(alpha = 0.7) +
  geom_line(data = preds_resp,
            aes(temp_c_value, .fitted, group = frag_ID),
            linewidth = 0.6, color = "blue") +
  facet_wrap(~ species, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(x = "Temperature (ºC)",
       y = "Respiration",
       title = "Thermal performance: respiration by species")

resp_pred_plot

ggsave(here("Output", "TPC", "Graphs","resp_predicted_plot_clean.pdf"),
       device = "pdf", height = 8, width = 6, resp_pred_plot)


#all species all rates plot
species_all_plot <- PnR_clean_iqr %>% ggplot(aes(temp_c_value, Values, color = species)) +
  geom_point(alpha = 0.7) +
  geom_line(data = preds_all,
            aes(temp_c_value, .fitted, group = frag_ID, color = species),
            linewidth = 0.6) +
  facet_wrap(~ PR, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(x = "Temperature (ºC)",
       y = "Metabolic rate (umol.cm2.hr)",
       title = "Thermal performance")

species_all_plot

ggsave(here("Output", "TPC", "Graphs", "all_rates_all_species_plot_clean.pdf"),
       device = "pdf", height = 8, width = 6, species_all_plot)

##########################################################
####stats to look at differences in thermal performance metrics####
#interested in the effect of species (species) on these parameters:
#rmax, topt, e, breadth
#include random effect of genotype (frag_ID)

#load libraries
library(here)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(DHARMa)

#read in data
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean.csv"))

#first visualize data - plots of Topt and other variables from TPCs######
rmax_plot <- topt_df %>% 
  ggplot(aes(species, rmax, fill = species)) +
  geom_violin(drop = FALSE)+
  geom_point(alpha = 0.7) +
  facet_wrap(~ PR, scales = "free_y") +
  theme_bw(base_size = 12)+
  labs(x = "Coral Species", 
       y = "Rmax")

topt_plot <- topt_df %>% 
  ggplot(aes(species, topt, fill = species)) +
  geom_violin(drop = FALSE)+
  geom_point(alpha = 0.7) +
  facet_wrap(~ PR, scales = "free_y") +
  theme_bw(base_size = 12)+
  labs(x = "Coral Species", 
       y = "Topt")

e_plot <- topt_df %>% 
  ggplot(aes(species, e, fill = species)) +
  geom_violin(drop = FALSE)+
  geom_point(alpha = 0.7) +
  facet_wrap(~ PR, scales = "free_y") +
  theme_bw(base_size = 12)+
  labs(x = "Coral Species", 
       y = "e")

breadth_plot <- topt_df %>% 
  ggplot(aes(species, breadth, fill = species)) +
  geom_violin(drop = FALSE)+
  geom_point(alpha = 0.7) +
  facet_wrap(~ PR, scales = "free_y") +
  theme_bw(base_size = 12)+
  labs(x = "Coral Species", 
       y = "breadth")

#Meeting with Nyssa
#plot as mean + SE
#filter out run 4 - run predictions and see how the data looks after that (seems like it's giving the most errors/usses)
#plot predicted with Topt with geom vline to check and make sure
#only remove points we think are not accurate
#scatterplot with dry weight and rmax or other parameter
#is it actually species or just the physiology that's driving topt?
#ancova style - with individual lines for each species
#see which variables are related to eachother
#later - see what depths
#net or gross whatever is cleaner - just focus on one
#HI = species have different physiology (anovas of all parameters)
#H2 = differences in phys translate to thermal performance (regresssions)
#H3 = morphology
#drop frags that we don't have data for
#H1 and H2, pretty pictures, maps are nice - simple star with data collection and 
#temperature data - changes by 3C, 27-30-27
#HOBO data - look for this
#yearly - NOAA data
#thermal variability is wild and that's why its cool to look at this
#tomorrow - meet again

#check distributions of data as well
ggplot(topt_df, aes(rmax)) + #look pretty ok
  geom_histogram(bins = 30) + 
  facet_wrap(~PR, scales = "free")
ggplot(topt_df, aes(topt)) + #slight left skew for GP and NP, resp won't work
  geom_histogram(bins = 30) + 
  facet_wrap(~PR, scales = "free")
ggplot(topt_df, aes(e)) + #a couple high variables for GP and NP, resp weird
  geom_histogram(bins = 30) + 
  facet_wrap(~PR, scales = "free")
ggplot(topt_df, aes(breadth)) + #left skewed for GP, NP good, resp weird
  geom_histogram(bins = 30) + 
  facet_wrap(~PR, scales = "free")

#all params
topt_gp <- topt_df %>% filter(PR == "GrossPhoto")
topt_np <- topt_df %>% filter(PR == "NetPhoto")
topt_resp <- topt_df %>% filter(PR == "Respiration")

#run lm for rmax
rmax_lm <- lm(rmax ~ species, data = topt_gp)
anova(rmax_lm)

#checks
simres <- simulateResiduals(rmax_lm)
plot(simres)

emm_rmax <- emmeans(rmax_lm, ~ species)
pairs(emm_rmax, adjust = "tukey")
emm_tbl_rmax <- as_tibble(summary(emm_rmax, infer = TRUE))

ggplot(emm_tbl_rmax, aes(emmean, species, color = species)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.2) +
  theme_bw(base_size = 12) +
  labs(x = "rmax est", y = "Species")


#run lm for topt
topt_lm <- lm(topt ~ species, data = topt_gp)
anova(topt_lm)

#checks
simres <- simulateResiduals(topt_lm)
plot(simres)

emm_topt <- emmeans(topt_lm, ~ species)
pairs(emm_topt, adjust = "tukey")
emm_tbl_topt <- as_tibble(summary(emm_topt, infer = TRUE))

ggplot(emm_tbl_topt, aes(emmean, species, color = species)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.2) +
  theme_bw(base_size = 12) +
  labs(x = "Topt est", y = "Species")

#run lm for e
e_lm <- lm(e ~ species, data = topt_gp)
anova(e_lm)

#checks
simres <- simulateResiduals(e_lm)
plot(simres)

emm_e <- emmeans(e_lm, ~ species)
pairs(emm_e, adjust = "tukey")
emm_tbl_e <- as_tibble(summary(emm_e, infer = TRUE))

ggplot(emm_tbl_e, aes(emmean, species, color = species)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.2) +
  theme_bw(base_size = 12) +
  labs(x = "e est", y = "Species")

#run lm for breadth
breadth_lm <- lm(breadth ~ species, data = topt_gp)
anova(breadth_lm)

#checks
simres <- simulateResiduals(breadth_lm)
plot(simres)

emm_breadth <- emmeans(breadth_lm, ~ species)
pairs(emm_breadth, adjust = "tukey")
emm_tbl_breadth <- as_tibble(summary(emm_breadth, infer = TRUE))

ggplot(emm_tbl_breadth, aes(emmean, species, color = species)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.2) +
  theme_bw(base_size = 12) +
  labs(x = "breadth est", y = "Species")

