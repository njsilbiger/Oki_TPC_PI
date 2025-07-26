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
path.p<-here("Data","RespoFiles","RawO2") #the location of all your respirometry files

# bring in all of the individual files
filenames_final<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders

#basename above removes the subdirectory name from the file, re-name as file.names.full
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) 

#empty chamber volume
ch.vol <- 650 #mL #of small chambers 

######### Load and tidy files ###############
############################################
#Load your respiration data file, with all the times, water volumes(mL), #not doing dry weight just SA
#RespoMeta <- read_csv(here("Data","RespoFiles","Respo_Metadata_SGDDilutions_Cabral_Varari.csv"))
BioData <- read_csv(here("Data","RespoFiles","Fragment_MeasurementSampling_Cabral_Varari.csv"))

RespoMeta <- read_csv(here("Data","RespoFiles","PI_meta_test.csv"))
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
RespoR <- tibble(.rows =length(filenames_final),
                 FileID_csv = NA,
                 Intercept = NA,
                 umol.L.sec = NA,
                 Temp.C = NA)

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
  
 read_respo<-function(filename){
   read_csv(skip=1,file.path(path.p, paste0(file.names.full))) %>% # reads in each file in list
     dplyr::select(Date, Time, Value, Temp) %>% # keep only what we need: Time stamp per 1sec, Raw O2 value per 1sec, in situ temp per 1sec
     unite(Date,Time,col="Time",remove=T, sep = " ") %>%
     drop_na() %>% 
     mutate(Time = mdy_hms(Time))
 }
 
 RespoData<-file.names.full %>%
   map_dfr(read_respo, .id = file.names.full)
 
  ## cut the data by start and stop times from metadata
  #Use start time of each light step from the metadata to separate data by light stop
  
  oxy_subsets <- Sample_Info[FRow,] %>%
      pmap(function(Light_level, start_time, stop_time, ...) {
      data <- Respo.Data1  %>%
        filter(Time >= start_time & Time <= stop_time) %>%
        mutate(sec = row_number()) %>%# add an id for each row to help remove the first few mins
        mutate(Light_level = Light_level,
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
  rename<- sub(".csv","", filenames_final[i])
  
    ### plot and export the thinned data ####
  p1<- ggplot(combined_oxy, aes(x = sec, y = Value)) +
    geom_point(color = "dodgerblue") +
    labs(
      x = 'Time (seconds)',
      y = expression(paste(' O'[2],' (',mu,'mol/L)')),
      title = "original"
    )+
    facet_wrap(~Light_level)
  
  
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
    select(sec2, Value, Light_level, Temp)%>%
    mutate(sec2 = as.numeric(sec2))%>%
    nest_by(Light_level) %>%
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
 
  for(j in 1:length(df$Light_level)){
   pdf(paste0(here("Output"),"/",rename,"_",j,".pdf" ))
    plot(df$regs[[j]])
    dev.off() 
  }
  
  df %>%
    select(Light_level,Temp.C, RegStats ) %>%
    unnest(RegStats) %>%
    mutate(FileID = rename)
    
  ################################
  # fill in all the O2 consumption and rate data
  # need clarity on what this is
  RespoR[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
  RespoR[i,1] <- paste0(rename,".csv") #stores the file name in the Date column
  RespoR[i,4] <- mean(Respo.Data1$Temp, na.rm=T)  #stores the Temperature from the incubation in the Temp.C column
}  


######### end of for loop - celebrate victory of getting through that ###############
############################################

#export raw data and read back in as a failsafe 
#this allows me to not have to run the for loop again !!!!!
write_csv(RespoR, here("Data","RespoFiles","Respo_R.csv"))  

##### 

RespoR <- read_csv(here("Data","RespoFiles","Respo_R.csv"))

######### Calculate Respiration rate ###############
############################################

RespoR2 <- RespoR %>%
  #drop_na(FileID_csv) %>% # drop NAs
  left_join(Sample_Info) %>% # Join the raw respo calculations with the metadata
  #mutate(Ch.Volume.ml = ifelse(is.na(volume_ml),ch.vol,ch.vol-volume_ml)) %>% # add 6 L for volume of all blanks and subtract org volume from chamber vol for all else
  mutate(Ch.Volume.mL = 650-volume_mL) %>% # hannah changed all this volume stuff to match my project
  mutate(Ch.Volume.L = Ch.Volume.mL * 0.001) %>% # mL to L conversion
  mutate(umol.sec = umol.L.sec*Ch.Volume.L) %>% #Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
  mutate_if(sapply(., is.character), as.factor) %>% #convert character columns to factors
  mutate(colony_number= as.factor(colony_number), 
         date_block= factor(date)) #make the blank column a factor

# Remove duplicates when assemblages were run multiple times
# anti <- RespoR2 %>% 
#   filter(Date == "2022-07-14") %>% 
#   distinct(SampleID) %>% 
#   left_join(RespoR2) %>% 
#   filter(Date == "2022-07-13")
# 
# RespoR2 <- RespoR2 %>% 
#   anti_join(anti)


#View(RespoR2)


#Account for blank rate by light/Dark and Block (if we do one blank per block)

#View(RespoR)

####### normalize the respo rates to the blanks ##### 
### need to somehow get blanks for the DAY to be used for runs in that day per dilution 
### for example: blanks for day 1 dils 1-9 should be used to subtract from colonies 7 and 8 from Cabral run in that same day 
## do this for both light and dark 


RespoR_Normalized <- RespoR2 %>% 
  group_by(SGD_dil, date_block, light_dark, colony_number, site) %>%
  filter(colony_number== "BLANK") %>%
  # also add block here if one blank per block --> Hannah's blanks are weird, so took that out of this grouping 
  #group_by(BLANK)%>% # also add block here if one blank per block
  #summarise(umol.sec = mean(umol.sec, na.rm=TRUE)) %>% # get mean value of blanks per run --> don't want because each their own blank only keep the actual blanks
  #group_by(date_block, SGD_dil, light_dark) %>% 
  dplyr::select(blank.rate = umol.sec) %>% ## rename the blank column 
  #mutate(blank.rate = umol.sec) %>%
  #dplyr::select(date, SGD_dil, blank.rate, light_dark) %>% #keep only what we need
  ungroup() %>% 
  dplyr::select(SGD_dil, light_dark, blank.rate, date_block, site) %>%
  #dplyr::select(blank.rate = umol.sec) %>% 
  #arrange(FileID) %>% # rename the blank rate column
  right_join(RespoR2) %>% # join with the respo data
  #arrange(FileID_csv) %>% 
  mutate(umol.sec.corr = umol.sec - blank.rate, # subtract the blank rates from the raw rates #### HOW IS THIS NORMALIZING THE DATA IS EVERYTHING IS JUST GOING TO BE 0 #####   
         mmol.cm2.hr = 0.001*(umol.sec.corr*3600)/SA_cm2, # convert to mmol g-1 hr-1
         mmol.cm2.hr_uncorr = 0.001*(umol.sec*3600)/SA_cm2) %>% 
  filter(colony_number!="BLANK") %>% # remove the Blank data
  #ungroup() %>%
  dplyr::select(date, sample_ID, light_dark, run_block, SA_cm2, run_block, mmol.cm2.hr, chamber_channel, 
                Temp.C, mmol.cm2.hr_uncorr, colony_number, SGD_dil, salinity, pH, site, SGD_number) #keep only what we need


#######################
### making a df for just blank data for future use in plots ### 
RespoR_Normalized_Blanks <- RespoR2 %>% 
  group_by(SGD_dil, date_block, light_dark, colony_number, SGD_number, site) %>%
  filter(colony_number== "BLANK") %>%
  dplyr::select(blank.rate = umol.sec)

RespoR_Normalized_AvgBlanks <- RespoR2 %>% 
  group_by(SGD_number, light_dark, date_block, colony_number, site) %>% 
  filter(colony_number== "BLANK") %>%
  summarise(umol.sec = mean(umol.sec, na.rm=TRUE)) %>% # get mean value of blanks per each dilution 
  dplyr::select(blank.rate = umol.sec, site) %>%  ## rename the blank column 
  ungroup() %>% 
  dplyr::select(SGD_number, light_dark, blank.rate, date_block, site) %>%
  right_join(RespoR2) %>% # join with the respo data
  #arrange(FileID_csv) %>% 
  mutate(umol.sec.corr = umol.sec - blank.rate, # subtract the blank rates from the raw rates #### HOW IS THIS NORMALIZING THE DATA IS EVERYTHING IS JUST GOING TO BE 0 #####   
         mmol.cm2.hr = 0.001*(umol.sec.corr*3600)/SA_cm2, # convert to mmol g-1 hr-1
         mmol.cm2.hr_uncorr = 0.001*(umol.sec*3600)/SA_cm2) %>% 
  filter(colony_number!="BLANK") %>% # remove the Blank data
  #ungroup() %>%
  dplyr::select(date, sample_ID, light_dark, run_block, SA_cm2, run_block, mmol.cm2.hr, chamber_channel, 
                Temp.C, mmol.cm2.hr_uncorr, colony_number, salinity, pH, site, SGD_number)

#############################

########## CALCULATING R AND GP ###############
###############################################

# make the respiration values positive (pull out data for dark treatments)
RespoR_Normalized_dark <- RespoR_Normalized %>% 
  filter(light_dark == "DARK") %>% 
  mutate(mmol.cm2.hr = mmol.cm2.hr*-1,
         mmol.cm2.hr_uncorr = mmol.cm2.hr_uncorr*-1) %>% 
  mutate(mmol.cm2.hr = ifelse(mmol.cm2.hr < 0, 0, mmol.cm2.hr), # for any values below 0, make 0
         mmol.cm2.hr_uncorr = ifelse(mmol.cm2.hr_uncorr < 0, 0, mmol.cm2.hr_uncorr)) %>% 
  mutate(P_R = "R") # all dark run rates get R for respiration

# all light run rates get NP for net photosynthesis
RespoR_Normalized_light <- RespoR_Normalized %>% 
  filter(light_dark == "LIGHT") %>% 
  mutate(mmol.cm2.hr = ifelse(mmol.cm2.hr < 0, 0, mmol.cm2.hr), # for any values below 0, make 0
         mmol.cm2.hr_uncorr = ifelse(mmol.cm2.hr_uncorr < 0, 0, mmol.cm2.hr_uncorr)) %>% 
  mutate(P_R = "NP")

# rejoin data into single df
RespoR_Normalized2 <- full_join(RespoR_Normalized_light, RespoR_Normalized_dark) #%>% 
#drop_na(mmol.gram.hr) # removes anticipated sampleID's that were not actually run


#make column for GP and group by fragment ID and temp to keep R and NP together
RespoR_NormalizedGP <- RespoR_Normalized2 %>% 
  group_by(colony_number, SGD_number, site, salinity, pH) %>% 
  summarize(mmol.cm2.hr = sum(mmol.cm2.hr),
            mmol.cm2.hr_uncorr = sum(mmol.cm2.hr_uncorr), # NP + R = GP
            #Temp.C = mean(Temp.C)
  ) %>% 
  mutate(P_R="GP") %>% # Label for Gross Photosynthesis
  mutate(light_dark = "LIGHT") %>% 
  mutate(mmol.cm2.hr = ifelse(mmol.cm2.hr < 0, 0, mmol.cm2.hr), # for any values below 0, make 0
         mmol.cm2.hr_uncorr = ifelse(mmol.cm2.hr_uncorr < 0, 0, mmol.cm2.hr_uncorr))

# rejoin for full df with NP, R, and GP rates
RespoR_Normalized_Full <- RespoR_Normalized2 %>% 
  dplyr::select(colony_number, SGD_dil, pH, SGD_number, site, salinity, light_dark, P_R, mmol.cm2.hr, mmol.cm2.hr_uncorr) %>% 
  full_join(RespoR_NormalizedGP)


write_csv(RespoR_Normalized_Full , here("Data","RespoFiles","Respo_RNormalized_AllRates.csv"))  

############################################################################
############################################################################
############################################################################
####### END OF CALCULATIONS -- NOW VISUALIZE AND PLOT ########## 
## notes from Nyssa: 
## plot raw blank - dilutions on X and blank rate on y for light and dark 
## should look really close to each other day to day 
## 5 points on each line -- variance within each one is low 
## if wonky, average across dilutions and them remove date block and join by averaged dils 

## when plotting 
### plots: log(diltuion) 
## transform scale so on natiral sclae -- so log transform the axes (trans (x="log"))
## make set of plots as diltuon on x and y is set of plots with salinity, nutrients, pH 
## correlation coefficient with pH TA nutr and salinity on both axes -- correlated with each other 
## jamie and Dani B correlogram -- git hub 
## prpobably log sclae both of them 
## relationship between the biological parameters 

########## PLOT IT! ###############
###############################################

### new plots for meeting on 08/02/2023 ### 

### blank rates ### 
## NOT AVERAGED ## 
Blank_Rates_R <- RespoR_Normalized_Blanks %>% 
  ggplot(aes(x=SGD_number, 
             y=blank.rate)) +
  scale_x_continuous(trans="log10") +
  geom_point() + 
  facet_wrap(site~light_dark)

ggsave(here("Outputs", "RespoOutput","BlankRates_NotAveraged.jpg"))


### blank rates AVERAGED ACROSS DAYS ### 

RespoR_Normalized_AvgBlanks <- RespoR2 %>% 
  filter(colony_number== "BLANK") %>%
  group_by(SGD_number, light_dark, site) %>%
  summarise(umol.sec = mean(umol.sec, na.rm=TRUE)) %>% # get mean value of blanks per each dilution 
  dplyr::select(blank.rate = umol.sec, site) #%>%  ## rename the blank column 
#ungroup() %>%
# dplyr::select(blank.rate, light_dark, colony_number, SGD_dil, site, SGD_number) #keep only what we need


Blank_Rates_Averaged <- RespoR_Normalized_AvgBlanks %>% 
  ggplot(aes(x=SGD_number, 
             y=blank.rate)) +
  scale_x_continuous(trans="log10") +
  geom_point() + 
  facet_wrap(site~light_dark)


#################################
############# New plots - for 8/2/23 ##########
###################################
RespoR_Normalized_Full <- read_csv(here("Data","RespoFiles","Respo_RNormalized_AllRates.csv"))

my_pal <- pnw_palette(name="Starfish",n=2,type="discrete")


#### R with dilutions and u.mol 
RatesPlot_R <- RespoR_Normalized_Full %>%
  filter(P_R == "R") %>%
  ggplot(aes(x=SGD_number, 
             y=mmol.cm2.hr, 
             color=colony_number)) +
  facet_wrap(colony_number~site, scales = "free_y") +
  geom_point() +
  scale_x_continuous(trans="log10") +
  labs(title="Respiration")
ggsave(here("Outputs", "RespoOutput","Respiration_plot_08022023.jpg"), 
       width = 10, height = 10)


### GP with dilutions and u.mol 
RatesPlot_GP <- RespoR_Normalized_Full %>%
  filter(P_R == "GP") %>%
  ggplot(aes(x=SGD_number, 
             y=mmol.cm2.hr, 
             color=colony_number)) +
  facet_wrap(colony_number~site, scales = "free_y") +
  geom_point() +
  scale_x_continuous(trans="log10") +
  labs(title="GP")
ggsave(here("Outputs", "RespoOutput","GP_plot_08022023.jpg"), 
       width = 10, height = 10)




#### R and GP with salinity and other env parameters ######
RatesPlot_GP_EnvParameters <- RespoR_Normalized_Full %>%
  filter(P_R == "GP") %>%
  ggplot(aes(x=salinity, 
             y=mmol.gram.hr, 
             color=colony_number)) +
  facet_wrap(colony_number~site, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm")+
  labs(title="Respiration with Environmental Parameters")



RatesPlot_R_EnvParameters <- RespoR_Normalized_Full %>%
  filter(P_R == "R") %>%
  ggplot(aes(x=salinity, 
             y=mmol.gram.hr, 
             color=colony_number)) +
  geom_smooth(method = "lm")+
  facet_wrap(colony_number~site, scales = "free") +
  geom_point() +
  labs(title="Respiration with Environmental Parameters")

RatesPlot_R_EnvParameters <- RespoR_Normalized_Full %>%
  filter(P_R == "GP") %>%
  ggplot(aes(x=SGD_number, 
             y=salinity, 
             color=colony_number)) +
  scale_x_continuous(trans="log10") +
  facet_wrap(colony_number~site, scales = "free") +
  geom_point() +
  geom_line()+
  labs(title="Respiration with Environmental Parameters")


RatesPlot_R_EnvParameters <- RespoR_Normalized_Full %>%
  filter(P_R == "GP") %>%
  ggplot(aes(x=SGD_number, 
             y=pH, 
             color=colony_number)) +
  scale_x_continuous(trans="log10") +
  facet_wrap(colony_number~site, scale = "free") +
  geom_line() +
  geom_point()
labs(title="Respiration with Environmental Parameters")




###################################################
######### correlation plots between environmental parameters #########
######################################################

## data 
corr_matrix_envparams <- BioData[-c(1:6,13:17)] 

corr_matrix_envparams2 <- RespoR_Normalized_Full %>% 
  select(salinity, new_pH)

##visualize it 
ggcorr(corr_matrix_envparams2, method = c("everything", "pearson")) 

###############################
### Dani B's methods #######
cor_mat <- round(cor(corr_matrix_envparams2),4)
head(cor_mat)

#melt the correlation matrix means it reassembles data frame to be more effective to complete corr matrix
#to long format
melted_cormat <- melt(cor_mat)
head(melted_cormat)

#visulaize the correlation matrix in general
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Get lower and upper triangle of the correlation matrix
#Note that, a correlation matrix has redundant information. We’ll use the functions below to set half of it to NA
get_lower_tri<-function(cor_mat){
  cormat[upper.tri(cor_mat)] <- NA
  return(cor_mat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cor_mat){
  cor_mat[lower.tri(cor_mat)]<- NA
  return(cor_mat)
}

#apply upper tri calculation to graphc
upper_tri <- get_upper_tri(cor_mat)
upper_tri

#melt the correlation matrix
#melt the correlation data and drop the rows with NA values 
melted_cormat <- melt(upper_tri, na.rm = TRUE)

#heatmap of correlation matrix
#negative correlations are in purple color and positive correlations in red
#scale_fill_gradient2 is used with the argument limit = c(-1,1) as correlation coefficients range from -1 to 1
#coord_fixed() : this function ensures that one unit on the x-axis is the same length as one unit on the y-axis
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "midnightblue", high = "firebrick4", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

melted_cormat <- melted_cormat %>% mutate_at(vars(starts_with("value")), funs(round(., 2)))

# Create a ggheatmap with basic characteristics, etc. 
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "firebrick3", mid = "white", high = "dodgerblue3", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


# Print the heatmap
print(ggheatmap)

#add correlation coefficients to the heatmap
#geom_text() to add the correlation coefficients on the graph
#guides() to change the position of the legend title
#if else statement in melted data frame to quotes of black and white to adjust text color 
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 6) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 18, face="bold", color="black"),
        axis.text.y = element_text(size = 18, face="bold", color="black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(0.8, 0),
        legend.title = element_text(size = 18, face="bold", color="black"),
        legend.text = element_text(size = 20, face="bold", color="black"),
        legend.position = c(0.48, 0.75),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 12, barheight = 2, 
                               title.position = "top", title.hjust = 0.5, title.vjust = 1.0))

##########################################################################################
##########################################################################################

##### old plots ##### 

RespoR_Normalized_Full <- read_csv(here("Data","RespoFiles","Respo_RNormalized_AllRates.csv"))

my_pal <- pnw_palette(name="Starfish",n=2,type="discrete")

# plot GP, NP and R

### R first ####
RatesPlot_R <- RespoR_Normalized_Full %>%
  filter(P_R == "R") %>%
  #mutate(SGD_dil = as.factor(SGD_dil)) %>%
  #filter(SGD_dil==0.00|SGD_dil==0.01|SGD_dil==0.03|SGD_dil==0.05|SGD_dil==0.1|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.0|SGD_dil==4.0) %>% 
  ggplot(aes(x=SGD_number, 
             y=mmol.gram.hr,
             color = colony_number)) +
  facet_wrap(~site) +
  geom_line() +
  scale_x_continuous(breaks=c(0.0, 0.01, 0.03, 0.05, 0.1, 0.5, 1.0, 2.0, 4.0)) +
  #scale_x_log10(breaks=c(0.0, 0.01, 0.03, 0.05, 0.1, 0.5, 1.0, 2.0, 4.0)) +
  #scale_x_continuous(breaks=c(0, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.5)) +
  # breaks=c(0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)) + 
  #facet_wrap(~ colony_number) +
  labs(title="Respiration")
ggsave(here("Outputs", "RespoOutput","Respiration_Dils1_9.jpg"), 
       width = 10, height = 10)


#### now GP ### 
RatesPlot_GP_Dils1_9 <- RespoR_Normalized_Full %>%
  filter(P_R == "GP") %>%
  #mutate(colony_number = as.factor(colony_number)) %>%
  #filter(SGD_dil==0|SGD_dil==0.03|SGD_dil==0.1|SGD_dil==0.2|SGD_dil==0.3|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.0|SGD_dil==3.5) %>%
  ggplot(aes(x=SGD_dil, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  scale_x_continuous(breaks=c(0.0, 0.01, 0.03, 0.05, 0.1, 0.5, 1.0, 2.0, 4.0)) + 
  #scale_x_log10(#limits=c(0,5),
  #breaks=c(0.00, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)) +
  facet_wrap(~ colony_number) + 
  labs(title="GP") #+ #, scales = "fixed")
#theme_bw()+
#theme(strip.background = element_rect(fill = "white"))+
#geom_smooth() +
#labs(x = "Environmental Treatment (High or Low)",
# color = "Assemblage \n Treatment",
# y = "Rate (mmol O2 gram-1 hr-1)",
#title = "Rate of O2 production or consumption") +
#scale_color_manual(values=my_pal) 
ggsave(here("Outputs", "RespoOutput","GP_Dils1_9.jpg"), 
       width = 10, height = 10)









######## THESE ARE FROM THE TEST TRIALS ###########

## all dilutions on one plot for Respiration 
RatesPlot_R_alldils <- RespoR_Normalized_Full %>%
  filter(P_R == "R") %>%
  mutate(colony_number = as.factor(colony_number)) %>%
  #filter(colony_number=="2") %>%
  #filter(SGD_dil==0.00|SGD_dil==0.01|SGD_dil==0.03|SGD_dil==0.1|SGD_dil==0.2|SGD_dil==0.3|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.0|SGD_dil==3.5) %>% 
  ggplot(aes(x=SGD_dil, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  #scale_x_continuous() +
  scale_x_log10(breaks=c(0.00, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)) +
  #scale_x_continuous(breaks=c(0, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.5)) +
  # scale_x_continuous(limits=c(0,5),
  # breaks=c(0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)) + 
  facet_wrap(~ colony_number) +
  labs(title="All Dilutions- Respiration")
ggsave(here("Outputs", "RespoOutput","AllDilutions_R.jpg"), 
       width = 10, height = 10)
#, scales = "fixed") +
#geom_smooth(method="lm") +
#theme_bw() +
#theme(strip.background = element_rect(fill = "white"))+
#labs(x = "Environmental Treatment (High or Low)",
# color = "Assemblage \n Treatment",
# y = "Rate (mmol O2 gram-1 hr-1)",
# title = "Rate of O2 production or consumption") +
#scale_color_manual(values=my_pal) 


#RatesPlot_R

## to find specific curve for dilutions for actual runs 
RatesPlot_R_proposeddils1 <- RespoR_Normalized_Full %>%
  filter(P_R == "R") %>%
  mutate(colony_number = as.factor(colony_number)) %>%
  filter(SGD_dil==0.00|SGD_dil==0.01|SGD_dil==0.03|SGD_dil==0.1|SGD_dil==0.2|SGD_dil==0.3|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.0|SGD_dil==3.5) %>% 
  ggplot(aes(x=SGD_dil, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  #scale_x_continuous() +
  scale_x_log10(breaks=c(0.0, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.5, 5.0)) +
  #scale_x_continuous(breaks=c(0, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.5)) +
  # scale_x_continuous(limits=c(0,5),
  # breaks=c(0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)) + 
  facet_wrap(~ colony_number) +
  labs(title="Proposed Dils 1- Respiration")
ggsave(here("Outputs", "RespoOutput","ProposedRDils1.jpg"), 
       width = 10, height = 10)

## second proposed dilutions 
RatesPlot_R_proposeddils2 <- RespoR_Normalized_Full %>%
  filter(P_R == "R") %>%
  mutate(colony_number = as.factor(colony_number)) %>%
  filter(SGD_dil==0.00|SGD_dil==0.01|SGD_dil==0.03|SGD_dil==0.05|SGD_dil==0.1|SGD_dil==0.3|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.5) %>% 
  ggplot(aes(x=SGD_dil, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  #scale_x_continuous() +
  scale_x_log10(breaks=c(0.0, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1.0, 2.5)) +
  #scale_x_continuous(breaks=c(0, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.5)) +
  # scale_x_continuous(limits=c(0,5),
  # breaks=c(0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)) + 
  facet_wrap(~ colony_number) +
  labs(title="Proposed Dils 2- Respiration")
ggsave(here("Outputs", "RespoOutput","ProposedRDils2.jpg"), 
       width = 10, height = 10)

##
##
## third proposed dilutions 
RatesPlot_R_proposeddils3 <- RespoR_Normalized_Full %>%
  filter(P_R == "R") %>%
  mutate(colony_number = as.factor(colony_number)) %>%
  filter(SGD_dil==0.00|SGD_dil==0.01|SGD_dil==0.03|SGD_dil==0.05|SGD_dil==0.1|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.0|SGD_dil==4.0) %>% 
  ggplot(aes(x=SGD_dil, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  #scale_x_continuous() +
  scale_x_log10(breaks=c(0.0, 0.01, 0.03, 0.05, 0.1, 0.5, 1.0, 2.0, 4.0)) +
  #scale_x_continuous(breaks=c(0, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.5)) +
  # scale_x_continuous(limits=c(0,5),
  # breaks=c(0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)) + 
  facet_wrap(~ colony_number) +
  labs(title="Proposed Dils 3- Respiration")
ggsave(here("Outputs", "RespoOutput","ProposedRDils3.jpg"), 
       width = 10, height = 10)



############### Now do GP ############ 

## GP (NP-R) 
RatesPlot_GP_alldils <- RespoR_Normalized_Full %>%
  filter(P_R == "GP") %>%
  mutate(colony_number = as.factor(colony_number)) %>%
  #filter(SGD_dil==0|SGD_dil==0.03|SGD_dil==0.1|SGD_dil==0.2|SGD_dil==0.3|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.0|SGD_dil==3.5) %>%
  ggplot(aes(x=SGD_dil, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  #scale_x_continuous() + 
  scale_x_log10(#limits=c(0,5),
    breaks=c(0.00, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)) +
  facet_wrap(~ colony_number) + 
  labs(title="All Dilutions - GP") + #, scales = "fixed")
  #theme_bw()+
  #theme(strip.background = element_rect(fill = "white"))+
  #geom_smooth() +
  #labs(x = "Environmental Treatment (High or Low)",
  # color = "Assemblage \n Treatment",
  # y = "Rate (mmol O2 gram-1 hr-1)",
  #title = "Rate of O2 production or consumption") +
  #scale_color_manual(values=my_pal) 
  ggsave(here("Outputs", "RespoOutput","AllGPDils.jpg"), 
         width = 10, height = 10)

##
##

##proposed dilution series 1 
RatesPlot_GP_proposeddils1 <- RespoR_Normalized_Full %>%
  filter(P_R == "GP") %>%
  mutate(colony_number = as.factor(colony_number)) %>%
  filter(SGD_dil==0|SGD_dil==0.03|SGD_dil==0.1|SGD_dil==0.2|SGD_dil==0.3|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.0|SGD_dil==3.5) %>%
  ggplot(aes(x=SGD_dil, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  #scale_x_continuous() + 
  scale_x_log10(#limits=c(0,5),
    breaks=c(0.0, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.5, 5.0)) +
  facet_wrap(~ colony_number) + 
  labs(title="Proposed Dilutions 1- GP")
ggsave(here("Outputs", "RespoOutput","ProposedGPDils1.jpg"), 
       width = 10, height = 10)

##
##

## proposed dilution series 2 
RatesPlot_GP_proposeddils2 <- RespoR_Normalized_Full %>%
  filter(P_R == "GP") %>%
  mutate(colony_number = as.factor(colony_number)) %>%
  filter(SGD_dil==0.00|SGD_dil==0.01|SGD_dil==0.03|SGD_dil==0.05|SGD_dil==0.1|SGD_dil==0.3|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.5) %>%
  ggplot(aes(x=SGD_dil, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  #scale_x_continuous() + 
  scale_x_log10(#limits=c(0,5),
    breaks=c(0.0, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1.0, 2.5)) +
  facet_wrap(~ colony_number)#, scales = "fixed")
ggsave(here("Outputs", "RespoOutput","ProposedGPDils2.jpg"), 
       width = 10, height = 10)

##
##

### proposed dilution series 3 
RatesPlot_GP_proposeddils3 <- RespoR_Normalized_Full %>%
  filter(P_R == "GP") %>%
  mutate(colony_number = as.factor(colony_number)) %>%
  filter(SGD_dil==0.00|SGD_dil==0.01|SGD_dil==0.03|SGD_dil==0.05|SGD_dil==0.1|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.0|SGD_dil==4.0) %>%
  ggplot(aes(x=SGD_dil, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  #scale_x_continuous() + 
  scale_x_log10(#limits=c(0,5),
    breaks=c(0.0, 0.01, 0.03, 0.05, 0.1, 0.5, 1.0, 2.0, 4.0)) +
  facet_wrap(~ colony_number) + 
  labs(title= "Proposed GP Dilutions 3")#, scales = "fixed")
ggsave(here("Outputs", "RespoOutput","ProposedGPDils3.jpg"), 
       width = 10, height = 10)


#ggsave(here("Outputs", "RespoOutput","AllRates.pdf"), RatesPlot, device = "pdf", width = 10, height = 10)

#RatesPlot_GP + RatesPlot_R

####### now make plots with PARAMETERS #### 
###############################################
# parameters for dilutions = TA, fDOM, nutrients, salinity, pH, DO, and temp 

RatesPlot_R_alldils_salinity <- RespoR_Normalized_Full %>%
  filter(P_R == "R") %>%
  mutate(colony_number = as.factor(colony_number)) %>%
  #filter(SGD_dil==0.00|SGD_dil==0.01|SGD_dil==0.03|SGD_dil==0.1|SGD_dil==0.2|SGD_dil==0.3|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.0|SGD_dil==3.5) %>% 
  ggplot(aes(x=salinity, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  scale_x_continuous(breaks=c(39.9, 39.8, 39.7, 39.6, 39.4, 39.2, 38.9, 38.8, 38.0, 37.7, 37.2, 36.4)) +
  #scale_x_log10(breaks=c(0.00, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)) +
  #scale_x_continuous(breaks=c(0, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.5)) +
  # scale_x_continuous(limits=c(0,5),
  # breaks=c(0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)) + 
  facet_wrap(~ colony_number) +
  labs(title="All Dilutions- salinity")
ggsave(here("Outputs", "RespoOutput","AllDilutionsR_salinity.jpg"), 
       width = 10, height = 10)

RatesPlot_GP_alldils_salinity <- RespoR_Normalized_Full %>%
  filter(P_R == "GP") %>%
  mutate(colony_number = as.factor(colony_number)) %>%
  #filter(SGD_dil==0|SGD_dil==0.03|SGD_dil==0.1|SGD_dil==0.2|SGD_dil==0.3|SGD_dil==0.5|SGD_dil==1.0|SGD_dil==2.0|SGD_dil==3.5) %>%
  ggplot(aes(x=salinity, 
             y=mmol.gram.hr,
             color = colony_number)) +
  geom_line() +
  #scale_x_continuous() + 
  scale_x_log10(#limits=c(0,5),
    breaks=c(39.9, 39.8, 39.7, 39.6, 39.4, 39.2, 38.9, 38.8, 38.0, 37.7, 37.2, 36.4)) +
  facet_wrap(~ colony_number) + 
  labs(title="All Dilutions - GP") #+ #, scales = "fixed")
#theme_bw()+
#theme(strip.background = element_rect(fill = "white"))+
#geom_smooth() +
#labs(x = "Environmental Treatment (High or Low)",
# color = "Assemblage \n Treatment",
# y = "Rate (mmol O2 gram-1 hr-1)",
#title = "Rate of O2 production or consumption") +
#scale_color_manual(values=my_pal) 
ggsave(here("Outputs", "RespoOutput","AllGPDils_salinity.jpg"), 
       width = 10, height = 10)





############ if want to run some stats ######## 
###############################################
# quick modeling
# check assumptions for all (referencing Biometry notes below)
# if not seeing growth trends look at my photos and see pigment changes / mortality

library(agricolae) # HSD.test()

# models
GPData <- RespoR_Normalized_Full %>% filter(P_R == "GP")
model1 <- lm(data = GPData, mmol.gram.hr ~ SGD_dil)
anova(model1)
HSD.test(model1, "SGD_dil", console=TRUE)

NPData <- RespoR_Normalized_Full %>% filter(P_R == "NP")
model2 <- lm(data = NPData, mmol.gram.hr ~ SGD_dil)
anova(model2)
HSD.test(model2, "SGD_dil", console=TRUE)

RData <- RespoR_Normalized_Full %>% filter(P_R == "R")
model3 <- lm(data = RData, mmol.gram.hr ~ SGD_dil)
anova(model3)
HSD.test(model3, "SGD-dil", console=TRUE)


#Now let's check our assumptions
plot(model1)
plot(model2)
plot(model3)

# first plot: are variances equal? look at the spread of points of both groups on the different sides of the plot
## anything more than 3x greater are probably unequal and may need to be transformed
## spread is about the same here
# second plot: normality
# third plot: same thing as the first plot, but the square root of the standardized residuals. tests the same thing. are variances equal
# fourth plot: leverage plot.  casey doesn't really use this for anova

# normality was fishy so do a qqp
library(car)
qqp(model1) # just a little off, but it's good enough
qqp(model2) # just a little off, but it's good enough
qqp(model3) # just a little off, but it's good enough

#I'm not 100% sure about the normal probability plot. Let's try it with confidence intervals
library(car)
resid1<-residuals(model1)
qqp(resid1, "norm") # numbers indicate outliers

resid2<-residuals(model2)
qqp(resid2, "norm")

resid3<-residuals(model3)
qqp(resid3, "norm")

## pete also suggested using gaussian model 
## Gaussian = a modeling technique that uses a probability distribution to 
# estimate the likelihood of a given point in a continuous set.



#If you wanted to do a post-hoc test
# library(emmeans)
# emmeans(model1, pairwise~Temp*Genotype, adjust="tukey")

#I wanted to use agricolae, but had to futz with this a little to make this happen
#Tukey tests
# tx<-with(mydata, interaction(Temp, Genotype))
# rmod<-lm(r~tx, data=mydata)
# library(agricolae)
# HSD.test(rmod, "tx", console=TRUE)