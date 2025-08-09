########################################################################
### This calculates the Chlorophyll from a 96 well plate in corals from teh SpectraMax iD5
### 
### Created by Dr. Nyssa Silbiger
### Edited on 07/10/2025
#########################################################################

## create a folder for your day of sampling (foldername below) and put your plate.csv files in there. All your data will be exported to that folder

#libraries
library(here)
library(tidyverse)


## File names -------------------
foldername<-'test_chl' # folder of the day
filename<-'testchl2.csv' # data
sampleID<-'Chl_template_test.csv' # template of sample IDs
platename<-'plate1_test' # this will be the name of your file
metadata_file<-"testchl2_metadata.csv" # file with slurry vols and SA

## What is the path length?
PL<-0.71 # pathlength for donahue lab plate


# DONT CHANGE ANYTHING BELOW HERE ----------------------------------

# Read in the Sample IDs
sampleIDNames<-read_csv(here('Data',foldername,sampleID))
sampleIDNames<-sampleIDNames %>%
  select(Well.Location, Sample.Name) # pull out the location and the IDs

# read in metadata
metadata<-read_csv(here("Data", foldername, metadata_file))

### format the 96 well plate data ################
#read in the rows w/o the dye
ChlPlate<-read.csv(here('Data',foldername,filename), row.names = c("A","B","C","D","E","F","G","H"),nrows = 8, fileEncoding="latin1", skip = 2)
# make the rownames a column
ChlPlate$Rows<-rownames(ChlPlate)

## Pull out each of the wavelengths to make own column
#630 wavelength
ChlPlate_630<-ChlPlate %>%
  select(Rows,X1:X12)%>%
  pivot_longer(cols = X1:X12, values_to =  "ChlPlate_630", names_to = "Column")

#663 wavelength
ChlPlate_663<-ChlPlate %>%
  select(Rows,X1.1:X12.1)%>% # pull out the 663 columns
  pivot_longer(cols = X1.1:X12.1, values_to =  "ChlPlate_663", names_to = "Column") %>%
  separate(col = Column, extra = "drop", into = 'Column') # remove the .1

#750 wavelength
ChlPlate_750<-ChlPlate %>%
  select(Rows,X1.2:X12.2)%>% # pull out the 750 columns
  pivot_longer(cols = X1.2:X12.2, values_to =  "ChlPlate_750", names_to = "Column") %>%
  separate(col = Column, extra = "drop", into = 'Column') # remove the .2


## bring everything to one dataframe (reduce allows me to use the left_join function multiple times at once)
AllData<-Reduce(function(...) left_join(...), list(ChlPlate_630,ChlPlate_663,ChlPlate_750))

# remove the X in the column name and add a 0 in front of single digits to keep the well name similar to how the instrument exports
AllData<-AllData %>%
  mutate(Column = as.numeric(str_remove(AllData$Column, pattern = "X")))%>% # remove the X
  mutate(Column = ifelse(Column<10, gsub("(\\d)+", "0\\1", Column),Column))%>% # only add a 0 in front of 1-9
  mutate(Well.Location = paste0(Rows,Column)) # paste with row name

### Run pH Analysis ##################
# chl from Jeffry and Humphreys
ChlData_raw<-AllData %>%
  mutate(chla = 11.43*(ChlPlate_663 - ChlPlate_750 / PL) - 0.64*(ChlPlate_630 - ChlPlate_750/PL),
         chlc = 27.09*(ChlPlate_630 - ChlPlate_750 / PL) - 3.63*(ChlPlate_663 - ChlPlate_750/PL),
         Totalchl = chla+chlc,
           daterun = Sys.Date()) %>%
  left_join(sampleIDNames) %>% # join with the sampleIDs
  drop_na(Sample.Name) %>%
  left_join(metadata) %>%
  mutate(chla_ug_cm2 = (chla*Slurry_vol_ml)/SA_cm2,  # normalize to coral slurry and SA
         chlc_ug_cm2 = (chlc*Slurry_vol_ml)/SA_cm2,
         Totalchl_ug_cm2 = (Totalchl*Slurry_vol_ml)/SA_cm2)
         

# averaged data by sample ID
ChlData_ave<-ChlData_raw %>%
  group_by(Sample.Name) %>%
  summarise_at(vars(chla:Totalchl, chla_ug_cm2:Totalchl_ug_cm2), .funs = list(function(x){mean(x, na.rm = TRUE)},
                                                 function(x){sd(x, na.rm = TRUE)}/sqrt(n()))) %>%
  rename(chla_mean = chla_fn1, chlc_mean=chlc_fn1, TotalChl_mean = Totalchl_fn1,
         chla_SE = chla_fn2, chlc_SE = chlc_fn2, Totalchl_SE = Totalchl_fn2,
         chla_ug_cm2_mean = chla_ug_cm2_fn1, chlc_ug_cm2_mean = chlc_ug_cm2_fn1,
         Totalchl_ug_cm2_mean = Totalchl_ug_cm2_fn1, 
         chla_ug_cm2_SE = chla_ug_cm2_fn2, chlc_ug_cm2_SE = chlc_ug_cm2_fn2,
         Totalchl_ug_cm2_SE = Totalchl_ug_cm2_fn2) %>%
  select(Sample.Name,chla_mean,chla_SE,chlc_mean,chlc_SE,TotalChl_mean,Totalchl_SE,
         chla_ug_cm2_mean,chla_ug_cm2_SE, chlc_ug_cm2_mean,chlc_ug_cm2_SE,Totalchl_ug_cm2_mean,Totalchl_ug_cm2_SE ) # put in a better order
  


### Export the data
## all the info
write_csv(ChlData_ave, here("Data", foldername, paste0(platename,"_summary.csv"))) # print the summary data
write_csv(ChlData_raw, here("Data", foldername, paste0(platename,"_raw.csv"))) # print the raw data

