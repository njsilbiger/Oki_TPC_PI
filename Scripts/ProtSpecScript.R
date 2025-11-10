########################################################################
### This calculates the soluble protein from a 96 well plate in corals from the SpectraMax iD5
### 
### Created by Maya Powell
### modified Dr. Nyssa Silbiger's chla script and Dr. Hollie Putnam's protein script
###
### Edited on 08/08/2025
#########################################################################

## create a folder for your day of sampling (foldername below) and put your plate.csv files in there. All your data will be exported to that folder

#libraries
library(here)
library(tidyverse)
library(broom)

## File names -------------------
#EDIT FOR EACH RUN!!
foldername<-'ProtSpec' # folder of the day
filename<-'Plate4_Proteins_08.12.25.csv' # data
sampleID<-'protein_platemap_run4.csv' # template of sample IDs
platename<-'plate4' # this will be the name of your file
metadata_file<-"protein_meta_run4.csv" # file with slurry vols and SA

# Only edit the run# for the blanks below this line!! ---------------------
#I just did case-sensitive find-replace for e.g. "Run4" -> "Run4"

# Read in the Sample IDs
sampleIDNames<-read_csv(here('Data',foldername,sampleID))
sampleIDNames<-sampleIDNames %>%
  select(Well.Location, Sample.Type, Sample.Name) # pull out the location and the IDs

# read in metadata
metadata<-read_csv(here("Data", foldername, metadata_file))

### format the 96 well plate data ################
#read in the rows w/o the dye
#protPlate<-read.csv(here('Data',foldername,filename), row.names = c("A","B","C","D","E","F","G","H"),nrows = 8, fileEncoding="latin1", skip = 2)
protPlate<-read.csv(here('Data',foldername,filename))
plate_start_row <- which(protPlate[[2]] == "A")[1]
protPlate <- protPlate[(plate_start_row):(plate_start_row + 7), 3:14]
protPlate[] <- lapply(protPlate, as.numeric) #was getting an issue with some of columns that were characters instead of numeric 
#adjust column names to fit plate format
colnames(protPlate) <- as.character(1:12)
# make the rownames a column
protPlate$Rows<-as.character(c("A","B","C","D","E","F","G","H"))

#make into long format data
protPlate <- protPlate %>%
  select(Rows,1:12) %>%
  pivot_longer(cols = -Rows, values_to = "abs", names_to = "Column")

#add a 0 in front of single digits to keep the well name similar to how the instrument exports
AllData<-protPlate %>%
  mutate(Column = str_pad(as.character(Column), width = 2, pad = "0"))%>% # only add a 0 in front of 1-9
  mutate(Well.Location = paste0(Rows,Column)) # paste with row name

#merge with sample data
AllData <- AllData %>%
  left_join(sampleIDNames) %>% # join with the sampleIDs
  drop_na(Sample.Name) %>%
  left_join(metadata) %>%
  filter(Sample.Type != "Symbiont") # remove any symbiont ones - just looking at host (only issue for plate 1)

# Plot standard curve
# Create standard curve following kit instructions
#YOU WILL HAVE TO CHANGE THIS FOR EVERY RUN SORRY!!!
standards <- tribble(
  ~std, ~BSA_ug.mL,
  "StandardA_Run4",        2000,
  "StandardB_Run4",        1500,
  "StandardC_Run4",        1000,
  "StandardD_Run4",         750,
  "StandardE_Run4",         500,
  "StandardF_Run4",         250,
  "StandardG_Run4",         125,
  "StandardH_Run4",          25,
  "StandardI_Run4",           0
)


std_curve <- AllData %>%
  #unnest(merged) %>%
  filter(Sample.Type == "Standard") %>%
  select(Run, Well.Location, Sample.Name, abs) %>%
  rename(std = Sample.Name) %>%
  #mutate(std = str_sub(std, 9, 9)) %>%
  group_by(std) %>%
  summarise(mean_abs = mean(abs)) %>%                       # calculate mean of standard duplicates
  mutate(mean_abs_adj = mean_abs - mean_abs[std == "StandardI_Run4"]) %>%       # subtract blank absorbance value from all
  #YOU WILL NEED TO CHANGE THE STANDARD # EACH ROUND WHEN YOU RUN THIS SORRY
  left_join(standards)

## Fit linear model for standard curve
# mod <- lm(BSA_ug.mL ~ abs562, data = std_curve)
# coef(mod)

## Fit nonlinear model for standard curve
#mod <- nls(formula = BSA_ug.mL ~ z + a * exp(b * abs562), start = list(z = 0, a = 1, b = 1), data = std_curve)
#fitted <- mod %>% broom::augment()

mod <- lm <- lm(BSA_ug.mL ~ mean_abs_adj, data = std_curve)

fitted <- mod %>% broom::augment()

# Plot standard curve
std_curve_plot <- std_curve %>%
  ggplot(aes(x = mean_abs_adj, y = BSA_ug.mL)) +
  geom_point(color = "red", size = 3) 

std_curve_plot + 
  geom_line(data = fitted, aes(x = mean_abs_adj, y = .fitted)) +
  labs(title = "Standard curve")

# Calculate protein concentration for all samples using standard curve
prot <- AllData %>%
  #unnest(merged) %>%
  select(Run, Well.Location, Sample.Name, abs, Slurry_vol_ml, SA_cm2) %>%        # Select only needed columns
  filter(!is.na(Sample.Name)) %>% 
  group_by(Sample.Name) %>%
  summarise(mean_abs = mean(abs)) %>% 
  #CHANGE BASED ON RUN NUMBER
  mutate(mean_abs_adj = mean_abs - mean_abs[Sample.Name == "StandardI_Run4"]) %>% #subtract abs of blank from ALL samples
  filter(!grepl("Standard", Sample.Name)) %>%   # Get just samples (not standards)
  mutate(prot_ug_mL = map_dbl(mean_abs_adj, ~ predict(mod, newdata = data.frame(mean_abs_adj = .)))) %>%   # Use standard curve to convert absorbance to protein
  left_join(metadata) %>% #merge with sample data again to get slurry vol and SA
  mutate(prot_ug_total = (prot_ug_mL*Slurry_vol_ml)) %>% # normalize to coral slurry
  mutate(prot_ug_cm2 = (prot_ug_total/SA_cm2))

#visually check that all points fall within the standard curve
std_curve_plot + 
  geom_point(data = prot, aes(x = mean_abs_adj, y = prot_ug_mL), pch = "X", cex = 5, alpha = 0.3) +
  labs(title = "All samples projected on standard curve")

### Export the data
write_csv(prot, here("Data", foldername, paste0(platename,"_summary.csv"))) # print the summary data
#yay done! now do it again!

## Combine all data!! ---------------------
p1 <- read.csv(here('Data', 'ProtSpec',"plate1_summary.csv"))
p2 <- read.csv(here('Data', 'ProtSpec',"plate2_summary.csv"))
p3 <- read.csv(here('Data', 'ProtSpec',"plate3_summary.csv"))
p4 <- read.csv(here('Data', 'ProtSpec',"plate4_summary.csv"))
p5 <- read.csv(here('Data', 'ProtSpec',"plate5_summary.csv"))
p6 <- read.csv(here('Data', 'ProtSpec',"plate6_summary.csv"))

prot_all_summ <- rbind(p1,p2,p3,p4,p5,p6)

write_csv(prot_all_summ, here("Data", 'ProtSpec', "protein_all_summary.csv")) # print the summary data


####Hollie's script below to use as needed##
#using Hollie's script for protein analysis from Mo'orea Pocillopora repo

# # Import data
# 
# # Define function to read in protein data
# read_pro <- function(file) {
#   prot_data <- read_csv(file, skip = 31, n_max = 31) %>%
#     magrittr::set_colnames(c("row", 1:12, "wavelength")) %>%
#     fill(row) %>%
#     gather("col", "absorbance", -wavelength, -row) %>%
#     unite("well", c(row, col), sep = "")
# }
# 
# # List protein data files
# prot_path = here('Data/test_prot/')                                        # Path to prot data directory
# all_prot_files <- list.files(path = prot_path, pattern = "*.csv")          # List all files in directory
# prot_platemaps <- list.files(path = prot_path, pattern = "platemap")       # List platemap files
# prot_data_files <- setdiff(all_prot_files, prot_platemaps)                 # List data files
# 
# # Read in all files into tibble
# df <- tibble(file = prot_data_files) %>%
#   mutate(platemap = map(file, ~ read_csv(paste0(prot_path, tools::file_path_sans_ext(.), "_platemap.csv"))),
#          prot_data = map(file, ~ read_pro(paste0(prot_path, .))))
# 
# # Merge platemap and data for each plate
# df <- df %>%
#   mutate(merged = map2(platemap, prot_data, ~ right_join(.x, .y)))
# 
# # Plot standard curve
# # Create standard curve following kit instructions
# standards <- tribble(
#   ~std, ~BSA_ug.mL,
#   "A",        2000,
#   "B",        1500,
#   "C",        1000,
#   "D",         750,
#   "E",         500,
#   "F",         250,
#   "G",         125,
#   "H",          25,
#   "I",           0
# )
# 
# std_curve <- df %>%
#   unnest(merged) %>%
#   filter(grepl("Standard", full_sample_id)) %>%
#   select(plate, well, full_sample_id, abs562 = `562:562`) %>%
#   rename(std = full_sample_id) %>%
#   mutate(std = str_sub(std, 9, 9)) %>%
#   #group_by(std) %>%
#   #summarise(abs562 = mean(abs562)) %>%                       # calculate mean of standard duplicates
#   #mutate(abs562.adj = abs562 - abs562[std == "I"]) %>%       # subtract blank absorbace value from all
#   left_join(standards)
# 
# ## Fit linear model for standard curve
# # mod <- lm(BSA_ug.mL ~ abs562, data = std_curve)
# # coef(mod)
# 
# ## Fit nonlinear model for standard curve
# #mod <- nls(formula = BSA_ug.mL ~ z + a * exp(b * abs562), start = list(z = 0, a = 1, b = 1), data = std_curve)
# #fitted <- mod %>% broom::augment()
# 
# mod <- lm <- lm(BSA_ug.mL ~ abs562, data = std_curve)
# 
# fitted <- mod %>% broom::augment()
# 
# # Plot standard curve
# std_curve_plot <- std_curve %>%
#   ggplot(aes(x = abs562, y = BSA_ug.mL)) +
#   geom_point(color = "red", size = 3) 
# 
# std_curve_plot + 
#   geom_line(data = fitted, aes(x = abs562, y = .fitted)) +
#   labs(title = "Standard curve")
# 
# # Calculate protein concentrations
# # Calculate protein concentration for all samples using standard curve
# prot <- df %>%
#   unnest(merged) %>%
#   filter(!grepl("Standard", full_sample_id)) %>%                     # Get just samples (not standards)
#   select(plate, well, full_sample_id, abs562 = `562:562`) %>%        # Select only needed columns
#   filter(!is.na(full_sample_id)) %>%                                 # Filter out empty wells
#   filter(full_sample_id != "BK") %>%                                 # Filter out blank wells
#   mutate(prot_ug.mL = map_dbl(abs562, ~ predict(mod, newdata = data.frame(abs562 = .))))    # Use standard curve to convert absorbance to protein
# 
# std_curve_plot + 
#   geom_point(data = prot, aes(x = abs562, y = prot_ug.mL), pch = "X", cex = 5, alpha = 0.3) +
#   labs(title = "All samples projected on standard curve")
# 
# prot <- prot %>%
#   group_by(full_sample_id) %>%
#   summarise(n = n(), prot_ug.mL = mean(abs562))
# 
# # Normalize to surface area
# # Surface area data
# sa <- read.csv("output/surface.area.calc.csv")
# # Tissue homogenate volume data
# homog_vols <- read_csv("data/homogenate_vols.csv") 
# 
# # Coral sample metadata
# metadata <- read_csv("data/corals_sampled.csv") 
# 
# # Join homogenate volumes and surface area with sample metadata
# metadata <- left_join(metadata, homog_vols, by="full_sample_id") %>%
#   left_join(sa, by="full_sample_id")
# 
# # Join prot data with metadata
# prot <- left_join(prot, metadata, by="full_sample_id") %>%
#   mutate(prot_ug = prot_ug.mL * homog_vol_ml,
#          prot_ug.cm2 = prot_ug / surface.area.cm2,
#          prot_mg.cm2 = prot_ug.cm2 / 1000)
# 
# # Plot results by species and site
# monthorder <- c("January", "February","March", "April","May",
#                 "June",  "July", "August", "September", 
#                 "October","November")
# 
# prot$timepoint.x <- factor(prot$timepoint.x,levels=monthorder)
# prot$group <- paste0(prot$timepoint.x, "_", prot$site.x)
# 
# # Plot 
# prot %>%
#   ggplot(aes(x = timepoint.x, y = prot_mg.cm2, group=group, color=site.x)) +
#   geom_boxplot(aes(group=group, color=site.x))+
#   geom_jitter(width = 0.1) + # Plot all points
#   scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#   #stat_summary(fun.data = mean_cl_normal, fun.args = list(mult = 1),    # Plot standard error
#   #             geom = "errorbar", color = c("black"), width = 0.5) +
#   #stat_summary(fun.y = mean, geom = "point", color = "black")+           # Plot mean
#   labs(x = "", y = "Host Soluble Protein (mg/cm2)") +
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# # Write data to output file
# # Write protein data to output file
# prot %>%
#   select(full_sample_id, prot_ug, prot_ug.cm2) %>%
#   mutate(partner="host")%>%
#   write_csv(., path = "output/host_protein.csv")
