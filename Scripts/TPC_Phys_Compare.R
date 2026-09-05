###### Code for TPC and Physiology comparisons ####### 
### Created by: Maya Powell
#### Last updated on: Oct 2 2025

### Install Packages #####
## if these packages are not yet installed, install them 
## great for updates or new users 
#if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented')

####Read in required libraries#####
##### Include Versions of libraries
#library(lubridate)
library(tidyverse)
library(here)
library(PNWColors)
library(viridis)
library(car)
library(dplyr)
library(ggplot2)
library(powerjoin)
library(forcats)
library(car)
library(emmeans)
library(performance)
library(purrr)
library(rlang)
library(tibble)
library(here)
library(psych)
library(corrplot)
library(ggpubr)
library(vegan)
#remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

###### Initial Data Read In ########
#TPC data with only the seven species that we have physio data for
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4.csv"))
topt_df <- topt_df %>% filter(sample_ID != "B08_TPC") %>% select(-species,-full_species,-SA_cm2)
#topt_df <- topt_df %>% dplyr::select(rmax, topt, e, PR, frag_ID)
phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta_all.csv"))
topt_df <- topt_df %>% left_join(phys_meta, by = "frag_ID")

#additional setup
sp_cols <- c(
  "Acropora hyacinthus" = '#d8aedd',
  "Echinopora lamellosa" = '#ba7999',
  "Favites complanata"   = '#dd4124',
  "Montipora aequituberculata" = '#ed8b00',
  "Montipora vietnamensis" = '#efbc82',
  "Pachyseris rugosa" = '#edd746',
  "Pocillopora eydouxi" = '#d0e2af',
  "Porites cylindrica" = '#45681e',
  "Porites rus" = '#7bbcd5',
  "Turbinaria frondens" = '#00496f'
)

se_fun <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

###Correlation coefficient plot####
#generated below - only need once
avg_AFDW <- read.csv(here("Data", "Physiology", "Average_Ash_Free_Dry_Weight.csv")) #afdw_mg_cm2
avg_AFDW <- avg_AFDW %>% mutate(afdw_log = log(afdw_mg_cm2)) %>%
  mutate(dw_log = log(dw_mg_cm2)) %>%
  dplyr::select(frag_ID, dw_mg_cm2, dw_log, afdw_mg_cm2,afdw_log)
chla_avg <- read.csv(here("Data", "Physiology", "Chla_avg.csv")) #chla_ug_cm2_mean
chla_avg <- chla_avg %>% mutate(chla_log = log(chla_ug_cm2_mean)) %>%
  mutate(chla_sym_log = log(chla_pg_sym)) %>%
  dplyr::select(frag_ID, chla_ug_cm2_mean,chla_log, chla_pg_sym, chla_sym_log)
avg_sym <- read.csv(here("Data", "Physiology", "Average_Sym_Density.csv"))
avg_sym <- avg_sym %>% filter(frag_ID != "C07") %>% filter(frag_ID != "D10") #remove crazy outliers!
avg_sym <- avg_sym %>% mutate(sym_log = log(sym_cm2)) %>% dplyr::select(frag_ID, sym_cm2,sym_log)
prot <- read.csv(here("Data", "Physiology", "protein_all_summary.csv")) #prot_ug_cm2
prot <- prot %>% mutate(prot_log = log(prot_ug_cm2)) %>% dplyr::select(frag_ID, prot_ug_cm2,prot_log)
physio_list <- list(avg_AFDW,chla_avg,avg_sym,prot)
all_physio <- physio_list %>% reduce(left_join)
write_csv(all_physio, here("Data", "Physiology", "all_physio_data.csv"))
all_physio <- read_csv(here("Data", "Physiology", "all_physio_data.csv"))

#generate full dataframe
result <- Reduce(function(x, y) merge(x, y, all = TRUE), list(topt_df,avg_AFDW,chla_avg,avg_sym,prot))
#drop NA data from respo data because we donʻt have all reps
all_data <- result %>% drop_na(rmax)
write_csv(result, here("Data", "Physiology", "all_data_concatenated.csv"))

##Read in concatenated data (made above)
all_data <- read_csv(here("Data", "Physiology", "all_data_concatenated.csv"))

gp_data <- all_data %>% filter(PR == "GrossPhoto")
np_data <- all_data %>% filter(PR == "NetPhoto")
r_data <- all_data %>% filter(PR == "Respiration")
write_csv(np_data, here("Data", "Physiology", "np_data.csv"))

#use these for correlation plots

##Additional correlations with specific temperature data
respo_select_temps <- read_csv(here("Data", "RespoFiles","TPC", "respo_select_temps.csv"))
respo_select_temps <- respo_select_temps %>% left_join(all_data, by = "frag_ID")
gp_data <- respo_select_temps %>% filter(PR == "GrossPhoto") %>%
  select(frag_ID,rmax,topt,e) %>% rename(gp_rmax = rmax,
                                 gp_topt = topt,
                                 gp_e = e)
np_data <- respo_select_temps %>% filter(PR == "NetPhoto") %>%
  select(-PR) %>% rename(np_rmax = rmax,
                         np_topt = topt,
                         np_e = e)
r_data <- respo_select_temps %>% filter(PR == "Respiration") %>%
  select(frag_ID,rmax,topt,e) %>%
  rename(r_rmax = rmax,
         r_topt = topt,
         r_e = e)
respo_list <- list(np_data, gp_data, r_data)
respo_phys_full <- respo_list %>% reduce(left_join)

write_csv(respo_phys_full, here("Data", "RespoFiles","TPC", "Respo_Physiology_GP_NP_R_AllRates.csv"))
respo_phys_full <- read_csv(here("Data", "RespoFiles","TPC", "Respo_Physiology_GP_NP_R_AllRates.csv"))

###Correlation plot
vars <- respo_phys_full |>
  dplyr::select(
    gp_rmax, gp_topt, gp_e,
    np_rmax, np_topt, np_e,
    r_rmax, r_topt, r_e,
    NP_26, NP_29, NP_32,
    GP_26, GP_29, GP_32,
    R_26, R_29, R_32,
    NPR_26, NPR_29, NPR_32,
    dw_mg_cm2,
    #dw_log,
    afdw_mg_cm2,
    #afdw_log,
    chla_ug_cm2_mean,
    #chla_log,
    chla_pg_sym,
    #chla_sym_log,
    sym_cm2,
    #sym_log,
    prot_ug_cm2)
    #prot_log)

cor_out <- psych::corr.test(
  vars,
  use    = "pairwise",   # handles NAs
  method = "pearson",    # change to "spearman" if you prefer
  adjust = "BH"          # Benjamini-Hochberg FDR correction
)

cor_mat <- cor_out$r      # correlation coefficients
p_mat   <- cor_out$p      # p-values (adjusted if adjust != "none")

cor_long <- cor_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("var1") %>%
  pivot_longer(
    cols      = -var1,
    names_to  = "var2",
    values_to = "r"
  )

p_long <- p_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("var1") %>%
  pivot_longer(
    cols      = -var1,
    names_to  = "var2",
    values_to = "p"
  )

cor_plot_df <- cor_long %>%
  left_join(p_long, by = c("var1", "var2")) %>%
  mutate(
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    label = ifelse(
      sig == "",
      sprintf("%.2f", r),                    # just r
      sprintf("%.2f\n%s", r, sig)            # r on first line, stars on second
    )
  )

var_levels <- colnames(cor_mat)
cor_plot_df <- cor_plot_df %>%
  mutate(
    var1 = factor(var1, levels = var_levels),
    var2 = factor(var2, levels = var_levels)
  )

cor_plot_df <- cor_plot_df %>%
  filter(as.numeric(var1) < as.numeric(var2))

corr_plot <- ggplot(cor_plot_df, aes(x = var1, y = var2, fill = r)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    limit = c(-1, 1),
    name  = "r"
  ) +
  geom_text(aes(label = label), size = 3.5) +
  coord_fixed() +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL,
    y = NULL,
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0),
    panel.grid  = element_blank()
  )
# corrplot(
#   cor_mat,
#   method      = "color",   # colored tiles
#   type        = "upper",   # upper triangle only
#   order       = "hclust",  # cluster similar variables
#   addCoef.col = "black",   # show r values on tiles
#   tl.col      = "black",   # text label color
#   tl.srt      = 45,        # rotate variable labels
#   p.mat       = p_mat,     # matrix of p-values
#   sig.level   = 0.05,      # significance cutoff
#   insig       = "blank"    # hide non-significant correlations
# )

ggsave(here("Output", "Physiology", "Correlations_Physio_GP_NP_R.pdf"), corr_plot, h = 15, w = 15)

####Topt and rmax graphs####
#np graphs

##Topt
#summary
np_topt_summary <- np_data %>%
  group_by(full_species) %>%
  summarise(
    n    = sum(!is.na(topt)),
    mean = mean(topt, na.rm = TRUE),
    se   = se_fun(topt),
    .groups = "drop")
#model
np_topt.mod.spp <- lm(topt~full_species, data = np_data)
Anova(np_topt.mod.spp) #ns
summary(np_topt.mod.spp)
check_model(np_topt.mod.spp)

#order data
means <- np_data %>% group_by(species) %>% summarize(m = mean(topt, na.rm = TRUE), .groups = "drop")
np_topt_ordered <- np_data %>% left_join(np_topt_summary, by = "full_species") %>% mutate(full_species = fct_reorder(full_species, mean))  # ascending by mean

#graph
np_topt_plot <- ggplot() +
  geom_jitter(data = np_topt_ordered, aes(x = full_species, y = topt, color = full_species, shape = life_history), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = np_topt_summary, aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = np_topt_summary, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  #stat_summary(data = np_topt_ordered, aes(x = full_species, y = np_topt), geom = "text", fun = max, vjust = -0.5, size = 8,
  #             label = c("a", "ab", "ab", "ab", "abc", "abc", "bc", "bc", "bc", "c"))+
  #a, ab, ab, ab, abc, abc, bc, bc, bc, c 
  theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #ylim(0,4100)+
  labs(title = "Net Photosynthesis",
       x = "Species", color = "Species",
       y = "Thermal Optimum (°C)")
np_topt_plot
ggsave(here("Output", "Physiology", "np_topt_species_jitter.pdf"), np_topt_plot, h = 5, w = 10)

##np_rmax
#summary
np_rmax_summary <- np_data %>%
  group_by(full_species) %>%
  summarise(
    n    = sum(!is.na(rmax)),
    mean = mean(rmax, na.rm = TRUE),
    se   = se_fun(rmax),
    .groups = "drop")
#model
np_rmax.mod.spp <- lm(rmax~full_species, data = np_data)
Anova(np_rmax.mod.spp) #0.01727 *
summary(np_rmax.mod.spp)
# Residual standard error: 0.2259 on 39 degrees of freedom
# Multiple R-squared:  0.3784,	Adjusted R-squared:  0.235 
# F-statistic: 2.638 on 9 and 39 DF,  p-value: 0.01727
check_model(np_rmax.mod.spp)
#pairwise comparisons
emm_obj <- emmeans::emmeans(np_rmax.mod.spp, ~ full_species)
emm_pairs <- pairs(emm_obj)
# Echinopora lamellosa - Montipora aequituberculata     0.4840 0.143 39   3.388  0.0458
# Favites complanata - Montipora aequituberculata       0.5440 0.143 39   3.808  0.0155
#almost:
# Montipora aequituberculata - Pachyseris rugosa       -0.4720 0.143 39  -3.304  0.0563

#order data
means <- np_data %>% group_by(species) %>% summarize(m = mean(rmax, na.rm = TRUE), .groups = "drop")
np_rmax_ordered <- np_data %>% left_join(np_rmax_summary, by = "full_species") %>% mutate(full_species = fct_reorder(full_species, mean))  # ascending by mean

#graph
np_rmax_plot <- ggplot() +
  geom_jitter(data = np_rmax_ordered, aes(x = full_species, y = rmax, color = full_species, shape = life_history), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = np_rmax_summary, aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = np_rmax_summary, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  stat_summary(data = np_rmax_ordered, aes(x = full_species, y = rmax), geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "b", "b"))+
  #a, ab, ab, ab, abc, abc, bc, bc, bc, c 
  theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #ylim(0,4100)+
  labs(title = "Net Photosynthesis",
       x = "Species", color = "Species",
       y = expression("Rate Maximum" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})))
np_rmax_plot
ggsave(here("Output", "Physiology", "np_rmax_species_jitter.pdf"), np_rmax_plot, h = 5, w = 10)


#gp graphs

##Topt
#summary
gp_topt_summary <- gp_data %>%
  group_by(full_species) %>%
  summarise(
    n    = sum(!is.na(topt)),
    mean = mean(topt, na.rm = TRUE),
    se   = se_fun(topt),
    .groups = "drop")
#model
gp_topt.mod.spp <- lm(topt~full_species, data = gp_data)
Anova(gp_topt.mod.spp) #ns p = 0.1014
summary(gp_topt.mod.spp)
check_model(gp_topt.mod.spp)

#order data
means <- gp_data %>% group_by(species) %>% summarize(m = mean(topt, na.rm = TRUE), .groups = "drop")
gp_topt_ordered <- gp_data %>% left_join(gp_topt_summary, by = "full_species") %>% mutate(full_species = fct_reorder(full_species, mean))  # ascending by mean

#graph
gp_topt_plot <- ggplot() +
  geom_jitter(data = gp_topt_ordered, aes(x = full_species, y = topt, color = full_species, shape = life_history), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = gp_topt_summary, aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = gp_topt_summary, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  #stat_summary(data = gp_topt_ordered, aes(x = full_species, y = gp_topt), geom = "text", fun = max, vjust = -0.5, size = 8,
  #             label = c("a", "ab", "ab", "ab", "abc", "abc", "bc", "bc", "bc", "c"))+
  #a, ab, ab, ab, abc, abc, bc, bc, bc, c 
  theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #ylim(0,4100)+
  labs(title = "Gross Photosynthesis",
       x = "Species", color = "Species",
       y = "Thermal Optimum (°C)")
gp_topt_plot
ggsave(here("Output", "Physiology", "gp_topt_species_jitter.pdf"), gp_topt_plot, h = 5, w = 10)

##gp_rmax
#summary
gp_rmax_summary <- gp_data %>%
  group_by(full_species) %>%
  summarise(
    n    = sum(!is.na(rmax)),
    mean = mean(rmax, na.rm = TRUE),
    se   = se_fun(rmax),
    .groups = "drop")
#model
gp_rmax.mod.spp <- lm(rmax~full_species, data = gp_data)
Anova(gp_rmax.mod.spp) #0.008926 **
summary(gp_rmax.mod.spp)
# Residual standard error: 0.2895 on 39 degrees of freedom
# Multiple R-squared:  0.4055,	Adjusted R-squared:  0.2683 
# F-statistic: 2.956 on 9 and 39 DF,  p-value: 0.008926
check_model(gp_rmax.mod.spp)
#pairwise comparisons
emm_obj <- emmeans::emmeans(gp_rmax.mod.spp, ~ full_species)
emm_pairs <- pairs(emm_obj)
# Montipora aequituberculata - Porites cylindrica      -0.6520 0.183 39  -3.560  0.0297
# Montipora aequituberculata - Porites rus             -0.7020 0.183 39  -3.833  0.0144
# Montipora aequituberculata - Pachyseris rugosa       -0.6160 0.183 39  -3.364  0.0486
# Favites complanata - Montipora aequituberculata       0.7260 0.183 39   3.965  0.0101
# Echinopora lamellosa - Montipora aequituberculata     0.6260 0.183 39   3.418  0.0425

#order data
means <- gp_data %>% group_by(species) %>% summarize(m = mean(rmax, na.rm = TRUE), .groups = "drop")
gp_rmax_ordered <- gp_data %>% left_join(gp_rmax_summary, by = "full_species") %>% mutate(full_species = fct_reorder(full_species, mean))  # ascending by mean

#graph
gp_rmax_plot <- ggplot() +
  geom_jitter(data = gp_rmax_ordered, aes(x = full_species, y = rmax, color = full_species, shape = life_history), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = gp_rmax_summary, aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = gp_rmax_summary, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  stat_summary(data = gp_rmax_ordered, aes(x = full_species, y = rmax), geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "ab", "ab", "ab", "ab", "b", "b", "b", "b", "b"))+
  #a, ab, ab, ab, abc, abc, bc, bc, bc, c 
  theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #ylim(0,4100)+
  labs(title = "Gross Photosynthesis",
       x = "Species", color = "Species",
       y = expression("Rate Maximum" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})))
gp_rmax_plot
ggsave(here("Output", "Physiology", "gp_rmax_species_jitter.pdf"), gp_rmax_plot, h = 5, w = 10)


#r graphs

##Topt
#summary
r_topt_summary <- r_data %>%
  group_by(full_species) %>%
  summarise(
    n    = sum(!is.na(topt)),
    mean = mean(topt, na.rm = TRUE),
    se   = se_fun(topt),
    .groups = "drop")
#model
r_topt.mod.spp <- lm(topt~full_species, data = r_data)
Anova(r_topt.mod.spp) #ns p = 0.4578
summary(r_topt.mod.spp)
check_model(r_topt.mod.spp)

#order data
means <- r_data %>% group_by(species) %>% summarize(m = mean(topt, na.rm = TRUE), .groups = "drop")
r_topt_ordered <- r_data %>% left_join(r_topt_summary, by = "full_species") %>% mutate(full_species = fct_reorder(full_species, mean))  # ascending by mean

#graph
r_topt_plot <- ggplot() +
  geom_jitter(data = r_topt_ordered, aes(x = full_species, y = topt, color = full_species, shape = life_history), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = r_topt_summary, aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = r_topt_summary, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  #stat_summary(data = r_topt_ordered, aes(x = full_species, y = r_topt), geom = "text", fun = max, vjust = -0.5, size = 8,
  #             label = c("a", "ab", "ab", "ab", "abc", "abc", "bc", "bc", "bc", "c"))+
  #a, ab, ab, ab, abc, abc, bc, bc, bc, c 
  theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #ylim(0,4100)+
  labs(title = "Respiration",
       x = "Species", color = "Species",
       y = "Thermal Optimum (°C)")
r_topt_plot
ggsave(here("Output", "Physiology", "r_topt_species_jitter.pdf"), r_topt_plot, h = 5, w = 10)

##r_rmax
#summary
r_rmax_summary <- r_data %>%
  group_by(full_species) %>%
  summarise(
    n    = sum(!is.na(rmax)),
    mean = mean(rmax, na.rm = TRUE),
    se   = se_fun(rmax),
    .groups = "drop")
#model
r_rmax.mod.spp <- lm(rmax~full_species, data = r_data)
Anova(r_rmax.mod.spp) #ns 0.7863
summary(r_rmax.mod.spp)
check_model(r_rmax.mod.spp) #bad

#order data
means <- r_data %>% group_by(species) %>% summarize(m = mean(rmax, na.rm = TRUE), .groups = "drop")
r_rmax_ordered <- r_data %>% left_join(r_rmax_summary, by = "full_species") %>% mutate(full_species = fct_reorder(full_species, mean))  # ascending by mean

#graph
r_rmax_plot <- ggplot() +
  geom_jitter(data = r_rmax_ordered, aes(x = full_species, y = rmax, color = full_species, shape = life_history), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = r_rmax_summary, aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = r_rmax_summary, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  #stat_summary(data = r_rmax_ordered, aes(x = full_species, y = rmax), geom = "text", fun = max, vjust = -0.5, size = 8,
  #             label = c("a", "ab", "ab", "ab", "ab", "b", "b", "b", "b", "b"))+
  #a, ab, ab, ab, abc, abc, bc, bc, bc, c 
  theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #ylim(0,4100)+
  labs(title = "Gross Photosynthesis",
       x = "Species", color = "Species",
       y = expression("Rate Maximum" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})))
r_rmax_plot
ggsave(here("Output", "Physiology", "r_rmax_species_jitter.pdf"), r_rmax_plot, h = 5, w = 10)

####Significant topt and rmax relationship plots based on corr plots####
#np rmax and chla pg sym
#np rmax and chla ug cm2
#np topt and chla pg sym
#np topt and chla ug cm2
#gp rmax and chla pg sym
#gp rmax and chla ug cm2
#r topt and prot ug cm2

#np rmax and chla pg sym
np_rmax_chla_sym_scatter <- ggplot(np_data) +
  geom_point(aes(y = rmax, x = chla_pg_sym, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~full_species, scales = "free")+
  geom_smooth(aes(y = rmax, x = chla_pg_sym, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Chlorophyll a" ~ (pg ~ symbiont^{-1})), 
       y = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
       color = "Species")
np_rmax_chla_sym_scatter
ggsave(here("Output", "Physiology", "np_rmax_chla_sym_scatter.pdf"), np_rmax_chla_sym_scatter, h = 5, w = 10)

#np rmax and chla ug cm
np_rmax_chla_scatter <- ggplot(np_data) +
  geom_point(aes(y = rmax, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
       color = "Species")
np_rmax_chla_scatter
ggsave(here("Output", "Physiology", "np_rmax_chla_scatter.pdf"), np_rmax_chla_scatter, h = 5, w = 10)

#np topt and chla pg sym
np_topt_chla_sym_scatter <- ggplot(np_data) +
  geom_point(aes(y = topt, x = chla_pg_sym, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~full_species, scales = "free")+
  geom_smooth(aes(y = topt, x = chla_pg_sym, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Chlorophyll a" ~ (pg ~ symbiont^{-1})), 
       y = "Thermal optimum (°C)",
       color = "Species")
np_topt_chla_sym_scatter
ggsave(here("Output", "Physiology", "np_topt_chla_sym_scatter.pdf"), np_topt_chla_sym_scatter, h = 5, w = 10)

#np topt and chla ug cm
np_topt_chla_scatter <- ggplot(np_data) +
  geom_point(aes(y = topt, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = "Thermal optimum (°C)",
       color = "Species")
np_topt_chla_scatter
ggsave(here("Output", "Physiology", "np_topt_chla_scatter.pdf"), np_topt_chla_scatter, h = 5, w = 10)

#np topt and chla ug cm by tissue biomass
np_topt_chla_afdw <- ggplot(np_data) +
  geom_point(aes(y = topt, x = chla_ug_cm2_mean, color = full_species, size = afdw_mg_cm2), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  scale_size_continuous(name = expression("Tissue biomass" ~ (mg ~ cm^{-2}))) +
  #xlim(2,25) +
  #ylim(25,32)+
  facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = "Thermal optimum (°C)",
       color = "Species")
np_topt_chla_afdw
ggsave(here("Output", "Physiology", "np_topt_chla_afdw.pdf"), np_topt_chla_afdw, h = 10, w = 20)

#np topt and chla ug cm by morphology
np_topt_chla_lifehx <- ggplot(np_data) +
  geom_point(aes(y = topt, x = chla_ug_cm2_mean, color = full_species), alpha = 1) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #scale_shape_discrete(name = "Life History Strategy") +
  #coord_transform(x = "log", y = "log")+
  #xlim(2,25) +
  #ylim(25,32)+
  facet_wrap(~full_species, scales = "free", ncol = 4)+
  theme(strip.text = element_text(face = "italic"), legend.position = "none") +
  #geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
  #            method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = "Thermal optimum (°C)",
       color = "Species")
np_topt_chla_lifehx
ggsave(here("Output", "Physiology", "np_topt_chla_spp_facet.pdf"), np_topt_chla_lifehx, h = 8, w = 16)

#np topt and chla ug cm by protein
np_topt_chla_prot <- ggplot(np_data) +
  geom_point(aes(y = topt, x = chla_ug_cm2_mean, color = full_species, size = prot_ug_cm2), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  scale_size_continuous(name = expression("Protein Content" ~ (mu*g ~ cm^{-2}))) +
  facet_wrap(~species, scales = "free")+
  xlim(2,25) +
  ylim(25,32)+
  #geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
  #            method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = "Thermal optimum (°C)",
       color = "Species")
np_topt_chla_prot
ggsave(here("Output", "Physiology", "np_topt_chla_prot.pdf"), np_topt_chla_prot, h = 10, w = 20)

#gp rmax and chla pg sym
gp_rmax_chla_sym_scatter <- ggplot(gp_data) +
  geom_point(aes(y = rmax, x = chla_pg_sym, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~full_species, scales = "free")+
  geom_smooth(aes(y = rmax, x = chla_pg_sym, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Chlorophyll a" ~ (pg ~ symbiont^{-1})), 
       y = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
       color = "Species")
gp_rmax_chla_sym_scatter
ggsave(here("Output", "Physiology", "gp_rmax_chla_sym_scatter.pdf"), gp_rmax_chla_sym_scatter, h = 5, w = 10)

#gp rmax and chla ug cm
gp_rmax_chla_scatter <- ggplot(gp_data) +
  geom_point(aes(y = rmax, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~full_species, scales = "free")+
  geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
       color = "Species")
gp_rmax_chla_scatter
ggsave(here("Output", "Physiology", "gp_rmax_chla_scatter.pdf"), gp_rmax_chla_scatter, h = 5, w = 10)

#r topt and prot ug cm
r_topt_prot_scatter <- ggplot(r_data) +
  geom_point(aes(y = topt, x = prot_ug_cm2, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log", y = "log")+
  #facet_wrap(~full_species, scales = "free")+
  geom_smooth(aes(y = topt, x = prot_ug_cm2, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Protein Content" ~ (mu*g ~ cm^{-2})),
       y = "Thermal optimum (°C)",
       color = "Species")
r_topt_prot_scatter
ggsave(here("Output", "Physiology", "r_topt_prot_scatter.pdf"), r_topt_prot_scatter, h = 5, w = 10)

###Stats
#maya trying to figure out effect size plot
sp_data <- np_data %>% filter(morphology == "branch/tabular")
np_topt_chla_mod <- lm(topt~chla_ug_cm2_mean, data = sp_data)
Anova(np_topt_chla_mod)
summary(np_topt_chla_mod) 
#"Pcyl" 0.2491
#"Mvie" 0.0151 *
#"Fcom" 0.2468
#"Prus" 0.009783 **
#"Peyd" 0.9732
#"Prug" 0.9157
#"Elam" 0.9743
#"Maeq" 0.07218 .
#"Tfro" 0.9265
#"Ahya" 0.7798

models<- np_data %>%
  nest(.by = species) %>% # nest all the data by species
  mutate(fit = map(data, ~lm(topt~chla_ug_cm2_mean, data = .)))
#column of dataframes is called data 
# "." is like "i" with the for loop
#
models

#effect size function:
get_effect_sizes <- function(data, response, predictor, group_var) {
  form <- as.formula(paste0("scale(", response, ") ~ scale(", predictor, ")"))
  
  data %>%
    filter(!is.na(.data[[response]]), !is.na(.data[[predictor]])) %>%
    nest(.by = all_of(group_var)) %>%
    mutate(
      n     = map_int(data, nrow),
      fit   = map(data, ~ lm(form, data = .x)),
      coeffs = map(fit, tidy, conf.int = TRUE)
    ) %>%
    select(-data, -fit) %>%
    unnest(coeffs) %>%
    filter(str_detect(term, "^scale\\(")) %>%   # keep only the slope term
    mutate(response = response, predictor = predictor)
}

morph_effects <- get_effect_sizes(np_data, response  = "topt", predictor = "chla_ug_cm2_mean", group_var = "morphology")
morph_effects

species_effects <- get_effect_sizes(np_data, response  = "topt", predictor = "chla_ug_cm2_mean", group_var = "full_species")
species_effects

effect_plot <- morph_effects %>%
  ggplot(aes(x = estimate, y = morphology, color = morphology)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), width = 0.2) +
  geom_point(size = 3) +
  #scale_color_manual(values = sp_cols)+
  scale_color_manual(values = morph_colors) +
  labs(x = "Standardized effect size\n(Chl a vs Thermal optimum)") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none", axis.title.y = element_blank())
  #theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_text(face = "italic"))
effect_plot

ggsave(here("Output", "Physiology", "effect_size_chla_topt_morphology.pdf"), effect_plot, h = 5, w = 5)

predictors_to_test <- c("chla_ug_cm2_mean", "prot_ug_cm2", "sym_cm2", "afdw_mg_cm2")

all_species_effects <- map_dfr(
  predictors_to_test,
  ~ get_effect_sizes(np_data, response = "topt", predictor = .x,
                     group_var = "full_species")
) %>%
  left_join(species_meta, by = "full_species")

all_morph_effects <- map_dfr(
  predictors_to_test,
  ~ get_effect_sizes(np_data, response = "topt", predictor = .x,
                     group_var = "morphology")
)

multi_effect_plot <- all_species_effects %>%
  ggplot(aes(x = estimate, y = morphology, color = morphology)) +
  #ggplot(aes(x = estimate, y = full_species, color = full_species)) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.3, size = 2) +
  scale_color_manual(values = morph_colors) +
  #scale_color_manual(values = sp_cols) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.7) +
  #geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), width = 0.2) +
  labs(x = "Standardized effect size on thermal optimum", y = "") +
  facet_wrap(~predictor, scales = "free_x", nrow = 1) +
  theme_bw(base_size = 22) +
  theme(legend.position = "none", axis.title.y = element_blank())
  #theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_text(face = "italic"))
multi_effect_plot

ggsave(here("Output","Physiology","physio_topt_effectsize_morphology_meanse.pdf"), multi_effect_plot, height = 5, width = 15)


#testing 
morph_interaction <- lm(scale(topt) ~ scale(chla_ug_cm2_mean) * morphology,
                        data = np_data)
tidy(morph_interaction, conf.int = TRUE)

species_interaction <- lm(scale(topt) ~ scale(chla_ug_cm2_mean) * full_species,
                          data = np_data)
tidy(species_interaction, conf.int = TRUE)

# A significant interaction term (e.g. "scale(chla_ug_cm2_mean):morphologyplating")
# means that group's slope differs significantly from the reference level.
# For an omnibus test of "does the slope differ across ALL groups" (not just
# vs. the reference), compare models with/without the interaction:
anova(lm(scale(topt) ~ scale(chla_ug_cm2_mean), data = np_data),
      species_interaction)



#####Model to look at inter vs intraspecific variation in physiology
library(lme4)

# 1. Null model: species as random intercept only
null_model <- lmer(topt ~ 1 + (1 | species/morphology), data = np_data)
summary(null_model)
# Groups             Name        Variance Std.Dev.
# morphology:species (Intercept) 0.05581  0.2362 species within morphology
# species            (Intercept) 0.04948  0.2224 species overall
# Residual                       0.69535  0.8339 majority of variance

vc_null <- as.data.frame(VarCorr(null_model))
var_species     <- vc_null$vcov[vc_null$grp == "full_species"]
var_resid_null  <- vc_null$vcov[vc_null$grp == "Residual"]
total_var       <- var_species + var_resid_null

# 2. Add chla as a fixed effect, same random-intercept structure
full_model <- lmer(rmax ~ scale(chla_ug_cm2_mean) + scale(sym_cm2) + scale(afdw_mg_cm2) + scale(prot_ug_cm2) + (1 | full_species), data = np_data)
summary(full_model)

vc_full        <- as.data.frame(VarCorr(full_model))
var_resid_full <- vc_full$vcov[vc_full$grp == "Residual"]

# 3. Variance "eaten up" by chla = the drop in residual variance
var_explained_by_physio <- var_resid_null - var_resid_full
var_unexplained        <- var_resid_full

# 4. Build the three-way variance table
variance_table <- data.frame(
  component  = c("Interspecific (species)",
                 "Intraspecific \u2013 explained by physio",
                 "Intraspecific \u2013 unexplained"),
  variance   = c(var_species, var_explained_by_physio, var_unexplained)
) %>%
  mutate(proportion = variance / sum(variance))

variance_table

#topt
#                             component  variance proportion
# 1             Interspecific (species) 0.1052947  0.1315124
# 2 Intraspecific – explained by physio 0.2647873  0.3307176
# 3         Intraspecific – unexplained 0.4305626  0.5377699

#rmax
#                             component   variance proportion
# 1             Interspecific (species) 0.01700589  0.2502145
# 2 Intraspecific – explained by physio 0.01897803  0.2792314
# 3         Intraspecific – unexplained 0.03198132  0.4705541


##np topt and chla
np_topt_chla_mod <- lm(topt~chla_ug_cm2_mean, data = np_data)
Anova(np_topt_chla_mod)
summary(np_topt_chla_mod) 
coeffs<-tidy(np_topt_chla_mod)
#Residual standard error: 0.8391 on 47 degrees of freedom
#Multiple R-squared:  0.1288,	Adjusted R-squared:  0.1103 
#F-statistic: 6.949 on 1 and 47 DF,  p-value: 0.01133 *
check_model(np_topt_chla_mod)

#np rmax and chla
np_rmax_chla_mod <- lm(rmax~chla_ug_cm2_mean, data = np_data)
Anova(np_rmax_chla_mod)
summary(np_rmax_chla_mod) 
# Residual standard error: 0.2411 on 47 degrees of freedom
# Multiple R-squared:  0.1463,	Adjusted R-squared:  0.1281 
# F-statistic: 8.055 on 1 and 47 DF,  p-value: 0.00668 **
check_model(np_rmax_chla_mod)

#gp rmax and chla
gp_rmax_chla_mod <- lm(rmax~chla_pg_sym, data = gp_data)
Anova(gp_rmax_chla_mod)
summary(gp_rmax_chla_mod) 
# Residual standard error: 0.3211 on 47 degrees of freedom
# Multiple R-squared:  0.1187,	Adjusted R-squared:  0.09992 
# F-statistic: 6.329 on 1 and 47 DF,  p-value: 0.01535
check_model(gp_rmax_chla_mod)

#r topt and prot
r_topt_prot_mod <- lm(topt~prot_ug_cm2, data = r_data)
Anova(r_topt_prot_mod)
summary(r_topt_prot_mod) 
# Residual standard error: 0.8712 on 14 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.3958,	Adjusted R-squared:  0.3526 
# F-statistic: 9.171 on 1 and 14 DF,  p-value: 0.009029 **
check_model(r_topt_prot_mod)

####P:R ratios

#read in data
respo_constant_temps <- read_csv(here("Data", "RespoFiles","TPC", "respo_constant_temps.csv"))
phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta_all.csv"))
respo_constant_temps <- respo_constant_temps %>% left_join(phys_meta, by = "frag_ID") %>%
  filter(frag_ID != "B03") #remove extreme outlier
respo_pared_temps <- respo_constant_temps %>% 
  filter(temp_c_value == "24.5" | temp_c_value == "28" | temp_c_value == "31" | temp_c_value == "34")

#quick plots to look at things
PR_plot <- ggplot(data = respo_constant_temps) +
  geom_jitter(aes(x = temp_c_value, y = NPR, color = full_species), width = 0.15, alpha = 0.8) +
  theme_bw(base_size = 22) +
  #theme(axis.text.x = element_blank())+
  #facet_wrap(~species, ncol = 5, scales = "free") +
  facet_wrap(~life_history) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  labs(color = "Species",
       y = "NP:R",
       #x = "Species")
       x = "Temperature (°C)")
PR_plot
ggsave(here("Output", "Physiology", "PR_temp_lifehx.pdf"), PR_plot, h = 8, w = 15)


####Morphology plots#####
gp_data <- all_data %>% filter(PR == "GrossPhoto")
np_data <- all_data %>% filter(PR == "NetPhoto")
r_data <- all_data %>% filter(PR == "Respiration")
physio_data_meta <- all_physio %>% left_join(phys_meta, by = "frag_ID")

rmax_summary <- np_data %>%
  group_by(morphology) %>%
  summarise(n = sum(!is.na(rmax)),
            mean = mean(rmax, na.rm = TRUE),
            se = se_fun(rmax),
            .groups = "drop")
rmax_morph_mod <- lm(rmax~morphology, data = np_data)
Anova(rmax_morph_mod) 

#chla pg sym 0.0007025 ***
#prot 0.02302 *

#0.0001898 ***

#species difs
#prot 0.01727 *
#e 0.0006691 ***

#morph difs
#prot 0.0001898 ***

summary(rmax_morph_mod)
check_model(rmax_morph_mod)
emm_obj <- emmeans::emmeans(rmax_morph_mod, ~ morphology)
emm_pairs <- pairs(emm_obj)

#prot
# contrast                         estimate  SE df t.ratio p.value
# (branch/tabular) - (massive/sub)     -246 152 92  -1.627  0.2397
# (branch/tabular) - plating           -495 180 92  -2.746  0.0197
# (massive/sub) - plating              -249 183 92  -1.360  0.3662

#rmax
# morphology         n  mean     se
# <chr>          <int> <dbl>  <dbl>
#   1 branch/tabular    19 0.864 0.0345
# 2 massive/sub       20 1.05  0.0548
# 3 plating           10 0.666 0.0855

#rmax
# contrast                         estimate     SE df t.ratio p.value
# (branch/tabular) - (massive/sub)   -0.182 0.0701 46  -2.592  0.0335
# (branch/tabular) - plating          0.198 0.0855 46   2.311  0.0642
# (massive/sub) - plating             0.380 0.0848 46   4.475  0.0001

#chla
# contrast                         estimate    SE df t.ratio p.value
# (branch/tabular) - (massive/sub)   -0.592 0.285 95  -2.076  0.1002
# (branch/tabular) - plating          0.761 0.348 95   2.188  0.0785
# (massive/sub) - plating             1.353 0.345 95   3.924  0.0005

means <- physio_data_meta %>% group_by(morphology) %>% summarize(m = mean(prot_ug_cm2, na.rm = TRUE), .groups = "drop")

physio_data_ordered <- physio_data_meta %>% 
  left_join(prot_summary, by = "morphology") %>% 
  mutate(morphology = fct_reorder(morphology, mean))  # ascending by mean

prot_plot <- ggplot() +
  geom_jitter(data = physio_data_ordered, aes(x = morphology, y = prot_ug_cm2, color = morphology), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = prot_summary, aes(x = morphology, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6)+
  geom_point(data = prot_summary, aes(x = morphology, y = mean), size = 2) +
  stat_summary(data = physio_data_ordered, aes(x = morphology, y = prot_ug_cm2), geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "b", "c"))+
  theme_bw(base_size = 22)+
  ylim(1,10)+
  theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x = element_blank()) +
  #theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1, face = "italic")) +
  scale_color_manual(values = morph_colors, name = "Morphology", labels = c("Encrusting/Plating", "Branching/Tabular", "Massive/Submassive")) +
  labs(y = expression("Chlorophyll a" ~ (pg ~ symbiont^{-2})))
prot_plot

ggsave(here("Output", "Physiology", "prot_sym_morphology_jitter.pdf"), prot_plot, h = 5, w = 10)

#####Ordination plots#####
#load physio data

#load metadata
phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta_all.csv"))
#load topt data
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4.csv"))
topt_df <- topt_df %>% filter(sample_ID != "B08_TPC") %>% select(-ctmin,-eh,-q10,-thermal_tolerance,-skewness,-thermal_safety_margin) %>%
  drop_na(e) #cleanup dataframe so things will run, all parameters taken out have too many NAs or infinity values
topt_matrix <- topt_df %>% select(rmax:frag_ID) #generate data for matrix
#separate out data, make sure to put frag_ID as rownames so you can re-join with metadata later
topt_r_data <- topt_matrix %>% filter(PR == "Respiration") %>% select(-PR) %>% column_to_rownames("frag_ID")
topt_np_data <- topt_matrix %>% filter(PR == "NetPhoto") %>% select(-PR) %>% column_to_rownames("frag_ID")
topt_gp_data <- topt_matrix %>% filter(PR == "GrossPhoto") %>% select(-PR) %>% column_to_rownames("frag_ID")
#make them matrices
topt_r <- as.matrix(topt_r_data)
topt_np <- as.matrix(topt_np_data)
topt_gp <- as.matrix(topt_gp_data)
#quick glance at data to look for strong patterns
pairs(x = topt_gp, gap = 0, cex.labels = 0.5) #look similar even between metrics
#scale across variable types
topt_r <- scale(topt_r)
topt_np <- scale(topt_np)
topt_g <- scale(topt_gp)
#generate pca data
pca_topt_r <- prcomp(topt_r)
pca_topt_np <- prcomp(topt_np)
pca_topt_gp <- prcomp(topt_gp)
#collapse PCA data
pc_axes_r <- as.data.frame(pca_topt_r$x)
pc_axes_np <- as.data.frame(pca_topt_np$x)
pc_axes_gp <- as.data.frame(pca_topt_gp$x)
#add frag_ID back
pc_axes_r$frag_ID <- rownames(pc_axes_r) 
pc_axes_np$frag_ID <- rownames(pc_axes_np) 
pc_axes_gp$frag_ID <- rownames(pc_axes_gp) 
#put it back with metadata
pca_r <- pc_axes_r %>% left_join(phys_meta, by = "frag_ID")
pca_np <- pc_axes_np %>% left_join(phys_meta, by = "frag_ID")
pca_gp <- pc_axes_gp %>% left_join(phys_meta, by = "frag_ID")

#plot
pca_spp <- ggplot(pca_r, aes(x = PC1, y = PC2, color = full_species, linetype = perf_imperf)) +
  geom_point(size = 3) +
  scale_color_manual(values = sp_cols) +
  stat_ellipse(aes(group = perf_imperf), level = 0.95, alpha = 0.5, color = "black", linewidth = 0.8) +
  theme_classic(base_size = 22)
pca_spp

#ggsave(here("Output", "Physiology", "r_topt_params_pc_spp_perf.pdf"), pca_spp, h = 5, w = 8)

topt_fit_r <- envfit(pca_topt_r$x[, c("PC1", "PC2")], topt_r_data, permutations = 999, na.rm = TRUE)
topt_fit_np <- envfit(pca_topt_np$x[, c("PC1", "PC2")], topt_np_data, permutations = 999, na.rm = TRUE)
topt_fit_gp <- envfit(pca_topt_gp$x[, c("PC1", "PC2")], topt_gp_data, permutations = 999, na.rm = TRUE)

#topt_fit #can look at p-values of data to see what is significant in driving differences

topt_scores_r <- as.data.frame(scores(topt_fit_r, display = "vectors")) %>% mutate(variable = rownames(.))
topt_scores_np <- as.data.frame(scores(topt_fit_np, display = "vectors")) %>% mutate(variable = rownames(.))
topt_scores_gp <- as.data.frame(scores(topt_fit_gp, display = "vectors")) %>% mutate(variable = rownames(.))

morph_colors <- c(
  "branch/tabular" = "#d6604d",
  "massive/sub" = "#4393c3",
  "plating" = "#74c476")

pca_spp_arrows <- ggplot(pca_np, aes(x = PC1, y = PC2, color = morphology, fill = morphology)) +
  geom_point(size = 3) +
  scale_color_manual(values = morph_colors, name = "Morphology", labels = c("Branching/Tabular", "Massive/Submassive", "Encrusting/Plating")) +
  scale_fill_manual(values = morph_colors, name = "Morphology", labels = c("Branching/Tabular", "Massive/Submassive", "Encrusting/Plating")) +
  #scale_color_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #scale_fill_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  stat_ellipse(geom = "polygon", alpha = 0.1)+
  #stat_ellipse(aes(group = perf_imperf), level = 0.95, alpha = 0.5, color = "black", linewidth = 0.8) +
  geom_segment(data = topt_scores_r,
               aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  #geom_text(data = topt_scores_r,
  #          aes(x = PC1 * 2.2, y = PC2 * 2.2, label = variable),
  #          color = "black", size = 5, inherit.aes = FALSE) +
  theme_classic(base_size = 22)

pca_spp_arrows

ggsave(here("Output", "Physiology", "np_topt_params_morph_arrows_nolabel.pdf"), pca_spp_arrows, h = 8, w = 12)

#ordination of physio data
#load physio data and generate dataframe with id as rownames
all_physio <- read_csv(here("Data", "Physiology", "all_physio_data.csv"))
all_physio_data <- all_physio %>% select(-dw_log,-afdw_log,-chla_log,-chla_sym_log,-sym_log,-prot_log) %>% 
  column_to_rownames("frag_ID") %>% drop_na()
all_physio_data_log <- all_physio %>% select(frag_ID,dw_log,afdw_log,chla_log,chla_sym_log,sym_log,prot_log) %>% 
  column_to_rownames("frag_ID") %>% drop_na()
#make them matrices
all_physio_mat <- as.matrix(all_physio_data)
all_physio_log_mat <- as.matrix(all_physio_data_log)
#quick glance at data to look for strong patterns
#pairs(x = all_physio_mat, gap = 0, cex.labels = 0.5) #look similar even between metrics
#scale across variable types
all_physio_mat <- scale(all_physio_mat)
all_physio_log_mat <- scale(all_physio_log_mat)
#generate pca data
pca_all_physio <- prcomp(all_physio_mat)
pca_all_physio_log <- prcomp(all_physio_log_mat)
#collapse PCA data
pc_axes_physio <- as.data.frame(pca_all_physio$x)
pc_axes_physio_log <- as.data.frame(pca_all_physio_log$x)
#add frag_ID back
pc_axes_physio$frag_ID <- rownames(pc_axes_physio) 
pc_axes_physio_log$frag_ID <- rownames(pc_axes_physio_log) 
#put it back with metadata
pca_physio <- pc_axes_physio %>% left_join(phys_meta, by = "frag_ID")
pca_physio_log <- pc_axes_physio_log %>% left_join(phys_meta, by = "frag_ID")

#plot
pca_physio_spp <- ggplot(pca_physio_log, aes(x = PC1, y = PC2, color = full_species, linetype = perf_imperf)) +
  geom_point(size = 3) +
  scale_color_manual(values = sp_cols) +
  #stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 0.8)
  stat_ellipse(aes(group = perf_imperf), level = 0.95, alpha = 0.5, color = "black", linewidth = 0.8) +
  theme_classic(base_size = 22)
pca_physio_spp

#ggsave(here("Output", "Physiology", "physio_params_log_pc_spp_perf.pdf"), pca_physio_spp, h = 5, w = 8)

physio_fit <- envfit(pca_all_physio$x[, c("PC1", "PC2")], all_physio_data, permutations = 999, na.rm = TRUE)
physio_log_fit <- envfit(pca_all_physio_log$x[, c("PC1", "PC2")], all_physio_data_log, permutations = 999, na.rm = TRUE)

physio_scores <- as.data.frame(scores(physio_fit, display = "vectors")) %>% mutate(variable = rownames(.))
physio_log_scores <- as.data.frame(scores(physio_log_fit, display = "vectors")) %>% mutate(variable = rownames(.))

pca_physio_spp_arrows <- ggplot(pca_physio, aes(x = PC1, y = PC2, color = morphology, fill = morphology)) +
  geom_point(size = 3) +
  scale_color_manual(values = morph_colors, name = "Morphology", labels = c("Branching/Tabular", "Massive/Submassive", "Encrusting/Plating")) +
  scale_fill_manual(values = morph_colors, name = "Morphology", labels = c("Branching/Tabular", "Massive/Submassive", "Encrusting/Plating")) +
  #scale_color_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #scale_fill_manual(values = sp_cols, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  stat_ellipse(geom = "polygon", alpha = 0.1)+
  #stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 0.8)+
  #stat_ellipse(aes(fill = group), level = 0.95, alpha = 0.5, color = "black", linewidth = 0.8) +
  geom_segment(data = physio_scores,
               aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  #geom_text(data = physio_scores,
  #          aes(x = PC1 * 3.2, y = PC2 * 3.2, label = variable),
  #          color = "black", size = 5, inherit.aes = FALSE) +
  theme_classic(base_size = 22)

pca_physio_spp_arrows

ggsave(here("Output", "Physiology", "physio_params_morph_arrows_nolabels.pdf"), pca_physio_spp_arrows, h = 8, w = 12)


#stats for beta dispersion and diversity

#generate distance matrix
physio_dist <- vegdist(pca_all_physio$x, method = "euclidean")

#beta dispersion comparisons
bet.phys <- betadisper(physio_dist,pca_physio$species)
anova(bet.phys) 
#species: p= 0.01511 *
#morphology: p = 0.1687
#lifehx: p = 0.001415 **
#perf: p = 0.196
#plot(bet.phys)
permutest(bet.phys, pairwise = TRUE, permutations = 999)
#significant different dispersion between species:
# Ahya               Elam      Fcom      Maeq      Mvie      Pcyl      Peyd      Prug      Prus  Tfro
# Ahya           0.2270000 0.0380000 0.0820000 0.0430000 0.0040000 0.6900000 0.0090000 0.0110000 0.059
# Elam 0.2243857           0.1480000 0.7190000 0.1660000 0.0130000 0.1690000 0.0630000 0.0520000 0.466
# Fcom 0.0460584 0.1683633           0.1960000 0.9210000 0.6300000 0.0570000 0.8020000 0.7730000 0.393
# Maeq 0.0899390 0.7288105 0.1947960           0.2060000 0.0140000 0.0690000 0.0930000 0.0700000 0.636
# Mvie 0.0526044 0.1693093 0.9214060 0.2063258           0.6970000 0.0620000 0.8920000 0.8520000 0.378
# Pcyl 0.0032141 0.0213470 0.6087902 0.0253814 0.6963382           0.0050000 0.7670000 0.8380000 0.078
# Peyd 0.6568682 0.1694244 0.0690225 0.0761317 0.0724571 0.0078358           0.0190000 0.0180000 0.048
# Prug 0.0154313 0.0737466 0.8142147 0.0913770 0.9001915 0.7673487 0.0266018           0.9400000 0.216
# Prus 0.0147479 0.0677470 0.7740440 0.0820575 0.8556820 0.8318837 0.0265848 0.9445994           0.180
# Tfro 0.0688079 0.4424315 0.3836850 0.5849779 0.3488034 0.0826819 0.0697288 0.2048564 0.1877570     

physio_perm <- adonis2(physio_dist ~ perf_imperf, data = pca_physio, permutations = 999)
physio_perm
#species: R2 = 0.41614, F= 6.4939,p=0.001 ***
#morphology: R2 =  0.06089,   F = 2.8852 p = 0.014 *
#lifehx R2 = 0.25432, F = 10.004, p = 0.001 ***
#perf R2 = 0.04086, F = 3.8345, p = 0.02 *

pairwise.adonis2(physio_dist ~ morphology, data=pca_physio, permutations = 999)
#branch/tabular_vs_massive/sub 0.272
#branch/tabular_vs_plating 0.035 *
#massive/sub_vs_plating 0.004 **

####TOPT
#generate distance matrix
np_dist <- vegdist(pca_topt_np$x, method = "euclidean")

#beta dispersion comparisons
bet.np <- betadisper(np_dist,pca_np$perf_imperf)
anova(bet.np) 
#species: p= 0.307
#morphology: p = 0.6865
#lifehx: p = 0.1694
#perf: p = 0.1142
#plot(bet.phys)
#permutest(bet.phys, pairwise = TRUE, permutations = 999)

np_perm <- adonis2(np_dist ~ morphology, data = pca_np, permutations = 999)
np_perm
#species: R2 = 0.36631, F = 2.3765, p = 0.001 ***
#morphology: R2 =  0.11404,   F = 2.8318 p = 0.003 **
#lifehx R2 = 0.11652, F = 1.8903, p = 0.044 *
#perf: R2 = 0.01904, F= 0.8736, p=0.503
pairwise.adonis2(np_dist ~ morphology, data=pca_np, permutations = 999)
#species differences
#Elam_vs_Fcom 0.018 *
#Elam_vs_Peyd 0.047 *
#Elam_vs_Mvie 0.021 *
#Elam_vs_Maeq 0.009 **
#Elam_vs_Prus 0.028 *
#Elam_vs_Pcyl 0.027 *
#Fcom_vs_Ahya 0.044 *
#Fcom_vs_Tfro 0.017 *
#Fcom_vs_Maeq 0.009 **
#Fcom_vs_Prus 0.049 *
#Prug_vs_Maeq 0.015 *
#Peyd_vs_Ahya 0.006 **
#Peyd_vs_Mvie 0.034 *
#Peyd_vs_Maeq 0.05 *
#Ahya_vs_Maeq 0.023 *
#Ahya_vs_Prus 0.029 *
#Ahya_vs_Pcyl 0.009 **
#Mvie_vs_Maeq 0.009 **
#Mvie_vs_Prus 0.034 *
#Maeq_vs_Pcyl 0.008 **

#massive/sub_vs_branch/tabular p = 0.239
#massive/sub_vs_plating p = 0.001 ***
#branch/tabular_vs_plating p = 0.029 *

#######DO QUICK PHYSIOLOGY PLOTS AND CHECK STATS FOR PAIRED DOWN DATASET####
#use NP dataset for now
all_data <- read_csv(here("Data", "Physiology", "all_data_concatenated.csv"))
gp_data <- all_data %>% filter(PR == "GrossPhoto")
np_data <- all_data %>% filter(PR == "NetPhoto")
r_data <- all_data %>% filter(PR == "Respiration")

#write for loop - for column, test species
columns <- c("rmax","topt","dw_mg_cm2", "afdw_mg_cm2","chla_ug_cm2_mean","chla_pg_sym","sym_cm2","prot_ug_cm2")

# Optional: nicer y-axis labels per variable
y_labels <- list(
  rmax              = "Maximum rate (Rmax)",
  topt              = "Thermal optimum (Topt, °C)",
  dw_mg_cm2         = expression("Dry Weight" ~ (mg ~ cm^{-2})),
  afdw_mg_cm2       = expression("Biomass" ~ (mg ~ cm^{-2})),
  chla_ug_cm2_mean  = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})),
  chla_pg_sym       = expression("Chlorophyll a content" ~ (pg ~ cm^{-2})),
  sym_cm2           = expression("Symbiont density" ~ (cells ~ cm^{-2})),
  prot_ug_cm2       = expression("Host protein content" ~ (mu*g ~ cm^{-2}))
)

analyze_var <- function(var, data) {
  var_sym <- sym(var)
  
  message("Processing variable: ", var)
  
  ## 1. Summarise by species for this variable
  summary_tbl <- data %>%
    group_by(full_species) %>%
    summarise(
      n    = sum(!is.na(!!var_sym)),
      mean = mean(!!var_sym, na.rm = TRUE),
      se   = se_fun(!!var_sym),
      .groups = "drop"
    )
  
  ## 2. Fit model: var ~ full_species
  form <- new_formula(lhs = expr(!!var_sym), rhs = expr(full_species))
  mod  <- lm(form, data = data)
  
  # ANOVA table (coerce to tibble with term column)
  anova_tbl <- car::Anova(mod) %>%
    as.data.frame() %>%
    rownames_to_column("term")
  
  # Grab p-value for full_species
  p_full <- anova_tbl %>%
    filter(term == "full_species") %>%
    pull(`Pr(>F)`) %>%
    first()
  
  # Run model checks (plots) if desired
  performance::check_model(mod)
  
  ## 3. emmeans if significant
  emm_obj   <- NULL
  emm_pairs <- NULL
  
  if (!is.na(p_full) && p_full < 0.05) {
    emm_obj   <- emmeans::emmeans(mod, ~ full_species)
    emm_pairs <- pairs(emm_obj)
  }
  
  ## 4. Plot jitter + means + SE for this variable
  df_ordered <- data %>%
    left_join(summary_tbl, by = "full_species") %>%
    mutate(full_species = fct_reorder(full_species, mean))
  
  y_lab <- y_labels[[var]]
  if (is.null(y_lab)) y_lab <- var
  
  p <- ggplot() +
    geom_jitter(
      data = df_ordered,
      aes(x = full_species, y = !!var_sym, color = full_species),
      width = 0.15, alpha = 0.8
    ) +
    geom_errorbar(
      data = summary_tbl,
      aes(x = full_species, ymin = mean - se, ymax = mean + se),
      width = 0.2, linewidth = 0.6
    ) +
    geom_point(
      data = summary_tbl,
      aes(x = full_species, y = mean),
      size = 2
    ) +
    theme_bw(base_size = 22) +
    theme(
      legend.position = "right",
      axis.text.x     = element_blank(),
      axis.title.x    = element_blank()
    ) +
    scale_color_manual(
      values = sp_cols,
      labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')"))
    ) +
    labs(
      x     = "Species",
      color = "Species",
      y     = y_lab
    )
  # (geom_text / stat_summary for letters intentionally omitted)
  
  # Save plot
  ggsave(
    filename = here("Output", "Physiology", paste0(var, "_species_jitter.pdf")),
    plot     = p,
    h        = 5,
    w        = 10
  )
  
  # Return a tidy list of results for this variable
  list(
    variable   = var,
    summary    = summary_tbl,
    model      = mod,
    anova      = anova_tbl,
    emm_obj    = emm_obj,
    emm_pairs  = emm_pairs,
    plot       = p
  )
}

# Run the pipeline for each variable using purrr
results <- columns %>%
  set_names() %>%
  map(~ analyze_var(.x, np_data))

# Examples of accessing outputs:
# results[["chla_ug_cm2_mean"]]$summary   # summary table
# results[["chla_ug_cm2_mean"]]$anova     # ANOVA table
# results[["chla_ug_cm2_mean"]]$emm_pairs # emmeans pairs (if significant)
# results[["chla_ug_cm2_mean"]]$plot      # ggplot object