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

###### Initial Data Read In ########
#TPC data with only the seven species that we have physio data for
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4.csv"))
topt_df <- topt_df %>% filter(sample_ID != "B08_TPC")
topt_df <- topt_df %>% dplyr::select(rmax, topt, e, PR, frag_ID)
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

#np topt and chla ug cm by life history
np_topt_chla_lifehx <- ggplot(np_data) +
  geom_point(aes(y = topt, x = chla_ug_cm2_mean, color = full_species, shape = life_history), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  scale_shape_discrete(name = "Life History Strategy") +
  #coord_transform(x = "log", y = "log")+
  #xlim(2,25) +
  #ylim(25,32)+
  facet_wrap(~species, scales = "free")+
  geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
  labs(x = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})), 
       y = "Thermal optimum (°C)",
       color = "Species")
np_topt_chla_lifehx
ggsave(here("Output", "Physiology", "np_topt_chla_lifehx.pdf"), np_topt_chla_lifehx, h = 8, w = 15)

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
  geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)+
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
#np topt and chla
np_topt_chla_mod <- lm(topt~chla_ug_cm2_mean, data = np_data)
Anova(np_topt_chla_mod)
summary(np_topt_chla_mod) 
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