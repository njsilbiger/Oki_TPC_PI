#Examine all metabolic data


# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(here)

#species colors
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

####Topt and other parameter graphs####
#read in dataframes and generate prediction dfs for each metric to graph
# PnR_clean <- read_csv(here("Data","RespoFiles","TPC","PnR_clean_no4.csv"))
# #add species names to predictions and topt dfs and pnr data
# BioData <- read_csv(here("Data","RespoFiles","TPC","Fragment_Measurements_TPC.csv"))
# BioSp <- BioData %>% dplyr::select(frag_ID, full_species)
# PnR_clean <-  PnR_clean %>% left_join(BioSp, by = "frag_ID")
# write_csv(PnR_clean, here("Data","RespoFiles","TPC","PnR_clean_no4.csv")) #P and R rate data cleaned
PnR_clean <- read_csv(here("Data","RespoFiles","TPC","PnR_clean_no4.csv"))
preds_all <- read_csv(here("Data","RespoFiles","TPC","Preds_data_clean_no4.csv"))

preds_gp <- preds_all %>% filter(PR == "GrossPhoto")
preds_np <- preds_all %>% filter(PR == "NetPhoto")
preds_resp <- preds_all %>% filter(PR == "Respiration")

preds_all_sp <- read_csv(here("Data","RespoFiles","TPC","Preds_data_clean_no4_species.csv"))

preds_gp_sp <- preds_all_sp %>% filter(PR == "GrossPhoto")
preds_np_sp <- preds_all_sp %>% filter(PR == "NetPhoto")
preds_resp_sp <- preds_all_sp %>% filter(PR == "Respiration")

#and topt data
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4.csv"))
topt_gp <- topt_df %>% filter(PR == "GrossPhoto")
topt_np <- topt_df %>% filter(PR == "NetPhoto")
topt_resp <- topt_df %>% filter(PR == "Respiration")

######plots of predicted TPCs with data#####

#gross photo
gp_pred_plot <- PnR_clean %>% filter(PR == "GrossPhoto") %>% 
  ggplot(aes(x = temp_c_value, y = Values, color = full_species)) +
  geom_point(alpha = 0.7, shape = 21) +
  geom_line(data = preds_gp,
            aes(temp_c_value, .fitted, group = frag_ID),
            linewidth = 0.6) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "italic"), legend.position = "none")+
  scale_color_manual(values = sp_cols)+
  #geom_vline(data = topt_gp,aes(xintercept = topt),linewidth = 0.3, color = "red") +
  #geom_hline(data = topt_gp,aes(yintercept = rmax),linewidth = 0.3, color = "darkgreen") +
  #facet_wrap(~ full_species, scales = "free_y") +
  #ylim(0.48,2.5) +
  facet_wrap(~ full_species, scales = "free", nrow = 2, ncol = 5) +
  labs(x = "Temperature (ºC)",
       y = expression("GP Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))

gp_pred_plot

#ggsave(here("Output", "TPC", "Graphs", "gp_predicted_plot.pdf"),device = "pdf", height = 8, width = 8, gp_pred_plot)

#net photo
np_pred_plot <- PnR_clean %>% filter(PR == "NetPhoto") %>% 
  ggplot(aes(x = temp_c_value, y = Values, color = full_species)) +
  geom_point(alpha = 0.7, shape = 21) +
  geom_line(data = preds_np,
            aes(temp_c_value, .fitted, group = frag_ID),
            linewidth = 0.6) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "italic"), legend.position = "none")+
  scale_color_manual(values = sp_cols)+
  #geom_vline(data = topt_np,aes(xintercept = topt),linewidth = 0.3, color = "red") +
  #geom_hline(data = topt_np,aes(yintercept = rmax),linewidth = 0.3, color = "darkgreen") +
  facet_wrap(~ full_species, scales = "free", nrow = 2, ncol = 5) +
  labs(x = "Temperature (ºC)",
       y = expression("NP Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))

np_pred_plot

#gsave(here("Output", "TPC", "Graphs","np_predicted_plot.pdf"),device = "pdf", height = 8, width = 8, np_pred_plot)

#respiration 
resp_pred_plot <- PnR_clean %>% filter(PR == "Respiration") %>% 
  ggplot(aes(x = temp_c_value, y = Values, color = full_species)) +
  geom_point(alpha = 0.7, shape = 21) +
  geom_line(data = preds_resp,
            aes(temp_c_value, .fitted, group = frag_ID),
            linewidth = 0.6) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "italic"), legend.position = "none")+
  scale_color_manual(values = sp_cols)+
  #geom_vline(data = topt_resp,aes(xintercept = topt),linewidth = 0.3, color = "red") +
  #geom_hline(data = topt_resp,aes(yintercept = rmax),linewidth = 0.3, color = "darkgreen") +
  facet_wrap(~ full_species, scales = "free", nrow = 2, ncol = 5) +
  labs(x = "Temperature (ºC)",
       y = expression("R Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))

resp_pred_plot
#ggsave(here("Output", "TPC", "Graphs","resp_predicted_plot.pdf"),device = "pdf", height = 8, width = 6, resp_pred_plot)

all_pred_plots <- ggarrange(np_pred_plot, gp_pred_plot, resp_pred_plot, 
                            nrow = 3, ncol = 1, labels = c("A", "B", "C"), 
                            font.label = list(size = 20, color = "black"))

ggsave(here("Output","TPC","Graphs","tpc_pred_all_gp_np_r.pdf"), all_pred_plots, h = 12, w = 12)

#TPC prediction plots by species

######plots of predicted TPCs with data#####

#gross photo
gp_pred_plot_sp <- PnR_clean %>% filter(PR == "GrossPhoto") %>% 
  ggplot(aes(x = temp_c_value, y = Values, color = full_species)) +
  geom_point(alpha = 0.7, shape = 21) +
  geom_line(data = preds_gp_sp,
            aes(temp_c_value, .fitted, group = full_species),
            linewidth = 0.6) +
  geom_ribbon(data = preds_gp_sp, aes(x = temp_c_value, ymin = conf_lower, ymax = conf_upper, 
                                      group = full_species, fill = full_species), alpha = 0.2, inherit.aes = F) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "italic"), legend.position = "none")+
  scale_color_manual(values = sp_cols)+
  scale_fill_manual(values = sp_cols) +
  #geom_vline(data = topt_gp,aes(xintercept = topt),linewidth = 0.3, color = "red") +
  #geom_hline(data = topt_gp,aes(yintercept = rmax),linewidth = 0.3, color = "darkgreen") +
  #facet_wrap(~ full_species, scales = "free_y") +
  #ylim(0.48,2.5) +
  facet_wrap(~ full_species, scales = "free", nrow = 2, ncol = 5) +
  labs(x = "Temperature (ºC)",
       y = expression("GP Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))

gp_pred_plot_sp

#ggsave(here("Output", "TPC", "Graphs", "gp_predicted_plot.pdf"),device = "pdf", height = 8, width = 8, gp_pred_plot_sp)

#net photo
np_pred_plot_sp <- PnR_clean %>% filter(PR == "NetPhoto") %>% 
  ggplot(aes(x = temp_c_value, y = Values, color = full_species)) +
  geom_point(alpha = 0.7, shape = 21) +
  geom_line(data = preds_np_sp,
            aes(temp_c_value, .fitted, group = full_species),
            linewidth = 0.6) +
  geom_ribbon(data = preds_np_sp, aes(x = temp_c_value, ymin = conf_lower, ymax = conf_upper, 
                                      group = full_species, fill = full_species), alpha = 0.2,inherit.aes = F) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "italic"), legend.position = "none")+
  scale_color_manual(values = sp_cols)+
  scale_fill_manual(values = sp_cols) +
  #geom_vline(data = topt_np,aes(xintercept = topt),linewidth = 0.3, color = "red") +
  #geom_hline(data = topt_np,aes(yintercept = rmax),linewidth = 0.3, color = "darkgreen") +
  facet_wrap(~ full_species, scales = "free", nrow = 2, ncol = 5) +
  labs(x = "Temperature (ºC)",
       y = expression("NP Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))

np_pred_plot_sp

#gsave(here("Output", "TPC", "Graphs","np_predicted_plot.pdf"),device = "pdf", height = 8, width = 8, np_pred_plot_sp)

#respiration 
resp_pred_plot_sp <- PnR_clean %>% filter(PR == "Respiration") %>% 
  ggplot(aes(x = temp_c_value, y = Values, color = full_species)) +
  geom_point(alpha = 0.7, shape = 21) +
  geom_line(data = preds_resp_sp,
            aes(temp_c_value, .fitted, group = ),
            linewidth = 0.6) +
  geom_ribbon(data = preds_resp_sp, aes(x = temp_c_value, ymin = conf_lower, ymax = conf_upper, 
                                        group = full_species, fill = full_species), alpha = 0.2, inherit.aes = F) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "italic"), legend.position = "none")+
  scale_color_manual(values = sp_cols)+
  scale_fill_manual(values = sp_cols) +
  #geom_vline(data = topt_resp,aes(xintercept = topt),linewidth = 0.3, color = "red") +
  #geom_hline(data = topt_resp,aes(yintercept = rmax),linewidth = 0.3, color = "darkgreen") +
  facet_wrap(~ full_species, scales = "free", nrow = 2, ncol = 5) +
  labs(x = "Temperature (ºC)",
       y = expression("R Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))

resp_pred_plot_sp

all_pred_plot_sp <- ggarrange(np_pred_plot_sp, gp_pred_plot_sp, resp_pred_plot_sp, 
                            nrow = 3, ncol = 1, labels = c("A", "B", "C"), 
                            font.label = list(size = 20, color = "black"))

ggsave(here("Output","TPC","Graphs","tpc_pred_all_gp_np_r_sp.pdf"), all_pred_plot_sp, h = 12, w = 12)



#all rates stacked plots
PnR_gp <- PnR_clean %>% filter(PR == "GrossPhoto")
species_gp_plot <- PnR_gp %>% ggplot(aes(temp_c_value, Values, color = full_species)) +
  geom_point(alpha = 0.7) +
  geom_line(data = preds_gp,
            aes(temp_c_value, .fitted, group = frag_ID, color = full_species),
            linewidth = 0.6) +
  #facet_wrap(~ PR, scales = "free_y") +
  scale_color_discrete(labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 12) +
  labs(x = "Temperature (ºC)",
       y = expression("Gross Photosynthesis" ~ (mu*mol ~ cm^{-2} ~ h^{-1})),
       title = "Thermal performance",
       color = "Species")

species_gp_plot

ggsave(here("Output", "TPC", "Graphs", "gp_stacked_plot_7sp.pdf"),
       device = "pdf", height = 8, width = 8, species_gp_plot)

PnR_np <- PnR_clean %>% filter(PR == "NetPhoto")
species_np_plot <- PnR_np %>% ggplot(aes(temp_c_value, Values, color = full_species)) +
  geom_point(alpha = 0.7) +
  geom_line(data = preds_np,
            aes(temp_c_value, .fitted, group = frag_ID, color = full_species),
            linewidth = 0.6) +
  #facet_wrap(~ PR, scales = "free_y") +
  scale_color_discrete(labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 12) +
  labs(x = "Temperature (ºC)",
       y = expression("Net Photosynthesis" ~ (mu*mol ~ cm^{-2} ~ h^{-1})),
       title = "Thermal performance",
       color = "Species")

species_np_plot

ggsave(here("Output", "TPC", "Graphs", "np_stacked_plot_7sp.pdf"),
       device = "pdf", height = 8, width = 8, species_np_plot)

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

#read in data
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4.csv"))
topt_df <- topt_df %>% filter(sample_ID != "B08_TPC")

#convert to long so you can plot everything at once
metrics <- c("topt", "rmax", "e", "breadth") #select metrics
grouping_vars <- c("species", "PR", "frag_ID", "full_species") #select grouping variables

topt_long <- topt_df %>% #make long dataframe
  dplyr::select(all_of(c(grouping_vars, metrics))) %>%
  pivot_longer(cols = all_of(metrics),
               names_to = "metric",
               values_to = "value")

#now generate mean and se dataframe
se_fun <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

topt_summary <- topt_long %>%
  group_by(full_species,PR, metric) %>%
  summarise(
    n    = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    se   = se_fun(value),
    .groups = "drop")

#first visualize data - plots of Topt and other variables from TPCs######
topt_plot <- ggplot() +
  geom_jitter(data = topt_long,
              aes(x = full_species, y = value, color = full_species),
              width = 0.15, alpha = 0.8) +
  geom_errorbar(data = topt_summary,
                aes(x = full_species,
                    ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.6) +
  geom_point(data = topt_summary,
             aes(x = full_species, y = mean),
             size = 2) +
  facet_wrap(PR ~ metric, scales = "free_y") +
  theme_bw(base_size = 12) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))
#labs(x = )
topt_plot

ggsave(here("Output", "TPC", "Graphs","Topt_allparams_clean_no4.pdf"),
       device = "pdf", height = 8, width = 8, topt_plot)

#individual dataframes to put in individual plots so they can be ordered (live laugh love i guess)
rmax_mean_gp <- topt_summary %>% filter(PR == "GrossPhoto") %>% filter(metric == "rmax")
topt_rmax_gp <- topt_long %>% 
  filter(PR == "GrossPhoto") %>% 
  filter(metric == "rmax") %>%
  left_join(rmax_mean_gp, by = "full_species") %>% 
  mutate(full_species = fct_reorder(full_species, mean))

rmax_mean_np <- topt_summary %>% filter(PR == "NetPhoto") %>% filter(metric == "rmax")
topt_rmax_np <- topt_long %>% 
  filter(PR == "NetPhoto") %>% 
  filter(metric == "rmax") %>%
  left_join(rmax_mean_np, by = "full_species") %>% 
  mutate(full_species = fct_reorder(full_species, mean))

topt_mean_gp <- topt_summary %>% filter(PR == "GrossPhoto") %>% filter(metric == "topt")
topt_topt_gp <- topt_long %>% 
  filter(PR == "GrossPhoto") %>% 
  filter(metric == "topt") %>%
  left_join(topt_mean_gp, by = "full_species") %>% 
  mutate(full_species = fct_reorder(full_species, mean))

topt_mean_np <- topt_summary %>% filter(PR == "NetPhoto") %>% filter(metric == "topt")
topt_topt_np <- topt_long %>% 
  filter(PR == "NetPhoto") %>% 
  filter(metric == "topt") %>%
  left_join(topt_mean_np, by = "full_species") %>% 
  mutate(full_species = fct_reorder(full_species, mean))

e_mean_gp <- topt_summary %>% filter(PR == "GrossPhoto") %>% filter(metric == "e")
topt_e_gp <- topt_long %>% 
  filter(PR == "GrossPhoto") %>% 
  filter(metric == "e") %>%
  left_join(e_mean_gp, by = "full_species") %>% 
  mutate(full_species = fct_reorder(full_species, mean))

e_mean_np <- topt_summary %>% filter(PR == "NetPhoto") %>% filter(metric == "e")
topt_e_np <- topt_long %>% 
  filter(PR == "NetPhoto") %>% 
  filter(metric == "e") %>%
  left_join(e_mean_np, by = "full_species") %>% 
  mutate(full_species = fct_reorder(full_species, mean))

breadth_mean_gp <- topt_summary %>% filter(PR == "GrossPhoto") %>% filter(metric == "breadth")
topt_breadth_gp <- topt_long %>% 
  filter(PR == "GrossPhoto") %>% 
  filter(metric == "breadth") %>%
  left_join(breadth_mean_gp, by = "full_species") %>% 
  mutate(full_species = fct_reorder(full_species, mean))

breadth_mean_np <- topt_summary %>% filter(PR == "NetPhoto") %>% filter(metric == "breadth")
topt_breadth_np <- topt_long %>% 
  filter(PR == "NetPhoto") %>% 
  filter(metric == "breadth") %>%
  left_join(breadth_mean_np, by = "full_species") %>% 
  mutate(full_species = fct_reorder(full_species, mean))

topt_breadth_plot <- ggplot() +
  geom_jitter(data = topt_breadth_np, aes(x = full_species, y = value, color = full_species), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = breadth_mean_np,aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = breadth_mean_np, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  stat_summary(data = topt_breadth_np, aes(x = full_species, y = value), geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "b"))+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "right")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), legend.position = "none")+
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  labs(x = "Species", y = "Breadth (°C)", color = "Species")

topt_e_plot <- ggplot() +
  geom_jitter(data = topt_e_np, aes(x = full_species, y = value, color = full_species), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = e_mean_np,aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = e_mean_np, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  ylim(0,1)+
  stat_summary(data = topt_e_np, aes(x = full_species, y = value), geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "b"))+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "right")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), , legend.position = "none", axis.title.y = element_text(face = "italic"))+
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  labs(x = "Species", y = "E (eV)", color = "Species")

topt_topt_plot <- ggplot() +
  geom_jitter(data = topt_topt_np, aes(x = full_species, y = value, color = full_species), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = topt_mean_np,aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = topt_mean_np, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "right")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), legend.position = "none")+
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  labs(x = "Species", y = "Thermal optimum (°C)", color = "Species")

topt_rmax_plot <- ggplot() +
  geom_jitter(data = topt_rmax_np, aes(x = full_species, y = value, color = full_species), width = 0.15, alpha = 0.8) +
  geom_errorbar(data = rmax_mean_np,aes(x = full_species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = rmax_mean_np, aes(x = full_species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  stat_summary(data = topt_rmax_np, aes(x = full_species, y = value), geom = "text", fun = max, vjust = -0.5, size = 8,
               label = c("a", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "b", "b"))+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),legend.position = "right")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), legend.position = "none")+
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  ylim(0,2)+
  labs(x = "Species", y = expression("Rmax" ~ (mu*mol ~ cm^{-2} ~ h^{-1})), , color = "Species")

library(ggpubr)
np_topt_plots <- ggarrange(topt_rmax_plot, topt_topt_plot, 
                           nrow = 1, ncol = 2, legend = "right", common.legend = TRUE)
np_topt_plots

ggsave(here("Output","TPC","Graphs","np_topt_plots_2.pdf"), np_topt_plots, h = 5, w = 15, dpi = 300)


#topt plot all together
topt_plot <- ggplot() +
  geom_jitter(data = topt_long,
              aes(x = full_species, y = value, color = full_species),
              width = 0.15, alpha = 0.8) +
  geom_errorbar(data = topt_summary,
                aes(x = full_species,
                    ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.6) +
  geom_point(data = topt_summary,
             aes(x = full_species, y = mean),
             size = 2) +
  facet_wrap(PR ~ metric, scales = "free_y") +
  theme_bw(base_size = 12) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))
#labs(x = )
topt_plot



#check distributions of data as well
ggplot(topt_long, aes(value)) +
  geom_histogram(bins = 30) + 
  facet_wrap(PR~metric, scales = "free")

#run a loop to model the data and get the output and the pairwise comparisons
library(broom)

PR_levels <- c("GrossPhoto", "NetPhoto", "Respiration")
resp_metrics <- c("rmax", "topt", "e", "breadth")

#data lists
emm_list   <- list()
emm_pair_list <- list()
anova_list <- list()

for (pr in PR_levels) {
  dat_pr <- topt_long %>% filter(PR == pr)
  
  for (resp in resp_metrics) {
    dat <- dat_pr %>% filter(metric == resp)
    
    fit <- lm(value ~ species, data = dat)
    
    # ANOVA table -> tidy tibble
    an_tbl <- broom::tidy(anova(fit)) %>%
      mutate(PR = pr, metric = resp, .before = 1)
    anova_list[[paste(pr, resp, sep = "__")]] <- an_tbl
    
    # EMMs per species (+ CIs) -> tibble
    emm_obj <- emmeans::emmeans(fit, ~ species)
    emm_pairs <- pairs(emm_obj)
    emm_pair_tbl <- as_tibble(summary(emm_pairs, infer = TRUE)) %>%
      mutate(PR = pr, metric = resp, .before = 1)
    emm_tbl <- as_tibble(summary(emm_obj, infer = TRUE)) %>%
      mutate(PR = pr, metric = resp, .before = 1)
    emm_list[[paste(pr, resp, sep = "__")]] <- emm_tbl
    emm_pair_list[[paste(pr, resp, sep = "__")]] <- emm_pair_tbl
  }
}

#add values to table
emm_all   <- bind_rows(emm_list)
anova_all <- bind_rows(anova_list)
emm_pair_set <- bind_rows(emm_pair_list)

#NP
write_csv(emm_all, here("Data","RespoFiles","TPC","emmeans_all_PR_metrics.csv"))
write_csv(anova_all, here("Data","RespoFiles","TPC","anova_all_PR_metrics.csv"))

emm_plot <- ggplot(emm_all, aes(x = emmean, y = species, color = species)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.2) +
  theme_bw(base_size = 12) +
  facet_wrap(PR ~ metric, scales = "free") +
  labs(title = "emmeans + confidence intervals")
emm_plot

ggsave(here("Output", "TPC", "Graphs","Emmeans_allparams_clean_no4.pdf"),
       device = "pdf", height = 8, width = 8, emm_plot)

### Compare with Danielle's Data ###
PnR_clean <- read_csv(here("Data","RespoFiles","TPC","PnR_clean_no4.csv"))

species_cols <- c(
  "Ahya" = '#d8aedd',
  "Elam" = '#ba7999',
  "Fcom"   = '#dd4124',
  "Maeq" = '#ed8b00',
  "Mvie" = '#efbc82',
  "Prug" = '#edd746',
  "Peyd" = '#d0e2af',
  "Pcyl" = '#45681e',
  "Prus" = '#7bbcd5',
  "Tfro" = '#00496f'
)

se_fun <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

#use PnR_clean to calculate NP:R for every single fragment at select temperatures
#do 26, 29, and 32
PnR_28 <- PnR_clean %>% filter(temp_c_value == 28)
PnR_28_summary <- PnR_28 %>%
  group_by(species, PR) %>%
  summarise(n = sum(!is.na(Values)),
            mean = mean(Values, na.rm = TRUE),
            se = se_fun(Values),
            .groups = "drop")

R_26 <- PnR_clean %>% filter(temp_c_value == 26) %>% filter(PR == "Respiration") %>%
  select(Values, frag_ID) %>% rename(R_26 = Values)
R_29 <- PnR_clean %>% filter(temp_c_value == 29) %>% filter(PR == "Respiration") %>%
  select(Values, frag_ID) %>% rename(R_29 = Values)
R_32 <- PnR_clean %>% filter(temp_c_value == 32) %>% filter(PR == "Respiration") %>%
  select(Values, frag_ID) %>% rename(R_32 = Values)

NP_26 <- PnR_clean %>% filter(temp_c_value == 26) %>% filter(PR == "NetPhoto") %>%
  select(Values, frag_ID) %>% rename(NP_26 = Values)
NP_29 <- PnR_clean %>% filter(temp_c_value == 29) %>% filter(PR == "NetPhoto") %>%
  select(Values, frag_ID) %>% rename(NP_29 = Values)
NP_32 <- PnR_clean %>% filter(temp_c_value == 32) %>% filter(PR == "NetPhoto") %>%
  select(Values, frag_ID) %>% rename(NP_32 = Values)

GP_26 <- PnR_clean %>% filter(temp_c_value == 26) %>% filter(PR == "GrossPhoto") %>%
  select(Values, frag_ID) %>% rename(GP_26 = Values)
GP_29 <- PnR_clean %>% filter(temp_c_value == 29) %>% filter(PR == "GrossPhoto") %>%
  select(Values, frag_ID) %>% rename(GP_29 = Values)
GP_32 <- PnR_clean %>% filter(temp_c_value == 32) %>% filter(PR == "GrossPhoto") %>%
  select(Values, frag_ID) %>% rename(GP_32 = Values)

df_full_list <- list(R_26, R_29, R_32, NP_26, NP_29, NP_32, GP_26, GP_29, GP_32)

respo_select_temps <- df_full_list %>% reduce(left_join)
respo_select_temps <- respo_select_temps %>% drop_na(NP_26) %>%
  mutate(NPR_26 = NP_26/R_26,
         NPR_29 = NP_29/R_29,
         NPR_32 = NP_32/R_32)
write_csv(respo_select_temps, here("Data", "RespoFiles","TPC", "respo_select_temps.csv"))


R <- PnR_clean %>% filter(PR == "Respiration") %>%
  select(frag_ID, temp_c_value, Values) %>% rename(R = Values)
NP <- PnR_clean %>% filter(PR == "NetPhoto") %>%
  select(frag_ID, temp_c_value, Values) %>% rename(NP = Values)
GP <- PnR_clean %>% filter(PR == "GrossPhoto") %>%
  select(frag_ID, temp_c_value, Values) %>% rename(GP = Values)

df_list <- list(R, NP, GP)

respo_constant_temps <- df_list %>% reduce(left_join)
respo_constant_temps <- respo_constant_temps %>% drop_na(NP) %>%
  mutate(NPR = NP/R,
         GPR = GP/R)
write_csv(respo_constant_temps, here("Data", "RespoFiles","TPC", "respo_constant_temps.csv"))

sp_mod <- lm(Values~species, data = GP_28)
Anova(sp_mod)
#summary(sp_mod)
#check_model(sp_mod)

emm_obj <- emmeans::emmeans(sp_mod, ~ species)
emm_pairs <- pairs(emm_obj)
emm_pairs

#GP 28 - Maeq lower than Elam, Fcom, Pcyl, and Prus
#Elam - Maeq  0.60595 0.174 40   3.482  0.0358
#Fcom - Maeq  0.65241 0.174 40   3.749  0.0178
#Maeq - Pcyl -0.59749 0.174 40  -3.433  0.0405
#Maeq - Prus -0.67170 0.174 40  -3.860  0.0132

#R 28 - Maeq lower than Prus and Pcyl
# Maeq - Prus -0.23000 0.0637 40  -3.612  0.0256
# Maeq - Pcyl -0.24705 0.0637 40  -3.880  0.0125

#GP 31 - 
# Fcom - Maeq  7.64e-01 0.184 40   4.150  0.0058
# Maeq - Pcyl -6.18e-01 0.184 40  -3.358  0.0488
# Maeq - Prug -6.18e-01 0.184 40  -3.358  0.0488
#Maeq - Prus -7.32e-01 0.184 40  -3.979  0.0095


#R 31 - Maeq only lower than Prus (almost Pcyl)
# Maeq - Prus -0.28185 0.0718 40  -3.927  0.0109
# Maeq - Pcyl -0.23840 0.0718 40  -3.322  0.0533

PnR_28_ordered <- PnR_28 %>% 
  left_join(PnR_28_summary) %>% 
  group_by(PR) %>%
  mutate(species = fct_reorder(species, mean))  # ascending by mean

PnR_28_plot <- ggplot() +
  geom_jitter(data = PnR_28_ordered, aes(x = species, y = Values, color = species), width = 0.15, alpha = 0.8) +
  facet_wrap(.~PR)+
  #stat_summary(data = PnR_28_ordered, aes(x = species, y = Values), geom = "text", fun = max, vjust = -0.5, size = 8,
  #             label = c("a", "a", "a", "a", "ab", "ab", "b"))+
  geom_errorbar(data = PnR_28_summary, aes(x = species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = PnR_28_summary, aes(x = species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x = element_blank()) +
  #theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1, face = "italic")) +
  scale_color_manual(values = species_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #ylim(25,325)+
  labs(x = "Species", color = "Species", y = expression("28 °C Metabolic Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))
PnR_28_plot

ggsave(here("Output", "TPC", "Graphs","PnR_means_28C.pdf"), device = "pdf", height = 8, width = 12, PnR_28_plot)

PnR_31 <- PnR_clean %>% filter(temp_c_value == 31)
PnR_31_summary <- PnR_31 %>%
  group_by(species, PR) %>%
  summarise(n = sum(!is.na(Values)),
            mean = mean(Values, na.rm = TRUE),
            se = se_fun(Values),
            .groups = "drop")

PnR_31_ordered <- PnR_31 %>% 
  left_join(PnR_31_summary) %>% 
  group_by(PR)%>%
  mutate(species = fct_reorder(species, mean))  # ascending by mean

PnR_31_plot <- ggplot() +
  geom_jitter(data = PnR_31_ordered, aes(x = species, y = Values, color = species), width = 0.15, alpha = 0.8) +
  facet_wrap(.~PR)+
  #stat_summary(data = PnR_31_ordered, aes(x = species, y = Values), geom = "text", fun = max, vjust = -0.5, size = 8,
  #             label = c("a", "a", "a", "a", "ab", "ab", "b"))+
  geom_errorbar(data = PnR_31_summary, aes(x = species, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.6) +
  geom_point(data = PnR_31_summary, aes(x = species, y = mean), size = 2) +
  theme_bw(base_size = 22) +
  theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x = element_blank()) +
  #theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1, face = "italic")) +
  scale_color_manual(values = species_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #ylim(25,325)+
  labs(x = "Species", color = "Species", y = expression("31 °C Metabolic Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))
PnR_31_plot

ggsave(here("Output", "TPC", "Graphs","PnR_means_31C.pdf"), device = "pdf", height = 8, width = 12, PnR_31_plot)
