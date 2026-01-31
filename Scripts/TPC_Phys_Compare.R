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

###### Initial Data Read In ########
#TPC data with only the seven species that we have physio data for
topt_df <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4.csv"))
topt_df <- topt_df %>% filter(sample_ID != "B08_TPC")
topt_df <- topt_df %>% dplyr::select(rmax, topt, e, PR, frag_ID, species, full_species)

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
avg_DW <- read.csv(here("Data", "Physiology", "Average_Dry_Weight.csv")) #dw_mg_cm2
avg_DW <- avg_DW %>% mutate(dw_log = log(dw_mg_cm2)) %>% dplyr::select(frag_ID, dw_mg_cm2,dw_log)
chla_avg <- read.csv(here("Data", "Physiology", "Chla_avg.csv")) #chla_ug_cm2_mean
chla_avg <- chla_avg %>% filter(frag_ID != "C07")#remove outlier
chla_avg <- chla_avg %>% mutate(chla_log = log(chla_ug_cm2_mean)) %>% dplyr::select(frag_ID, chla_ug_cm2_mean,chla_log)
avg_sym <- read.csv(here("Data", "Physiology", "Average_Sym_Density.csv"))
avg_sym <- avg_sym %>% filter(frag_ID != "C07") %>% filter(frag_ID != "D10") #remove crazy outliers!
avg_sym <- avg_sym %>% mutate(sym_log = log(sym_cm2)) %>% dplyr::select(frag_ID, sym_cm2,sym_log)
prot <- read.csv(here("Data", "Physiology", "protein_all_summary.csv")) #prot_ug_cm2
prot <- prot %>% mutate(prot_log = log(prot_ug_cm2)) %>% dplyr::select(frag_ID, prot_ug_cm2,prot_log)

#generate full dataframe
result <- Reduce(function(x, y) merge(x, y, all = TRUE), list(topt_df,avg_DW,chla_avg,avg_sym,prot))
#drop NA data from respo data because we donʻt have all reps
all_data <- result %>% drop_na(rmax)
write_csv(result, here("Data", "Physiology", "all_data_concatenated.csv"))
all_data <- read_csv(here("Data", "Physiology", "all_data_concatenated.csv"))

sp_keep <- c("Fcom","Prus", "Peyd", "Elam", "Maeq", "Tfro", "Ahya")
data_7sp <- all_data %>% filter(species %in% sp_keep)

gp_data <- all_data %>% filter(PR == "GrossPhoto")
np_data <- all_data %>% filter(PR == "NetPhoto")
r_data <- all_data %>% filter(PR == "Respiration")

library(psych)
library(corrplot)
vars <- np_data |>
  dplyr::select(
    rmax,
    topt,
    dw_mg_cm2,
    #dw_log,
    chla_ug_cm2_mean,
    #chla_log,
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

ggsave(here("Output", "Physiology", "corr_np.pdf"), corr_plot, h = 6, w = 6)



####DO QUICK PHYSIOLOGY PLOTS AND CHECK STATS FOR PAIRED DOWN DATASET####
#use NP dataset for now
all_data <- read_csv(here("Data", "Physiology", "all_data_concatenated.csv"))
np_data <- all_data %>% filter(PR == "NetPhoto")

#write for loop - for column, test species

library(dplyr)
library(ggplot2)
library(forcats)
library(car)
library(emmeans)
library(performance)
library(purrr)
library(rlang)
library(tibble)
library(here)

columns <- c("rmax","topt","dw_mg_cm2","chla_ug_cm2_mean","sym_cm2","prot_ug_cm2")

# Optional: nicer y-axis labels per variable
y_labels <- list(
  rmax              = "Maximum rate (Rmax)",
  topt              = "Thermal optimum (Topt, °C)",
  dw_mg_cm2         = expression("Tissue biomass" ~ (mg ~ cm^{-2})),
  chla_ug_cm2_mean  = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})),
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


#DW chla scatter
avg_DW <- read.csv(here("Data", "Physiology", "Average_Dry_Weight.csv"))
avg_DW <- avg_DW %>% rename_at('species_long', ~'full_species')
avg_DW <- avg_DW %>% mutate(dw_log = log(dw_mg_cm2))

chla_avg <- read.csv(here("Data", "Physiology", "Chla_avg.csv"))
chla_avg <- chla_avg %>% mutate(chla_log = log(chla_ug_cm2_mean))

chla_dw <- power_full_join(chla_avg, avg_DW, by = "frag_ID", conflict = coalesce_xy) %>% 
  drop_na()

dw_chla_scat <- ggplot(chla_dw) +
  geom_point(aes(x = dw_mg_cm2, y = chla_ug_cm2_mean, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  labs(x = expression("Tissue biomass" ~ (mg ~ cm^{-2})), 
       y = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})), 
       color = "Species") +
  geom_smooth(aes(x = dw_mg_cm2, y = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black", linewidth = 1.1)
dw_chla_scat

ggsave(here("Output", "Physiology", "dryweight_chla_scatter.pdf"), dw_chla_scat, h = 6, w = 10)

##### dry weight #####
#avg_DW <- read.csv(here("Data", "Physiology", "Average_Dry_Weight.csv"))
#avg_DW <- avg_DW %>% rename_at('species_long', ~'full_species')

library(powerjoin)
DW_topt <- power_full_join(avg_DW, topt_df, by = "frag_ID", conflict = coalesce_xy) %>% 
  drop_na() %>%
  filter(PR != "Respiration")

#scatterplot of dw and topt params

#topt
dw_topt_scatter <- ggplot(filter(all_data, PR == "NetPhoto")) +
  geom_point(aes(y = topt, x = dw_mg_cm2, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  #coord_transform(x = "log")+
  #facet_wrap(~PR, scales = "free")+
  geom_smooth(aes(y = topt, x = dw_mg_cm2, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(x = expression("Tissue biomass" ~ (mg ~ cm^{-2})), 
       y = "Thermal optimum (°C)",
       color = "Species")
dw_topt_scatter

by(DW_topt, DW_topt$PR, function(d) summary(lm(topt ~ dw_mg_cm2, data = d))) #*
library(performance)
check_model(lm(topt ~ dw_mg_cm2, data = filter(DW_topt, PR == "NetPhoto")))
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(topt ~ dw_mg_cm2, data = d)))
dw_topt_mod <- lm(topt ~ dw_mg_cm2, data = filter(DW_topt, PR == "NetPhoto"))
summary(dw_topt_mod)
anova(dw_topt_mod)

#topt net photo
#dw_log        0.5526     0.2563   2.156   0.0395 * 
#topt GP
#dw_log        0.1981     0.1701   1.165    0.252

#topt
#DW_topt$PR: GrossPhoto
#dw_mg_cm2    0.006035   0.002387   2.528   0.0172 *  significant
#none individual sig dif

ggsave(here("Output", "Physiology", "dryweight_topt_np.pdf"), dw_topt_scatter, h = 5, w = 10)

#rmax

dw_rmax_scatter <- ggplot(filter(DW_topt, PR == "NetPhoto")) +
  geom_point(aes(y = rmax, x = dw_mg_cm2, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  coord_transform(x = "log")+
  #facet_wrap(~PR, scales = "free")+
  #geom_smooth(aes(y = rmax, x = dw_mg_cm2, group = 1),
  #            method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(x = expression("Tissue biomass" ~ (mg ~ cm^{-2})), 
       y = expression("Rmax" ~ (mu*mol ~ cm^{-2} ~ h^{-1})), 
       color = "Species")
dw_rmax_scatter

ggsave(here("Output", "Physiology", "dryweight_rmax_log.pdf"), dw_rmax_scatter, h = 4, w = 12)

by(DW_topt, DW_topt$PR, function(d) summary(lm(rmax ~ dw_mg_cm2, data = d))) #ns
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(rmax ~ dw_mg_cm2, data = d)))
#none sig

#e
dw_e_scatter <- ggplot(filter(DW_topt, PR == "GrossPhoto")) +
  geom_point(aes(y = e, x = dw_mg_cm2, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  coord_transform(x = "log")+
  #facet_wrap(~PR, scales = "free")+
  #geom_smooth(aes(y = e, x = dw_mg_cm2, group = 1),
  #            method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(x = expression("Tissue biomass" ~ (mg ~ cm^{-2})), 
       y = ("E (eV)"),  
       color = "Species")
dw_e_scatter

ggsave(here("Output", "Physiology", "dryweight_e_log.pdf"), dw_e_scatter, h = 4, w = 12)

by(DW_topt, DW_topt$PR, function(d) summary(lm(e ~ dw_log, data = d)))
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(e ~ dw_mg_cm2, data = d)))

#e DW_topt$PR: GrossPhoto
#dw_log      -0.16421    0.07198  -2.281  0.03307 * 

#e
#DW_topt$PR: GrossPhoto
#dw_mg_cm2   -0.0019061  0.0007072  -2.695   0.0136 *  
#none significant


#breadth
dw_breadth_scatter <- ggplot(filter(DW_topt, PR == "GrossPhoto")) +
  geom_point(aes(y = breadth, x = dw_mg_cm2, color = full_species), alpha = 1) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  coord_transform(x = "log")+
  #facet_wrap(~PR, scales = "free")+
  geom_smooth(aes(y = breadth, x = dw_mg_cm2, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(x = expression("Tissue biomass" ~ (mg ~ cm^{-2})), 
       y = expression("Breadth (°C)"), 
       color = "Species")
dw_breadth_scatter

dw_breadth_scatter <- dw_breadth_scatter +
  #species-specific regression lines
  #geom_smooth(aes(y = breadth, x = dw_mg_cm2, color = full_species),
  #            method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = breadth, x = dw_log, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_breadth_log.pdf"), dw_breadth_scatter, h = 4, w = 12, dpi = 300)

by(DW_topt, DW_topt$PR, function(d) summary(lm(breadth ~ dw_log, data = d))) #ns
by(DW_topt, list(DW_topt$PR, DW_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(breadth ~ dw_mg_cm2, data = d)))

#breadth DW_topt$PR: GrossPhoto
#dw_log        0.6504     0.2504   2.598   0.0146 *  

#breadth
#DW_topt$PR: GrossPhoto
#dw_mg_cm2   0.005614   0.002442   2.299   0.0289 * 

#breadth
#GrossPhoto
#Turbinaria frondens
#dw_mg_cm2   0.035475   0.008652    4.10   0.0262 *


library(ggpubr)
dw_topt_log <- ggarrange(dw_rmax_scatter, dw_topt_scatter, dw_e_scatter, dw_breadth_scatter, 
                           nrow = 2, ncol = 2, legend = "right", common.legend = TRUE)
dw_topt_log

ggsave(here("Output","TPC","Graphs","dw_topt_log_np.pdf"), dw_topt_log, h = 8, w = 14, dpi = 300)


##### chlorophyll a #####
chla_avg <- read.csv(here("Data", "Physiology", "Chla_avg.csv"))

chla_topt <- power_full_join(topt_df, chla_avg, by = "frag_ID", conflict = coalesce_xy) %>% 
  drop_na() %>%
  filter(PR != "Respiration")

#scatterplot of chla and topt params

#topt
chla_topt_scatter <- ggplot(filter(all_data, PR == "GrossPhoto"))+
  geom_point(aes(y = topt, x = chla_ug_cm2_mean, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free_y")+
  #ylim(25,250)+
  labs(x = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})), y = "Thermal optimum (°C)", color = "Species")
chla_topt_scatter

chla_topt_scatter <- chla_topt_scatter +
  #species-specific regression lines
  #geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, color = full_species),
  #            method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "chla_topt_scatter_reg.pdf"), chla_topt_scatter, h = 4, w = 12)

by(all_data, all_data$PR, function(d) summary(lm(topt ~ chla_log, data = d))) #ns
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(topt ~ chla_ug_cm2_mean, data = d))) #ns

#rmax
chla_rmax_scatter <- ggplot(filter(all_data, PR == "GrossPhoto"))+
  geom_point(aes(y = rmax, x = chla_ug_cm2_mean, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  labs(x = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})), 
       y = expression("Rate Maximum" ~ (mu*mol ~ cm^{-2} ~ h^{-1})), 
       color = "Species")
chla_rmax_scatter

chla_rmax_scatter <- chla_rmax_scatter +
  #species-specific regression lines
  geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "chla_rmax_scatter_reglines.pdf"), chla_rmax_scatter, h = 4, w = 12)

by(all_data, all_data$PR, function(d) summary(lm(rmax ~ chla_ug_cm2_mean, data = d))) #ns
by(all_data, list(all_data$PR, all_data$full_species),
   function(d) if(nrow(d) > 1) summary(lm(rmax ~ chla_ug_cm2_mean, data = d)))

#rmax chla
#GrossPhoto Echinopora lamellosa
#chla_ug_cm2_mean  0.13595    0.01915   7.098  0.00575 **

#e
chla_e_scatter <- ggplot(chla_topt) +
  geom_point(aes(y = e, x = chla_ug_cm2_mean, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  labs(x = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})), 
       y = ("Activation energy (eV)"), 
       color = "Species")
chla_e_scatter

chla_e_scatter <- chla_e_scatter +
  #species-specific regression lines
  geom_smooth(aes(y = e, x = chla_ug_cm2_mean, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = e, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "chla_e_scatter_reglines.pdf"), chla_e_scatter, h = 4, w = 12)

by(all_data, all_data$PR, function(d) summary(lm(e ~ chla_log, data = d))) #ns
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(e ~ chla_ug_cm2_mean, data = d))) #ns


#breadth
chla_breadth_scatter <- ggplot(chla_topt) +
  geom_point(aes(y = breadth, x = chla_ug_cm2_mean, color = full_species)) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  facet_wrap(~PR, scales = "free")+
  labs(x = expression("Chlorophyll a content" ~ (mu*g ~ cm^{-2})), 
       y = expression("Thermal performance breadth (°C)"), 
       color = "Species")
chla_breadth_scatter


chla_breadth_scatter <- chla_breadth_scatter +
  #species-specific regression lines
  geom_smooth(aes(y = breadth, x = chla_ug_cm2_mean, color = full_species),
              method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = breadth, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black",
              linetype = "longdash", linewidth = 1.1)

ggsave(here("Output", "Physiology", "chla_breadth_scatter_reglines.pdf"), chla_breadth_scatter, h = 4, w = 12)

by(chla_topt, chla_topt$PR, function(d) summary(lm(breadth ~ chla_ug_cm2_mean, data = d))) #ns
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(breadth ~ chla_ug_cm2_mean, data = d))) #ns

#sym density
by(all_data, all_data$PR, function(d) summary(lm(topt ~ sym_log, data = d))) #ns
by(all_data, all_data$PR, function(d) summary(lm(rmax ~ sym_log, data = d))) #ns
by(all_data, all_data$PR, function(d) summary(lm(e ~ sym_log, data = d))) #ns

#protein
by(all_data, all_data$PR, function(d) summary(lm(topt ~ prot_log, data = d))) #ns
by(all_data, all_data$PR, function(d) summary(lm(rmax ~ prot_log, data = d))) #ns
by(all_data, all_data$PR, function(d) summary(lm(e ~ prot_log, data = d))) #ns

#####MAYA DELETE THIS BS


chla_topt <- power_full_join(avg_chla, topt_df, by = "frag_ID", conflict = coalesce_xy) %>% 
  drop_na() %>%
  filter(PR != "Respiration")

#scatterplot of chla and topt params

#topt
chla_topt_scatter <- ggplot(filter(chla_topt, PR == "NetPhoto")) +
  geom_point(aes(y = topt, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  coord_transform(x = "log")+
  #facet_wrap(~PR, scales = "free")+
  geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(x = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       y = "Thermal optimum (°C)",
       color = "Species")
chla_topt_scatter

chla_topt_scatter <- chla_topt_scatter +
  #species-specific regression lines
  #geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, color = full_species),
  #            method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = topt, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black", linewidth = 1.1)

by(chla_topt, chla_topt$PR, function(d) summary(lm(topt ~ chla_log, data = d))) #*
library(performance)
check_model(lm(topt ~ chla_ug_cm2_mean, data = filter(chla_topt, PR == "NetPhoto")))
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(topt ~ chla_log, data = d)))

#topt net photo
#chla_ug_cm2_mean        0.5526     0.2563   2.156   0.0395 * 

#topt
#chla_topt$PR: NetPhoto
#chla_ug_cm2_mean    0.006035   0.002387   2.528   0.0172 *  significant
#none individual sig dif

ggsave(here("Output", "Physiology", "dryweight_topt_log_np.pdf"), chla_topt_scatter, h = 6, w = 10)

#rmax

chla_rmax_scatter <- ggplot(filter(chla_topt, PR == "NetPhoto")) +
  geom_point(aes(y = rmax, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  coord_transform(x = "log")+
  #facet_wrap(~PR, scales = "free")+
  geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(x = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       y = expression("Rmax" ~ (mu*mol ~ cm^{-2} ~ h^{-1})), 
       color = "Species")
chla_rmax_scatter

chla_rmax_scatter <- chla_rmax_scatter +
  #species-specific regression lines
  # geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, color = full_species),
  #             method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = rmax, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_rmax_log.pdf"), chla_rmax_scatter, h = 4, w = 12)

by(chla_topt, chla_topt$PR, function(d) summary(lm(rmax ~ chla_log, data = d))) #ns
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(rmax ~ chla_ug_cm2_mean, data = d)))
#none sig

#e
chla_e_scatter <- ggplot(filter(chla_topt, PR == "NetPhoto")) +
  geom_point(aes(y = e, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  coord_transform(x = "log")+
  #facet_wrap(~PR, scales = "free")+
  #geom_smooth(aes(y = e, x = chla_ug_cm2_mean, group = 1),
  #            method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(x = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       y = ("E (eV)"),  
       color = "Species")
chla_e_scatter

chla_e_scatter <- chla_e_scatter +
  #species-specific regression lines
  #geom_smooth(aes(y = e, x = chla_ug_cm2_mean, color = full_species),
  #            method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = e, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = FALSE, color = "black", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_e_log.pdf"), chla_e_scatter, h = 4, w = 12)

by(chla_topt, chla_topt$PR, function(d) summary(lm(e ~ chla_log, data = d)))
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(e ~ chla_ug_cm2_mean, data = d)))

#e chla_topt$PR: GrossPhoto
#chla_ug_cm2_mean      -0.16421    0.07198  -2.281  0.03307 * 

#e
#chla_topt$PR: GrossPhoto
#chla_ug_cm2_mean   -0.0019061  0.0007072  -2.695   0.0136 *  
#none significant


#breadth
chla_breadth_scatter <- ggplot(filter(chla_topt, PR == "NetPhoto")) +
  geom_point(aes(y = breadth, x = chla_ug_cm2_mean, color = full_species), alpha = 0.5) +
  scale_color_manual(values = sp_cols, labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  theme_bw(base_size = 22) +
  coord_transform(x = "log")+
  #facet_wrap(~PR, scales = "free")+
  geom_smooth(aes(y = breadth, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1) +
  labs(x = expression("Dry weight" ~ (mg ~ cm^{-2})), 
       y = expression("Breadth (°C)"), 
       color = "Species")
chla_breadth_scatter

chla_breadth_scatter <- chla_breadth_scatter +
  #species-specific regression lines
  #geom_smooth(aes(y = breadth, x = chla_ug_cm2_mean, color = full_species),
  #            method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line within each facet
  geom_smooth(aes(y = breadth, x = chla_ug_cm2_mean, group = 1),
              method = "lm", se = TRUE, color = "black", linewidth = 1.1)

ggsave(here("Output", "Physiology", "dryweight_breadth_log.pdf"), chla_breadth_scatter, h = 4, w = 12)

by(chla_topt, chla_topt$PR, function(d) summary(lm(breadth ~ chla_log, data = d))) #ns
by(chla_topt, list(chla_topt$PR, chla_topt$full_species),
   function(d) if(nrow(d) > 1) summary(lm(breadth ~ chla_ug_cm2_mean, data = d)))

#breadth chla_topt$PR: NetPhoto
#chla_ug_cm2_mean        0.6504     0.2504   2.598   0.0146 *  

#breadth
#chla_topt$PR: NetPhoto
#chla_ug_cm2_mean   0.005614   0.002442   2.299   0.0289 * 

#breadth
#NetPhoto
#Turbinaria frondens
#chla_ug_cm2_mean   0.035475   0.008652    4.10   0.0262 *


library(ggpubr)
chla_topt_log <- ggarrange(chla_rmax_scatter, chla_topt_scatter, chla_e_scatter, chla_breadth_scatter, 
                         nrow = 2, ncol = 2, legend = "right", common.legend = TRUE)
chla_topt_log

ggsave(here("Output","TPC","Graphs","chla_topt_log_np.pdf"), chla_topt_log, h = 8, w = 12)

