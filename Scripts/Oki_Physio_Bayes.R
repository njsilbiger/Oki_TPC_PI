#####Okinawa Bulk Coral Physiology 2025
###Bayesian analysis
##Maya Powell
#Created November 3rd 2026

######Setup

#Load in libraries

#install packages as needed
library(tidyverse)
library(here)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(brms)
library(tidybayes)
library(vegan)
library(pairwiseAdonis)

#### Initial data processing ####

# Initial tissue biomass calculations

#code to calculate average dry weight per fragment
#needed initially but not after writing average dry weight datasheet
#read in metadata
# phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta_all.csv"))
# 
# DW <- read.csv(here("Data", "Physiology", "afdw_sesoko.csv"))
# 
# DW <- DW %>% left_join(phys_meta, by = "frag_ID")
# 
# #calculate dry weight & ash free dry weight
# DW <- DW %>%
#   mutate(dw_g_mL = dry_pan - pan_weight,
#          dw_total_g = dw_g_mL*blastate_vol_mL,
#          dw_g_cm2 = dw_total_g/SA_cm2,
#          dw_mg_cm2 = dw_g_cm2*1000,
#          afdw_g_mL = dry_pan - ash_pan,
#          afdw_total_g = afdw_g_mL*blastate_vol_mL,
#          afdw_g_cm2 = afdw_total_g/SA_cm2,
#          afdw_mg_cm2 = afdw_g_cm2*1000)
# 
# #average across the reps
# #use all 6 reps across two years for dry weight
# avg_DW <- DW %>%
#   group_by(frag_ID) %>%
#   dplyr::select(-ash_pan,-AFDW,-afdw_g_mL,-afdw_total_g,-afdw_g_cm2,-afdw_mg_cm2) %>%
#   summarise(dw_g_mL = mean(dw_g_mL),
#             dw_total_g = mean(dw_total_g),
#             dw_g_cm2 = mean(dw_g_cm2),
#             dw_mg_cm2 = mean(dw_mg_cm2))
# 
# #use only the 3 reps over one year for AFDW
# avg_AFDW <- DW %>%
#   drop_na(afdw_g_cm2) %>%
#   group_by(frag_ID) %>%
#   summarise(afdw_g_mL = mean(afdw_g_mL),
#             afdw_total_g = mean(afdw_total_g),
#             afdw_g_cm2 = mean(afdw_g_cm2),
#             afdw_mg_cm2 = mean(afdw_mg_cm2)) %>%
#   left_join(avg_DW, by = "frag_ID") %>%
#   left_join(phys_meta, by = "frag_ID")
# 
# write.csv(avg_AFDW, here("Data", "Physiology", "Average_Ash_Free_Dry_Weight.csv"), row.names = FALSE)

## Initial chlorophyll a processing

#read in data
#data is corrected for blank rates and surface area and blanks are removed
# #uncomment this section if you need to redo prep for chla dataframe
# chla_plate1 <- read.csv(here("Data", "ChlaSpec", "oki_chla_plate1_summary.csv"))
# chla_plate2 <- read.csv(here("Data", "ChlaSpec", "oki_chla_plate2_summary.csv"))
# chla_plate3 <- read.csv(here("Data", "ChlaSpec", "oki_chla_plate3_summary.csv"))
# chla_plate4 <- read.csv(here("Data", "ChlaSpec", "oki_chla_plate4_summary.csv"))
# 
# # concatenate all plates
# chla_all_corr <- bind_rows(chla_plate1,chla_plate2,chla_plate3,chla_plate4)
# chla_all_corr <- chla_all_corr %>% rename(frag_ID = Sample.Name)
# 
# #add metadata, SA, and symbiont density to calc chla/sym
# SA_frag_final <- read.csv(here("Data","SurfaceArea","SA_fragments_final.csv"))
# chla_meta <- SA_frag_final %>% filter(Frag_type == "physio") %>%
#   select(frag_ID, full_species, species, calc_SA)
# avg_sym <- read.csv(here("Data", "Physiology", "Average_Sym_Density.csv"))
# sym_select <- avg_sym %>%
#   select(frag_ID, sym_total)
# 
# chla_data <- left_join(chla_all_corr, chla_meta, by = "frag_ID") %>%
#   rename(SA_cm2 = calc_SA)
# chla_data <- left_join(chla_data, sym_select, by = "frag_ID") #only have 7 species here for sym so will have to subset
# 
# #save data
# write_csv(chla_data, here("Data", "Physiology", "Chla_avg.csv"))
#
#just did this once to calculate and remove outliers, can always re-calculate above
# chla_avg <- chla_avg %>%
#   #filter(full_species %in% sp_keep) %>%
#   filter(frag_ID != "E04", frag_ID != "E05") %>% #remove outlier frag
#   mutate(chla_ug_total = (chla_mean*50)) %>% #multiply by slurry volume to get total chla
#   mutate(chla_pg_sym = (chla_ug_total/sym_total)*1000000) #divide by total symbionts and multiply by 1000000 to get pg
# write_csv(chla_avg, here("Data", "Physiology", "Chla_avg.csv"))

## Needed initially to calculate avg per frag sym density

#code to calculate average Symbiont Density per fragment
#needed initially but not after writing average Symbiont Density datasheet
# phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta_all.csv"))
# 
# #read in sym photo id sheet and sym counts
# counts <- read.csv(here("Data", "Physiology", "sym_counts_all_sept9.csv"))
# photos <- read.csv(here("Data", "Physiology", "sym_photos_all.csv"))
# 
# #combine
# sym_raw <- photos %>% left_join(counts, by = "image_id")
# 
# #add together counts from photos from 2 focus settings for photos
# #img_rm <- c("P8184974","P8184975","P8184976","P8174396","P8174397","P8174662","P8174663","P8174668","P8174669")
# sym_add <- sym_raw %>%
#   group_by(frag_ID,rep) %>% mutate(count_2 = sum(count)) %>%
#   filter(image_no == 1) #keep only 1 row of data for each
#   #filter(!image_id %in% img_rm) #filter out images to remove based on SE >15% of mean (as calculated below)
# 
# #average sym counts across 6 reps for each id
# #and calculate SE
# sym_avg <- sym_add %>%
#   group_by(frag_ID) %>% mutate(count_avg = mean(count_2)) %>%
#   mutate(count_se = sd(count_2)/sqrt(length(count_2))) %>%
#   #check to see which have highest error and look at add data for those >15%
#   mutate(count_var = count_se/count_avg) %>%
#   filter(rep == 6) #keep only one row of data for each
# 
# #error >0.15 = H01, A04, E04, F07, H08
# #all error less than 0.2, just keeping all!
# #H01 - two counts ~200, others ~250
# #A04 - one count 133, others ~80-100
# #E04 - one count 5, others ~15-20
# #F07 - two counts ~80, all others ~30
# #H08 - two counts in the 300s while all others in the 100s
# 
# #now calculate actual symbiont density
# #the 1x1mm grid should be 0.1uL
# #we want cells/mL so 10,000 x our counts (1000uL in a mL/0.1uL = x 10,000)
# #homogenate volume for all is 50mL
# #SA is calculated as well here
# #combine with metadata for calculations
# 
# sym <- sym_avg %>%
#   select(frag_ID, count_avg) %>% #keep only relevant columns
#   left_join(phys_meta, by = "frag_ID") %>% #add metadata
#   mutate(sym_mL = count_avg*10000) %>%
#   mutate(sym_total = sym_mL*blastate_vol_mL) %>%
#   mutate(sym_cm2 = sym_total/SA_cm2)
# 
# write.csv(sym, here("Data", "Physiology", "Average_Sym_Density.csv"), row.names = FALSE)

## Concatenate and save all data

#generated below - only need once
# avg_AFDW <- read.csv(here("Data", "Physiology", "Average_Ash_Free_Dry_Weight.csv")) #afdw_mg_cm2
# avg_AFDW <- avg_AFDW %>% mutate(afdw_log = log(afdw_mg_cm2)) %>%
#   mutate(dw_log = log(dw_mg_cm2)) %>%
#   dplyr::select(frag_ID, dw_mg_cm2, dw_log, afdw_mg_cm2,afdw_log)
# chla_avg <- read.csv(here("Data", "Physiology", "Chla_avg.csv")) #chla_ug_cm2_mean
# chla_avg <- chla_avg %>% mutate(chla_log = log(chla_ug_cm2_mean)) %>%
#   mutate(chla_sym_log = log(chla_pg_sym)) %>%
#   dplyr::select(frag_ID, chla_ug_cm2_mean,chla_log, chla_pg_sym, chla_sym_log)
# avg_sym <- read.csv(here("Data", "Physiology", "Average_Sym_Density.csv"))
# avg_sym <- avg_sym %>% filter(frag_ID != "C07") %>% filter(frag_ID != "D10") #remove crazy outliers!
# avg_sym <- avg_sym %>% mutate(sym_log = log(sym_cm2)) %>% dplyr::select(frag_ID, sym_cm2,sym_log)
# prot <- read.csv(here("Data", "Physiology", "protein_all_summary.csv")) #prot_ug_cm2
# prot <- prot %>% mutate(prot_log = log(prot_ug_cm2)) %>% dplyr::select(frag_ID, prot_ug_cm2,prot_log)
# physio_list <- list(avg_AFDW,chla_avg,avg_sym,prot)
# all_physio <- physio_list %>% reduce(left_join)
# write_csv(all_physio, here("Data", "Physiology", "all_physio_data.csv"))
#all_physio <- read_csv(here("Data", "Physiology", "all_physio_data.csv"))


#### Read in data ####

all_physio <- read.csv(here("Data", "Physiology", "all_physio_data.csv"))
all_physio <- all_physio |> mutate(sym_cm2_106 = sym_cm2/1000000,
                                   sym_cm2_106_log = log(sym_cm2_106)) 
phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta_all.csv"))

physio <- phys_meta %>% left_join(all_physio, by = "frag_ID")

### Initial data checks ###

#names(all_physio)
metrics_all <- physio |> 
  select(afdw_mg_cm2,chla_ug_cm2_mean,chla_pg_sym,sym_cm2_106,prot_ug_cm2,
         afdw_log,chla_log,chla_sym_log,sym_cm2_106_log,prot_log) |> 
  names()

hist <- physio |>
  pivot_longer(all_of(metrics_all)) |>
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free")
hist

qq <- physio |>
  pivot_longer(all_of(metrics_all)) |>
  ggplot(aes(sample = value)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~name, scales = "free")
qq

norm_tests <- physio |>
  pivot_longer(all_of(metrics_all)) |>
  group_by(name) |>
  summarise(p = shapiro.test(value)$p.value)
norm_tests

#use these as the final metrics to assess for physio
metrics_final <- physio |> 
  select(sym_cm2_106,afdw_log,chla_log,chla_sym_log,prot_log) |> 
  names()
#"sym_cm2_106","afdw_log","chla_log","chla_sym_log","sym_cm2_106_log","prot_log" 

### Bayesian model setup ####
seed <- 8
out_dir <- here("Data/Bayes_Physio")
mod_dir <- file.path(out_dir, "models")
dir.create(mod_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic(base_size = 12))

species_colors <- c(
  "Acropora hyacinthus" = '#d8aedd',
  "Echinopora lamellosa" = '#ba7999',
  "Favites complanata"   = '#dd4124',
  "Montipora aequituberculata" = '#ed8b00',
  "Montipora vietnamensis" = '#efbc82',
  "Pachyseris rugosa" = '#edd746',
  "Pocillopora grandis" = '#d0e2af',
  "Porites cylindrica" = '#45681e',
  "Porites rus" = '#7bbcd5',
  "Turbinaria frondens" = '#00496f'
)

#### Bayesian models ####

#create bayesian linear models
#priors for gaussian:
#Normal(0, 2.5 * SD(y))   (differences in original units)
d <- physio |> mutate(full_species = factor(full_species))

metrics_final <- physio |> 
  select(sym_cm2_106,afdw_log,chla_log,chla_sym_log,prot_log) |> 
  names()

metric_data <- function(m) {
  d |>
    transmute(full_species, y = .data[[m]]) |>
    filter(!is.na(y)) |>
    droplevels()
}

# Model formula (same for every metric)
species_mod <- bf(y ~ full_species, family = gaussian())

fit_model <- function(m, prior_only = FALSE) {
  dm <- metric_data(m)
  b_prior <- prior_string(sprintf("normal(0, %s)", signif(2.5 * sd(dm$y), 3)),
                          class = "b")
  brm(
    species_mod,
    data         = dm,
    prior        = b_prior,
    sample_prior = if (prior_only) "only" else "no",
    chains = 4, iter = 2000, warmup = 1000, cores = 4,
    control = list(adapt_delta = 0.95),
    seed = seed,
    file = file.path(mod_dir, paste0(m, "_gaussian", if (prior_only) "_prioronly" else "")),
    file_refit = "on_change"
  )
}

#fit all metrics
##THIS STEP RUNS THE MODELS AND TAKES A WHILE (~5min)
fits <- map(set_names(metrics_final), fit_model)

#model summaries
map(fits, summary)

# Posterior draws of each species' expected value
get_epred <- function(fit, m) {
  dm   <- metric_data(m)
  newd <- tibble(full_species = factor(levels(dm$full_species),
                                       levels = levels(dm$full_species)))
  fit |>
    epred_draws(newdata = newd) |>
    ungroup() |>
    mutate(metric = m)
}

# All pairwise contrasts between species
get_contrasts <- function(fit, m) {
  out <- get_epred(fit, m) |>
    compare_levels(.epred, by = full_species, comparison = "pairwise") |>
    ungroup() |>
    group_by(contrast = full_species) |>
    summarise(
      median = median(.epred),
      lower  = quantile(.epred, 0.025),
      upper  = quantile(.epred, 0.975),
      pd     = max(mean(.epred > 0), mean(.epred < 0)),  # probability of direction
      .groups = "drop"
    ) |>
    mutate(metric = m,
           excludes_zero = sign(lower) == sign(upper),
           .before = 1)
  
  # For log-transformed metrics, also give the ratio between species
  # (exp of a log difference = ratio of geometric means)
  if (str_detect(m, "_log$")) {
    out <- out |>
      mutate(ratio_median = exp(median),
             ratio_lower  = exp(lower),
             ratio_upper  = exp(upper))
  }
  out
}

AllContrasts <- imap_dfr(fits, get_contrasts)
#AllContrasts

write_csv(AllContrasts, here("Output/Physio_Bayes", "species_contrasts.csv"))

# Species estimates (posterior median + 95% CrI) - useful for plotting later
AllEpred <- imap_dfr(fits, get_epred)
SpeciesEst <- AllEpred |>
  group_by(metric, full_species) |>
  median_qi(.epred)

#### Posterior predictive (pp!) checks by species ####
pp_species <- function(fit, m) {
  p1 <- pp_check(fit, ndraws = 100) +
    ggtitle(paste(m))
  p2 <- pp_check(fit, type = "dens_overlay_grouped",
                 group = "full_species", ndraws = 50) +
    ggtitle("Density by species")
  p <- p1 / p2 + plot_layout(heights = c(1, 2))
  ggsave(here("Output/Physio_Bayes", paste0(m, "_ppcheck.png")), p, width = 12, height = 10)
  p
}

pp_plots <- imap(fits, pp_species)
#pp_plots$afdw_log

#### Species plots - raw data + posterior estimates ####

#labels
metric_labels <- c(
  afdw_mg_cm2      = expression("Biomass" ~ (mg ~ cm^{-2})),
  chla_ug_cm2_mean = expression("Chlorophyll a" ~ (mu*g ~ cm^{-2})),
  chla_pg_sym      = expression("Chlorophyll a" ~ (pg ~ symbiont^{-2})),
  prot_ug_cm2      = expression("Protein" ~ (mu*g ~ cm^{-2})),
  sym_cm2_106      = expression("Symbionts" ~ (cells ~ cm^{-2} ~ x ~ 10^{6}))
)

# Raw (untransformed) column for each modelled metric, used for points + labels
raw_cols <- c(
  afdw_log     = "afdw_mg_cm2",
  chla_log     = "chla_ug_cm2_mean",
  chla_sym_log = "chla_pg_sym",
  prot_log     = "prot_ug_cm2",
  sym_cm2_106  = "sym_cm2_106"
)

plot_species <- function(m) {
  is_log <- str_detect(m, "_log$")
  raw    <- raw_cols[[m]]
  
  # Posterior draws of species means, back-transformed for log metrics
  # (exp of the mean log = geometric mean, in original units)
  ep <- AllEpred |>
    filter(metric == m) |>
    mutate(est = if (is_log) exp(.epred) else .epred)
  
  # Order species by posterior medians
  sp_order <- ep |>
    group_by(full_species) |>
    summarise(med = median(est), .groups = "drop") |>
    arrange(med) |>
    pull(full_species) |>
    as.character()
  
  ep <- ep |> mutate(full_species = factor(full_species, levels = sp_order))
  
  # Raw data in original units (same rows that went into the model)
  raw_d <- d |>
    filter(!is.na(.data[[m]])) |>
    transmute(full_species = factor(full_species, levels = sp_order),
              y = .data[[raw]])
  
  p <- ggplot() +
    geom_jitter(data = raw_d,
                aes(x = full_species, y = y, color = full_species),
                width = 0.12, height = 0, size = 2.5, alpha = 0.8) +
    stat_halfeye(data = ep,
                 aes(x = full_species, y = est, fill = full_species),
                 position = position_nudge(x = 0.25), width = 0.5,
                 .width = c(0.50, 0.95), point_interval = median_qi,
                 slab_alpha = 0.5, color = "black", point_size = 1.8,
                 show.legend = FALSE) +
    scale_fill_manual(values = species_colors, guide = "none") +
    scale_color_manual(values = species_colors, breaks = sp_order, name = "Species") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = NULL, y = metric_labels[[raw]]) +
    theme_bw(base_size = 14) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(face = "italic"))
  
  if (is_log) {
    p <- p + scale_y_log10(breaks = scales::breaks_log(n = 6),
                           labels = scales::label_comma())
  }
  p
}

species_plots <- map(set_names(metrics_final), plot_species)

# Save each plot
iwalk(species_plots, \(p, m)
      ggsave(here("Output/Physio_Bayes", paste0(m, "_species_bayes.png")), p,
             width = 8, height = 4.5))

#### Ordination plot of data ####

#ordination of physio data
#load physio data and generate dataframe with id as rownames
all_physio_data <- all_physio %>% #read in - chunk 1 if needed
  dplyr::select(frag_ID,afdw_log,chla_log,chla_sym_log,sym_cm2_106,prot_log) %>% 
  column_to_rownames("frag_ID") %>% drop_na()
#make them matrices
all_physio_mat <- as.matrix(all_physio_data)
#quick glance at data to look for strong patterns
#pairs(x = all_physio_mat, gap = 0, cex.labels = 0.5) #look similar even between metrics
#scale across variable types
all_physio_mat <- scale(all_physio_mat)
#generate pca data
pca_all_physio <- prcomp(all_physio_mat)
#collapse PCA data
pc_axes_physio <- as.data.frame(pca_all_physio$x)
#add frag_ID back
pc_axes_physio$frag_ID <- rownames(pc_axes_physio) 
#put it back with metadata
pca_physio <- pc_axes_physio %>% left_join(phys_meta, by = "frag_ID")

#add env fit arrows
physio_fit <- envfit(pca_all_physio$x[, c("PC1", "PC2")], all_physio_data, permutations = 999, na.rm = TRUE)

physio_scores <- as.data.frame(scores(physio_fit, display = "vectors")) %>% mutate(variable = rownames(.))

pca_physio_spp_arrows <- ggplot(pca_physio, aes(x = PC1, y = PC2, color = full_species, fill = full_species)) +
  #geom_point(size = 3, alpha = 0.1) +
  scale_color_manual(values = species_colors, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  scale_fill_manual(values = species_colors, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  #stat_ellipse(geom = "polygon", alpha = 0.1, linewidth = 0.2)+
  stat_ellipse(level = 0.0001, geom = "point", shape = 21) +
  #stat_ellipse(level = 0.95, alpha = 0, linewidth = 0)+
  #stat_ellipse(aes(fill = group), level = 0.95, alpha = 0.5, color = "black", linewidth = 0.8) +
  geom_segment(data = physio_scores,
               aes(x = 0, y = 0, xend = PC1*1, yend = PC2*1),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  #labs(title = "Physiology") +
  geom_text(data = physio_scores,
            aes(x = PC1 * 1.2, y = PC2 * 1.2, label = variable),
            color = "black", size = 5, inherit.aes = FALSE) +
  theme_classic(base_size = 22)

pca_physio_spp_arrows

#ggsave(here("Output", "Physiology", "physio_params_spp_arrows_nolabels.pdf"), pca_physio_spp_arrows, h = 8, w = 12)

### Dispersion and PERMANOVA stats ###

# generate distance matrix
physio_dist <- vegdist(pca_all_physio$x, method = "euclidean")

#beta dispersion comparisons
bet.phys <- betadisper(physio_dist,pca_physio$species)
anova(bet.phys) 
#species: p=0.02419 *

#plot(bet.phys)
#permutest(bet.phys, pairwise = TRUE, permutations = 999)

physio_perm <- adonis2(physio_dist ~ species, data = pca_physio, permutations = 999)
physio_perm
#species: R2 = 0.41009, F= 6.3339,p=0.001 ***

pw_comps_sp <- pairwise.adonis2(physio_dist ~ species, data=pca_physio, permutations = 999)

# table of pairwise comparison output
pw_sp_list <- pw_comps_sp[names(pw_comps_sp) != "parent_call"]
pw_sp_table <- do.call(rbind, lapply(names(pw_sp_list), function(nm) {
  df <- as.data.frame(pw_sp_list[[nm]])
  df$Term <- rownames(df)
  df$comparison <- nm
  df[, c("comparison", "Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
}))

write.csv(pw_sp_table, here("Output/Physiology/PERMANOVA_ord_pairwise_comparison_physio.csv"))

#### All plots together ###

physio_plots_bayes <- ggarrange(species_plots$chla_log, species_plots$chla_sym_log, species_plots$sym_cm2_106, species_plots$afdw_log, species_plots$prot_log,pca_physio_spp_arrows,
                                ncol = 3, nrow = 2, legend = "right", common.legend = TRUE,
                                labels = c("A", "B", "C", "D", "E", "F"), 
                                font.label = list(size = 20, color = "black"))
physio_plots_bayes

ggsave(here("Output/Physio_Bayes","physio_plots_bayes.png"), physio_plots_bayes, h = 8, w = 16)
