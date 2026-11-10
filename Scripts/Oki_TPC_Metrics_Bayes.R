#####Okinawa Bulk Coral Physiology 2025
###Bayesian analysis
##Maya Powell
#Created November 3rd 2026

###### Setup ####

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

#### Read in data ####

all_tpc <- read.csv(here("Data/RespoFiles/TPC/bayes_params_no4.csv"))
all_tpc$full_species <- gsub('eydouxi', 'grandis', all_tpc$full_species)
all_tpc <- all_tpc |> filter(frag_ID != "B08")
#in looking through some of the models, rhat is much larger than 1 for quite a few - issue?? need to check with Nyssa

# Select data we want to examine and make separate info for each resp metric
r_tpc <- all_tpc |> filter(PR == "Respiration") |> 
  select(frag_ID,T_opt_median,e_median,rmax_median, r_tref_median) |> 
  rename(topt_r = T_opt_median, e_r = e_median, rmax_r = rmax_median, r_tref_r = r_tref_median) 
np_tpc <- all_tpc |> filter(PR == "NetPhoto") |> 
  select(frag_ID,T_opt_median,e_median,rmax_median, r_tref_median) |> 
  rename(topt_np = T_opt_median, e_np = e_median, rmax_np = rmax_median, r_tref_np = r_tref_median)
gp_tpc <- all_tpc |> filter(PR == "GrossPhoto") |> 
  select(frag_ID,T_opt_median,e_median,rmax_median, r_tref_median) |> 
  rename(topt_gp = T_opt_median, e_gp = e_median, rmax_gp = rmax_median, r_tref_gp = r_tref_median)

# Join all together with physio metadata
phys_meta <- read.csv(here("Data", "Physiology", "Physio_meta_all.csv"))

tpc <- gp_tpc %>%
  left_join(np_tpc, by = "frag_ID") |> 
  left_join(r_tpc, by = "frag_ID") |> 
  left_join(phys_meta, by = "frag_ID")

#tref is 27°C as set in Bayes_TPC script
#use this to get P:R ratios
tpc <- tpc |> 
  mutate(pr_gp = r_tref_gp/r_tref_r,
         pr_np = r_tref_np/r_tref_r)
tpc$pr_gp <- replace(tpc$pr_gp, tpc$pr_gp>5, NA) #replace crazy outliers
tpc$topt_np <- replace(tpc$topt_np, tpc$topt_np<25, NA)

### Initial data checks ###

#names(all_physio)
metrics_all <- tpc |> 
  select(topt_gp, topt_np, topt_r, e_gp, e_np, e_r, rmax_gp, rmax_np, rmax_r, pr_gp, pr_np) |> 
  names()

hist <- tpc |>
  pivot_longer(all_of(metrics_all)) |>
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free")
hist

qq <- tpc |>
  pivot_longer(all_of(metrics_all)) |>
  ggplot(aes(sample = value)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~name, scales = "free")
qq

norm_tests <- tpc |>
  pivot_longer(all_of(metrics_all)) |>
  group_by(name) |>
  summarise(p = shapiro.test(value)$p.value)
norm_tests
#all normal except r topt is very not normal!!
#topt_r  0.00492
#TRANSFORM R MAYBE SQRT

#use these as the final metrics to assess for physio
#"topt_gp" "topt_np" "topt_r"  "e_gp"    "e_np"    "e_r"     "rmax_gp" "rmax_np" "rmax_r" "pr_gp" "pr_np"

### Bayesian model setup ####
seed <- 8
out_dir <- here("Data/Bayes_TPC_Metrics")
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
d <- tpc |> mutate(full_species = factor(full_species))

metrics_final <- tpc |> 
  select(topt_gp, topt_np, topt_r, e_gp, e_np, e_r, rmax_gp, rmax_np, rmax_r, pr_gp, pr_np) |> 
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

write_csv(AllContrasts, here("Output/TPC_Metrics_Bayes", "species_contrasts.csv"))

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
  ggsave(here("Output/TPC_Metrics_Bayes", paste0(m, "_ppcheck.png")), p, width = 12, height = 10)
  p
}

pp_plots <- imap(fits, pp_species)
#pp_plots$afdw_log

#### Species plots - raw data + posterior estimates ####

#labels
metric_labels <- c(
  topt_gp = "Thermal optimum (°C)",
  topt_np = "Thermal optimum (°C)",
  topt_r  = "Thermal optimum (°C)",
  e_gp = expression("Energy" ~ (joules ~ mol^{-1})),
  e_np = expression("Energy" ~ (joules ~ mol^{-1})),
  e_r  = expression("Energy" ~ (joules ~ mol^{-1})),
  rmax_gp = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
  rmax_np = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
  rmax_r  = expression("Rate Max" ~ (mu*mol ~ cm^{-2} ~ hr^{-1})),
  pr_gp = "GP:R Ratio at 27°C",
  pr_np = "NP:R Ratio at 27°C"
)

# Raw (untransformed) column for each modelled metric, used for points + labels
raw_cols <- c(
  topt_gp = "topt_gp",
  topt_np = "topt_np",
  topt_r  = "topt_r",
  e_gp = "e_gp",
  e_np = "e_np",
  e_r  = "e_r",
  rmax_gp = "rmax_gp",
  rmax_np = "rmax_np",
  rmax_r  = "rmax_r",
  pr_gp = "pr_gp",
  pr_np = "pr_np")

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
                 .width = c(0.5, 0.95), point_interval = median_qi,
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
# iwalk(species_plots, \(p, m)
#       ggsave(here("Output/TPC_Metrics_Bayes", paste0(m, "_species_bayes.png")), p,
#              width = 8, height = 4.5))

#### Ordination plot of data ####

#load metadata

#separate out data, make sure to put frag_ID as rownames so you can re-join with metadata later
topt_r_data <- tpc %>% select(frag_ID, topt_r, e_r, rmax_r) %>% column_to_rownames("frag_ID")
topt_np_data <- tpc %>% select(frag_ID, topt_np, e_np, rmax_np, pr_np) %>% 
  drop_na() |> 
  column_to_rownames("frag_ID")
topt_gp_data <- tpc %>% select(frag_ID, topt_gp, e_gp, rmax_gp, pr_gp) %>% 
  drop_na() |> 
  column_to_rownames("frag_ID")
#make them matrices
topt_r <- as.matrix(topt_r_data)
topt_np <- as.matrix(topt_np_data)
topt_gp <- as.matrix(topt_gp_data)
#quick glance at data to look for strong patterns
pairs(x = topt_r, gap = 0, cex.labels = 0.5) #look similar even between metrics
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

topt_fit_r <- envfit(pca_topt_r$x[, c("PC1", "PC2")], topt_r_data, permutations = 999, na.rm = TRUE)
topt_fit_np <- envfit(pca_topt_np$x[, c("PC1", "PC2")], topt_np_data, permutations = 999, na.rm = TRUE)
topt_fit_gp <- envfit(pca_topt_gp$x[, c("PC1", "PC2")], topt_gp_data, permutations = 999, na.rm = TRUE)

#topt_fit #can look at p-values of data to see what is significant in driving differences

topt_scores_r <- as.data.frame(scores(topt_fit_r, display = "vectors")) %>% mutate(variable = rownames(.))
topt_scores_np <- as.data.frame(scores(topt_fit_np, display = "vectors")) %>% mutate(variable = rownames(.))
topt_scores_gp <- as.data.frame(scores(topt_fit_gp, display = "vectors")) %>% mutate(variable = rownames(.))

gp_pca_spp_arrows <- ggplot(pca_gp, aes(x = PC1, y = PC2, color = full_species, fill = full_species)) +
  scale_color_manual(values = species_colors, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  scale_fill_manual(values = species_colors, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  stat_ellipse(level = 0.00001, geom = "point", shape = 21) +
  #stat_ellipse(level = 0.95, alpha = 0, linewidth = 0)+
  geom_segment(data = topt_scores_gp,
               aes(x = 0, y = 0, xend = PC1*1, yend = PC2*1),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = topt_scores_gp,
            aes(x = PC1 * 1.2, y = PC2 * 1.2, label = variable),
            color = "black", size = 5, inherit.aes = FALSE) +
  theme_classic(base_size = 18)

gp_pca_spp_arrows

#ggsave(here("Output", "Physiology", "gp_pca_spp_arrows.png"), gp_pca_spp_arrows, h = 8, w = 12)

np_pca_spp_arrows <- ggplot(pca_np, aes(x = PC1, y = PC2, color = full_species, fill = full_species)) +
  scale_color_manual(values = species_colors, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  scale_fill_manual(values = species_colors, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  stat_ellipse(level = 0.00001, geom = "point", shape = 21) +
  #stat_ellipse(level = 0.95, alpha = 0, linewidth = 0)+
  geom_segment(data = topt_scores_np,
               aes(x = 0, y = 0, xend = PC1*1, yend = PC2*1),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = topt_scores_np,
            aes(x = PC1 * 1.2, y = PC2 * 1.2, label = variable),
            color = "black", size = 5, inherit.aes = FALSE) +
  theme_classic(base_size = 18)

np_pca_spp_arrows
#ggsave(here("Output", "Physiology", "np_pca_spp_arrows.png"), np_pca_spp_arrows, h = 8, w = 12)


r_pca_spp_arrows <- ggplot(pca_r, aes(x = PC1, y = PC2, color = full_species, fill = full_species)) +
  scale_color_manual(values = species_colors, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  scale_fill_manual(values = species_colors, name = "Species", labels = function(x) parse(text = paste0("italic('", gsub("'", "\\\\'", x), "')")))+
  stat_ellipse(level = 0.00001, geom = "point", shape = 21) +
  #stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 1)+
  geom_segment(data = topt_scores_r,
               aes(x = 0, y = 0, xend = PC1*1, yend = PC2*1),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = topt_scores_r,
            aes(x = PC1 * 1.2, y = PC2 * 1.2, label = variable),
            color = "black", size = 5, inherit.aes = FALSE) +
  theme_classic(base_size = 18)

r_pca_spp_arrows

#ggsave(here("Output", "Physiology", "r_pca_spp_arrows.png"), r_pca_spp_arrows, h = 8, w = 12)

### Dispersion and PERMANOVA stats ###

#generate distance matrix
np_dist <- vegdist(pca_topt_np$x, method = "euclidean")

set.seed(8)
#beta dispersion comparisons
bet.np <- betadisper(np_dist,pca_np$species)
anova(bet.np) 
#species: p= 0.7557
#plot(bet.phys)
#permutest(bet.phys, pairwise = TRUE, permutations = 999)

np_perm <- adonis2(np_dist ~ species, data = pca_np, permutations = 999)
np_perm
#species: R2 = 0.37422, F = 2.4585, p = 0.015 *
pw_comps_np <- pairwise.adonis2(np_dist ~ species, data=pca_np, permutations = 999)

#table of pairwise comparison output
pw_list_np <- pw_comps_np[names(pw_comps_np) != "parent_call"]
pw_table_np <- do.call(rbind, lapply(names(pw_list_np), function(nm) {
  df <- as.data.frame(pw_list_np[[nm]])
  df$Term <- rownames(df)
  df$comparison <- nm
  df[, c("comparison", "Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
}))

write.csv(pw_table_np, here("Output/Physiology/PERMANOVA_ord_pairwise_comparisons_np_tpc_params.csv"))

#gp
#generate distance matrix
gp_dist <- vegdist(pca_topt_gp$x, method = "euclidean")

#beta dispersion comparisons
bet.gp <- betadisper(gp_dist,pca_gp$species)
anova(bet.gp) 
#species: p= 0.3966

#plot(bet.phys)
#permutest(bet.phys, pairwise = TRUE, permutations = 999)

gp_perm <- adonis2(gp_dist ~ species, data = pca_gp, permutations = 999)
gp_perm
#species: 0.27322 1.5873  0.132

#pairwise.adonis2(gp_dist ~ morphology, data=pca_gp, permutations = 999)

#r
#generate distance matrix
r_dist <- vegdist(pca_topt_r$x, method = "euclidean")

#beta dispersion comparisons
bet.r <- betadisper(r_dist,pca_r$species)
anova(bet.r) 
#species: p= 0.6735

#plot(bet.phys)
#permutest(bet.phys, pairwise = TRUE, permutations = 999)

r_perm <- adonis2(r_dist ~ species, data = pca_r, permutations = 999)
r_perm
#0.21271 1.1708  0.302

#ggsave(here("Output", "Physiology", "physio_params_spp_arrows_nolabels.pdf"), pca_physio_spp_arrows, h = 8, w = 12)

#### All plots together ###
###EDIT THIS!!
tpc_schematic <- readRDS(here("Output/Okinawa_Map/tpc_schematic.rds"))

np_plots_bayes <- ggarrange(species_plots$topt_np, species_plots$rmax_np, species_plots$e_np, 
                            species_plots$pr_np, tpc_schematic,
                                ncol = 3, nrow = 2, legend = "right", common.legend = TRUE,
                                labels = c("A", "B", "C", "D", "E"), 
                                font.label = list(size = 20, color = "black"))
np_plots_bayes

ggsave(here("Output/TPC_Metrics_Bayes","np_plots_bayes.pdf"), np_plots_bayes, h = 8, w = 16)


gp_r_plots_bayes <- ggarrange(species_plots$topt_gp, species_plots$rmax_gp, species_plots$e_gp,gp_pca_spp_arrows,
                              species_plots$topt_r, species_plots$rmax_r, species_plots$e_r, r_pca_spp_arrows,
                            ncol = 4, nrow = 2, legend = "right", common.legend = TRUE,
                            labels = c("A", "B", "C", "D", "E","F","G","H"), 
                            font.label = list(size = 20, color = "black"))
gp_r_plots_bayes

ggsave(here("Output/TPC_Metrics_Bayes","gp_r_plots_bayes.pdf"), gp_r_plots_bayes, h = 8, w = 16)

all_tpc_metrics_bayes <- ggarrange(species_plots$topt_np, species_plots$rmax_np, species_plots$e_np, species_plots$pr_np,
                                 species_plots$topt_gp, species_plots$rmax_gp, species_plots$e_gp, species_plots$pr_gp,
                                 species_plots$topt_r, species_plots$rmax_r, species_plots$e_r, tpc_schematic,
                                 np_pca_spp_arrows, gp_pca_spp_arrows, r_pca_spp_arrows,
                                 ncol = 4, nrow = 4, legend = "right", common.legend = TRUE,
                                 labels = c("A", "B", "C", "D", "E","F","G","H","I","J","K","L","M","N","O"), 
                                 font.label = list(size = 20, color = "black"))
all_tpc_metrics_bayes

ggsave(here("Output/TPC_Metrics_Bayes","all_tpc_metrics_bayes.png"), all_tpc_metrics_bayes, h = 14, w = 20)


