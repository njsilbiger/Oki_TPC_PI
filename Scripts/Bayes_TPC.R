#Bayesian TPC analysis
#Maya Powell

#load libraries
library(here)
library(tidyverse)
library(broom)
library(purrr)
#install.packages("nimble")
library(nimble)
#library(remotes)
#remotes::install_github("johnwilliamsmithjr/bayesTPC")
library(bayesTPC)
library(coda)
library(ggpubr)

#Bayesian TPC package
#bayesTPC: Bayesian inference for thermal performance curves in R
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70004
#https://github.com/johnwilliamsmithjr/bayesTPC

# #look at info for shsch model
# get_default_model_specification("pawar_shsch")
# # bayesTPC Model Specification of Type: pawar_shsch
# # Model Formula:
# #   m[i] <- ( (e_h > e) * r_tref * exp((e/(8.62e-05)) * ((1/(T_ref + 273.15)) - (1/(Temp + 273.15))))/(1 + (e/(e_h - e)) * exp((e_h/(8.62e-05)) * (1/(T_opt + 273.15) - 1/(Temp + 273.15)))) )
# # Model Distribution:
# #   Trait[i] ~ T(dnorm(mean = m[i], tau = 1/sigma.sq), 0, )
# # Model Parameters and Priors:
# #   e ~ dunif(0, 1)
# # e_h ~ dunif(0, 30)
# # r_tref ~ dunif(0, 10)
# # T_opt ~ dunif(0, 50)
# # sigma.sq ~ dexp(1)
# # Model Constants:
# #   T_ref = 20

#read in data
df_clean <- read_csv(here("Data", "RespoFiles", "TPC", "PnR_clean_no4.csv"))

# set up values for the model
T_REF     <- 27 #use same as nls tref
NITER     <- 20000
BURN      <- 5000
NCHAINS   <- 3 # >1 allows gelman.diag()
new_temps <- c(24.5, 26, 27, 28, 29, 30, 31, 32, 34) # prediction temps

# Priors per rate type (anything not listed uses the bayesTPC default which is 0,50)
# literature values:
#30.5-33.5: Nyssa paper
#29.5-33.7 :Assessment of temperature optimum signatures of corals at both latitudinal extremes of the Red Sea
#10-35: Upper-mesophotic and shallow reef corals exhibit similar thermal tolerance, sensitivity and optima
#30.5 ± 1.8: Thermal extremes likely trigger metabolic imbalance in coral holobionts (meta analysis)
#32.56 - 36.91: Nutrient and sediment loading affect multiple facets of functionality in a tropical branching coral 
#34.2 - 34.9: Chronic low-level nutrient enrichment benefits coral thermal performance in a fore reef habitat
#33.8 ± 3.5: Latitudinal variation in thermal performance of the common coral Pocillopora spp

priors_by_PR <- list(
  NetPhoto    = list(T_opt = "dunif(22, 36)"), #using 22 to 36 since that's just outside the range of our TPCs
  GrossPhoto  = list(T_opt = "dunif(22, 36)"),
  Respiration = list(T_opt = "dunif(24, 40)") #look at literature and see what 
)

# Use the package's own model formula for predictions
tpc_formula <- get_formula("pawar_shsch")
eval_tpc <- function(draws, Temp, T_ref = T_REF) {
  eval(tpc_formula, envir = list(e = draws$e, e_h = draws$e_h,
                                 r_tref = draws$r_tref, T_opt = draws$T_opt,
                                 Temp = Temp, T_ref = T_ref))
}

# Median + 95% HPD of a vector of draws
summ_draws <- function(x) {
  hpd <- HPDinterval(as.mcmc(x), prob = 0.95)
  tibble(median = median(x), lower = hpd[1], upper = hpd[2])
}

fits_bayes   <- list() #fitted btpc objects, named "PR_species"
waic_bayes   <- list() #aic values
preds_bayes  <- list() #predicted values
params_bayes <- list() #parameters
draws_bayes  <- list() #full posterior draws, for between-species comparisons

#run for loop!!!!
#across different PR and species
#nicely keeps you updated with progress with outputs
#took about 10mins with 20,000 iterations across 3 parameters and 10 frag_ID
for (pr in unique(df_clean$PR)) {
  for (sp in unique(df_clean$frag_ID[df_clean$PR == pr])) {
    
    key <- paste(pr, sp, sep = "_")
    
    my_df <- df_clean |>
      filter(PR == pr, frag_ID == sp) |>
      mutate(Values = as.numeric(Values),
             temp_c_value = as.numeric(temp_c_value)) |>
      drop_na(Values, temp_c_value)
    
    # Likelihood is truncated at 0, so negative values can't be fit
    stopifnot("Negative values found" = all(my_df$Values >= 0)) #keep this and use below as needed
    # n_neg <- sum(my_df$Values < 0)
    # if (n_neg > 0) {
    #   message(sprintf("%s: dropping %d negative value(s)", key, n_neg))
    #   my_df <- my_df |> filter(Values >= 0)
    # }
    
    dat <- list(Trait = my_df$Values, Temp = my_df$temp_c_value)
    
    fit <- tryCatch(
      b_TPC(
        data      = dat,
        model     = "pawar_shsch",
        niter     = NITER,
        burn      = BURN,
        nchains   = NCHAINS,
        constants = list(T_ref = T_REF),
        inits     = list(e = 0.6, e_h = 3,
                         r_tref = min(median(dat$Trait), 9.9),  # inside dunif(0, 10)
                         T_opt = 30),
        priors    = priors_by_PR[[pr]],
        verbose   = FALSE
      ),
      error = function(err) {
        message(sprintf("Fit failed for %s: %s", key, conditionMessage(err)))
        NULL
      }
    )
    if (is.null(fit)) next
    
    fits_bayes[[key]] <- fit
    waic_bayes[[key]] <- tryCatch(get_WAIC(fit), error = function(err) NA)
    
    # Posterior draws (burn-in already removed by b_TPC; chains combined)
    d <- as_tibble(as.matrix(fit$samples))
    
    # Convergence: max Gelman-Rubin R-hat and min effective sample size
    rhat_max <- tryCatch(max(gelman.diag(fit$samples, multivariate = FALSE)$psrf[, 1]),
                         error = function(err) NA_real_)
    ess_min  <- min(effectiveSize(fit$samples))
    
    # Predictions with 95% HPD at your measurement temps
    preds_bayes[[key]] <- map_dfr(new_temps, function(t)
      summ_draws(eval_tpc(d, t)) |> mutate(temp_c_value = t)) |>
      transmute(PR = pr, frag_ID = sp, temp_c_value,
                .fitted = median, conf_lower = lower, conf_upper = upper)
    
    # T_opt is the peak temperature in this parameterization;
    # rmax = curve evaluated at each draw's own T_opt
    d <- d |> mutate(rmax = eval_tpc(d, d$T_opt))
    
    draws_bayes[[key]] <- d |> mutate(PR = pr, frag_ID = sp)
    
    params_bayes[[key]] <- d |>
      select(rmax, T_opt, e, e_h, r_tref, sigma.sq) |>
      pivot_longer(everything(), names_to = "param") |>
      group_by(param) |>
      summarise(summ_draws(value), .groups = "drop") |>
      mutate(PR = pr, frag_ID = sp, n = nrow(my_df),
             rhat_max = rhat_max, ess_min = ess_min)
  }
}

#issue I've run into: each call of b_TPC uses nimble to compile the model into C++
#then loads in r session as shared libraries (DLLs)
#R can only hold so many (~600) so you need to restart R before doing this
#takes ~40min-1hr to run through 150 samples (3 x 50 samples for np, gp, and r)

preds_bayes_df  <- bind_rows(preds_bayes)
params_bayes_df <- bind_rows(params_bayes)
draws_bayes_df  <- bind_rows(draws_bayes)

#Tpc parameters dataframe - wide version, one row per PR x species (like topt_df_sp)
params_bayes_wide <- params_bayes_df |>
  pivot_wider(names_from = param,
              values_from = c(median, lower, upper),
              names_glue = "{param}_{.value}")

# Flag fits that need attention
flagged_fits <- params_bayes_wide |> filter(rhat_max > 1.1 | ess_min < 400)

BioData <- read_csv(here("Data","RespoFiles","TPC","Fragment_Measurements_TPC.csv"))
species_meta <- BioData |> group_by(frag_ID) |> distinct(frag_ID, .keep_all = TRUE) |> select(frag_ID, full_species, species_ID)
params_bayes_wide <- params_bayes_wide |> left_join(species_meta, by = "frag_ID")

#save data
#thin draws so that the file is small enough for github (checked fits against larger dataset, still looks great)
draws_thin <- draws_bayes_df |>
  group_by(PR, frag_ID) |>
  slice(seq(1, n(), by = 5)) |>
  ungroup()
#save
saveRDS(draws_thin, here("Data","RespoFiles","TPC","bayes_draws_thin_no4.rds")) #using RDS bc it's huge and better for the MCMC output

#can also save as csv as needed
write_csv(params_bayes_wide, here("Data","RespoFiles","TPC","bayes_params_no4.csv"))

#diagnostics for all fits
# plot(fits_bayes[["GrossPhoto_Pcyl"]])
# traceplot(fits_bayes[["GrossPhoto_Pcyl"]])
# ppo_plot(fits_bayes[["GrossPhoto_Pcyl"]])   # is T_opt just reflecting the prior?
# summary(fits_bayes[["GrossPhoto_Pcyl"]])
# plot_prediction(fits_bayes[["GrossPhoto_Pcyl"]])

#read in data
PnR_clean <- read_csv(here("Data","RespoFiles","TPC","PnR_clean_no4.csv"))
PnR_clean <- PnR_clean |> left_join(species_meta, by = "frag_ID")

#draws_bayes_df <- readRDS(here("Data","RespoFiles","TPC","bayes_draws_no4.rds")) #too large
draws_thin <- readRDS(here("Data","RespoFiles","TPC","bayes_draws_thin_no4.rds"))
params_bayes_wide <- read_csv(here("Data","RespoFiles","TPC","bayes_params_no4.csv"))

#little metadata action to join
frag_lookup <- PnR_clean |> distinct(frag_ID, species, full_species)

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

#create smooth curves
T_REF       <- 27
tpc_formula <- get_formula("pawar_shsch")
eval_tpc <- function(draws, Temp, T_ref = T_REF) {
  eval(tpc_formula, envir = list(e = draws$e, e_h = draws$e_h,
                                 r_tref = draws$r_tref, T_opt = draws$T_opt,
                                 Temp = Temp, T_ref = T_ref))
}
summ_draws <- function(x) {
  hpd <- HPDinterval(as.mcmc(x), prob = 0.95)
  tibble(median = median(x), lower = hpd[1], upper = hpd[2])
}

plot_temps <- seq(min(PnR_clean$temp_c_value), max(PnR_clean$temp_c_value), by = 0.1)

curves_bayes <- draws_thin |>
  group_by(PR, frag_ID) |>
  #slice(seq(1, n(), by = 10)) |> # already thinned above - just use this
  group_modify(~ map_dfr(plot_temps, function(t)
    summ_draws(eval_tpc(.x, t)) |> mutate(temp_c_value = t))) |>
  ungroup() |>
  rename(.fitted = median, conf_lower = lower, conf_upper = upper) |>
  left_join(frag_lookup, by = "frag_ID")

params_plot <- params_bayes_wide

#plot all rates - function
plot_bayes_tpc <- function(pr, ylab, show_rmax = FALSE) {
  dat   <- PnR_clean   |> filter(PR == pr)
  preds <- curves_bayes |> filter(PR == pr)
  topt  <- params_plot  |> filter(PR == pr)
  
  p <- ggplot(dat, aes(x = temp_c_value, y = Values, color = full_species)) +
    geom_line(data = preds, aes(temp_c_value, .fitted, group = frag_ID),linewidth = 0.6) +
    # geom_ribbon(data = preds,
    #             aes(x = temp_c_value, ymin = conf_lower, ymax = conf_upper,
    #                 fill = full_species),
    #             inherit.aes = FALSE, alpha = 0.2) +
    # geom_line(data = preds, aes(temp_c_value, .fitted), linewidth = 0.7) +
    geom_point(alpha = 0.7, shape = 21) +
    scale_color_manual(values = sp_cols) +
    #scale_fill_manual(values = sp_cols) +
    facet_wrap(~ full_species, scales = "free", nrow = 2, ncol = 5) +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(face = "italic"), legend.position = "none") +
    labs(x = "Temperature (ºC)", y = ylab)
  
  if (show_rmax) {
    p <- p + geom_hline(data = topt, aes(yintercept = rmax_median),
                        linewidth = 0.3, color = "darkgreen")
  }
  p
}

gp_pred_plot <- plot_bayes_tpc("GrossPhoto", expression("GP Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))
np_pred_plot <- plot_bayes_tpc("NetPhoto", expression("NP Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))
resp_pred_plot <- plot_bayes_tpc("Respiration", expression("R Rate" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))

gp_pred_plot
np_pred_plot
resp_pred_plot

all_pred_plots <- ggarrange(np_pred_plot, gp_pred_plot, resp_pred_plot, 
                            nrow = 3, ncol = 1, labels = c("A", "B", "C"), 
                            font.label = list(size = 20, color = "black"))

ggsave(here("Output","TPC","Graphs","tpc_pred_all_gp_np_r_bayesian.pdf"), all_pred_plots, h = 12, w = 12)


##look at parameters

# species averages of the fragment-level estimates
params_sp_avg <- params_plot |>
  group_by(PR, full_species) |>
  summarise(across(c(T_opt_median, T_opt_lower, T_opt_upper,
                     rmax_median, rmax_lower, rmax_upper),
                   \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

# species order by NetPhoto species average, smallest to largest
topt_order <- params_sp_avg |> filter(PR == "NetPhoto") |> arrange(T_opt_median) |> pull(full_species)
rmax_order <- params_sp_avg |> filter(PR == "NetPhoto") |> arrange(rmax_median)  |> pull(full_species)

topt_compare_plot <- params_sp_avg |>
  mutate(full_species = factor(full_species, levels = topt_order)) |>
  ggplot(aes(x = full_species, y = T_opt_median, color = full_species)) +
  geom_pointrange(aes(ymin = T_opt_lower, ymax = T_opt_upper)) +
  scale_color_manual(values = sp_cols) +
  facet_wrap(~ PR, ncol = 1) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = NULL, y = "Topt (ºC), mean median and mean 95% HPD")

rmax_compare_plot <- params_sp_avg |>
  mutate(full_species = factor(full_species, levels = rmax_order)) |>
  ggplot(aes(x = full_species, y = rmax_median, color = full_species)) +
  geom_pointrange(aes(ymin = rmax_lower, ymax = rmax_upper)) +
  scale_color_manual(values = sp_cols) +
  facet_wrap(~ PR, ncol = 1, scales = "free") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = NULL, y = expression("Rmax" ~ (mu*mol ~ cm^{-2} ~ h^{-1})))

topt_compare_plot
rmax_compare_plot

###compare values across bayesian and nls estimates

PnR_clean <- read_csv(here("Data","RespoFiles","TPC","PnR_clean_no4.csv"))
PnR_clean <- PnR_clean |> left_join(species_meta, by = "frag_ID")

# 1. fragment Bayesian (from the bayesTPC loop)
params_bayes_wide <- read_csv(here("Data","RespoFiles","TPC","bayes_params_no4.csv"))

# 2. species-level nls (topt_df_sp from the nls species loop)
topt_sp <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4_species.csv")) |>
  select(rmax, topt, ctmin, ctmax, e, eh, q10, thermal_safety_margin,
         thermal_tolerance, breadth, skewness, PR, species) |>
  distinct(PR, species, .keep_all = TRUE)

# 3. fragment-level nls
topt_frag <- read_csv(here("Data","RespoFiles","TPC","Topt_data_clean_no4.csv"))

# Make sure every fragment has its species (avoids .x/.y columns from earlier joins)
frag_lookup <- PnR_clean |> distinct(frag_ID, species, full_species)
sp_lookup   <- PnR_clean |> distinct(species, full_species)

topt_frag <- topt_frag |>
  select(-any_of(c("species", "full_species"))) |>
  left_join(frag_lookup, by = "frag_ID")

# Parameters shared by all three approaches
# (Bayes names -> nls names: T_opt = topt, e_h = eh)
params_to_compare <- c("topt", "rmax", "e", "eh")

# columns: method, PR, species, param, estimate, lower, upper, n

# Bayesian: posterior median and 95% HPD
# fragment-level Bayesian results with species attached
frag_lookup <- PnR_clean |> distinct(frag_ID, species, full_species)

bayes_frag <- params_bayes_wide |>
  select(-any_of(c("species", "full_species"))) |>
  left_join(frag_lookup, by = "frag_ID") |>
  transmute(PR, species, frag_ID,
            topt = T_opt_median,
            rmax = rmax_median,
            e    = e_median,
            eh   = e_h_median)

# Bayesian fragment level: species mean of fragment posterior medians, 95% CI of the mean
frag_summary_bayes <- bayes_frag |>
  pivot_longer(all_of(params_to_compare), names_to = "param", values_to = "value") |>
  filter(is.finite(value)) |>
  group_by(PR, species, param) |>
  summarise(n        = n(),
            estimate = mean(value),
            sd       = sd(value),
            se       = sd / sqrt(n),
            median   = median(value),
            lower    = estimate - qt(0.975, df = pmax(n - 1, 1)) * se,
            upper    = estimate + qt(0.975, df = pmax(n - 1, 1)) * se,
            .groups  = "drop")

bayes_frag_long <- frag_summary_bayes |>
  select(PR, species, param, estimate, lower, upper, n) |>
  mutate(method = "Bayesian (fragment mean)")

# nls species: point estimates only (no parameter CIs were bootstrapped)
nls_sp_long <- topt_sp |>
  select(PR, species, all_of(params_to_compare)) |>
  pivot_longer(all_of(params_to_compare), names_to = "param", values_to = "estimate") |>
  mutate(lower = NA_real_, upper = NA_real_,
         method = "nls (species)", n = NA_integer_)

# nls fragment: species mean with 95% CI of the mean across fragments
frag_summary <- topt_frag |>
  select(PR, species, frag_ID, all_of(params_to_compare)) |>
  pivot_longer(all_of(params_to_compare), names_to = "param", values_to = "value") |>
  filter(is.finite(value)) |>
  group_by(PR, species, param) |>
  summarise(n        = n(),
            estimate = mean(value),
            sd       = sd(value),
            se       = sd / sqrt(n),
            median   = median(value),
            lower    = estimate - qt(0.975, df = pmax(n - 1, 1)) * se,
            upper    = estimate + qt(0.975, df = pmax(n - 1, 1)) * se,
            .groups  = "drop")

nls_frag_long <- frag_summary |>
  select(PR, species, param, estimate, lower, upper, n) |>
  mutate(method = "nls (fragment mean)")

compare_long <- bind_rows(bayes_frag_long, nls_sp_long, nls_frag_long) |>
  left_join(sp_lookup, by = "species") |>
  mutate(method = factor(method, levels = c("Bayesian (fragment mean)",
                                            "nls (species)",
                                            "nls (fragment mean)")),
         param  = factor(param, levels = params_to_compare))

#pivot!! 

compare_wide <- compare_long |>
  mutate(method_short = recode(method,
                               "Bayesian (fragment mean)"  = "bayes",
                               "nls (species)"       = "nls_sp",
                               "nls (fragment mean)" = "nls_frag")) |>
  select(PR, species, full_species, param, method_short, estimate) |>
  pivot_wider(names_from = c(param, method_short),
              values_from = estimate,
              names_glue = "{param}_{method_short}") |>
  mutate(topt_bayes_minus_nls_sp   = topt_bayes - topt_nls_sp,
         topt_bayes_minus_nls_frag = topt_bayes - topt_nls_frag) |>
  arrange(PR, species)

# Readable version for a report: "estimate [lower, upper]", rounded
compare_table <- compare_long |>
  mutate(value = case_when(
    is.na(lower) ~ sprintf("%.2f", estimate),
    TRUE         ~ sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper))) |>
  select(PR, full_species, param, method, value) |>
  pivot_wider(names_from = method, values_from = value) |>
  arrange(PR, param, full_species)

#compare_wide
#compare_table

write_csv(compare_wide,  here("Data","RespoFiles","TPC","TPC_method_comparison_wide.csv"))
write_csv(compare_table, here("Data","RespoFiles","TPC","TPC_method_comparison_table.csv"))
write_csv(frag_summary,  here("Data","RespoFiles","TPC","TPC_fragment_species_means.csv"))

#plot2compare

compare_plot <- compare_long |>
  ggplot(aes(x = estimate, y = full_species, color = method)) +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  position = position_dodge(width = 0.6), size = 0.3) +
  facet_grid(PR ~ param, scales = "free") +
  scale_color_manual(values = c("Bayesian (fragment mean)"  = "#1b9e77",
                                "nls (species)"       = "#d95f02",
                                "nls (fragment mean)" = "#7570b3")) +
  theme_bw(base_size = 22) +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "bottom") +
  labs(x = "Estimate (Bayesian: median & 95% HPD; nls: mean & 95% CI)",
       y = NULL, color = NULL)

compare_plot
ggsave(here("Output/TPC/Graphs/tpc_param_comparison_bayes_nls.pdf"), 
       height = 10, width = 20, compare_plot)

# Topt only
topt_plot <- compare_long |>
  filter(param == "topt") |>
  ggplot(aes(x = estimate, y = full_species, color = method)) +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  position = position_dodge(width = 0.6)) +
  facet_wrap(~ PR, scales = "free") +
  scale_color_manual(values = c("Bayesian (fragment mean)"  = "#1b9e77",
                                "nls (species)"       = "#d95f02",
                                "nls (fragment mean)" = "#7570b3")) +
  theme_bw(base_size = 22) +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "bottom") +
  labs(x = "Topt (ºC)", y = NULL, color = NULL)

topt_plot

ggsave(here("Output/TPC/Graphs/topt_comparison_bayes_nls.pdf"), height = 8, width = 20, topt_plot)

#### Learning about Bayesian TPC package ###
# #Bayesian TPC package
# #bayesTPC: Bayesian inference for thermal performance curves in R
# #https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70004
# #https://github.com/johnwilliamsmithjr/bayesTPC
# 
# #learn about bayesTPC
# get_models()
# # [1] "poisson_glm_lin"    "poisson_glm_quad"   "binomial_glm_lin"   "binomial_glm_quad" 
# # [5] "bernoulli_glm_lin"  "bernoulli_glm_quad" "briere"             "gaussian"          
# # [9] "kamykowski"         "pawar_shsch"        "quadratic"          "ratkowsky"         
# # [13] "stinner"            "weibull" 
# 
# #we want to use:
# #pawar_shsch which is the sharpe-schoolfield model
# 
# #fit models using b_TPC()
# #things you can adjust:
# #number of iterations (niter)
# #burn-in period length (burn)
# #number of MCMC chains (nchains). 
# #Any of four sampling methods implemented in nimble can be specified using the samplerType argument.
# 
# #model info functions:
# # get_formula()	Gets the formula for a model
# # get_model_params()	Gets all parameters to be fitted
# # get_default_priors()	Gets all priors for fitted parameters
# # get_model_constants()	Gets model constants, if they exist
# # get_default_constants()	Gets default values for model constants, if they exist
# # get_models()	Lists the implemented TPC models in the package
# # get_default_model_specification()	Obtain the default specification and priors for implemented TPC models
# 
# #example
# # simulate data
# N <- 16
# q <- .75
# T_min <- 10
# T_max <- 35
# Temps <- rep(c(15, 20, 25, 30), N / 4)
# Traits <- rep(0, N)
# for (i in 1:N) {
#   while (Traits[i] <= 0) {
#     Traits[i] <- rnorm(
#       1,
#       -1 * q * (Temps[i] - T_max) * (Temps[i] - T_min) * (Temps[i] > T_min) * (Temps[i] < T_max), 2
#     )
#   }
# }
# trait_list <- list(Trait = Traits, Temp = Temps)
# 
# # create model
# ## Not run: 
# quadratic_model <- b_TPC(data = trait_list, model = "quadratic")
# quadratic_model <- b_TPC(
#   data = trait_list, model = "quadratic", niter = 8000,
#   inits = list(T_min = 15, T_max = 30),
#   priors = list(q = "dunif(0, .5)", sigma.sq = "dexp(1)")
# )
# 
# #MCMC diagnostics for the object:
# #sample portion in the mcmc.list format from the coda package
# #MCMC diagnostic plots
# 
# #look at traceplot:
# traceplot(quadratic_model)
# #shows the sampled values by sequential iteration
# #If burn has been specified, traceplot() only shows the samples after the burn-in period.
# #here burn hasnt been specified
# 
# #examine the relationship between the prior information and the posterior samples
# ppo_plot(quadratic_model)
# #degree of overlap between the priors specified versus a kernel density estimation of the posterior sample
# 
# #other visualizations can be done with code
# #e.g. gelman.diag() and gelman.plot()
# 
# #can pass the model to any diagnostic tool that accepts the coda::mcmc.list object type.
# 
# #Summaries and visualization
# #print(), summary(), and predict()
# print(quadratic_model) #quick overview of the fitted model
# summary(quadratic_model) #detailed summary of the MCMC results and returns summary statistics of the sample including means and credible intervals for all model parameters
# quad_predictions <- predict(quadratic_model) #fitted centre (mean or median) and bounding (95% quantiles or HPD intervals) of the TPC based on the MCMC samples.
# #We additionally provide the sample maximum a posteriori (MAP) estimator—the MCMC sample with the highest posterior probability among the obtained samples, similar to a classical maximum likelihood estimator
# 
# #visualize the (pairwise) joint posterior distribution of parameters to represent the relative density of samples.
# ipairs(quadratic_model)
# 
# #plots the median and 95% Highest Posterior Density (HPD) interval of the fitted function 
# #i.e. plugging the samples into the TPC function and calculating the median and HPD interval at all evaluated temperatures
# plot(quadratic_model)
# 
# posterior_predictive(quadratic_model)
# #simulates draws from the posterior predictive distribution - both the samples describing the TPC function and the observational model
# #uses these samples to calculate the mean/median and the HPD interval of those simulated points
# 
# plot_prediction(quadratic_model)
# #simulates from the posterior and then plots them
# #or can be fed input from posterior_predictive() 
# #to visualize a previously calculated prediction.
# 
# #comparing models:
# #bayesTPC allows access to the Widely Applicable Information Criterion (WAIC, Gelman et al., 2013) 
# get_WAIC(quadratic_model)
# #This can be used to compare our fitted models
# #The preferred model will be the one with the lowest wAIC value.
# 
# #Read in data from Respo_process_TPC script
# 
# # read in data
# df_clean <- read_csv(here("Data","RespoFiles","TPC","PnR_clean_no4.csv"))
# #removes 289 data points (makes sense, one whole run plus quite a few)
# 
# #look at priors for shsch model
# get_default_model_specification("pawar_shsch")
# # bayesTPC Model Specification of Type: pawar_shsch
# # Model Formula:
# #   m[i] <- ( (e_h > e) * r_tref * exp((e/(8.62e-05)) * ((1/(T_ref + 273.15)) - (1/(Temp + 273.15))))/(1 + (e/(e_h - e)) * exp((e_h/(8.62e-05)) * (1/(T_opt + 273.15) - 1/(Temp + 273.15)))) )
# # Model Distribution:
# #   Trait[i] ~ T(dnorm(mean = m[i], tau = 1/sigma.sq), 0, )
# # Model Parameters and Priors:
# #   e ~ dunif(0, 1)
# # e_h ~ dunif(0, 30)
# # r_tref ~ dunif(0, 10)
# # T_opt ~ dunif(0, 50)
# # sigma.sq ~ dexp(1)
# # Model Constants:
#   T_ref = 20