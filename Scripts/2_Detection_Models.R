# -------------------------------------------------------
#
#                    Load/install Packages
#
# -------------------------------------------------------


library(tidyverse)
library(spOccupancy)
library(terra)
library(sp)
library(sf)
library(raster)

set.seed(123)
options(scipen = 9999)
setwd(".")


# -------------------------------------------------------
#
#                    Read In data
#
# -------------------------------------------------------

# Occupancy data
leo_dat <- read.csv("./Data/LeopardOccu.csv")

# Camera locations
cam_loc <- vect("D:/Persian_Leopard_GIS/CameraLocationsNew/CameraLocationsNew.shp")

# Initial detection covariates
det_covs <- read.csv("./Data/Det_Covs_Nima.csv")

# -------------------------------------------------------
#
#                    Data Wrangling
#
# -------------------------------------------------------

# Camera locations to dataframe
cam_loc_df <- as.data.frame(cam_loc)

# Keeping only non-detection/detection
y <- leo_dat[,-c(1,2)]
# View(y)

# Removing first 48 days since only one camera was operational for a total of 120 days
y120 <- y[,-c(1:48)]
colnames(y120) <- 1:ncol(y120)

# Removing the first 48 and last 22 days to a period when all cameras were active
y88 <- y[,-c(1:48, 137:168)]
y[,-c(1:48, 146:168)]
colnames(y88) <- 1:ncol(y88)

# Creating a spOccupancy data
spOccu_dat120 <- list(y = y120,
                   coords = cam_loc_df[,c(2,3)]   
)

spOccu_dat88 <- list(y = y88,
                      coords = cam_loc_df[,c(2,3)]   
)



# -------------------------------------------------------
#
#         Cumulative Detection Probability
#
# -------------------------------------------------------

# MCMC
n.samples <- 80000
n.burn <- 20000
n.thin <- 5
n.chains <- 3

n.report = 20000

# Initial values
inits120 <- list(beta = 0, 
              alpha = 0,
              z = apply(spOccu_dat120$y, 1, max, na.rm = TRUE)
)

inits88 <- list(beta = 0, 
                 alpha = 0,
                 z = apply(spOccu_dat88$y, 1, max, na.rm = TRUE)
)



# Priors
priors <- list(beta.normal = list(mean = 0, var = 2.72),
               alpha.normal = list(mean = 0, var = 2.72) 
)

# Fit null detection model for 120 days
p_null120 <- PGOcc(occ.formula = ~ 1, 
                  det.formula = ~1, 
                  data = spOccu_dat120, 
                  inits = inits120, 
                  n.samples = n.samples, 
                  priors = priors, 
                  n.omp.threads = 1, 
                  verbose = TRUE, 
                  n.report = n.report, 
                  n.burn = n.burn,
                  n.thin = n.thin, 
                  n.chains = n.chains)


# Check for convergence
p_null120_rhat<- p_null120$rhat
p_null120_rhat_vals <- unlist(p_null120_rhat)
which(p_null120_rhat_vals > 1.1)

# Fit null detection model for 88 days
p_null88 <- PGOcc(occ.formula = ~ 1, 
                   det.formula = ~1, 
                   data = spOccu_dat88, 
                   inits = inits88, 
                   n.samples = n.samples, 
                   priors = priors, 
                   n.omp.threads = 1, 
                   verbose = TRUE, 
                   n.report = n.report, 
                   n.burn = n.burn,
                   n.thin = n.thin, 
                   n.chains = n.chains)


# Check for convergence
p_null88_rhat<- p_null88$rhat
p_null88_rhat_vals <- unlist(p_null88_rhat)
which(p_null88_rhat_vals > 1.1)

# Extract posterior samples
det_post120 <- p_null120$alpha.samples 
det_post88 <- p_null88$alpha.samples 

# Transform from logit to probability scale
det_prob120 <- plogis(det_post120[, 1]) 
det_prob88 <- plogis(det_post88[, 1]) 

# Number of survey days
n_days120 <- ncol(y120)
n_days88 <- ncol(y88)

# Create df
plot_data120 <- data.frame(
  Day = 1:n_days120,
  Mean = sapply(1:n_days120, function(d) mean(1 - (1 - det_prob120)^d)),
  Lower = sapply(1:n_days120, function(d) quantile(1 - (1 - det_prob120)^d, 0.025)),
  Upper = sapply(1:n_days120, function(d) quantile(1 - (1 - det_prob120)^d, 0.975)),
  Group = "120 days"
)

plot_data88 <- data.frame(
  Day = 1:n_days88,
  Mean = sapply(1:n_days88, function(d) mean(1 - (1 - det_prob88)^d)),
  Lower = sapply(1:n_days88, function(d) quantile(1 - (1 - det_prob88)^d, 0.025)),
  Upper = sapply(1:n_days88, function(d) quantile(1 - (1 - det_prob88)^d, 0.975)),
  Group = "88 days"
)

# Combine data
plot_data <- rbind(plot_data120, plot_data88)


# Plot
det_prob_plot <- ggplot(plot_data, aes(x = Day, y = Mean, color = Group, fill = Group)) +
                    geom_line() +
                    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
                    # geom_vline(xintercept = 10, linetype = "dashed", color = "black") +
                    scale_color_manual(values = c("120 days" = "blue", "88 days" = "red")) +
                    scale_fill_manual(values = c("120 days" = "blue", "88 days" = "red")) +
                    scale_x_continuous(breaks = seq(0, max(plot_data$Day), by = 20)) +
                    labs(
                      y = "Cumulative Detection Probability",
                      x = "Survey Day",
                      title = "Cumulative Detection Probability by Duration"
                    ) +
                    theme_minimal(base_size = 14) +
                    theme(
                      panel.background = element_rect(fill = "white", color = NA),
                      plot.background = element_rect(fill = "white", color = NA),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.title.x = element_text(margin = margin(t = 10)),
                      axis.title.y = element_text(margin = margin(r = 10))
)

# Print
print(det_prob_plot)

# little to no difference in detection
# 120 days retains more information since there is only 120 sites

# Round detection probs
plot_data120$Mean <- round(plot_data120$Mean, 2)
plot_data120$Lower <- round(plot_data120$Lower, 2)
plot_data120$Upper <- round(plot_data120$Upper, 2)
print(plot_data120)

# At 10 days detection probability is 0.45 and 
# at 12 days det prob is 0.51 and would result in fewer NAs than 10 day occasion lengths

# Save
ggsave("./Outputs/Cumulative_Detection_Prob_Plot.png", plot = det_prob_plot, width = 8, height = 6, dpi = 300)
write.csv(plot_data120, "./Outputs/Cumulative_Detection_Prob_120days.csv")


# -------------------------------------------------------
#
#              Occasions
#
# -------------------------------------------------------

# Creating a matrix
n_periods <- 10
days_per_period <- 12
n_sites <- nrow(y120)
y_occ <- matrix(NA, nrow = n_sites, ncol = n_periods)

for (i in 1:n_periods) {
  # Define column indices for the current period
  col_start <- (i - 1) * days_per_period + 1
  col_end <- min(i * days_per_period, ncol(y))
  cols <- col_start:col_end
  
  # Subset the period's columns
  period_data <- y120[, cols]
  
  # Apply row-wise logic
  y_occ[, i] <- apply(period_data, 1, function(row) {
    if (all(is.na(row))) {
      return(NA)           # No survey effort
    } else if (any(row == 1, na.rm = TRUE)) {
      return(1)            # At least one detection
    } else if (any(row == 0, na.rm = TRUE)) {
      return(0)            # At least one non-detection
    } else {
      return(NA)           # All NA, again as a fallback
    }
  })
}

# Check output
print(y_occ)

saveRDS(y_occ, "./Data/y_occ.rds")

# -------------------------------------------------------
#
#              Detection Covariates
#
# -------------------------------------------------------

# Days a camera was active during each occasion
periods <- split(1:ncol(y120), ceiling(seq_along(1:ncol(y120))/days_per_period))
Days_Active <- sapply(periods, function(cols) rowSums(!is.na(y120[, cols, drop = FALSE])))
Days_Active <- as.data.frame(Days_Active)
colnames(Days_Active) <- paste0("Days_Active_", 1:n_periods)
print(Days_Active)
saveRDS(Days_Active, "./Data/Days_Active.rds")

# Survey period
Survey_Period <- as.data.frame(sapply(1:n_periods, function(x) rep(x, n_sites)))


# -------------------------------------------------------
#
#             spOccupancy Data
#
# -------------------------------------------------------

# Creating a spOccupancy data
spOccu_dat <- list(y = y_occ,
                   det.covs = list(Days_Active = Days_Active,
                                   Survey_Period = Survey_Period,
                                   Camera_Type = det_covs$camera.type,
                                   Picture_Background = det_covs$picture.background
                   ),
                   coords = cam_loc_df[,c(2,3)]
)
str(spOccu_dat)

# -------------------------------------------------------
#
#             Detection Models
#
# -------------------------------------------------------

# MCMC specifications
n.samples <- 150000
n.burn <- 50000
n.thin <- 10
n.chains <- 3

# Initial values
inits <- list(beta = 0, 
              alpha = 0,
              z = apply(spOccu_dat$y, 1, max, na.rm = TRUE)
)

# Priors
priors <- list(beta.normal = list(mean = 0, var = 2.72),
               alpha.normal = list(mean = 0, var = 2.72)
)

# ---------------
# Formulas
# ---------------

# List of variable names
vars <- c("scale(Days_Active)", "(1|Survey_Period)", "as.factor(Camera_Type)", "as.factor(Picture_Background)")

# Generate all additive combinations
formulas <- lapply(1:length(vars), function(k) {
  combn(vars, k, FUN = function(x) paste("~", paste(x, collapse = " + ")), simplify = TRUE)
})

# Flatten into a character vector
formulas <- unlist(formulas)

# Add null model
formulas <- c("~ 1", formulas)

print(formulas)

# ---------------
# Fit Models
# ---------------

# Model storage
fit_list <- list()

# WAIC storage
waic_df <- data.frame(
                      Model = integer(),
                      Formula = character(), 
                      elpd = numeric(),
                      pD = numeric(), 
                      WAIC = numeric()
)


# Rhat storage
rhat_list <- list()

# Fit a model for each formula
for (i in seq_along(formulas)) {
  f <- formulas[[i]]
  message(sprintf("Fitting model %d of %d: formula = %s",
                  i, length(formulas), f))
  
  # Convert string to formula
  det_formula <- as.formula(f)
  
  # Fit model 
  fit <- PGOcc(occ.formula = ~ 1,
               det.formula = det_formula,
               data = spOccu_dat,
               inits = inits,
               priors = priors,
               n.samples = n.samples,
               n.burn = n.burn,
               n.thin = n.thin,
               n.chains = n.chains,
               verbose = FALSE
               # ,
               # n.report = 50000
  )

  # Save model
  fit_list[[i]] <- fit
  
  # Extract and save WAIC
  waic_val <- waicOcc(fit)
  waic_df <- rbind(waic_df, data.frame(
    Model = i,
    Formula = f,
    elpd = unname(waic_val["elpd"]),
    pD = unname(waic_val["pD"]),
    WAIC = unname(waic_val["WAIC"]),
    stringsAsFactors = FALSE
  ))
  
  # Save Rhat
  rhat_list[[f]] <- fit$rhat
  
}

# Check for convergence
rhat_vals <- unlist(rhat_list)
which(rhat_vals > 1.1)

# Rank by WAIC
waic_df <- waic_df[order(waic_df$WAIC), ]
print(waic_df)

# Check model fit of best detection model
ppc.out <- ppcOcc(fit_list[[waic_df[1,"Model"]]], fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out)
 
# Best detection model
summary(fit_list[[waic_df[1,"Model"]]])

# Save detection model ranking
write.csv(waic_df, "./Outputs/Detection_Model_Ranking.csv")


# Extract posterior samples
det_post_bm <- fit_list[[waic_df[1,"Model"]]]$alpha.samples 
 
# Transform from logit to probability scale
det_prob_bm <- plogis(det_post_bm[, 1]) 

# Number of occasions
n_days <- 10
 
# Create df
plot_data_bm <- data.frame(
  Day = 1:n_days,
  Mean = sapply(1:n_days, function(d) mean(1 - (1 - det_prob_bm)^d)),
  Lower = sapply(1:n_days, function(d) quantile(1 - (1 - det_prob_bm)^d, 0.025)),
  Upper = sapply(1:n_days, function(d) quantile(1 - (1 - det_prob_bm)^d, 0.975))
)

# Plot
ggplot(plot_data_bm, aes(x = Day, y = Mean)) +
      geom_line(color = "blue") +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.2) +
      scale_x_continuous(breaks = seq(0, max(plot_data_bm$Day), by = 1)) +
      labs(
        y = "Cumulative Detection Probability",
        x = "Occasion",
        title = "Cumulative Detection Probability"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))
      )

print(plot_data_bm)


# ------------------- End of Script -------------------