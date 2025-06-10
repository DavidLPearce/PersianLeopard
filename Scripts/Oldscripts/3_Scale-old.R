# -------------------------------------------------------
#
#                    Load/install Packages
#
# -------------------------------------------------------


library(tidyverse)
library(spOccupancy)

set.seed(123)
options(scipen = 9999)
setwd(".")


# -------------------------------------------------------
#
#                    Read In data
#
# -------------------------------------------------------

spOccu_dat <- readRDS("spOccu_dat.rds")
str(spOccu_dat)

# -------------------------------------------------------
#
#                       Models
#
# -------------------------------------------------------

# ---------------------
# Model specifications
# ---------------------

# MCMC
batch.length <- 50
n.batch <- 1500
batch.length * n.batch # Total number of MCMC samples per chain
n.burn <- 20000
n.thin <- 20
n.chains <- 3

total_iterations <- batch.length * n.batch
retained_per_chain <- (total_iterations - n.burn) / n.thin
posterior_samples <- retained_per_chain * n.chains
posterior_samples


# Pair-wise distances between all sites
dist.mat <- dist(spOccu_dat$coords)

# Exponential covariance model
cov.model <- "matern"

# Initial values
inits <- list(alpha = 0, 
              beta = 0, 
              z = apply(spOccu_dat$y, 1, max, na.rm = TRUE), 
              sigma.sq = 2, 
              phi = 3 / mean(dist.mat),
              nu = 1.5,
              w = rep(0, nrow(spOccu_dat$y))
)

# Tuning
tuning <- list(phi = 0.5,
               nu = 0.5
)

# Priors
min.dist <- min(dist.mat)
max.dist <- max(dist.mat)
priors <- list(beta.normal = list(mean = 0, var = 2.72), 
               alpha.normal = list(mean = 0, var = 2.72), 
               sigma.sq.ig = c(2, 1), 
               phi.unif = c(3/max.dist, 3/min.dist),
               nu.unif = c(0.5, 2.5)
)

# ---------------------
# Loop specifications
# ---------------------

# Define buffer distances 
buffers <- seq(100, 15000, by = 100)

# Create an empty list to store models or WAICs
waic_df <- matrix(NA, nrow = length(buffers), ncol = 4)
colnames(waic_df) <- c("Buffer", "elpd", "pD", "WAIC")
waic_df <- as.data.frame(waic_df)

# Rhat list
rhat_list <- list()

# ---------------------
# Elevation
# ---------------------
pb <- txtProgressBar(min = 0, max = length(buffers), style = 3)# Progress bar
for (i in seq_along(buffers)) {
  
  # Create covariate name
  covariate_name <- paste0("elev_mn_", buffers[i], "m")
  
  # Create the full formula string
  covariate_name <- paste("~ scale(", covariate_name, ")")
  
  # Convert to formula object
  occ_formula <- as.formula(covariate_name)

  # Fit the model
  model_fit <-spPGOcc(occ.formula = occ_formula, 
                      det.formula = ~ as.factor(Camera_Type), 
                      data = spOccu_dat,
                      inits = inits, 
                      priors = priors,
                      NNGP = TRUE,
                      cov.model = cov.model,
                      n.neighbors = 15,
                      n.batch = n.batch, 
                      batch.length = batch.length, 
                      tuning = tuning, 
                      n.burn = n.burn, 
                      n.thin = n.thin, 
                      n.chains = n.chains,
                      verbose = FALSE)

  # Save WAIC
  waic_value <- waicOcc(model_fit)
  waic_df[i, 1] <- buffers[i]
  waic_df[i, 2] <- waic_value["elpd"]
  waic_df[i, 3] <- waic_value["pD"]
  waic_df[i, 4] <- waic_value["WAIC"]
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
  
  # Update progress bar
  setTxtProgressBar(pb, i)
} # ----- End Loop -----

# Check Rhat
rhat_vals <- unlist(rhat_list)
which(rhat_vals > 1.1)

# Order by WAIC (ascending)
waic_df <- as.data.frame(waic_df)
waic_df <- waic_df[order(waic_df$WAIC), ]
print(waic_df)
plot(waic_df$Buffer, waic_df$WAIC)

# Save WAIC table
write.csv(waic_df,"./Outputs/ElevationScale_WAIC.csv")

# ---------------------
# Slope
# ---------------------

# reset waic_df and rhat_list
waic_df[,] <- NA
waic_df <- waic_df[order(as.numeric(rownames(waic_df))), ]
rhat_list <- list()

pb <- txtProgressBar(min = 0, max = length(buffers), style = 3)# Progress bar
for (i in seq_along(buffers)) {
  
  # Create covariate name
  covariate_name <- paste0("slope_mn_", buffers[i], "m")
  
  # Create the full formula string
  covariate_name <- paste("~ scale(", covariate_name, ")")
  
  # Convert to formula object
  occ_formula <- as.formula(covariate_name)
  
  # Fit the model
  model_fit <-spPGOcc(occ.formula = occ_formula, 
                      det.formula = ~ as.factor(Camera_Type), 
                      data = spOccu_dat,
                      inits = inits, 
                      priors = priors,
                      NNGP = TRUE,
                      cov.model = cov.model,
                      n.neighbors = 15,
                      n.batch = n.batch, 
                      batch.length = batch.length, 
                      tuning = tuning, 
                      n.burn = n.burn, 
                      n.thin = n.thin, 
                      n.chains = n.chains,
                      verbose = FALSE)
  
  # Save WAIC
  waic_value <- waicOcc(model_fit)
  waic_df[i, 1] <- buffers[i]
  waic_df[i, 2] <- waic_value["elpd"]
  waic_df[i, 3] <- waic_value["pD"]
  waic_df[i, 4] <- waic_value["WAIC"]
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
  
  # Update progress bar
  setTxtProgressBar(pb, i)
} # ----- End Loop -----

# Check Rhat
rhat_vals <- unlist(rhat_list)
which(rhat_vals > 1.1)

# Order by WAIC (ascending)
waic_df <- as.data.frame(waic_df)
waic_df <- waic_df[order(waic_df$WAIC), ]
print(waic_df)
plot(waic_df$Buffer, waic_df$WAIC)

# Save WAIC table
write.csv(waic_df,"./Outputs/SlopeScale_WAIC.csv")

# ---------------------
# VRM
# ---------------------

# reset waic_df and rhat_list
waic_df[,] <- NA
waic_df <- waic_df[order(as.numeric(rownames(waic_df))), ]
rhat_list <- list()

pb <- txtProgressBar(min = 0, max = length(buffers), style = 3)# Progress bar
for (i in seq_along(buffers)) {
  
  # Create covariate name
  covariate_name <- paste0("VRM_mn_", buffers[i], "m")
  
  # Create the full formula string
  covariate_name <- paste("~ scale(", covariate_name, ")")
  
  # Convert to formula object
  occ_formula <- as.formula(covariate_name)
  
  # Fit the model
  model_fit <-spPGOcc(occ.formula = occ_formula, 
                      det.formula = ~ as.factor(Camera_Type), 
                      data = spOccu_dat,
                      inits = inits, 
                      priors = priors,
                      NNGP = TRUE,
                      cov.model = cov.model,
                      n.neighbors = 15,
                      n.batch = n.batch, 
                      batch.length = batch.length, 
                      tuning = tuning, 
                      n.burn = n.burn, 
                      n.thin = n.thin, 
                      n.chains = n.chains,
                      verbose = FALSE)
  
  # Save WAIC
  waic_value <- waicOcc(model_fit)
  waic_df[i, 1] <- buffers[i]
  waic_df[i, 2] <- waic_value["elpd"]
  waic_df[i, 3] <- waic_value["pD"]
  waic_df[i, 4] <- waic_value["WAIC"]
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
  
  # Update progress bar
  setTxtProgressBar(pb, i)
} # ----- End Loop -----

# Check Rhat
rhat_vals <- unlist(rhat_list)
which(rhat_vals > 1.1)

# Order by WAIC (ascending)
waic_df <- as.data.frame(waic_df)
waic_df <- waic_df[order(waic_df$WAIC), ]
print(waic_df)
plot(waic_df$Buffer, waic_df$WAIC)

# Save WAIC table
waic_df<-read.csv("./Outputs/VRMScale_WAIC.csv")
write.csv(waic_df,"./Outputs/VRMScale_WAIC.csv")

# ---------------------
# N Pastures
# ---------------------

# reset waic_df and rhat_list
waic_df[,] <- NA
waic_df <- waic_df[order(as.numeric(rownames(waic_df))), ]
rhat_list <- list()

pb <- txtProgressBar(min = 0, max = length(buffers), style = 3)# Progress bar
for (i in seq_along(buffers)) {
  
  # Create covariate name
  covariate_name <- paste0("Npastures_", buffers[i], "m")
  
  # Create the full formula string
  covariate_name <- paste("~ ", covariate_name, "")
  
  # Convert to formula object
  occ_formula <- as.formula(covariate_name)
  
  # Fit the model
  model_fit <-spPGOcc(occ.formula = occ_formula, 
                      det.formula = ~ as.factor(Camera_Type), 
                      data = spOccu_dat,
                      inits = inits, 
                      priors = priors,
                      NNGP = TRUE,
                      cov.model = cov.model,
                      n.neighbors = 15,
                      n.batch = n.batch, 
                      batch.length = batch.length, 
                      tuning = tuning, 
                      n.burn = n.burn, 
                      n.thin = n.thin, 
                      n.chains = n.chains,
                      verbose = FALSE)
  
  # Save WAIC
  waic_value <- waicOcc(model_fit)
  waic_df[i, 1] <- buffers[i]
  waic_df[i, 2] <- waic_value["elpd"]
  waic_df[i, 3] <- waic_value["pD"]
  waic_df[i, 4] <- waic_value["WAIC"]
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
  
  # Update progress bar
  setTxtProgressBar(pb, i)
} # ----- End Loop -----

# Check Rhat
rhat_vals <- unlist(rhat_list)
which(rhat_vals > 1.1)

# Order by WAIC (ascending)
waic_df <- as.data.frame(waic_df)
waic_df <- waic_df[order(waic_df$WAIC), ]
print(waic_df)

# Save WAIC table
write.csv(waic_df,"./Outputs/NpasturesScale_WAIC.csv")

# ---------------------
# N Villages
# ---------------------

# reset waic_df and rhat_list
waic_df[,] <- NA
waic_df <- waic_df[order(as.numeric(rownames(waic_df))), ]
rhat_list <- list()

pb <- txtProgressBar(min = 0, max = length(buffers), style = 3)# Progress bar
for (i in seq_along(buffers)) {
  
  # Create covariate name
  covariate_name <- paste0("Nvillages_", buffers[i], "m")
  
  # Create the full formula string
  covariate_name <- paste("~ ", covariate_name, "")
  
  # Convert to formula object
  occ_formula <- as.formula(covariate_name)
  
  # Fit the model
  model_fit <-spPGOcc(occ.formula = occ_formula, 
                      det.formula = ~ as.factor(Camera_Type), 
                      data = spOccu_dat,
                      inits = inits, 
                      priors = priors,
                      NNGP = TRUE,
                      cov.model = cov.model,
                      n.neighbors = 15,
                      n.batch = n.batch, 
                      batch.length = batch.length, 
                      tuning = tuning, 
                      n.burn = n.burn, 
                      n.thin = n.thin, 
                      n.chains = n.chains,
                      verbose = FALSE)
  
  # Save WAIC
  waic_value <- waicOcc(model_fit)
  waic_df[i, 1] <- buffers[i]
  waic_df[i, 2] <- waic_value["elpd"]
  waic_df[i, 3] <- waic_value["pD"]
  waic_df[i, 4] <- waic_value["WAIC"]
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
  
  # Update progress bar
  setTxtProgressBar(pb, i)
} # ----- End Loop -----

# Check Rhat
rhat_vals <- unlist(rhat_list)
which(rhat_vals > 1.1)

# Order by WAIC (ascending)
waic_df <- as.data.frame(waic_df)
waic_df <- waic_df[order(waic_df$WAIC), ]
print(waic_df)

# Save WAIC table
write.csv(waic_df,"./Outputs/NvillagesScale_WAIC.csv")


# ---------------------
# N River
# ---------------------

# reset waic_df and rhat_list
waic_df[,] <- NA
waic_df <- waic_df[order(as.numeric(rownames(waic_df))), ]
rhat_list <- list()

pb <- txtProgressBar(min = 0, max = length(buffers), style = 3)# Progress bar
for (i in seq_along(buffers)) {
  
  # Create covariate name
  covariate_name <- paste0("Nriver_", buffers[i], "m")
  
  # Create the full formula string
  covariate_name <- paste("~ ", covariate_name, "")
  
  # Convert to formula object
  occ_formula <- as.formula(covariate_name)
  
  # Fit the model
  model_fit <-spPGOcc(occ.formula = occ_formula, 
                      det.formula = ~ as.factor(Camera_Type), 
                      data = spOccu_dat,
                      inits = inits, 
                      priors = priors,
                      NNGP = TRUE,
                      cov.model = cov.model,
                      n.neighbors = 15,
                      n.batch = n.batch, 
                      batch.length = batch.length, 
                      tuning = tuning, 
                      n.burn = n.burn, 
                      n.thin = n.thin, 
                      n.chains = n.chains,
                      verbose = FALSE)
  
  # Save WAIC
  waic_value <- waicOcc(model_fit)
  waic_df[i, 1] <- buffers[i]
  waic_df[i, 2] <- waic_value["elpd"]
  waic_df[i, 3] <- waic_value["pD"]
  waic_df[i, 4] <- waic_value["WAIC"]
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
  
  # Update progress bar
  setTxtProgressBar(pb, i)
} # ----- End Loop -----

# Check Rhat
rhat_vals <- unlist(rhat_list)
which(rhat_vals > 1.1)

# Order by WAIC (ascending)
waic_df <- as.data.frame(waic_df)
waic_df <- waic_df[order(waic_df$WAIC), ]
print(waic_df)

# Save WAIC table
write.csv(waic_df,"./Outputs/NriverScale_WAIC.csv")

# ---------------------
# VRM Agg Indx
# ---------------------

# VRM aggregation index and patch only have values for all sites
# starting at 1000 meters
buffers <- seq(1000, 15000, by = 100)

# Create an empty list to store models or WAICs
waic_df <- matrix(NA, nrow = length(buffers), ncol = 4)
colnames(waic_df) <- c("Buffer", "elpd", "pD", "WAIC")
as.data.frame(waic_df)

pb <- txtProgressBar(min = 0, max = length(buffers), style = 3)# Progress bar
for (i in seq_along(buffers)) {
  
  # Create covariate name
  covariate_name <- paste0("VRM_ai_", buffers[i], "m")
  
  # Create the full formula string
  covariate_name <- paste("~ scale(", covariate_name, ")")
  
  # Convert to formula object
  occ_formula <- as.formula(covariate_name)
  
  # Fit the model
  model_fit <-spPGOcc(occ.formula = occ_formula, 
                      det.formula = ~ as.factor(Camera_Type), 
                      data = spOccu_dat,
                      inits = inits, 
                      priors = priors,
                      NNGP = TRUE,
                      cov.model = cov.model,
                      n.neighbors = 15,
                      n.batch = n.batch, 
                      batch.length = batch.length, 
                      tuning = tuning, 
                      n.burn = n.burn, 
                      n.thin = n.thin, 
                      n.chains = n.chains,
                      verbose = FALSE)
  
  # Save WAIC
  waic_value <- waicOcc(model_fit)
  waic_df[i, 1] <- buffers[i]
  waic_df[i, 2] <- waic_value["elpd"]
  waic_df[i, 3] <- waic_value["pD"]
  waic_df[i, 4] <- waic_value["WAIC"]
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
  
  # Update progress bar
  setTxtProgressBar(pb, i)
} # ----- End Loop -----

# Check Rhat
rhat_vals <- unlist(rhat_list)
which(rhat_vals > 1.1)

# Order by WAIC (ascending)
waic_df <- as.data.frame(waic_df)
waic_df <- waic_df[order(waic_df$WAIC), ]
print(waic_df)

# Save WAIC table
write.csv(waic_df,"./Outputs/VRMAggIndxiScale_WAIC.csv")

# ---------------------
# VRM Patch Area
# ---------------------

# reset waic_df and rhat_list
waic_df[,] <- NA
waic_df <- waic_df[order(as.numeric(rownames(waic_df))), ]
rhat_list <- list()

pb <- txtProgressBar(min = 0, max = length(buffers), style = 3)# Progress bar
for (i in seq_along(buffers)) {
  
  # Create covariate name
  covariate_name <- paste0("VRM_mnParea_", buffers[i], "m")
  
  # Create the full formula string
  covariate_name <- paste("~ scale(", covariate_name, ")")
  
  # Convert to formula object
  occ_formula <- as.formula(covariate_name)
  
  # Fit the model
  model_fit <-spPGOcc(occ.formula = occ_formula, 
                      det.formula = ~ as.factor(Camera_Type), 
                      data = spOccu_dat,
                      inits = inits, 
                      priors = priors,
                      NNGP = TRUE,
                      cov.model = cov.model,
                      n.neighbors = 15,
                      n.batch = n.batch, 
                      batch.length = batch.length, 
                      tuning = tuning, 
                      n.burn = n.burn, 
                      n.thin = n.thin, 
                      n.chains = n.chains,
                      verbose = FALSE)
  
  # Save WAIC
  waic_value <- waicOcc(model_fit)
  waic_df[i, 1] <- buffers[i]
  waic_df[i, 2] <- waic_value["elpd"]
  waic_df[i, 3] <- waic_value["pD"]
  waic_df[i, 4] <- waic_value["WAIC"]
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
  
  # Update progress bar
  setTxtProgressBar(pb, i)
} # ----- End Loop -----

# Check Rhat
rhat_vals <- unlist(rhat_list)
which(rhat_vals > 1.1)

# Order by WAIC (ascending)
waic_df <- as.data.frame(waic_df)
waic_df <- waic_df[order(waic_df$WAIC), ]
print(waic_df)

# Save WAIC table
write.csv(waic_df,"./Outputs/VRMPatchAreaScale_WAIC.csv")

 