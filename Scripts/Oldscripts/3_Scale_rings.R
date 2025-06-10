
# ------------------------------------------------------------------------------
#
#                         Load/Install Packages
#
# ------------------------------------------------------------------------------


library(tidyverse)
library(spOccupancy)


set.seed(123)
options(scipen = 9999)
setwd(".")

# ------------------------------------------------------------------------------
#
#                             Read in data
#
# ------------------------------------------------------------------------------

# Covariates
dist_df <- readRDS("Data/OccuCovData/dist_df.rds")
elev_df <- readRDS("Data/OccuCovData/elev_df.rds")
slope_df <- readRDS("Data/OccuCovData/slope_df.rds")
aspect_df <- readRDS("Data/OccuCovData/aspect_df.rds")
VRM_df <- readRDS("Data/OccuCovData/VRM_df.rds")
VRM_AI_df <- readRDS("Data/OccuCovData/VRM_AI_df.rds")
VRM_Patch_df <- readRDS("Data/OccuCovData/Vrm_Patch_df.rds")
Npast_df <- readRDS("Data/OccuCovData/Npast_df.rds")
Nvil_df <- readRDS("Data/OccuCovData/Nvil_df.rds")
Nriver_df <- readRDS("Data/OccuCovData/Nriver_df.rds")

# Distance and number of pastures, villages, and rivers will
# not be used with distance-weighted scale of effect.
# Number of pastures, villages, and rivers will will just be
# in buffers

# Removing site, lat, and long
elev_df <- elev_df[,-c(1:3)]
slope_df <- slope_df[,-c(1:3)]
aspect_df <- aspect_df[,-c(1:3)]
VRM_df <- VRM_df[,-c(1:3)]
VRM_AI_df <- VRM_AI_df[,-c(1:3)] 
VRM_Patch_df <- VRM_Patch_df[,-c(1:3)]  
Npast_df <- Npast_df[,-c(1:3)]
Nvil_df <- Nvil_df[,-c(1:3)]
Nriver_df <- Nriver_df[,-c(1:3)]

# Initial detection covariates
Days_Active <- readRDS("./Data/Days_Active.rds")

# Detection matrix
y_occ <- readRDS("./Data/y_occ.rds")
 
# Camera locations
cam_loc_df <- dist_df[,c(1:3)]


# ------------------------------------------------------------------------------
#
#                 Distance-weighted Scale of Effect
#
# ------------------------------------------------------------------------------

# ---------------------------
# Function
# ---------------------------
# buffer_radii = buffer_vec
# sigma = 1500

# Distance weighted function
weight_fun <- function(buffer_radii, sigma) {
  # Calculate ring areas
  areas <- pi * (buffer_radii^2 - c(0, head(buffer_radii, -1))^2)
  
  # Gaussian weights before normalization
  raw_weights <- exp(- (buffer_radii^2) / (2 * sigma^2)) * areas
  
  # Normalize to sum to 1
  weights <- raw_weights / sum(raw_weights)
  
  return(weights)
}

# Example 
# buffer_vec <- seq(250, 15000, by = 250)
# weights <- weight_fun(buffer_vec, sigma = 500)
# max_idx <- which.max(weights)
# max_buffer <- buffer_vec[max_idx]
# 
# # Example
# plot(buffer_vec, weights, type = "l", main = "", ylab = "Distance-weight", xlab = "Distance (m)")
# abline(v = max_buffer, lty = 2, col = "red")


# ---------------------------
# Calculate Distance-weighted
# scale of effect
# ---------------------------

# Buffers and sigmas
buffer_vec <- seq(500, 15000, by = 250) 
sigma_vec <-  seq(500, 15000, by = 250) 

# Replicate data for each sigma 
elev_list <- replicate(length(sigma_vec), elev_df, simplify = FALSE)
slope_list <- replicate(length(sigma_vec), slope_df, simplify = FALSE)
aspect_list <- replicate(length(sigma_vec), aspect_df, simplify = FALSE)
VRM_list <- replicate(length(sigma_vec), VRM_df, simplify = FALSE)
VRM_AI_list <- replicate(length(sigma_vec), VRM_AI_df, simplify = FALSE)
VRM_Patch_list <- replicate(length(sigma_vec), VRM_Patch_df, simplify = FALSE)

# list to store weighted elevation at varying sigmas
weighted_list <- vector("list", length = length(sigma_vec))

# Elevation
cov_list <- elev_list
for (i in seq_along(sigma_vec)) {
  sigma <- sigma_vec[i]
  weights <- weight_fun(buffer_vec, sigma)
  cov_matrix <- cov_list[[i]]
  weighted_matrix <- sweep(cov_matrix, 2, weights, FUN = "*")
  colnames(weighted_matrix) <- paste0(colnames(cov_matrix), "_sigma", sigma)
  weighted_list[[i]] <- weighted_matrix
}

elev_weighted_df <- do.call(cbind, weighted_list)
head(elev_weighted_df[,1:10])
saveRDS(elev_weighted_df, "Data/OccuCovData/elev_weighted_df.rds")

# Slope 
cov_list <- slope_list
weighted_list <- vector("list", length = length(sigma_vec))
for (i in seq_along(sigma_vec)) {
  sigma <- sigma_vec[i]
  weights <- weight_fun(buffer_vec, sigma)
  cov_matrix <- cov_list[[i]]
  weighted_matrix <- sweep(cov_matrix, 2, weights, FUN = "*")
  colnames(weighted_matrix) <- paste0(colnames(cov_matrix), "_sigma", sigma)
  weighted_list[[i]] <- weighted_matrix
} 

slope_weighted_df <- do.call(cbind, weighted_list)
head(slope_weighted_df[,1:10])
saveRDS(slope_weighted_df, "Data/OccuCovData/slope_weighted_df.rds")

# Aspect
cov_list <- aspect_list
weighted_list <- vector("list", length = length(sigma_vec))
for (i in seq_along(sigma_vec)) {
  sigma <- sigma_vec[i]
  weights <- weight_fun(buffer_vec, sigma)
  cov_matrix <- cov_list[[i]]
  weighted_matrix <- sweep(cov_matrix, 2, weights, FUN = "*")
  colnames(weighted_matrix) <- paste0(colnames(cov_matrix), "_sigma", sigma)
  weighted_list[[i]] <- weighted_matrix
} 

aspect_weighted_df <- do.call(cbind, weighted_list)
head(aspect_weighted_df[,1:10])
saveRDS(aspect_weighted_df, "Data/OccuCovData/aspect_weighted_df.rds")

# VRM
cov_list <- VRM_list
weighted_list <- vector("list", length = length(sigma_vec))
for (i in seq_along(sigma_vec)) {
  sigma <- sigma_vec[i]
  weights <- weight_fun(buffer_vec, sigma)
  cov_matrix <- cov_list[[i]]
  weighted_matrix <- sweep(cov_matrix, 2, weights, FUN = "*")
  colnames(weighted_matrix) <- paste0(colnames(cov_matrix), "_sigma", sigma)
  weighted_list[[i]] <- weighted_matrix
} 

VRM_weighted_df <- do.call(cbind, weighted_list)
head(VRM_weighted_df[,1:10])
saveRDS(VRM_weighted_df, "Data/OccuCovData/VRM_weighted_df.rds")

# VRM AI
cov_list <- VRM_AI_list
weighted_list <- vector("list", length = length(sigma_vec))
for (i in seq_along(sigma_vec)) {
  sigma <- sigma_vec[i]
  weights <- weight_fun(buffer_vec, sigma)
  cov_matrix <- cov_list[[i]]
  weighted_matrix <- sweep(cov_matrix, 2, weights, FUN = "*")
  colnames(weighted_matrix) <- paste0(colnames(cov_matrix), "_sigma", sigma)
  weighted_list[[i]] <- weighted_matrix
}

VRM_AI_weighted_df <- do.call(cbind, weighted_list)
head(VRM_AI_weighted_df[,1:10])
saveRDS(VRM_AI_weighted_df, "Data/OccuCovData/VRM_AI_weighted_df.rds")


# VRM Patch
cov_list <- VRM_Patch_list
weighted_list <- vector("list", length = length(sigma_vec))
for (i in seq_along(sigma_vec)) {
  sigma <- sigma_vec[i]
  weights <- weight_fun(buffer_vec, sigma)
  cov_matrix <- cov_list[[i]]
  weighted_matrix <- sweep(cov_matrix, 2, weights, FUN = "*")
  colnames(weighted_matrix) <- paste0(colnames(cov_matrix), "_sigma", sigma)
  weighted_list[[i]] <- weighted_matrix
} 

VRM_Patch_weighted_df <- do.call(cbind, weighted_list)
head(VRM_Patch_weighted_df[,1:10])
saveRDS(VRM_Patch_weighted_df, "Data/OccuCovData/VRM_Patch_weighted_df.rds")

# ------------------------------------------------------------------------------
#                             Determinig Scale
# ------------------------------------------------------------------------------

# Elevation
spOccu_elev <- list(y = y_occ,
                   occ.covs = elev_weighted_df,
                   det.covs = list(Days_Active = Days_Active
                                   ),
                   coords = cam_loc_df[,c(2,3)])
# Slope
spOccu_slope <- list(y = y_occ,
                    occ.covs = slope_weighted_df,
                    det.covs = list(Days_Active = Days_Active
                    ),
                    coords = cam_loc_df[,c(2,3)])
# Aspect
spOccu_aspect <- list(y = y_occ,
                     occ.covs = aspect_weighted_df,
                     det.covs = list(Days_Active = Days_Active
                     ),
                     coords = cam_loc_df[,c(2,3)])
# VRM
spOccu_VRM <- list(y = y_occ,
                      occ.covs = VRM_weighted_df,
                      det.covs = list(Days_Active = Days_Active
                      ),
                      coords = cam_loc_df[,c(2,3)])
# VRM AI
spOccu_VRM_AI <- list(y = y_occ,
                   occ.covs = VRM_AI_weighted_df,
                   det.covs = list(Days_Active = Days_Active
                   ),
                   coords = cam_loc_df[,c(2,3)])
# VRM Patch
spOccu_VRM_Patch <- list(y = y_occ,
                      occ.covs = VRM_Patch_weighted_df,
                      det.covs = list(Days_Active = Days_Active
                      ),
                      coords = cam_loc_df[,c(2,3)])
# N Pasture
spOccu_Npast <- list(y = y_occ,
                         occ.covs = Npast_df,
                         det.covs = list(Days_Active = Days_Active
                         ),
                         coords = cam_loc_df[,c(2,3)])
# N Village
spOccu_Nvil <- list(y = y_occ,
                     occ.covs = Nvil_df,
                     det.covs = list(Days_Active = Days_Active
                     ),
                     coords = cam_loc_df[,c(2,3)])
# N River
spOccu_Nriver <- list(y = y_occ,
                    occ.covs = Nriver_df,
                    det.covs = list(Days_Active = Days_Active
                    ),
                    coords = cam_loc_df[,c(2,3)])

# ---------------------
# Model Specifications
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
dist.mat <- dist(spOccu_elev$coords)

# Exponential covariance model
cov.model <- "matern"

# Initial values
inits <- list(alpha = 0, 
              beta = 0, 
              z = apply(spOccu_elev$y, 1, max, na.rm = TRUE), 
              sigma.sq = 2, 
              phi = 3 / mean(dist.mat),
              nu = 1.5,
              w = rep(0, nrow(spOccu_elev$y))
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

# Create a grid of buffer-sigma combos
combo_grid <- expand.grid(Buffer = buffer_vec, Sigma = sigma_vec)

# Create storage for WAIC and Rhat
waic_df <- data.frame(Ring = NA, Sigma = NA, elpd = NA, pD = NA, WAIC = NA)
rhat_list <- list()

# -----------------------------
# Distance-weighted      
# -----------------------------

# ---------------------
# Elevation
# ---------------------
for (i in seq_len(nrow(combo_grid))) {
  
  buffer_val <- combo_grid$Buffer[i]
  sigma_val  <- combo_grid$Sigma[i]
  
  # Print current buffer and sigma
  message(sprintf("Fitting model %d of %d: buffer = %dm, sigma = %d",
                  i, nrow(combo_grid), buffer_val, sigma_val))
  
  # Construct covariate name 
  cov_name <- paste0("elev_ring", buffer_val, "_sigma", sigma_val)

  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula, 
                       det.formula = ~ Days_Active, 
                       data = spOccu_elev,
                       inits = inits, 
                       priors = priors,
                       NNGP = FALSE,
                       cov.model = cov.model,
                       n.neighbors = 5,
                       n.batch = n.batch, 
                       batch.length = batch.length, 
                       tuning = tuning, 
                       n.burn = n.burn, 
                       n.thin = n.thin, 
                       n.chains = n.chains,
                       n.omp.threads = 35,
                       verbose = FALSE)
  
  # Extract WAIC
  waic_val <- waicOcc(model_fit)
  
  # Save to waic_df
  waic_df[i, ] <- c(buffer_val, 
                    sigma_val,
                    waic_val["elpd"],
                    waic_val["pD"],
                    waic_val["WAIC"])
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
  
}

# Check for convergence
elev_rhat_vals <- unlist(rhat_list)
which(elev_rhat_vals > 1.1)

# Sort by WAIC
elev_waic <- waic_df[order(waic_df$WAIC), ]
head(elev_waic, 25)

# Save WAIC 
write.csv(elev_waic, "./Outputs/Elevation_BufferSigmaCombo_WAIC.csv")


# ---------------------
# Slope
# ---------------------

# Clear WAIC and Rhat
waic_df <- data.frame(Buffer = NA, Sigma = NA, elpd = NA, pD = NA, WAIC = NA)
rhat_list <- list()

for (i in seq_len(nrow(combo_grid))) {
  
  buffer_val <- combo_grid$Buffer[i]
  sigma_val  <- combo_grid$Sigma[i]
  
  # Print current buffer and sigma
  message(sprintf("Fitting model %d of %d: buffer = %dm, sigma = %d",
                  i, nrow(combo_grid), buffer_val, sigma_val))
  
  # Construct covariate name 
  cov_name <- paste0("slope_ring", buffer_val, "_sigma", sigma_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula, 
                       det.formula = ~ Days_Active, 
                       data = spOccu_slope,
                       inits = inits, 
                       priors = priors,
                       NNGP = FALSE,
                       cov.model = cov.model,
                       n.neighbors = 5,
                       n.batch = n.batch, 
                       batch.length = batch.length, 
                       tuning = tuning, 
                       n.burn = n.burn, 
                       n.thin = n.thin, 
                       n.chains = n.chains,
                       n.omp.threads = 35,
                       verbose = FALSE)
  
  # Extract WAIC
  waic_val <- waicOcc(model_fit)
  
  # Save to waic_df
  waic_df[i, ] <- c(buffer_val, 
                    sigma_val,
                    waic_val["elpd"],
                    waic_val["pD"],
                    waic_val["WAIC"])
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat

}

# Check for convergence
slope_rhat_vals <- unlist(rhat_list)
which(slope_rhat_vals > 1.1)

# Sort by WAIC
slope_waic <- waic_df[order(waic_df$WAIC), ]
head(slope_waic, 25)

# Save WAIC 
write.csv(slope_waic, "./Outputs/Slope_BufferSigmaCombo_WAIC.csv")


# ---------------------
# Aspect
# ---------------------

# Clear WAIC and Rhat
waic_df <- data.frame(Buffer = NA, Sigma = NA, elpd = NA, pD = NA, WAIC = NA)
rhat_list <- list()

for (i in seq_len(nrow(combo_grid))) {

  buffer_val <- combo_grid$Buffer[i]
  sigma_val  <- combo_grid$Sigma[i]

  # Print current buffer and sigma
  message(sprintf("Fitting model %d of %d: buffer = %dm, sigma = %d",
                  i, nrow(combo_grid), buffer_val, sigma_val))

  # Construct covariate name 
  cov_name <- paste0("aspect_ring", buffer_val, "_sigma", sigma_val)

  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))

  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_aspect,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov.model = cov.model,
                       n.neighbors = 5,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 35,
                       verbose = FALSE)

  # Extract WAIC
  waic_val <- waicOcc(model_fit)

  # Save to waic_df
  waic_df[i, ] <- c(buffer_val,
                    sigma_val,
                    waic_val["elpd"],
                    waic_val["pD"],
                    waic_val["WAIC"])

  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat

}

# Check for convergence
aspect_rhat_vals <- unlist(rhat_list)
which(aspect_rhat_vals > 1.1)

# Sort by WAIC
aspect_waic <- waic_df[order(waic_df$WAIC), ]
head(aspect_waic, 25)

# Save WAIC
write.csv(aspect_waic, "./Outputs/Aspect_BufferSigmaCombo_WAIC.csv")

# ---------------------
# VRM
# ---------------------

# Clear WAIC and Rhat
waic_df <- data.frame(Buffer = NA, Sigma = NA, elpd = NA, pD = NA, WAIC = NA)
rhat_list <- list()

for (i in seq_len(nrow(combo_grid))) {
  
  buffer_val <- combo_grid$Buffer[i]
  sigma_val  <- combo_grid$Sigma[i]
  
  # Print current buffer and sigma
  message(sprintf("Fitting model %d of %d: buffer = %dm, sigma = %d",
                  i, nrow(combo_grid), buffer_val, sigma_val))
  
  # Construct covariate name 
  cov_name <- paste0("VRM_ring", buffer_val, "_sigma", sigma_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula, 
                       det.formula = ~ Days_Active, 
                       data = spOccu_VRM,
                       inits = inits, 
                       priors = priors,
                       NNGP = FALSE,
                       cov.model = cov.model,
                       n.neighbors = 5,
                       n.batch = n.batch, 
                       batch.length = batch.length, 
                       tuning = tuning, 
                       n.burn = n.burn, 
                       n.thin = n.thin, 
                       n.chains = n.chains,
                       n.omp.threads = 35,
                       verbose = FALSE)
  
  # Extract WAIC
  waic_val <- waicOcc(model_fit)
  
  # Save to waic_df
  waic_df[i, ] <- c(buffer_val, 
                    sigma_val,
                    waic_val["elpd"],
                    waic_val["pD"],
                    waic_val["WAIC"])
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat

}

# Check for convergence
VRM_rhat_vals <- unlist(rhat_list)
which(VRM_rhat_vals > 1.1)

# Sort by WAIC
VRM_waic <- waic_df[order(waic_df$WAIC), ]
head(VRM_waic, 25)

# Save WAIC 
write.csv(VRM_waic, "./Outputs/VRM_BufferSigmaCombo_WAIC.csv")


# ---------------------
# VRM AI
# ---------------------

# Clear WAIC and Rhat
waic_df <- data.frame(Buffer = NA, Sigma = NA, elpd = NA, pD = NA, WAIC = NA)
rhat_list <- list()

for (i in seq_len(nrow(combo_grid))) {
  
  buffer_val <- combo_grid$Buffer[i]
  sigma_val  <- combo_grid$Sigma[i]
  
  # Print current buffer and sigma
  message(sprintf("Fitting model %d of %d: buffer = %dm, sigma = %d",
                  i, nrow(combo_grid), buffer_val, sigma_val))
  
  # Construct covariate name 
  cov_name <- paste0("VRM_AI_ring", buffer_val, "_sigma", sigma_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula, 
                       det.formula = ~ Days_Active, 
                       data = spOccu_VRM_AI,
                       inits = inits, 
                       priors = priors,
                       NNGP = FALSE,
                       cov.model = cov.model,
                       n.neighbors = 5,
                       n.batch = n.batch, 
                       batch.length = batch.length, 
                       tuning = tuning, 
                       n.burn = n.burn, 
                       n.thin = n.thin, 
                       n.chains = n.chains,
                       n.omp.threads = 35,
                       verbose = FALSE)
  
  # Extract WAIC
  waic_val <- waicOcc(model_fit)
  
  # Save to waic_df
  waic_df[i, ] <- c(buffer_val, 
                    sigma_val,
                    waic_val["elpd"],
                    waic_val["pD"],
                    waic_val["WAIC"])
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat

}

# Check for convergence
VRM_AI_rhat_vals <- unlist(rhat_list)
which(VRM_AI_rhat_vals > 1.1)

# Sort by WAIC
VRM_AI_waic <- waic_df[order(waic_df$WAIC), ]
head(VRM_AI_waic, 25)

# Save WAIC 
write.csv(VRM_AI_waic, "./Outputs/VRM_AI_BufferSigmaCombo_WAIC.csv")

# ---------------------
# VRM Patch
# ---------------------

# Clear WAIC and Rhat
waic_df <- data.frame(Buffer = NA, Sigma = NA, elpd = NA, pD = NA, WAIC = NA)
rhat_list <- list()

for (i in seq_len(nrow(combo_grid))) {
  
  buffer_val <- combo_grid$Buffer[i]
  sigma_val  <- combo_grid$Sigma[i]
  
  # Print current buffer and sigma
  message(sprintf("Fitting model %d of %d: buffer = %dm, sigma = %d",
                  i, nrow(combo_grid), buffer_val, sigma_val))
  
  # Construct covariate name 
  cov_name <- paste0("VRM_Parea_ring", buffer_val, "_sigma", sigma_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula, 
                       det.formula = ~ Days_Active, 
                       data = spOccu_VRM_Patch,
                       inits = inits, 
                       priors = priors,
                       NNGP = FALSE,
                       cov.model = cov.model,
                       n.neighbors = 5,
                       n.batch = n.batch, 
                       batch.length = batch.length, 
                       tuning = tuning, 
                       n.burn = n.burn, 
                       n.thin = n.thin, 
                       n.chains = n.chains,
                       n.omp.threads = 35,
                       verbose = FALSE)
  
  # Extract WAIC
  waic_val <- waicOcc(model_fit)
  
  # Save to waic_df
  waic_df[i, ] <- c(buffer_val, 
                    sigma_val,
                    waic_val["elpd"],
                    waic_val["pD"],
                    waic_val["WAIC"])
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat

}

# Check for convergence
VRM_Patch_rhat_vals <- unlist(rhat_list)
which(VRM_Patch_rhat_vals > 1.1)

# Sort by WAIC
VRM_Patch_waic <- waic_df[order(waic_df$WAIC), ]
head(VRM_Patch_waic, 25)

# Save WAIC 
write.csv(VRM_Patch_waic, "./Outputs/VRM_Patch_BufferSigmaCombo_WAIC.csv")


# -----------------------------
# Buffer      
# -----------------------------


# ---------------------
# N Pasture
# ---------------------

# Clear WAIC and Rhat
waic_df <- data.frame(Buffer = NA, elpd = NA, pD = NA, WAIC = NA)
rhat_list <- list()
i =2
for (i in 1:length(buffer_vec)) {
  
  buffer_val <- buffer_vec[i]

  # Print current buffer and sigma
  message(sprintf("Fitting model %d of %d: buffer = %dm",
                  i, length(buffer_vec), buffer_val))
  
  # Construct covariate name 
  cov_name <- paste0("Npastures_buf", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula, 
                       det.formula = ~ Days_Active, 
                       data = spOccu_Npast,
                       inits = inits, 
                       priors = priors,
                       NNGP = FALSE,
                       cov.model = cov.model,
                       n.neighbors = 5,
                       n.batch = n.batch, 
                       batch.length = batch.length, 
                       tuning = tuning, 
                       n.burn = n.burn, 
                       n.thin = n.thin, 
                       n.chains = n.chains,
                       n.omp.threads = 35,
                       verbose = FALSE)
  
  # Extract WAIC
  waic_val <- waicOcc(model_fit)
  
  # Save to waic_df
  waic_df[i, ] <- c(buffer_val, 
                    waic_val["elpd"],
                    waic_val["pD"],
                    waic_val["WAIC"])
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
  
}

# Check for convergence
Npast_rhat_vals <- unlist(rhat_list)
which(Npast_rhat_vals > 1.1)

# Sort by WAIC
Npast_waic <- waic_df[order(waic_df$WAIC), ]
head(Npast_waic, 25)

# Save WAIC 
write.csv(Npast_waic, "./Outputs/Npast_Buffer_WAIC.csv")

# ---------------------
# N Villages
# ---------------------

# Clear WAIC and Rhat
waic_df <- data.frame(Buffer = NA, elpd = NA, pD = NA, WAIC = NA)
rhat_list <- list()

for (i in 1:length(buffer_vec)) {
  
  buffer_val <- buffer_vec[i]
  
  # Print current buffer and sigma
  message(sprintf("Fitting model %d of %d: buffer = %dm",
                  i, length(buffer_vec), buffer_val))
  
  # Construct covariate name 
  cov_name <- paste0("Nvillages_buf", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula, 
                       det.formula = ~ Days_Active, 
                       data = spOccu_Nvil,
                       inits = inits, 
                       priors = priors,
                       NNGP = FALSE,
                       cov.model = cov.model,
                       n.neighbors = 5,
                       n.batch = n.batch, 
                       batch.length = batch.length, 
                       tuning = tuning, 
                       n.burn = n.burn, 
                       n.thin = n.thin, 
                       n.chains = n.chains,
                       n.omp.threads = 35,
                       verbose = FALSE)
  
  # Extract WAIC
  waic_val <- waicOcc(model_fit)
  
  # Save to waic_df
  waic_df[i, ] <- c(buffer_val, 
                    waic_val["elpd"],
                    waic_val["pD"],
                    waic_val["WAIC"])
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat

}

# Check for convergence
Nvil_rhat_vals <- unlist(rhat_list)
which(Nvil_rhat_vals > 1.1)

# Sort by WAIC
Nvil_waic <- waic_df[order(waic_df$WAIC), ]
head(Nvil_waic, 25)

# Save WAIC 
write.csv(Nvil_waic, "./Outputs/Nvil_Buffer_WAIC.csv")


# ---------------------
# N Rivers
# ---------------------

# Clear WAIC and Rhat
waic_df <- data.frame(Buffer = NA, elpd = NA, pD = NA, WAIC = NA)
rhat_list <- list()

for (i in 1:length(buffer_vec)) {
  
  buffer_val <- buffer_vec[i]
  
  # Print current buffer and sigma
  message(sprintf("Fitting model %d of %d: buffer = %dm",
                  i, length(buffer_vec), buffer_val))
  
  # Construct covariate name 
  cov_name <- paste0("Nriver_buf", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula, 
                       det.formula = ~ Days_Active, 
                       data = spOccu_Nriver,
                       inits = inits, 
                       priors = priors,
                       NNGP = FALSE,
                       cov.model = cov.model,
                       n.neighbors = 5,
                       n.batch = n.batch, 
                       batch.length = batch.length, 
                       tuning = tuning, 
                       n.burn = n.burn, 
                       n.thin = n.thin, 
                       n.chains = n.chains,
                       n.omp.threads = 35,
                       verbose = FALSE)
  
  # Extract WAIC
  waic_val <- waicOcc(model_fit)
  
  # Save to waic_df
  waic_df[i, ] <- c(buffer_val, 
                    waic_val["elpd"],
                    waic_val["pD"],
                    waic_val["WAIC"])
  
  # Save Rhat
  rhat_list[[i]] <- model_fit$rhat
}

# Check for convergence
Nriver_rhat_vals <- unlist(rhat_list)
which(Nriver_rhat_vals > 1.1)

# Sort by WAIC
Nriver_waic <- waic_df[order(waic_df$WAIC), ]
head(Nriver_waic, 25)

# Save WAIC 
write.csv(Nriver_waic, "./Outputs/Nriver_Buffer_WAIC.csv")


# ------------- end of script -------------