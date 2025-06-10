
# ------------------------------------------------------------------------------
#
#                         Load/Install Packages
#
# ------------------------------------------------------------------------------


library(tidyverse)
library(spOccupancy)
library(progress)

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
elev_MN_df <- readRDS("Data/OccuCovData/elev_MN_df.rds")
elev_SD_df <- readRDS("Data/OccuCovData/elev_SD_df.rds")
slope_MN_df <- readRDS("Data/OccuCovData/slope_MN_df.rds")
slope_SD_df <- readRDS("Data/OccuCovData/slope_SD_df.rds")
aspect_SINE_df <- readRDS("Data/OccuCovData/aspect_SINE_df.rds")
aspect_COS_df <- readRDS("Data/OccuCovData/aspect_COS_df.rds")
VRML_MN_df <- readRDS("Data/OccuCovData/VRML_MN_df.rds")
VRML_PRP_df <- readRDS("Data/OccuCovData/VRML_PRP_df.rds")
VRML_PD_df <- readRDS("Data/OccuCovData/VRML_PD_df.rds")
VRML_PA_df <- readRDS("Data/OccuCovData/VRML_PA_df.rds")
VRML_COH_df <- readRDS("Data/OccuCovData/VRML_COH_df.rds")
Npast_df <- readRDS("Data/OccuCovData/Npast_df.rds")
Nvil_df <- readRDS("Data/OccuCovData/Nvil_df.rds")
Nriver_df <- readRDS("Data/OccuCovData/Nriver_df.rds")


# Initial detection covariates
Days_Active <- readRDS("./Data/Days_Active.rds")

# Detection matrix
y_occ <- readRDS("./Data/y_occ.rds")

# Camera locations
cam_loc_df <- dist_df[,c(1:3)]


# ------------------------------------------------------------------------------
#
#                             Format Data
#
# ------------------------------------------------------------------------------

# Removing site, lat, and long
elev_MN_df <- elev_MN_df[,-c(1:3)] 
elev_SD_df <- elev_SD_df[,-c(1:3)] 
slope_MN_df <- slope_MN_df[,-c(1:3)] 
slope_SD_df <- slope_SD_df[,-c(1:3)] 
aspect_SINE_df <- aspect_SINE_df[,-c(1:3)] 
aspect_COS_df <- aspect_COS_df[,-c(1:3)] 
VRML_MN_df <- VRML_MN_df[,-c(1:3)]
VRML_PRP_df <- VRML_PRP_df[,-c(1:3)] 
VRML_PD_df <- VRML_PD_df[,-c(1:3)]
VRML_PA_df <- VRML_PA_df[,-c(1:3)]
VRML_COH_df <- VRML_COH_df[,-c(1:3)]
Npast_df <- Npast_df[,-c(1:3)]
Nvil_df <- Nvil_df[,-c(1:3)] 
Nriver_df <- Nriver_df[,-c(1:3)]


# Scale
elev_MN_df <- scale(elev_MN_df, center = TRUE, scale = TRUE)
elev_SD_df <- scale(elev_SD_df, center = TRUE, scale = TRUE)
slope_MN_df <- scale(slope_MN_df, center = TRUE, scale = TRUE)
slope_SD_df <- scale(slope_SD_df, center = TRUE, scale = TRUE)
aspect_SINE_df <- scale(aspect_SINE_df, center = TRUE, scale = TRUE)
aspect_COS_df <- scale(aspect_COS_df, center = TRUE, scale = TRUE)
VRML_MN_df <- scale(VRML_MN_df, center = TRUE, scale = TRUE)
VRML_PRP_df <- scale(VRML_PRP_df, center = TRUE, scale = TRUE)
VRML_PD_df <- scale(VRML_PD_df, center = TRUE, scale = TRUE)
VRML_PA_df <- scale(VRML_PA_df, center = TRUE, scale = TRUE)
VRML_COH_df <- scale(VRML_COH_df, center = TRUE, scale = TRUE)
# Npast_df <- scale(Npast_df, center = TRUE, scale = TRUE)
# Nvil_df <- scale(Nvil_df, center = TRUE, scale = TRUE)
# Nriver_df <- scale(Nriver_df, center = TRUE, scale = TRUE)

# ------------------------------------------------------------------------------
#
#                           Covariate scale
#
# ------------------------------------------------------------------------------


# ----------------
# Data objects
# ----------------

# Elevation Mean
spOccu_elev_MN <- list(y = y_occ,
                   occ.covs = elev_MN_df,
                   det.covs = list(Days_Active = Days_Active
                                   ),
                   coords = cam_loc_df[,c(2,3)])

# Elevation Standard Deviation
spOccu_elev_SD <- list(y = y_occ,
                       occ.covs = elev_SD_df,
                       det.covs = list(Days_Active = Days_Active
                       ),
                       coords = cam_loc_df[,c(2,3)])
# Slope Mean
spOccu_slope_MN <- list(y = y_occ,
                    occ.covs = slope_MN_df,
                    det.covs = list(Days_Active = Days_Active
                    ),
                    coords = cam_loc_df[,c(2,3)])

# Slope Standard Deviation
spOccu_slope_SD <- list(y = y_occ,
                        occ.covs = slope_SD_df,
                        det.covs = list(Days_Active = Days_Active
                        ),
                        coords = cam_loc_df[,c(2,3)])
# Aspect Sine
spOccu_aspect_SINE <- list(y = y_occ,
                     occ.covs = aspect_SINE_df,
                     det.covs = list(Days_Active = Days_Active
                     ),
                     coords = cam_loc_df[,c(2,3)])

# Aspect Cos
spOccu_aspect_COS <- list(y = y_occ,
                          occ.covs = aspect_COS_df,
                          det.covs = list(Days_Active = Days_Active
                          ),
                          coords = cam_loc_df[,c(2,3)])
# VRML Mean
spOccu_VRML_MN <- list(y = y_occ,
                      occ.covs = VRML_MN_df,
                      det.covs = list(Days_Active = Days_Active
                      ),
                      coords = cam_loc_df[,c(2,3)])
# VRML Proportion
spOccu_VRML_PRP <- list(y = y_occ,
                   occ.covs = VRML_PRP_df,
                   det.covs = list(Days_Active = Days_Active
                   ),
                   coords = cam_loc_df[,c(2,3)])

# VRML Patch Density
spOccu_VRML_PD <- list(y = y_occ,
                       occ.covs = VRML_PD_df,
                       det.covs = list(Days_Active = Days_Active
                       ),
                       coords = cam_loc_df[,c(2,3)])

# VRML Patch Area
spOccu_VRML_PA <- list(y = y_occ,
                      occ.covs = VRML_PA_df,
                      det.covs = list(Days_Active = Days_Active
                      ),
                      coords = cam_loc_df[,c(2,3)])

# VRML Cohesion Index
spOccu_VRML_COH <- list(y = y_occ,
                      occ.covs = VRML_COH_df,
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
batch.length <- 25
n.batch <- 10000
batch.length * n.batch # Total number of MCMC samples per chain
n.burn <- 100000
n.thin <- 20
n.chains <- 3


total_iterations <- batch.length * n.batch
retained_per_chain <- (total_iterations - n.burn) / n.thin
posterior_samples <- retained_per_chain * n.chains
posterior_samples


# Pair-wise distances between all sites
dist_mat <- dist(spOccu_elev_MN$coords)

# Exponential covariance model
cov_model <- "gaussian"

# Initial values
inits <- list(alpha = 0,
              beta = 0,
              z = apply(spOccu_elev_MN$y, 1, max, na.rm = TRUE),
              sigma.sq = 2,
              phi = 3 / median(dist_mat), #3 / mean(dist_mat),
              w = rep(0, nrow(spOccu_elev_MN$y))
)

# Tuning
tuning <- list(phi = 0.5
)

# Priors
min.dist <- min(dist_mat)
max.dist <- max(dist_mat)

priors <- list(beta.normal = list(mean = 0, var = 2.72),
               alpha.normal = list(mean = 0, var = 2.72),
               sigma.sq.ig = c(5, 1),
                c(3/max.dist, 3/min.dist)
)

# ---------------------
# Loop specifications
# ---------------------

# Define buffers
buffer_vec <- seq(250, 15000, by = 250)

# ---------------------
# Elevation MN
# ---------------------

# Create storage for WAIC 
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Create storage for Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "Elevation MN | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {

  # Get buffer size
  buffer_val <- buffer_vec[i]

  # Construct covariate name
  cov_name <- paste0("elev_MN_buffer", buffer_val)

  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))

  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_elev_MN,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))

}

# Check for convergence
elev_MN_rhat <- unlist(rhat_list)
which(elev_MN_rhat > 1.1)

# Sort by WAIC
elev_MN_waic <- waic_df[order(waic_df$WAIC), ]
head(elev_MN_waic, 25)

# Save WAIC
write.csv(elev_MN_waic, "./Outputs/Scale_WAIC/Elevation_MN_WAIC.csv")

# ---------------------
# Elevation SD
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "Elevation SD | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

colnames(spOccu_elev_SD$occ.covs) <- gsub("elev_MN", "elev_SD", colnames(spOccu_elev_SD$occ.covs))

# Start model fitting loop
for (i in 1:length(buffer_vec)) {
  
  # Get buffer size
  buffer_val <- buffer_vec[i]
  
  # Construct covariate name
  cov_name <- paste0("elev_SD_buffer", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_elev_SD,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
  
}

# Check for convergence
elev_SD_rhat <- unlist(rhat_list)
which(elev_SD_rhat > 1.1)

# Sort by WAIC
elev_SD_waic <- waic_df[order(waic_df$WAIC), ]
head(elev_SD_waic, 25)

# Save WAIC
write.csv(elev_SD_waic, "./Outputs/Scale_WAIC/Elevation_SD_WAIC.csv")

# ---------------------
# Slope MN
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()


# Progress bar
pb <- progress_bar$new(
  format = "Slope MN | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {
  
  # Get buffer size
  buffer_val <- buffer_vec[i]
  
  # Construct covariate name
  cov_name <- paste0("slope_MN_buffer", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_slope_MN,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
  
}

# Check for convergence
slope_MN_rhat <- unlist(rhat_list)
which(slope_MN_rhat > 1.1)

# Sort by WAIC
slope_MN_waic <- waic_df[order(waic_df$WAIC), ]
head(slope_MN_waic, 25)

# Save WAIC
write.csv(slope_MN_waic, "./Outputs/Scale_WAIC/Slope_MN_WAIC.csv")


# ---------------------
# Slope SD
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "Slope SD | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {
  
  # Get buffer size
  buffer_val <- buffer_vec[i]

  # Construct covariate name
  cov_name <- paste0("slope_SD_buffer", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_slope_SD,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
  
}

# Check for convergence
slope_SD_rhat <- unlist(rhat_list)
which(slope_SD_rhat > 1.1)

# Sort by WAIC
slope_SD_waic <- waic_df[order(waic_df$WAIC), ]
head(slope_SD_waic, 25)

# Save WAIC
write.csv(slope_SD_waic, "./Outputs/Scale_WAIC/Slope_SD_WAIC.csv")

# ---------------------
# Aspect Sine
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "Aspect Sine | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {
  
  # Get buffer size
  buffer_val <- buffer_vec[i]

  # Construct covariate name
  cov_name <- paste0("aspect_SINE_buffer", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_aspect_SINE,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
  
}

# Check for convergence
aspect_SINE_rhat <- unlist(rhat_list)
which(aspect_SINE_rhat > 1.1)

# Sort by WAIC
aspect_SINE_waic <- waic_df[order(waic_df$WAIC), ]
head(aspect_SINE_waic, 25)

# Save WAIC
write.csv(aspect_SINE_waic, "./Outputs/Scale_WAIC/Aspect_SINE_WAIC.csv")

# ---------------------
# Aspect Cosine
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "Aspect Cosine | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {
  
  # Get buffer size
  buffer_val <- buffer_vec[i]

  # Construct covariate name
  cov_name <- paste0("aspect_COS_buffer", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_aspect_COS,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
  
}

# Check for convergence
aspect_COS_rhat <- unlist(rhat_list)
which(aspect_COS_rhat > 1.1)

# Sort by WAIC
aspect_COS_waic <- waic_df[order(waic_df$WAIC), ]
head(aspect_COS_waic, 25)

# Save WAIC
write.csv(aspect_COS_waic, "./Outputs/Scale_WAIC/Aspect_COS_WAIC.csv")

# ---------------------
# VRML MN 
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "VRML MN | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {
  
  # Get buffer size
  buffer_val <- buffer_vec[i]

  # Construct covariate name
  cov_name <- paste0("VRML_MN_buffer", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_VRML_MN,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
  
}

# Check for convergence
VRML_MN_rhat <- unlist(rhat_list)
which(VRML_MN_rhat > 1.1)

# Sort by WAIC
VRML_MN_waic <- waic_df[order(waic_df$WAIC), ]
head(VRML_MN_waic, 25)

# Save WAIC
write.csv(VRML_MN_waic, "./Outputs/Scale_WAIC/VRML_MN_WAIC.csv")

# ---------------------
# VRML PRP
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "VRML PRP | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {
  
  # Get buffer size
  buffer_val <- buffer_vec[i]

  # Construct covariate name
  cov_name <- paste0("VRML_PRP_buffer", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_VRML_PRP,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
  
}

# Check for convergence
VRML_PRP_rhat <- unlist(rhat_list)
which(VRML_PRP_rhat > 1.1)

# Sort by WAIC
VRML_PRP_waic <- waic_df[order(waic_df$WAIC), ]
head(VRML_PRP_waic, 25)

# Save WAIC
write.csv(VRML_PRP_waic, "./Outputs/Scale_WAIC/VRML_PRP_WAIC.csv")

# ---------------------
# VRML PD
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "VRML PD | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {
  
  # Get buffer size
  buffer_val <- buffer_vec[i]
  
  # Construct covariate name
  cov_name <- paste0("VRML_PD_buffer", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_VRML_PD,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
  
}

# Check for convergence
VRML_PD_rhat <- unlist(rhat_list)
which(VRML_PD_rhat > 1.1)

# Sort by WAIC
VRML_PD_waic <- waic_df[order(waic_df$WAIC), ]
head(VRML_PD_waic, 25)

# Save WAIC
write.csv(VRML_PD_waic, "./Outputs/Scale_WAIC/VRML_PD_WAIC.csv")

# ---------------------
# VRML PA
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "VRML PA | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {
  
  # Get buffer size
  buffer_val <- buffer_vec[i]
  
  # Construct covariate name
  cov_name <- paste0("VRML_PA_buffer", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_VRML_PA,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
  
}

# Check for convergence
VRML_PA_rhat <- unlist(rhat_list)
which(VRML_PA_rhat > 1.1)

# Sort by WAIC
VRML_PA_waic <- waic_df[order(waic_df$WAIC), ]
head(VRML_PA_waic, 25)

# Save WAIC
write.csv(VRML_PA_waic, "./Outputs/Scale_WAIC/VRML_PA_WAIC.csv")

# ---------------------
# VRML COH
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "VRML COH | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {
  
  # Get buffer size
  buffer_val <- buffer_vec[i]
  
  # Construct covariate name
  cov_name <- paste0("VRML_COH_buffer", buffer_val)
  
  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))
  
  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_VRML_COH,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
  
}

# Check for convergence
VRML_COH_rhat <- unlist(rhat_list)
which(VRML_COH_rhat > 1.1)

# Sort by WAIC
VRML_COH_waic <- waic_df[order(waic_df$WAIC), ]
head(VRML_COH_waic, 25)

# Save WAIC
write.csv(VRML_COH_waic, "./Outputs/Scale_WAIC/VRML_COH_WAIC.csv")


# ---------------------
# N Pasture
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "N Pasture | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {

  # Get buffer size
  buffer_val <- buffer_vec[i]

  # Construct covariate name
  cov_name <- paste0("Npastures_buffer", buffer_val)

  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))

  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_Npast,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))

}

# Check for convergence
Npast_rhat_vals <- unlist(rhat_list)
which(Npast_rhat_vals > 1.1)

# Sort by WAIC
Npast_waic <- waic_df[order(waic_df$WAIC), ]
head(Npast_waic, 25)

# Save WAIC
write.csv(Npast_waic, "./Outputs/Scale_WAIC/Npast_Buffer_WAIC.csv")

# ---------------------
# N Villages
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "N Villages | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {

  # Get buffer size
  buffer_val <- buffer_vec[i]

  # Construct covariate name
  cov_name <- paste0("Nvillages_buffer", buffer_val)

  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))

  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_Nvil,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))

}

# Check for convergence
Nvil_rhat_vals <- unlist(rhat_list)
which(Nvil_rhat_vals > 1.1)

# Sort by WAIC
Nvil_waic <- waic_df[order(waic_df$WAIC), ]
head(Nvil_waic, 25)

# Save WAIC
write.csv(Nvil_waic, "./Outputs/Scale_WAIC/Nvil_Buffer_WAIC.csv")


# ---------------------
# N Rivers
# ---------------------

# Clear WAIC
waic_df <- data.frame(
  Buffer = rep(NA, length(buffer_vec)),
  elpd = rep(NA, length(buffer_vec)),
  pD = rep(NA, length(buffer_vec)),
  WAIC = rep(NA, length(buffer_vec))
)

# Clear Rhat
rhat_list <- list()

# Progress bar
pb <- progress_bar$new(
  format = "N Rivers | Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(buffer_vec),
  clear = FALSE,
  width = 100
)

# Start model fitting loop
for (i in 1:length(buffer_vec)) {

  # Get buffer size
  buffer_val <- buffer_vec[i]

  # Construct covariate name
  cov_name <- paste0("Nriver_buffer", buffer_val)

  # Create formula
  occ_formula <- as.formula(paste("~", cov_name))

  # Fit model
  model_fit <- spPGOcc(occ.formula = occ_formula,
                       det.formula = ~ Days_Active,
                       data = spOccu_Nriver,
                       inits = inits,
                       priors = priors,
                       NNGP = FALSE,
                       cov_model = cov_model,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       tuning = tuning,
                       n.burn = n.burn,
                       n.thin = n.thin,
                       n.chains = n.chains,
                       n.omp.threads = 6,
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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))
}

# Check for convergence
Nriver_rhat_vals <- unlist(rhat_list)
which(Nriver_rhat_vals > 1.1)

# Sort by WAIC
Nriver_waic <- waic_df[order(waic_df$WAIC), ]
head(Nriver_waic, 25)

# Save WAIC
write.csv(Nriver_waic, "./Outputs/Scale_WAIC/Nriver_Buffer_WAIC.csv")


# ------------- End of Script -------------