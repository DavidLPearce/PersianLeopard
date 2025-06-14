
# ------------------------------------------------------------------------------
#
#                         Load/Install Packages
#
# ------------------------------------------------------------------------------

library(tidyverse)
library(spOccupancy)
library(progress)
library(reshape2) # used for heatmap
library(psych) # used for pairs panel plotting

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

# WAIC Tables
elev_MN_waic <- read.csv("./Outputs/Scale_WAIC/Elevation_MN_WAIC.csv")
elev_SD_waic <- read.csv("./Outputs/Scale_WAIC/Elevation_SD_WAIC.csv")
slope_MN_waic <- read.csv("./Outputs/Scale_WAIC/Slope_MN_WAIC.csv")
slope_SD_waic <- read.csv("./Outputs/Scale_WAIC/Slope_SD_WAIC.csv")
aspect_SINE_waic <- read.csv("./Outputs/Scale_WAIC/Aspect_SINE_WAIC.csv")
aspect_COS_waic <- read.csv("./Outputs/Scale_WAIC/Aspect_COS_WAIC.csv")
VRML_MN_waic <- read.csv("./Outputs/Scale_WAIC/VRML_MN_WAIC.csv")
VRML_PRP_waic <- read.csv("./Outputs/Scale_WAIC/VRML_PRP_WAIC.csv")
VRML_PD_waic <- read.csv("./Outputs/Scale_WAIC/VRML_PD_WAIC.csv")
VRML_PA_waic <- read.csv("./Outputs/Scale_WAIC/VRML_PA_WAIC.csv")
VRML_COH_waic <- read.csv("./Outputs/Scale_WAIC/VRML_COH_WAIC.csv")
Npast_waic <- read.csv("./Outputs/Scale_WAIC/Npast_Buffer_WAIC.csv")
Nvil_waic <- read.csv("./Outputs/Scale_WAIC/Nvil_Buffer_WAIC.csv")
Nriver_waic <- read.csv("./Outputs/Scale_WAIC/Nriver_Buffer_WAIC.csv")

# Initial detection covariates
Days_Active <- readRDS("./Data/Days_Active.rds")

# Detection matrix
y_occ <- readRDS("./Data/y_occ.rds")

# Camera locations
cam_loc_df <- dist_df[,c(1:3)]


# ------------------------------------------------------------------------------
#
#                               Scale and Correlations
#
# ------------------------------------------------------------------------------

# Removing site, lat, and long
dist_df <- dist_df[,-c(1:3)] 
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

# ---------------------------
# "Best" Buffer Sizes
# ---------------------------
head(elev_MN_waic, 5)
head(elev_SD_waic, 5)
head(slope_MN_waic, 5)
head(slope_SD_waic, 5)
head(aspect_SINE_waic, 5)
head(aspect_COS_waic, 5)
head(VRML_MN_waic, 5)
head(VRML_PRP_waic, 5)
head(VRML_PD_waic, 5)
head(VRML_PA_waic, 5)
head(VRML_COH_waic, 5)
head(Npast_waic, 5)
head(Nvil_waic, 5)
head(Nriver_waic, 5)

# Buffer column index
elev_MN_idx <- elev_MN_waic[1,"X"]
elev_SD_idx <- elev_SD_waic[1,"X"]
slope_MN_idx <- slope_MN_waic[1,"X"]
slope_SD_idx <- slope_SD_waic[1,"X"]
aspect_SINE_idx <- aspect_SINE_waic[1,"X"]
aspect_COS_idx <- aspect_COS_waic[1,"X"]
VRML_MN_idx <- VRML_MN_waic[1,"X"]
VRML_PRP_idx <- VRML_PRP_waic[1,"X"]
VRML_PD_idx <- VRML_PD_waic[1,"X"]
VRML_PA_idx <- VRML_PA_waic[1,"X"]
VRML_COH_idx <- VRML_COH_waic[1,"X"]
Npast_idx <- Npast_waic[1,"X"]
Nvil_idx <- Nvil_waic[1,"X"]
Nriver_idx <- Nriver_waic[1,"X"]

# Extract "best" buffer size
best_elev_MN <- as.data.frame(elev_MN_df[ , elev_MN_idx]) 
best_elev_SD <- as.data.frame(elev_SD_df[ , elev_SD_idx])  
best_slope_MN <- as.data.frame(slope_MN_df[ , slope_MN_idx])  
best_slope_SD <- as.data.frame(slope_SD_df[ , slope_SD_idx])  
best_aspect_SINE <- as.data.frame(aspect_SINE_df[ , aspect_SINE_idx])  
best_aspect_COS <- as.data.frame(aspect_COS_df[ , aspect_COS_idx])  
best_VRML_MN <- as.data.frame(VRML_MN_df[ , VRML_MN_idx])  
best_VRML_PRP <- as.data.frame(VRML_PRP_df[ , VRML_PRP_idx])  
best_VRML_PD <- as.data.frame(VRML_PD_df[ , VRML_PD_idx])  
best_VRML_PA <- as.data.frame(VRML_PA_df[ , VRML_PA_idx]) 
best_VRML_COH <- as.data.frame(VRML_COH_df[ , VRML_COH_idx])  
best_Npast <- as.data.frame(Npast_df[ , Npast_idx])  
best_Nvil <- as.data.frame(Nvil_df[ , Nvil_idx])
best_Nriver <- as.data.frame(Nriver_df[ , Nriver_idx])  

# Renaming columns
colnames(best_elev_MN) <- colnames(elev_MN_df)[elev_MN_idx]
colnames(best_elev_SD) <- colnames(elev_SD_df)[elev_SD_idx]  
colnames(best_slope_MN) <- colnames(slope_MN_df)[slope_MN_idx]  
colnames(best_slope_SD) <- colnames(slope_SD_df)[slope_SD_idx] 
colnames(best_aspect_SINE) <- colnames(aspect_SINE_df)[aspect_SINE_idx]  
colnames(best_aspect_COS) <- colnames(aspect_COS_df)[aspect_COS_idx] 
colnames(best_VRML_MN) <- colnames(VRML_MN_df)[VRML_MN_idx]  
colnames(best_VRML_PRP) <- colnames(VRML_PRP_df)[VRML_PRP_idx]  
colnames(best_VRML_PD) <- colnames(VRML_PD_df)[VRML_PD_idx]  
colnames(best_VRML_PA) <- colnames(VRML_PA_df)[VRML_PA_idx]
colnames(best_VRML_COH) <- colnames(VRML_COH_df)[VRML_COH_idx]  
colnames(best_Npast) <- colnames(Npast_df)[Npast_idx] 
colnames(best_Nvil) <- colnames(Nvil_df)[Nvil_idx]
colnames(best_Nriver) <- colnames(Nriver_df)[Nriver_idx] 


# Combine all covariates into one dataframe
cov_df <- cbind(dist_df,
                best_elev_MN,
                best_elev_SD,
                best_slope_MN,
                best_slope_SD,  
                best_aspect_SINE,
                best_aspect_COS,  
                best_VRML_MN,  
                best_VRML_PRP, 
                best_VRML_PD, 
                best_VRML_PA, 
                best_VRML_COH,  
                best_Npast, 
                best_Nvil,
                best_Nriver
)

# Scale
cov_df <- scale(cov_df, center = TRUE, scale = TRUE)

# ---------------------------
# Check Correlations
# ---------------------------

# Create correlation matrix
cov_corr_mat <- cor(cov_df)  

# Correlations with values
pairs.panels(cov_corr_mat, gap = 0, bg = c("blue", "red"), pch = 21, main = "") # view
png("./Outputs/Correlation_Mat.png", width = 1800, height = 1200, res = 150) # export
pairs.panels(cov_corr_mat, gap = 0, bg = c("blue", "red"), pch = 21, main = "")
dev.off() # for export
dev.off() # for plot panel

# Format for heat map
cov_corr_mat_melted <- melt(cov_corr_mat)

# Heat Map
cov_corr_heatmap <- ggplot(cov_corr_mat_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(hjust = 1, size = 14)) +
  labs(title = "", x = "", y = "")

# View
print(cov_corr_heatmap)

# Save
ggsave("./Outputs/Correlation_Heatmap.jpg",
       plot = cov_corr_heatmap, width = 14, height = 12, dpi = 300)
dev.off()  


# Diagional matrix with strength of correlation
corrplot(cov_corr_mat, 
         type = "lower", 
         diag = FALSE, 
         mar = c(1,0.5,5,1), 
         tl.col = 'black',
         tl.pos = 'ld', 
         tl.srt = 45, 
         xpd = TRUE, 
         main = "Pearson's correlations")



# ---------------------------
# Formulas
# ---------------------------

# Spliting the data into two data sets
# 1 for nearest distance to
# 1 for number of X in buffer
View(cov_df)
dist_cov_df <- cov_df[,-c(15:17)]
num_cov_df <- cov_df[,-c(1:3)]

# Run correlations on both df
dist_corr_mat <- cor(dist_cov_df)
num_corr_mat <- cor(num_cov_df)

# Covariate names
dist_covs <- colnames(dist_cov_df)
num_covs <- colnames(num_corr_mat)

# Individual covariats
dist_1formulas <- paste("~", dist_covs)
num_1formulas <- paste("~", num_covs)

# Get all 2-way combinations
dist_2combos <- combn(seq_along(dist_covs), 2, simplify = FALSE)
num_2combos <- combn(seq_along(num_covs), 2, simplify = FALSE)

# Filter combinations with correlation ≤ 0.40
dist_low_2combos <- Filter(function(idx) {
  abs(dist_corr_mat[idx[1], idx[2]]) <= 0.40
}, dist_2combos)

num_low_2combos <- Filter(function(idx) {
  abs(num_corr_mat[idx[1], idx[2]]) <= 0.40
}, num_2combos)

# Convert to formulas
dist_2formulas <- sapply(dist_low_2combos, function(idx) {
  paste("~", paste(dist_covs[idx], collapse = " + "))
})

num_2formulas <- sapply(num_low_2combos, function(idx) {
  paste("~", paste(num_covs[idx], collapse = " + "))
})



# Get all 3-way combinations
dist_3combos <- combn(seq_along(dist_covs), 3, simplify = FALSE)
num_3combos <- combn(seq_along(num_covs), 3, simplify = FALSE)

dist_low_2combos <- Filter(function(idx) {
  abs(dist_corr_mat[idx[1], idx[2]]) <= 0.40
}, dist_2combos)


# Filter if all pairwise correlations are ≤ 0.40
dist_low_3combos <- Filter(function(idx) {
  i <- idx[1]; j <- idx[2]; k <- idx[3]
  all(abs(dist_corr_mat[c(i,j,k), c(i,j,k)][lower.tri(matrix(0, 3, 3))]) <= 0.40)
}, dist_3combos)

num_low_3combos <- Filter(function(idx) {
  i <- idx[1]; j <- idx[2]; k <- idx[3]
  all(abs(num_corr_mat[c(i,j,k), c(i,j,k)][lower.tri(matrix(0, 3, 3))]) <= 0.40)
}, num_3combos)


# Convert to formulas
dist_3formulas <- sapply(dist_low_3combos, function(idx) {
  paste("~", paste(dist_covs[idx], collapse = " + "))
})

num_3formulas <- sapply(num_low_3combos, function(idx) {
  paste("~", paste(num_covs[idx], collapse = " + "))
})


# Combining formulas
all_1formulas <- union(dist_1formulas, num_1formulas)
all_2formulas <- union(dist_2formulas, num_2formulas)
all_3formulas <- union(dist_3formulas, num_3formulas)

# Combine into one object
formulas <- c(all_1formulas, all_2formulas, all_3formulas)

# Add null model
formulas <- c("~ 1", formulas)

# Take a look
print(formulas)

# ------------------------------------------------------------------------------
#
#                               Candidate Models
#
# ------------------------------------------------------------------------------

# Create spOccupancy data
spOccu_dat <- list(y = y_occ,
                    occ.covs = cov_df,
                    det.covs = list(Days_Active = Days_Active
                    ),
                    coords = cam_loc_df[,c(2,3)])

# Export
saveRDS(spOccu_dat, "Data/spOccu_dat.rds")

# ---------------------------
# MCMC Specifications
# ---------------------------

# MCMC
batch.length <- 25
n.batch <- 30000
batch.length * n.batch # Total number of MCMC samples per chain
n.burn <- 200000
n.thin <- 20
n.chains <- 3


total_iterations <- batch.length * n.batch
retained_per_chain <- (total_iterations - n.burn) / n.thin
posterior_samples <- retained_per_chain * n.chains
posterior_samples


# Pair-wise distances between all sites
dist_mat <- dist(spOccu_dat$coords)

# Exponential covariance model
cov_model <- "gaussian"

# Initial values
inits <- list(alpha = 0,
              beta = 0,
              z = apply(spOccu_dat$y, 1, max, na.rm = TRUE),
              sigma.sq = 2,
              phi = 3 / median(dist_mat), #3 / mean(dist_mat),
              w = rep(0, nrow(spOccu_dat$y))
)

# Tuning
tuning <- list(phi = 0.5
)

# Priors
min_dist <- min(dist_mat)
max_dist <- max(dist_mat)

priors <- list(beta.normal = list(mean = 0, var = 2.72),
               alpha.normal = list(mean = 0, var = 2.72),
               sigma.sq.ig = c(5, 1),
               c(3/max_dist, 3/min_dist)
)


# ---------------------------
# Loop Objects
# ---------------------------

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

# ---------------------------
# Candidate Model Loop
# ---------------------------

# Progress bar
pb <- progress_bar$new(
  format = "Model :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = length(formulas),
  clear = FALSE,
  width = 100
)

# Fit a model for each formula
for (i in seq_along(formulas)) {

  # Extract formula
  f <- formulas[[i]]

  # Format from string to formula
  occ_formula <- as.formula(f)

  # Fit model
  fit <- spPGOcc(occ.formula = occ_formula,
                 det.formula = ~ scale(Days_Active),
                 data = spOccu_dat,
                 inits = inits,
                 priors = priors,
                 NNGP = FALSE,
                 cov.model = cov_model,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 tuning = tuning,
                 n.burn = n.burn,
                 n.thin = n.thin,
                 n.chains = n.chains,
                 verbose = FALSE)

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
  
  # Update progress bar
  pb$tick(tokens = list(current = i))

} # End Loop

# Check for convergence
rhat_vals <- unlist(rhat_list)
which(rhat_vals > 1.1)

# Rank by WAIC
waic_df <- waic_df[order(waic_df$WAIC), ]
head(waic_df,5)

# Check model fit of best occupancy model
ppc_out_ft <- ppcOcc(fit_list[[waic_df[1,"Model"]]], fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out_ft) # 0.7228

ppc_out_ch2 <- ppcOcc(fit_list[[waic_df[1,"Model"]]], fit.stat = 'chi-squared', group = 1)
summary(ppc_out_ch2) # 0.64

# Best occupancy model
bm_fit <- fit_list[[waic_df[1,"Model"]]]
summary(bm_fit)

# Save occupancy WAIC
write.csv(waic_df, "./Outputs/Candidate_Model_Ranking.csv")

# Save best occupancy model
saveRDS(bm_fit, "./Outputs/Best_Candidate_Model.rds")

# ------------------- End of Script -------------------
