
# ------------------------------------------------------------------------------
#
#                         Load/Install Packages
#
# ------------------------------------------------------------------------------

library(tidyverse)
library(gridExtra)
library(spOccupancy)
library(progress)
library(terra)
library(coda)
library(stars)

set.seed(123)
options(scipen = 9999)
setwd(".")

# ------------------------------------------------------------------------------
#
#                             Read in data
#
# ------------------------------------------------------------------------------


# Best model
model_fit <- readRDS("./Outputs/Best_Candidate_Model.rds")


# Prediction covs
pred_cov_df <- readRDS("Data/OccuCovData/pred_cov_df.rds")
nrow(pred_cov_df)


# ------------------------------------------------------------------------------
#
#                               Prediction
#
# ------------------------------------------------------------------------------

# ---------------------------
# Coefficient Estimates
# ---------------------------

# Model summary
summary(model_fit)
str(model_fit)

# Extract samples
beta_df <- as.data.frame(model_fit$beta.samples)
alpha_df <- as.data.frame(model_fit$alpha.samples)

# Rename columns
colnames(alpha_df) <- c("Detection Intercept", "Days Active") # Detection

colnames(beta_df) <- c("Spatial Use Intercept",                 # Occupancy
                       "Aspect Sine", 
                       "Aspect Cosine", 
                       "VRML Cohesion Index")

# Reshape 
alpha_long <- alpha_df %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value") %>%
  mutate(Type = "Detection")

beta_long <- beta_df %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value") %>%
  mutate(Type = "Spatial Use")

# Combine and set factor levels 
param_data <- bind_rows(alpha_long, 
                        beta_long) %>%
  mutate(
    Parameter = factor(Parameter, levels = rev(c(
      "Detection Intercept", "Days Active", 
      "Spatial Use Intercept", "Aspect Sine", 
      "Aspect Cosine", "VRML Cohesion Index"
    ))),
    Type = factor(Type, levels = c("Detection", 
                                   "Spatial Use"))
  )

# Filter to 95% central interval for each parameter
param_data <- param_data %>%
  group_by(Parameter) %>%
  filter(Value > quantile(Value, 0.025) & Value < quantile(Value, 0.975)) %>%
  ungroup()


# Define fill colors
fill_colors <- c("Detection" = "red", 
                 "Spatial Use" = "blue")

# Coefficient Plot
coef_plot <- ggplot(param_data, aes(x = Value, y = Parameter, fill = Type)) +
                geom_violin(scale = "width", trim = TRUE, alpha = 0.5) +
                geom_boxplot(width = 0.1, fill = "white", outlier.shape = "black", alpha = 0.7) +  # Boxplot layer
                geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
                scale_x_continuous(
                  limits = c(-5, 5), 
                  breaks = seq(-5, 5, 1), 
                  expand = expansion(mult = c(0.05, 0.05))
                ) +
                scale_y_discrete(
                  expand = expansion(mult = c(0.15, 0.15))
                ) +
                scale_fill_manual(values = fill_colors, guide = "none") +
                labs(x = "Estimate", 
                     y = "Parameter") +
                theme_minimal(base_size = 14) +
                theme(
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "white", color = NA),
                  plot.background = element_rect(fill = "white", color = NA),
                  axis.line = element_line(color = "black")
                )

# View
print(coef_plot)

# Save
ggsave("./Outputs/Coefficient_Estimate_Plot.jpg", 
       plot = coef_plot, width = 10, height = 6, dpi = 300)

# ---------------------------
# Spatial Parameters
# ---------------------------

# Extract samples
theta_df <- as.data.frame(model_fit$theta.samples)

# Rename columns
colnames(theta_df) <- c("Spatial Variance", "Spatial Decay")

# Pivot longer
theta_long <- theta_df %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")

# Order
theta_long$Parameter <- factor(theta_long$Parameter, levels = c("Spatial Variance", "Spatial Decay"))

# 95% credible interval
theta_long <- theta_long %>%
  group_by(Parameter) %>%
  filter(
    Value >= quantile(Value, 0.025),
    Value <= quantile(Value, 0.975)
  ) %>%
  ungroup()

# Separate data for each parameter
spatial_variance <- theta_long %>% filter(Parameter == "Spatial Variance")
spatial_decay <- theta_long %>% filter(Parameter == "Spatial Decay")

# Spatial Variance
sp_plot1 <- ggplot(spatial_variance, aes(x = Value, fill = Parameter)) +
              geom_density(alpha = 0.5, color = NA) +
              scale_fill_manual(values = c("Spatial Variance" = "orange")) +
              scale_x_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
              labs(x = "", 
                   y = "Density",
                   title = "A)") +
              theme_minimal(base_size = 12) +
              theme(
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white", color = NA),
                axis.line = element_line(color = "black")
              )

# Spatial Decay
sp_plot2 <- ggplot(spatial_decay, aes(x = Value, fill = Parameter)) +
              geom_density(alpha = 0.5, color = NA) +
              scale_fill_manual(values = c("Spatial Decay" = "purple")) +
              scale_x_continuous(limits = c(0, 0.002), breaks = seq(0, 0.002, 0.0004)) +
              labs(x = "Estimate", 
                   y = "Density", 
                   title = "B)") +
              theme_minimal(base_size = 12) +
              theme(
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white", color = NA),
                axis.line = element_line(color = "black")
              )

# Multipanel figure
grid.arrange(sp_plot1, sp_plot2, ncol = 1)


# Save
jpeg("./Outputs/Spatial_Param_Plot.jpg", width = 6, height = 8, units = "in", res = 300)
grid.arrange(sp_plot1, sp_plot2, ncol = 1)
dev.off()


# ---------------------------
# Parameter/Coefficient Table
# ---------------------------

# Combine coefficient and spatial parameter estimates

# Add type label to theta_long to match others
theta_summary <- theta_long %>%
  mutate(Type = "Spatial")

# Combine all parameters
all_params <- bind_rows(param_data, theta_summary)

# Summarize
summary_table <- all_params %>%
  group_by(Type, Parameter) %>%
  summarise(
    Mean = mean(Value),
    Median = median(Value),
    LCI = quantile(Value, 0.025),
    UCI = quantile(Value, 0.975),
    .groups = "drop"
  ) %>%
  arrange(factor(Type, levels = c("Occupancy", "Detection", "Spatial")))

print(summary_table)

# Export
write.csv(summary_table , "Outputs/Parameter_Estimate_Table.csv")

 
# ---------------------------
# Predicted Covariate Effects
# ---------------------------

# Scale predicted values by mean and sd used to fit the model
pred_cov_df[,4:6] <- scale(pred_cov_df[,4:6] ,center = TRUE, scale = TRUE)

# Extract occupancy beta estimates
beta0_samples <- beta_df[,"Spatial Use Intercept"]       # Intercept
beta1_samples <- beta_df[,"Aspect Sine"]   # aspect_SINE_buffer4000 
beta2_samples <- beta_df[,"Aspect Cosine"] # aspect_COS_buffer5000
beta3_samples <- beta_df[,"VRML Cohesion Index"]       # VRML_COH_buffer750

# Create a prediction of aspect_SINE_buffer4000 values
cov1_pred_vals <- seq(min(pred_cov_df[, "aspect_SINE_buffer4000"]), max(pred_cov_df[, "aspect_SINE_buffer4000"]), length.out = 1000)
cov2_pred_vals <- seq(min(pred_cov_df[, "aspect_COS_buffer5000"]), max(pred_cov_df[, "aspect_COS_buffer5000"]), length.out = 1000)
cov3_pred_vals <- seq(min(pred_cov_df[, "VRML_COH_buffer750"]), max(pred_cov_df[, "VRML_COH_buffer750"]), length.out = 1000)

# Matrices for storing predictions
cov1_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov1_pred_vals))
cov2_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov2_pred_vals))
cov3_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov3_pred_vals))

# Generate predictions
for (i in 1:length(beta0_samples)) {
  cov1_preds[i, ] <- beta0_samples[i] + beta1_samples[i] * cov1_pred_vals  
  cov2_preds[i, ] <- beta0_samples[i] + beta2_samples[i] * cov2_pred_vals
  cov3_preds[i, ] <- beta0_samples[i] + beta3_samples[i] * cov3_pred_vals
  
}


# Calculate mean predictions
cov1_preds_mean <- apply(cov1_preds, 2, mean)# mean
cov1_preds_LCI <- apply(cov1_preds, 2, quantile, probs = 0.025) # LCI
cov1_preds_HCI <- apply(cov1_preds, 2, quantile, probs = 0.975) # HCI

cov2_preds_mean <- apply(cov2_preds, 2, mean)
cov2_preds_LCI <- apply(cov2_preds, 2, quantile, probs = 0.025)
cov2_preds_HCI <- apply(cov2_preds, 2, quantile, probs = 0.975)

cov3_preds_mean <- apply(cov3_preds, 2, mean)
cov3_preds_LCI <- apply(cov3_preds, 2, quantile, probs = 0.025)
cov3_preds_HCI <- apply(cov3_preds, 2, quantile, probs = 0.975)

# Combine into a single data frame
cov1_pred_df <- data.frame(
  cov1_scaled = cov1_pred_vals,
  cov1_preds_mean = cov1_preds_mean,
  cov1_preds_LCI = cov1_preds_LCI,
  cov1_preds_HCI = cov1_preds_HCI)

cov2_pred_df <- data.frame(
  cov2_scaled = cov2_pred_vals,
  cov2_preds_mean = cov2_preds_mean,
  cov2_preds_LCI = cov2_preds_LCI,
  cov2_preds_HCI = cov2_preds_HCI)

cov3_pred_df <- data.frame(
  cov3_scaled = cov3_pred_vals,
  cov3_preds_mean = cov3_preds_mean,
  cov3_preds_LCI = cov3_preds_LCI,
  cov3_preds_HCI = cov3_preds_HCI)

# Aspect Sine 
cov1_pred_plot <- ggplot(cov1_pred_df, aes(x = cov1_scaled, y = cov1_preds_mean)) +
                    geom_ribbon(aes(ymin = cov1_preds_LCI, ymax = cov1_preds_HCI), 
                                fill = "gray70", alpha = 0.4) +
                    geom_line(color = "black", linewidth = 1.2) +
                    scale_x_continuous(limits = c(-2, 3), breaks = seq(-2, 3, 1)) +
                    scale_y_continuous(limits = c(-12, 12), breaks = seq(-12, 12, 4)) +
                    labs(
                      x = "", 
                      y = "Predicted Effect", 
                      title = "A)"
                    ) +
                    theme_minimal(base_size = 12) +
                    theme(
                      panel.grid = element_blank(),
                      axis.line = element_line(color = "black"),
                      plot.title = element_text(face = "bold", hjust = 0),
                      plot.background = element_rect(fill = "white", color = NA),
                      panel.background = element_rect(fill = "white", color = NA)
                    )
print(cov1_pred_plot)

# Aspect Cosine 
cov2_pred_plot <- ggplot(cov2_pred_df, aes(x = cov2_scaled, y = cov2_preds_mean)) +
  geom_ribbon(aes(ymin = cov2_preds_LCI, ymax = cov2_preds_HCI), 
              fill = "gray70", alpha = 0.4) +
  geom_line(color = "black", linewidth = 1.2) +
  scale_x_continuous(limits = c(-3.5, 3), breaks = seq(-3, 3, 1)) +
  scale_y_continuous(limits = c(-10, 14), breaks = seq(-10, 14, 4)) +
  labs(
    x = "", 
    y = "Predicted Effect", 
    title = "B)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
print(cov2_pred_plot)

# VRML Cohesion Index
cov3_pred_plot <- ggplot(cov3_pred_df, aes(x = cov3_scaled, y = cov3_preds_mean)) +
  geom_ribbon(aes(ymin = cov3_preds_LCI, ymax = cov3_preds_HCI), 
              fill = "gray70", alpha = 0.4) +
  geom_line(color = "black", linewidth = 1.2) +
  scale_x_continuous(limits = c(-2.5, 3.5), breaks = seq(-3, 3, 1)) +
  scale_y_continuous(limits = c(-10, 15), breaks = seq(-10, 15, 4)) +
  labs(
    x = "Covariate (scaled)", 
    y = "Predicted Effect", 
    title = "C)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
print(cov3_pred_plot)

# Multipanel figure
grid.arrange(cov1_pred_plot, cov2_pred_plot, cov3_pred_plot, ncol = 1)


# Save
jpeg("./Outputs/OccuCov_Prediction_Plot.jpg", width = 6, height = 8, units = "in", res = 300)
grid.arrange(cov1_pred_plot, cov2_pred_plot, cov3_pred_plot, ncol = 1)
dev.off()


# ---------------------------
# Prediction Maps
# ---------------------------

## *** Memory allocation not enough on laptop, ran on TAMU HPRC Faster Cluster ***

# # Matrix of prediction
# X.0 <- as.matrix(data.frame(intercept = median(model_fit$beta.samples[,1]), # intercept median is 1.173018
#                             OccuCov1 = pred_cov_df$aspect_SINE_buffer4000, 
#                             OccuCov2 = pred_cov_df$aspect_COS_buffer5000,
#                             OccuCov3 = pred_cov_df$VRML_COH_buffer750
# ))
# 
# # Predict
# pred_out <- predict(model_fit,
#                     X.0,
#                     coords.0 = pred_cov_df[,2:3],
#                     ignore.RE = FALSE,
#                     type = 'occupancy',
#                     verbose = TRUE,
#                     n.report = 1000,
#                     n.omp.threads = 5
# )
# 
# # Combine into dataframe
# occu_map_dat <- data.frame(psi_mean = apply(pred_out$psi.0.samples, 2, mean),
#                             psi_sd = apply(pred_out$psi.0.samples, 2, sd),
#                             x = pred_cov_df[,'x'],
#                            y = pred_cov_df[,'y']
# )
# 
# # Convert to stars
# dat.stars <- st_as_stars(occu_map_dat, dims = c('x', 'y'))

# Read in dat.stars that was ran on TAMU HPRC
dat.stars <- readRDS("Data/dat.stars_500x500_study_area.rds")


# # Convert and export as a raster
# occu_rast <- rast(dat.stars)   # As raster
# crs(occu_rast) <- "EPSG:32640" # Correct CRS to UTM zone 40N
# writeRaster(occu_rast,         # Export
#             "E:/Persian_Leopard_GIS/LeopardOccupancyRaster/occu_rast.tif",
#             filetype = "GTiff", overwrite = TRUE)


# Mean Occupancy
MNoccu_plot <- ggplot() + 
  geom_stars(data = dat.stars, aes(x = x, y = y, fill = psi_mean)) +
  scale_fill_viridis_c(
    na.value = 'transparent',
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    name = "Mean"
  ) +
  labs(
    x = 'Easting',
    y = 'Northing',
    title = 'A)'
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    panel.border = element_blank(),  # Remove all borders
    axis.line = element_line(color = "black")  # Add bottom & left axis lines only
  )

print(MNoccu_plot)

# Standard Deviation Occupancy
SDoccu_plot <- ggplot() + 
  geom_stars(data = dat.stars, aes(x = x, y = y, fill = psi_sd)) +
  scale_fill_viridis_c(
    na.value = 'transparent',
    name = "Standard \nDeviation"
  ) +
  labs(
    x = 'Easting',
    y = '',
    title = 'B)'
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    panel.border = element_blank(),  # Remove all borders
    axis.line = element_line(color = "black")  # Add bottom & left axis lines only
  )



print(SDoccu_plot)


# Multipanel figure
grid.arrange(MNoccu_plot, SDoccu_plot, ncol = 2)


# Save
jpeg("./Outputs/Occupancy_Map.jpg", width = 12, height = 6, units = "in", res = 300)
grid.arrange(MNoccu_plot, SDoccu_plot, ncol = 2)
dev.off()


# Contoured Prediction Maps ----------------

MNoccu_contour_plot <- ggplot() +
  geom_contour_filled(
    data = as.data.frame(dat.stars), 
    aes(x = x, y = y, z = psi_mean), 
    bins = 4
  ) +
  scale_fill_viridis_d(
    name = "Mean",
    labels = function(x) gsub(" ", "", gsub(",", "–", gsub("\\(|\\[|\\]", "", x)))
  ) +
  coord_equal() +
  guides(fill = guide_legend(reverse = TRUE)) +   
  labs(
    x = "Easting",
    y = "Northing",
    title = "A)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),   
    axis.title = element_text(face = "bold"), 
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",   
    panel.border = element_blank(),  
    axis.line = element_line(color = "black")   
  )

print(MNoccu_contour_plot)

# SD
SDoccu_contour_plot <- ggplot() +
  geom_contour_filled(
    data = as.data.frame(dat.stars), 
    aes(x = x, y = y, z = psi_sd), 
    bins = 4
  ) +
  scale_fill_viridis_d(
    name = "Standard\nDeviation",
    labels = function(x) gsub(" ", "", gsub(",", "–", gsub("\\(|\\[|\\]", "", x))) 
  ) +
  coord_equal() +
  guides(fill = guide_legend(reverse = TRUE)) +   
  labs(
    x = "Easting",
    y = "",
    title = "B)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),   
    axis.title = element_text(face = "bold"), 
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",   
    panel.border = element_blank(),  
    axis.line = element_line(color = "black")   
  )

print(SDoccu_contour_plot)

# Multipanel figure
grid.arrange(MNoccu_contour_plot, SDoccu_contour_plot, ncol = 2)


# Save
jpeg("./Outputs/Occupancy_Contour_Map.jpg", width = 12, height = 6, units = "in", res = 300)
grid.arrange(MNoccu_contour_plot, SDoccu_contour_plot, ncol = 2)
dev.off()






# ------------------- End of Script -------------------