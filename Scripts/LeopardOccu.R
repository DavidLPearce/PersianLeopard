# -------------------------------------------------------
#
#                    Load/install Packages
#
# -------------------------------------------------------


library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(unmarked)
library(AICcmodavg)


set.seed(123)
options(scipen = 9999)
setwd(".")


# -------------------------------------------------------
#
#                    Read In data
#
# -------------------------------------------------------

# reading in occupancy data
leo_dat <- read.csv("./LeopardOccu.csv")

# subsetting out site name and ID
y <- leo_dat[,-c(1,2)]

# creating an unmarked frame
leo_umf <- unmarkedFrameOccu(y, 
                             siteCovs = NULL, 
                             obsCovs = NULL)


# -------------------------------------------------------
#
#         Naive Cumulative Detection Probability
#
# -------------------------------------------------------



fit1 <- occu( data = leo_umf,
              ~1 # Detection
              ~1,# Occupancy
              )


# Predicting detection using best detection model
# Whitetail deer
leo_daily_p <- predict(fit1, type = "det", species = "Leopard")


# Extract the predicted probabilities
leo_p_det <- leo_daily_p$Predicted

# Number of survey days
n_days <- 181



# Calculate cumulative detection probability over 28 days
# Cumulative detection probability for a single species over n days is:
# P_cumulative = 1 - (1 - p)^n_days
leo_cumulative_p <- 1 - cumprod(1 - leo_p_det[1:n_days])


# Create a data frame for plotting
plot_data <- data.frame(
  Day = 1:n_days,
  Leo_Cumulative = leo_cumulative_p
)

# Plot the cumulative detection probabilities
ggplot(plot_data, aes(x = Day)) +
  geom_smooth(aes(y = Leo_Cumulative, color = "Leopard"), 
              method = "loess", se = FALSE, linewidth = 1) +
  labs(y = "Cumulative Detection Probability\n", 
       x = "\n Survey Day", title = "Persian Leopard Cumulative Detection Probability p(.)") +
  scale_color_manual(name = "Legend",
                     values = c("Persian Leopard" = "black")) +
  scale_x_continuous(breaks = seq(0, max(plot_data$Day), by = 30)) +  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA) 
  )


leo_cum_det_df <- as.data.frame(leo_cumulative_p)
leo_cum_det_df[,2] <- 1:NROW(leo_cum_det_df)
colnames(leo_cum_det_df) <- c("Cumulative Detection Probability", "Survey Day")

write.csv(leo_cum_det_df, "Leopard_CumulativeDetProb.csv")
