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
