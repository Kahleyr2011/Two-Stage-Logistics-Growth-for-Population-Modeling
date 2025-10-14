#########################################
# Two-Stage Discrete Logistic Growth Model in R
# Stage 1: Fit r and K without death term
# Stage 2: Project forward with user-specified death rate
#########################################

# ---- 1. Load libraries ----
# Install if needed:
# install.packages(c("tidyverse"))
library(tidyverse)

# ---- 2. User Inputs ----
inputs <- list(
  P0 = 146,         # Initial population
  r_guess = 0.5,   # Initial guess for growth rate
  K_guess = 4000,  # Initial guess for carrying capacity
  r_min = 0,       # Lower bound for r
  r_max = 1,       # Upper bound for r
  K_min = 600,     # Lower bound for K
  K_max = 5000,    # Upper bound for K
  extend_years = 25,  # Forecast beyond dataset
  death_rate = 0   # STUDENT INPUT: Death rate for projections (try 0, 5, 10, 20, etc.)
)

# ---- 3. Import Data ----
# CSV format: Year,Population
# OPTION 1: Use file.choose() - opens file dialog for user to select data
# data <- read.csv(file.choose())

# OPTION 2: Place CSV file in same folder as this R script and use filename only
data <- read.csv("/Users/reginajordan/Library/Mobile Documents/com~apple~CloudDocs/Courses/Math Modeling/Whooping_Crane_Population_Data/whooping_crane_population_for_R.csv")

# OPTION 3: Use sample data if no file available (for testing)
# data <- data.frame(
#   Year = 2000:2020,
#   P_obs = c(18, 25, 32, 42, 55, 68, 82, 98, 115, 133, 152, 172, 193, 215, 238, 262, 287, 313, 340, 368, 397)
# )

colnames(data) <- c("Year", "P_obs")

# ---- 4. STAGE 1: Standard Logistic Model (No Death) ----
logistic_model <- function(params, years, P0) {
  r <- params[1]; K <- params[2]
  P <- numeric(length(years))
  P[1] <- P0
  for (i in 2:length(years)) {
    P[i] <- P[i-1] + r * P[i-1] * (1 - P[i-1] / K)
  }
  return(P)
}

# ---- 5. SSE Function for Stage 1 ----
sse_fn <- function(par, years, P_obs, P0) {
  P_pred <- logistic_model(par, years, P0)
  sum((P_obs - P_pred)^2, na.rm = TRUE)
}

# ---- 6. STAGE 1: Fit r and K (No Death Term) ----
inputs$K_guess <- max(min(inputs$K_guess, inputs$K_max), inputs$K_min)
inputs$r_guess <- max(min(inputs$r_guess, inputs$r_max), inputs$r_min)

fit_opt <- optim(
  par   = c(inputs$r_guess, inputs$K_guess),
  fn    = function(par) sse_fn(par, data$Year, data$P_obs, inputs$P0),
  method= "L-BFGS-B",
  lower = c(inputs$r_min, inputs$K_min),
  upper = c(inputs$r_max, inputs$K_max),
  control = list(maxit = 2000)
)

if (fit_opt$convergence != 0) {
  stop(paste("Optimizer did not converge. Code =", fit_opt$convergence))
}

r_est <- fit_opt$par[1]  # FIXED for all scenarios
K_est <- fit_opt$par[2]  # FIXED for all scenarios

cat("STAGE 1 - Parameter estimates (fitted to observed data):\n")
cat("r (growth rate) =", round(r_est, 4), "\n")
cat("K (carrying capacity) =", round(K_est, 2), "\n")
cat("These values are FIXED for all death scenarios.\n\n")

# ---- 7. STAGE 2: Projection Model with Death ----
logistic_model_with_death <- function(r, K, d, years, P0, transition_year = NULL) {
  P <- numeric(length(years))
  P[1] <- P0
  
  for (i in 2:length(years)) {
    # Apply death term only after transition year (if specified)
    if (!is.null(transition_year) && years[i] <= transition_year) {
      # Before transition: no death term
      P[i] <- P[i-1] + r * P[i-1] * (1 - P[i-1] / K)
    } else {
      # After transition: include death term
      P[i] <- P[i-1] + r * P[i-1] * (1 - P[i-1] / K) - d
    }
    # Prevent negative populations
    P[i] <- max(P[i], 0)
  }
  return(P)
}

# ---- 8. Generate Projections ----
last_year <- max(data$Year)
future_years <- (last_year+1):(last_year + inputs$extend_years)
all_years <- c(data$Year, future_years)

# Model with death term applied only to future projections
extended_P <- logistic_model_with_death(
  r = r_est, 
  K = K_est, 
  d = inputs$death_rate,
  years = all_years,
  P0 = inputs$P0,
  transition_year = last_year  # Death term starts after observed data
)

cat("STAGE 2 - Projection with death rate =", inputs$death_rate, "\n")
cat("Model equation for projections: P(n+1) = P(n) + r*P(n)*(1 - P(n)/K) - d\n\n")

# ---- 9. Yearly Results ----
yearly_results <- tibble(
  Year = as.character(all_years),
  P_obs = c(data$P_obs, rep(NA, length(future_years))),
  P_model = extended_P,
  Death_Applied = c(rep("No", length(data$Year)), rep("Yes", length(future_years))),
  Squared_Error = c((data$P_obs - extended_P[1:length(data$Year)])^2, 
                    rep(NA, length(future_years))),
  r_fixed = r_est,
  K_fixed = K_est,
  death_rate = inputs$death_rate
)

# ---- 10. Parameter Summary Row ----
param_row <- tibble(
  Year = "PARAMETERS",
  P_obs = NA,
  P_model = NA,
  Death_Applied = NA,
  Squared_Error = NA,
  r_fixed = r_est,
  K_fixed = K_est,
  death_rate = inputs$death_rate
)

# ---- 11. Final Results ----
final_results <- bind_rows(param_row, yearly_results)

# View in RStudio (spreadsheet-like)
View(final_results)

# Export to CSV (saves to current working directory)
filename <- paste0("logistic_model_death_", inputs$death_rate, ".csv")
write.csv(final_results, filename, row.names = FALSE)
cat("Results saved as:", filename, "\n")
cat("File location:", file.path(getwd(), filename), "\n")
cat("Working directory:", getwd(), "\n")

# Optional: Export to Excel (uncomment if needed)
# library(openxlsx)
# write.xlsx(final_results, paste0("logistic_model_death_", inputs$death_rate, ".xlsx"))

# Optional: Export just the yearly data (without parameter row)
# write.csv(yearly_results, paste0("yearly_data_death_", inputs$death_rate, ".csv"), row.names = FALSE)

# ---- 12. Plot Results ----
# Create a visual separator at the transition point
transition_line_data <- data.frame(x = as.numeric(as.character(last_year)) + 0.5)

ggplot(yearly_results, aes(x = as.numeric(Year))) +
  # Observed data
  geom_point(aes(y = P_obs), color = "blue", size = 3, na.rm = TRUE) +
  geom_line(aes(y = P_obs), color = "blue", linetype = "dashed",
            linewidth = 1, na.rm = TRUE) +
  # Model predictions
  geom_line(aes(y = P_model), color = "red", linewidth = 1) +
  # Transition line
  geom_vline(data = transition_line_data, aes(xintercept = x), 
             linetype = "dotted", color = "gray50", linewidth = 1) +
  # Labels and formatting
  labs(
    title = paste0("Two-Stage Logistic Model (Death Rate = ", inputs$death_rate, ")"),
    subtitle = paste0("Stage 1: Fit r=", round(r_est,3), ", K=", round(K_est,0), 
                      " | Stage 2: Fixed r,K + Death Term"),
    x = "Year",
    y = "Population",
    caption = "Blue = observed | Red = model | Dotted line = transition to death term"
  ) +
  annotate("text", x = min(as.numeric(yearly_results$Year)), 
           y = max(yearly_results$P_model, na.rm = TRUE) * 0.9,
           label = "No Death\n(Fitting Phase)", hjust = 0, color = "darkgreen") +
  annotate("text", x = max(as.numeric(yearly_results$Year)), 
           y = max(yearly_results$P_model, na.rm = TRUE) * 0.9,
           label = paste0("Death = ", inputs$death_rate, "\n(Projection Phase)"), 
           hjust = 1, color = "darkred") +
  theme_minimal()

# ---- 13. Instructions for Students ----
cat("\n" , rep("=", 80), "\n")
cat("SETUP INSTRUCTIONS FOR NEW USERS:\n")
cat(rep("=", 80), "\n")
cat("1. DATA FILE: Choose one option for importing data:\n")
cat("   Option A: Uncomment line with file.choose() - opens file dialog\n")
cat("   Option B: Place your CSV file in same folder as this R script\n")
cat("   Option C: Use the sample data (uncomment those lines)\n")
cat("2. CSV FORMAT: Your data file should have columns: Year, Population\n")
cat("3. OUTPUT: Results will be saved in the same folder as this R script\n")
cat(rep("=", 80), "\n\n")

cat(rep("=", 80), "\n")
cat("INSTRUCTIONS FOR STUDENTS:\n")
cat(rep("=", 80), "\n")
cat("1. The model first fits r and K to observed data (no death term)\n")
cat("2. These values are FIXED: r =", round(r_est,4), ", K =", round(K_est,2), "\n")
cat("3. To experiment with different death scenarios:\n")
cat("   - Change 'death_rate = X' in the inputs section (line 21)\n")
cat("   - Try values like: 0, 5, 10, 15, 20, 30, etc.\n")
cat("   - Re-run the code to see how death affects projections\n")
cat("4. The death term only applies to future projections, not the fitting period\n")
cat("5. Compare results with different death rates to understand the impact\n")
cat("6. Results automatically save as CSV files in your current folder\n")
cat(rep("=", 80), "\n")
