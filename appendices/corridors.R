# This code assumes your caribou GPS data is in a data frame called caribou_data
# with columns for animal ID (id), timestamp (timestamp), and coordinates (x, y)
# See Migration mapping section in caribou_methods.pdf
# Please note that this script is a template. You will need to customize the file paths, column names, environmental variables, and candidate model structures to fit your specific research questions and dataset for caribou. The source emphasizes the biological relevance of variable selection, so ensure your chosen covariates are appropriate for caribou ecology.

library(adehabitatLT)
library(adehabitatHR)
library(amt)
library(sf)
library(lme4)
library(terra)
library(ggplot2)

# --- 1. Data Preparation & NSD Calculation ---
# Assuming 'caribou_data' is a dataframe with id, timestamp, x, y
# Convert to a trajectory object for analysis
caribou_traj <- as.ltraj(xy = caribou_data[,c("x","y")], date = caribou_data$timestamp, id = caribou_data$id)

# Calculate Net Squared Displacement (NSD)
# The study used March 15 as the start date
# Here we calculate NSD from the first point of each individual's track
nsd_data <- ld(caribou_traj)

# Example plot for one individual (e.g., ID "C01")
c01_nsd <- nsd_data[nsd_data$id == "C01", ]
ggplot(c01_nsd, aes(x = date, y = R2n)) +
  geom_line() +
  labs(title = "NSD Plot for Caribou C01", x = "Date", y = "Net Squared Displacement (m^2)") +
  theme_minimal()
# Visual inspection of this plot helps identify migration timing


# --- 2. Define Seasons and Filter Data ---
# Based on your NSD plots, define start/end dates for each individual/season
# This part is manual and requires interpretation
# Example for one individual:
migrant_C01_data <- caribou_data[caribou_data$id == "C01",]

# Define seasons based on dates you identified from NSD plots
migrant_C01_data$season <- NA
migrant_C01_data$season[migrant_C01_data$timestamp >= as.POSIXct("2023-11-01") & migrant_C01_data$timestamp < as.POSIXct("2024-05-15")] <- "winter"
migrant_C01_data$season[migrant_C01_data$timestamp >= as.POSIXct("2024-06-01") & migrant_C01_data$timestamp < as.POSIXct("2024-09-15")] <- "summer"
# ... and so on for migration periods


# --- 3. Delineate Seasonal Ranges with BBMM ---
# Filter data for one season for one individual
winter_data_C01 <- migrant_C01_data[migrant_C01_data$season == "winter" & !is.na(migrant_C01_data$season),]

# Convert to a trajectory object for BBMM
winter_traj_C01 <- as.ltraj(xy = winter_data_C01[,c("x","y")], date = winter_data_C01$timestamp, id = winter_data_C01$id)

# Estimate BBMM parameters
# The study uses liker() for sigma1 and a fixed value for sigma2
# Note: The liker function can be computationally intensive
bb_params <- liker(winter_traj_C01, sig2 = 30, rangesig1 = c(1, 10))
# Let's assume liker outputted an optimal sigma1 of 5 for this example
sigma1_est <- 5

# Calculate the BBMM Utilization Distribution (UD)
bbmm_winter_C01 <- kernelbb(winter_traj_C01, sig1 = sigma1_est, sig2 = 30, grid = 50)

# Get the 95% UD contour to represent the winter range
winter_range_C01 <- getverticeshr(bbmm_winter_C01, percent = 95)

# Calculate area of the winter range
winter_range_C01_sf <- st_as_sf(winter_range_C01)
winter_area_ha <- st_area(winter_range_C01_sf) / 10000 # Convert m^2 to hectares
print(paste("Winter Range Area for C01:", round(winter_area_ha, 2), "ha"))


# --- 4. Identify Migration Routes and Stopover Sites ---
# Follow a similar BBMM process using data from the migration period
# For stopovers, identify the top 10% of the migration UD
# migration_ud <- kernelbb(...)
stopover_contour <- getverticeshr(migration_ud, percent = 10)
# Further analysis would be needed to check duration within these sites (>12h)


# --- 5. Analyze Altitudinal Variation ---
# Load your Digital Elevation Model (DEM)
dem <- raster("path/to/your/dem.tif")

# Convert your data to a spatial object
caribou_sf <- st_as_sf(caribou_data, coords = c("x", "y"), crs = st_crs(dem))

# Extract elevation for each point
caribou_data$elevation <- raster::extract(dem, caribou_sf)

# Calculate and plot 14-day moving average of elevation for one animal
library(zoo)
c01_elevation_data <- caribou_data[caribou_data$id == "C01", ]
c01_elevation_data <- c01_elevation_data[order(c01_elevation_data$timestamp),] # ensure chronological order

# Aggregate to daily mean elevation
daily_mean_elev <- aggregate(elevation ~ as.Date(timestamp), data = c01_elevation_data, FUN = mean)
names(daily_mean_elev) <- c("date", "mean_elevation")

# Calculate 14-day moving average
daily_mean_elev$moving_avg_elev <- rollmean(daily_mean_elev$mean_elevation, k = 14, fill = NA, align = "right")

# Plot the results
ggplot(daily_mean_elev, aes(x = date, y = moving_avg_elev)) +
  geom_line() +
  labs(title = "14-Day Moving Average Elevation for Caribou C01", x = "Date", y = "Mean Elevation (m)") +
  theme_minimal()
# This plot helps identify traditional vs. abbreviated altitudinal migration


#------------------------------------------------------------------
# Step 1: Prepare Caribou GPS Data and Environmental Rasters
#------------------------------------------------------------------

# Load your caribou GPS data (assuming a data frame named 'caribou_data')
# It should have columns like: animal_id, x, y, timestamp
# caribou_data <- read.csv("your_caribou_data.csv")
# caribou_data$timestamp <- as.POSIXct(caribou_data$timestamp) # Ensure timestamp is correct format

# Load your environmental raster layers (e.g., elevation, slope, ruggedness, roads, NDVI)
# Ensure all rasters have the same coordinate reference system (CRS) and resolution
elevation_raster <- raster("path/to/elevation.tif")
slope_raster <- raster("path/to/slope.tif")
ruggedness_raster <- raster("path/to/ruggedness.tif")
# Add other rasters as needed (e.g., distance_to_roads, NDVI)

# Create a raster stack of your covariates
covariate_stack <- stack(elevation_raster, slope_raster, ruggedness_raster)
names(covariate_stack) <- c("elevation", "slope", "ruggedness")

#------------------------------------------------------------------
# Step 2: Home-Range Scale Analysis (RSF)
#------------------------------------------------------------------

# Create a track object using amt
trk <- make_track(caribou_data, .x = x, .y = y, .t = timestamp, id = animal_id, crs = crs(covariate_stack))

# Calculate home ranges (95% KDE) for each caribou
# Note: This is one way to do it; the source used KDEs.
sp_points <- SpatialPointsDataFrame(coords = trk[, c("x_", "y_")], data = trk, proj4string = crs(covariate_stack))
kdes <- kernelUD(sp_points[, "id"], h = "href")
home_ranges <- getverticeshr(kdes, percent = 95)

# Generate used and available points
# 10 available points per used point
rsf_data <- data.frame()
for (id in unique(trk$id)) {
  used_points <- filter(trk, id == !!id)

  # Get the specific home range polygon for this individual
  individual_hr <- home_ranges[home_ranges@data$id == id, ]

  # Generate random available points within the home range
  available_points <- spsample(individual_hr, n = nrow(used_points) * 10, type = "random")

  # Combine used (case=1) and available (case=0) points
  used_df <- data.frame(x = used_points$x_, y = used_points$y_, case = 1, id = id)
  available_df <- data.frame(x = available_points@coords[,1], y = available_points@coords[,2], case = 0, id = id)

  rsf_data <- rbind(rsf_data, used_df, available_df)
}

# Extract covariate values for all points
rsf_data <- cbind(rsf_data, raster::extract(covariate_stack, rsf_data[, c("x", "y")]))

# Center and scale covariates
rsf_data_scaled <- rsf_data %>%
  mutate(across(c(elevation, slope, ruggedness), scale))

# Add squared terms if hypothesized (e.g., for elevation)
rsf_data_scaled$elevation2 <- rsf_data_scaled$elevation^2

# Fit the RSF model using mixed-effects logistic regression
# Individual ID is the random effect
rsf_model <- glmer(case ~ elevation + elevation2 + slope + ruggedness + (1 | id),
                   data = rsf_data_scaled, family = binomial(link = "logit"))

summary(rsf_model)

# Note: For model selection, you would build multiple candidate models and compare them using AICc.

#------------------------------------------------------------------
# Step 3: Fine-Scale Analysis (iSSA)
#------------------------------------------------------------------

# Resample track to a regular time interval (e.g., 2 hours, as in the source)
# This is crucial for creating steps of consistent duration.
trk_resampled <- trk %>% nest(data = -"id") %>%
  mutate(steps = map(data, ~ track_resample(., rate = hours(2), tolerance = minutes(15)))) %>%
  select(id, steps) %>% unnest(cols = steps)

# Create steps and generate 10 random steps for each used step
steps <- trk_resampled %>% steps_by_burst() %>%
  random_steps(n_control = 10)

# Extract covariates at the end of each step
steps_covariates <- steps %>%
  extract_covariates(covariate_stack) %>%
  mutate(log_sl_ = log(sl_), cos_ta_ = cos(ta_)) # Include movement parameters

# Center and scale covariates
steps_covariates_scaled <- steps_covariates %>%
  mutate(across(c(elevation, slope, ruggedness), scale))

# Add squared terms if hypothesized
steps_covariates_scaled$elevation2 <- steps_covariates_scaled$elevation^2

# Fit the iSSA model using conditional logistic regression (clogit) with a random effect
# 'case_' indicates used (TRUE) vs. random (FALSE) steps
# 'step_id_' is the stratum
# Use lme4::glmer for mixed-effects iSSA
issa_model <- glmer(case_ ~ elevation + elevation2 + slope + ruggedness +
                    sl_ + log_sl_ + cos_ta_ + (1 | id),
                    data = steps_covariates_scaled, family = binomial,
                    # This is a trick to approximate clogit with random effects
                    # by adding a fixed effect for each stratum
                    # It can be computationally intensive
                    # control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
                    )

# The source used lme4, but amt's `fit_clogit` is standard for iSSA.
# A more typical (non-mixed-effect) iSSA would be:
# issa_model_clogit <- steps_covariates_scaled %>%
#   fit_clogit(case_ ~ elevation + elevation2 + slope + ruggedness +
#              sl_ + log_sl_ + cos_ta_ + strata(step_id_))
# summary(issa_model_clogit)

summary(issa_model)

#------------------------------------------------------------------
# Step 4: Model Interpretation and Prediction
#------------------------------------------------------------------

# Interpret coefficients: Positive values indicate selection, negative values indicate avoidance.
# Check if 95% confidence intervals overlap zero to determine significance.

# Create predictive habitat suitability maps from the RSF model
# Use the 'predict' function with the raster stack
suitability_map <- predict(covariate_stack, rsf_model, type = "response", re.form = NA)
plot(suitability_map)
