# Load and structure weather data

weather_combined <- readRDS("weather_koeln-bonn.rds")

process_weather_data <- function(weather_combined) {
  
  # Calculate Average Temperature (Tavg)
  weather_combined$Tavg <- ( weather_combined$Tmax + weather_combined$Tmin ) / 2
  
  # Calculate Consecutive Dry Days
  
  # Initialize the new column with zeros 
  weather_combined$day_consec_dry <- 0 
  
  dry_days <- weather_combined$Prec == 0
  rle_dry <- rle(dry_days)
  consec_dry <- unlist(lapply(seq_along(rle_dry$lengths), function(i) {
    if (rle_dry$values[i]) {
      seq_len(rle_dry$lengths[i])
    } else {
      rep(0, rle_dry$lengths[i])
    }
  }))
  weather_combined$day_consec_dry <- consec_dry
  
  
  # Calculate Consecutive Wet Days
  
  # Initialize the new column with zeros 
  
  # Logical vector: TRUE when wet (Rain), FALSE otherwise
  wet_days <- weather_combined$Prec > 0
  
  # run length encoding on wet_days
  rle_wet <- rle(wet_days)
  
  consec_wet <- unlist(lapply(seq_along(rle_wet$lengths), function(i) {
    if (rle_wet$values[i]) {
      seq_len(rle_wet$lengths[i])
    } else {
      rep(0, rle_wet$lengths[i])
    }
  }))
  
  weather_combined$day_consec_wet <- consec_wet

  weather_list <<- split(x = weather_combined, f = weather_combined$id_seaon)
  
  scenarios <<- c("historical", "ssp126", "ssp245", "ssp370", "ssp585")

  return(weather_list)
}

# Apply function to receive structured weather file
weather_list <- process_weather_data(weather_combined)
 