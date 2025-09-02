# library(ggtext)
# library(png)
# library(grid)
# library(patchwork)
# library(tidyverse)
# library(decisionSupport)
# library(dplyr)
# library(ggh4x)
#@Paul: local pacakage??
library(introdataviz)
#library(lubridate) # Ensure lubridate is loaded for yday()

flist <- list.files('future_weather/', full.names = TRUE)

future_weather <- purrr::map(flist, function(file_path) {
  
  # Extract filename
  filename <- basename(file_path) # just the file name, no directory
  
  # Split by '.' and extract parts
  name_parts <- str_split(filename, '\\.')[[1]]
  
  # Check if name_parts is long enough
  if (length(name_parts) < 5) {
    stop(paste("Filename", filename, "does not have enough parts."))
  }
  
  # Read and process the file
  read.csv(file_path) %>%
    select(-X) %>%
    mutate(
      ssp = name_parts[3],
      gcm = name_parts[4],
      scenario_year = name_parts[5],
      yday = lubridate::yday(DATE),
      season = ifelse(yday >= 0, Year + 1, Year)
    )
})

future_weather <- do.call('rbind', future_weather)
future_weather$id = paste(future_weather$ssp, future_weather$gcm, future_weather$scenario_year, sep = '--')


#weed out incomplete seasons
drop_list <- future_weather %>%
  group_by(id, season) %>%
  summarise(n = n()) %>%
  filter(n < 365)

future_weather <- future_weather %>%
  filter(!(paste(id, season) %in% paste(drop_list$id, drop_list$season)))

hist_weather <- read.csv('weather_2020_koeln-bonn.csv') %>%
  mutate(scenario_year = 2020,
         ssp = 'historical',
         gcm = 'historical',
         yday = lubridate::yday(DATE),
         season = ifelse(yday >= 0,
                         yes = Year  +1,
                         no = Year),
         id = paste(ssp, gcm, scenario_year, sep = '--'))
drop_list <- hist_weather %>%
  group_by(id, season) %>%
  summarise(n = n()) %>%
  filter(n < 365)
hist_weather <- hist_weather %>%
  filter(!(paste(id, season) %in% paste(drop_list$id, drop_list$season)))

#combine hist and future weather
weather_combined <- future_weather %>%
  rbind(hist_weather)

#add unique id to each seaoson
weather_combined$id_seaon <- paste(weather_combined$id, weather_combined$season, sep = '--')

#write.csv(weather_combined, 'weather_onion_koeln-bonn.csv')

weather_combined <- read.csv("weather_onion_koeln-bonn.csv")

weather_combined$Tavg <- ( weather_combined$Tmax + weather_combined$Tmin ) / 2

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


weather_combined$day_consec_wet <- 0

# Logical vector: TRUE wenn wet (regen), FALSE sonst
wet_days <- weather_combined$Prec > 0

# run length encoding auf wet_days
rle_wet <- rle(wet_days)


consec_wet <- inverse.rle(list(
  lengths = rle_wet$lengths,
  values = ifelse(rle_wet$values, seq_len(max(rle_wet$lengths)), 0)
))


consec_wet <- unlist(lapply(seq_along(rle_wet$lengths), function(i) {
  if (rle_wet$values[i]) {
    seq_len(rle_wet$lengths[i])
  } else {
    rep(0, rle_wet$lengths[i])
  }
}))

weather_combined$day_consec_wet <- consec_wet



#split by each season
weather_list <- split(x = weather_combined, f = weather_combined$id_seaon)

# Define the scenario prefixes and corresponding variable names
scenarios <- c("historical", "ssp126", "ssp245", "ssp370", "ssp585")



# Weather Lists justs needs to be calculated once, model picks different random season each run.



##### Input Table Onion model


input_variables <- read.csv("input_table_onion.csv", header = TRUE, sep = ";")

make_variables <- function(est,n=1)
{ x<-random(rho=est, n=n)
for(i in colnames(x)) assign(i,as.numeric(x[1,i]),envir=.GlobalEnv)
}

make_variables(as.estimate(input_variables))

onion_climate_impact <- function(){ # Start of onion_climate_impact function
  
  # 1st step randomly selecet a season
  # Loop through each scenario and assign to List
  
  weather_scenario_list <- list()
  
  for (scenario in scenarios) {
    # Find all matching names
    matching_names <- names(weather_list)[grepl(paste0("^", scenario, "--"), names(weather_list))]
    
    # Randomly select one
    selected_name <- sample(matching_names, 1)
    
    # Store in the list with name like "weather_historical"
    weather_scenario_list[[paste0("weather_", scenario)]] <- weather_list[[selected_name]]
    
    # Optional: print the selected name
    cat("Selected for", scenario, ":", selected_name, "\n")
  }
  
  ## just for name storage
  
  weather_scenario_names <- weather_scenario_list
  
  
  ### Several steps are needed in order to calculate PAR, these are explained in the following
  
  # Calculate Extraterrestrial Radiation (Ra)
  #
  # This function calculates Ra (in MJ/m²/day) based on latitude and day of the year,
  # following FAO-56 methodology.
  # This is to calcualte extra terrestrial radiation, which is required in order to calculate PAR
  
  # Functions
  calc_Ra <- function(latitude_koeln_bonn_c, yday) {
    Gsc <- 0.0820
    phi <- latitude_koeln_bonn_c * pi / 180
    dr <- 1 + 0.033 * cos(2 * pi * yday / 365)
    delta <- 0.409 * sin((2 * pi * yday / 365) - 1.39)
    omega_s <- acos(-tan(phi) * tan(delta))
    Ra <- (24 * 60 / pi) * Gsc * dr *
      (omega_s * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(omega_s))
    return(Ra)
  }
  
  calc_Rs <- function(Ra, deltaT) {
    return(Krs_c * Ra * deltaT)
  }
  
  
  # Process each scenario
  weather_scenario_list <- lapply(names(weather_scenario_list), function(scenario) {
    df <- weather_scenario_list[[scenario]]
    
    # --- Calculate PAR ---
    Ra <- calc_Ra(latitude_koeln_bonn_c, df$yday)
    deltaT <- df$Tmax - df$Tmin
    Rs <- calc_Rs(Ra, deltaT)
    df$PAR <- Rs * 0.45
    
    # --- Calculate daily GDD ---
    df$GDD_daily <- pmax(0, ((df$Tmax + df$Tmin) / 2) - base_temp_p)
    
    # Sort and filter days after planting
    df_sorted <- df[order(df$yday), ]
    df_after_planting <- df_sorted[df_sorted$yday >= planting_yday_p, ]
    
    
    # PHASE 1: Emergence
    
    gdd_emergence_cumsum <- cumsum(df_after_planting$GDD_daily)
    emergence_index <- which(gdd_emergence_cumsum >= GDD_field_emergence_required_p)[1]
    
    if (!is.na(emergence_index)) {
      emergence_yday <- df_after_planting$yday[emergence_index]
      df$in_emergence_phase <- df$yday %in% df_after_planting$yday[1:emergence_index]
    } else {
      emergence_yday <- NA
      df$in_emergence_phase <- FALSE
    }
    
    
    # PHASE 2: Vegetative
    
    if (!is.na(emergence_yday)) {
      df_after_emergence <- df[df$yday >= emergence_yday, ]
      gdd_veg_cumsum <- cumsum(df_after_emergence$GDD_daily)
      bulbing_index <- which(gdd_veg_cumsum >= GDD_vegetative_required_p)[1]
      
      if (!is.na(bulbing_index)) {
        bulbing_yday <- df_after_emergence$yday[bulbing_index]
        df$in_vegetative_phase <- df$yday >= emergence_yday & df$yday < bulbing_yday
      } else {
        bulbing_yday <- NA
        df$in_vegetative_phase <- FALSE
      }
    } else {
      bulbing_yday <- NA
      df$in_vegetative_phase <- FALSE
    }
    
    
    # PHASE 3: Bulbing
    
    if (!is.na(bulbing_yday)) {
      df_after_bulbing <- df[df$yday >= bulbing_yday, ]
      gdd_bulbing_cumsum <- cumsum(df_after_bulbing$GDD_daily)
      maturation_index <- which(gdd_bulbing_cumsum >= GDD_bulbing_required_p)[1]
      
      if (!is.na(maturation_index)) {
        maturation_yday <- df_after_bulbing$yday[maturation_index]
        df$in_bulbing_phase <- df$yday >= bulbing_yday & df$yday < maturation_yday
      } else {
        maturation_yday <- NA
        df$in_bulbing_phase <- FALSE
      }
    } else {
      maturation_yday <- NA
      df$in_bulbing_phase <- FALSE
    }
    
    
    # PHASE 4: Maturation + Harvest Window
    
    if (!is.na(maturation_yday)) {
      df_after_maturation <- df[df$yday >= maturation_yday, ]
      gdd_mat_cumsum <- cumsum(df_after_maturation$GDD_daily)
      harvest_index <- which(gdd_mat_cumsum >= GDD_maturation_required_p)[1]
      
      if (!is.na(harvest_index)) {
        harvest_yday <- df_after_maturation$yday[harvest_index]
        df$in_maturation_phase <- df$yday >= maturation_yday & df$yday < harvest_yday
        df$in_harvest_window <- df$yday >= harvest_yday
      } else {
        harvest_yday <- NA
        df$in_maturation_phase <- FALSE
        df$in_harvest_window <- FALSE
      }
    } else {
      harvest_yday <- NA
      df$in_maturation_phase <- FALSE
      df$in_harvest_window <- FALSE
    }
    
    
    # Final annotation
    
    df$emergence_yday <- emergence_yday
    df$bulbing_yday <- bulbing_yday
    df$maturation_yday <- maturation_yday
    df$harvest_yday <- harvest_yday
    
    return(df)
  })
  
  # Reapply names (optional but safe)
  names(weather_scenario_list) <- names(weather_scenario_names)
  
  
  
  ###### Functions for yield reducing factors #########
  
  ## Abiotic stress factors
  
  # Seedbed preparation
  
  get_seedbed_stress <- function(Prec, Tavg, day_consec_wet) {
    # Trockenstress: zu wenig Regen + heißen Tagen
    if (day_consec_wet < day_consec_wet_seedbed_upper_limit_risk_high_p & Prec < Prec_seedbed_upper_limit_risk_high_p & Tavg > Tavg_seedbed_lower_limit_risk_high_p) {
      return(chance_event(
        chance = risk_seedbed_high_p,
        value_if = yield_reduction_seedbed_p,
        value_if_not = 0))
    } else if (day_consec_wet < day_consec_wet_seedbed_upper_limit_risk_medium_p & Prec < Prec_seedbed_upper_limit_risk_medium) {
      return(chance_event(
        chance = risk_seedbed_medium_p,
        value_if = yield_reduction_seedbed_p,
        value_if_not = 0))
    } else {
      return(0)
    }
  }
 
  
  #Drought Stress
  
  
  get_drought_stress <- function(Prec, Tavg, day_consec_wet) {
    if (Prec < Prec_drought_upper_limit_risk_high_p & Tavg > Tavg_drought_lower_limit_risk_high_p) {
      return(chance_event(
        chance =  risk_drought_high_p,
        value_if = yield_reduction_drought_p,
        value_if_not = 0))
    } else if (Prec < Prec_drought_upper_limit_risk_medium_p & day_consec_wet < day_consec_wet_drought_upper_limit_risk_medium_p) {
      return(chance_event(
        chance = risk_drought_medium_p,
        value_if = yield_reduction_drought_p,
        value_if_not = 0))
    } else {
      return(0)
    }
  }
  
  
  # Extreme Rainfall
  
  get_extreme_rain_stress <- function(Prec) {
    if (Prec > Prec_extreme_rain_lower_limit_risk_high_p) {
      return(chance_event(
        chance = risk_extreme_rain_high_p,
        value_if = yield_reduction_extreme_rain_p,
        value_if_not = 0))
    } else if (Prec > Prec_extreme_rain_lower_limit_risk_medium_p) {
      return(chance_event(
        chance = risk_extreme_rain_medium_p,
        value_if = yield_reduction_extreme_rain_p,
        value_if_not = 0))
    } else {
      return(0)
    }
  }
  
  ## Harvest Risk
  
  
  
  ### Biotic Stress Factors
  
  ## Weed Pressure
  # Low PAR , little bit of rain helps emergence, this needs to be updated any maybe different conditions for different phases
  
  get_weed_pressure_stress <- function(PAR, Prec, day_consec_wet) {
    if (PAR < PAR_weed_upper_limit_risk_high_p) {
      if (Prec >= Prec_weed_lower_limit_risk_high_p & day_consec_wet >= day_consec_wet_weed_lower_limit_risk_high_p) {
        return(chance_event(
          chance = risk_weed_high_p,
          value_if = yield_reduction_weed_p,
          value_if_not = 0,
        ))
      } else if (Prec >= Prec_weed_lower_limit_risk_medium_p) {
        return(chance_event(
          chance =  risk_weed_medium_p,
          value_if = yield_reduction_weed_p,
          value_if_not = 0,
        ))
      } else {
        return(0)
      }
    } else { # Added this else block
      return(0)
    }
  }
  
  ## Fungal Pathogens
  
  # Risk for Botrytis
  
  get_botrytis_stress <- function(Tavg, Prec, day_consec_wet) {
    if (Prec >= Prec_botrytis_lower_limit_risk_high_p & day_consec_wet >= day_consec_wet_botrytis_lower_limit_risk_high_p) {
      return(chance_event(
        chance = risk_neck_rot_high_p,
        value_if = yield_reduction_neck_rot_p,
        value_if_not = 0,
      ))
      
    } else if (Prec >= Prec_botrytis_lower_limit_risk_medium_p) {
      return(chance_event(
        chance = risk_neck_rot_medium_p,
        value_if = yield_reduction_neck_rot_p,
        value_if_not = 0,
      ))
    } else {
      return(0)
    }
  }
  
  
  
  # Fusarium (Fungi, soil, pathogen, wet & warm increases risk of infection)
  get_fusarium_stress <- function(Tavg, Prec, day_consec_wet) {
    if (Tavg >= Tavg_fusarium_lower_limit_all_risks_p & Tavg <= Tavg_fusarium_upper_limit_all_risks_p) {
      if (Prec >= Prec_fusarium_lower_limit_risk_high_p & day_consec_wet >= day_consec_wet_fusarium_lower_limit_risk_high_p) {
        return(chance_event(
          chance =  risk_fusarium_high_p,
          value_if = yield_reduction_fusarium_p,
          value_if_not = 0,
        ))
      } else if (Prec >= Prec_fusarium_lower_limit_risk_medium_p) {
        return(chance_event(
          chance = risk_fusarium_medium_p,
          value_if = yield_reduction_fusarium_p,
          value_if_not = 0,
        ))
      } else {
        return(0)
      }
    } else {
      return(0)
    }
  }
  
  # Downy Mildew (Peronospora destructor) – cold nights , wet periods
  
  get_downy_mildew_stress <- function(Tavg, Prec, day_consec_wet) {
    if (Tavg >= Tavg_downy_mildew_lower_limit_all_risks_p & Tavg <= Tavg_downy_mildew_upper_limit_all_risks_p) {
      if (Prec >= Prec_downy_mildew_lower_limit_risk_high_p & day_consec_wet >= day_consec_wet_downy_mildew_lower_limit_risk_high_p) {
        return(chance_event(
          chance = risk_downy_mildew_high_p,
          value_if = yield_reduction_downy_mildew_p,
          value_if_not = 0,
        ))
      } else if (Prec >= Prec_downy_mildew_lower_limit_risk_medium_p) {
        return(chance_event(
          chance = risk_downy_mildew_medium_p,
          value_if = yield_reduction_downy_mildew_p,
          value_if_not = 0,
        ))
      } else {
        return(0)
      }
    } else {
      return(0)
    }
  }
  
  ## Animal Pressure
  
  # Thripse – heiß & trocken
  get_thrips_stress <- function(Tavg, Prec) {
    if (Tavg > Tavg_thrips_lower_limit_risk_high_p & Prec < Prec_thrips_upper_limit_risk_high_p) {
      return(chance_event(
        chance = risk_thrips_high_p,
        value_if = yield_reduction_thrips_p,
        value_if_not = 0,
      ))
    } else {
      return(0)
    }
  }
  
  # Zwiebelfliege – warm, leicht feucht fördert Larvenaktivität
  get_onion_fly_stress <- function(Tavg, Prec) {
    if (Tavg >= Tavg_onion_fly_lower_limit_all_risks_p) {
      if (Prec >= Prec_onion_fly_lower_limit_risk_high_p) {
        return(chance_event(
          chance = risk_onion_fly_high_p,
          value_if = yield_reduction_onion_fly_p,
          value_if_not = 0
        ))
      } else { # This else belongs to 'if (Prec >= 1)'
        return(chance_event(
          chance = risk_onion_fly_medium_p,
          value_if = yield_reduction_onion_fly_p,
          value_if_not = 0
        ))
      }
    } else {
      return(0)
    }
  }
  
  # Glasflügelzikade – trocken & warm bevorzugt
  get_leafhopper_stress <- function(Tavg, Prec) {
    if (Tavg >= Tavg_leafhopper_lower_limit_risk_medium_p & Prec < Prec_leafhopper_upper_limit_risk_medium_p) {
      return(chance_event(
        chance  = risk_leafhopper_medium_p,
        value_if = yield_reduction_leafhopper_p,
        value_if_not = 0,
      ))
    } else {
      return(0)
    }
  }
  
  # Drahtwurm – moderate Temperaturen, feuchte Böden (Regen als Proxy)
  get_wireworm_stress <- function(Tavg, Prec) {
    if (Tavg >= Tavg_wireworm_lower_limit_all_risks_p & Tavg <= Tavg_wireworm_upper_limit_all_risks_p & Prec >= Prec_wireworm_lower_limit_all_risks_p) {
      return(chance_event(
        chance  = risk_wireworm_medium_p,
        value_if = yield_reduction_wireworm_p,
        value_if_not = 0,
      ))
    } else {
      return(0)
    }
  }
  
  #####  Apply stresss functions
  
  season_risks <- lapply(weather_scenario_list, function(df) {
    
    # Helper to calculate mean values safely
    safe_mean <- function(x, condition) { # MODIFIED
      if (any(condition, na.rm = TRUE)) {
        mean(x[condition], na.rm = TRUE)
      } else {
        0 # Return 0 if no valid elements for averaging in this phase
      }
    }
    safe_max  <- function(x, condition) { # MODIFIED
      if (any(condition, na.rm = TRUE)) {
        max(x[condition], na.rm = TRUE)
      } else {
        0 # Return 0 if no valid elements for max in this phase
      }
    }
    
    # --- EMERGENCE PHASE STRESSORS
    emergence_filter <- df$in_emergence_phase == TRUE
    
    emergence_Tavg <- safe_mean(df$Tavg, emergence_filter)
    emergence_Prec <- safe_mean(df$Prec, emergence_filter)
    emergence_PAR  <- safe_mean(df$PAR, emergence_filter)
    emergence_consec_wet <- safe_max(df$day_consec_wet, emergence_filter)
    
    emergence_stress <- 0
    if (any(emergence_filter, na.rm = TRUE)) {
      emergence_stress <- emergence_stress +
        get_weed_pressure_stress(PAR = emergence_PAR, Prec = emergence_Prec, day_consec_wet = emergence_consec_wet) +
        get_fusarium_stress(Prec = emergence_Prec, Tavg = emergence_Tavg, day_consec_wet = emergence_consec_wet) +
        get_onion_fly_stress(Tavg = emergence_Tavg, Prec = emergence_Prec) +
        get_wireworm_stress(Tavg = emergence_Tavg, Prec = emergence_Prec)
    }
    
    # --- VEGETATIVE PHASE STRESSORS
    vegetative_filter <- df$in_vegetative_phase == TRUE
    
    vegetative_Tavg <- safe_mean(df$Tavg, vegetative_filter)
    vegetative_Prec <- safe_mean(df$Prec, vegetative_filter)
    vegetative_PAR  <- safe_mean(df$PAR, vegetative_filter)
    vegetative_consec_wet <- safe_max(df$day_consec_wet, vegetative_filter)
    
    vegetative_stress <- 0
    if (any(vegetative_filter, na.rm = TRUE)) {
      vegetative_stress <- vegetative_stress +
        get_fusarium_stress(Prec = vegetative_Prec, Tavg = vegetative_Tavg, day_consec_wet = vegetative_consec_wet) +
        get_downy_mildew_stress(Prec = vegetative_Prec, Tavg = vegetative_Tavg, day_consec_wet = vegetative_consec_wet) +
        get_thrips_stress(Tavg = vegetative_Tavg, Prec = vegetative_Prec) +
        get_leafhopper_stress(Tavg = vegetative_Tavg, Prec = vegetative_Prec) +
        get_onion_fly_stress(Tavg = vegetative_Tavg, Prec = vegetative_Prec)
    }
    
    # --- BULBING PHASE STRESSORS
    bulbing_filter <- df$in_bulbing_phase == TRUE
    
    bulbing_Tavg <- safe_mean(df$Tavg, bulbing_filter)
    bulbing_Prec <- safe_mean(df$Prec, bulbing_filter)
    bulbing_PAR  <- safe_mean(df$PAR, bulbing_filter)
    bulbing_consec_wet <- safe_max(df$day_consec_wet, bulbing_filter)
    
    bulbing_stress <- 0
    if (any(bulbing_filter, na.rm = TRUE)) {
      bulbing_stress <- bulbing_stress +
        get_botrytis_stress(Prec = bulbing_Prec, Tavg = bulbing_Tavg, day_consec_wet = bulbing_consec_wet) +
        get_fusarium_stress(Prec = bulbing_Prec, Tavg = bulbing_Tavg, day_consec_wet = bulbing_consec_wet) +
        get_downy_mildew_stress(Prec = bulbing_Prec, Tavg = bulbing_Tavg, day_consec_wet = bulbing_consec_wet) +
        get_onion_fly_stress(Tavg = bulbing_Tavg, Prec = bulbing_Prec)
    }
    
    # --- MATURATION PHASE STRESSORS
    maturation_filter <- df$in_maturation_phase == TRUE
    
    maturation_Tavg <- safe_mean(df$Tavg, maturation_filter)
    maturation_Prec <- safe_mean(df$Prec, maturation_filter)
    maturation_PAR  <- safe_mean(df$PAR, maturation_filter)
    maturation_consec_wet <- safe_max(df$day_consec_wet, maturation_filter)
    
    maturation_stress <- 0
    if (any(maturation_filter, na.rm = TRUE)) {
      maturation_stress <- maturation_stress +
        get_botrytis_stress(Prec = maturation_Prec, Tavg = maturation_Tavg, day_consec_wet = maturation_consec_wet) +
        get_fusarium_stress(Prec = maturation_Prec, Tavg = maturation_Tavg, day_consec_wet = maturation_consec_wet) +
        get_downy_mildew_stress(Prec = maturation_Prec, Tavg = maturation_Tavg, day_consec_wet = maturation_consec_wet) +
        get_onion_fly_stress(Tavg = maturation_Tavg, Prec = maturation_Prec)
    }
    
    # --- RETURN: Named list per phase
    return(list(
      emergence_phase_stress = emergence_stress,
      vegetative_phase_stress = vegetative_stress,
      bulbing_phase_stress = bulbing_stress,
      maturation_phase_stress = maturation_stress
    ))
  })
  
  
  # Combine weather data and seasonal risk info
  weather_scenario_list <- Map(function(weather_df, risks) {
    # Add the risk variables as new columns (repeated to match number of rows)
    for (name in names(risks)) {
      weather_df[[name]] <- risks[[name]]
    }
    return(weather_df)
  }, weather_scenario_list, season_risks)
  
  
  
  ## Calculate Biomass Growth per Phase
  
  calculate_biomass_daily_emergence <- function(PAR, LAI_emergence_p, Tavg, Prec, emergence_phase_stress, LUE_onion_p) {
    
    
    # Temperature effect
    f_T <- ifelse(Tavg >= f_T_1_lower_p & Tavg <= f_T_1_upper_p, 1,
                  ifelse(Tavg < f_T_0_lower_p | Tavg > f_T_0_upper_p, 0, 0.5))
    
    # Water availability
    f_W <- ifelse(Prec >= f_W_1_lower_p & Prec <= f_W_1_upper_p, 1,
                  ifelse(Prec < f_W_0.5_p, 0.5, 0.7))
    
    # Stress (same scalar applied to all days in a phase)
    f_S <- 1 - emergence_phase_stress
    
    # Vectorized light interception & biomass increment
    delta_B <- LUE_onion_p * PAR * (1 - exp(-k * LAI_emergence_p)) * f_T * f_W * f_S
    return(delta_B)
  }
  
  calculate_biomass_daily_vegetative <- function(PAR, LAI_veg_p, Tavg,
                                                 Prec, vegetative_phase_stress, LUE_onion_p) {
    
    # Temperature effect
    f_T <- ifelse(Tavg >= f_T_1_lower_p & Tavg <= f_T_1_upper_p, 1,
                  ifelse(Tavg < f_T_0_lower_p | Tavg > f_T_0_upper_p, 0, 0.5))
    
    # Water availability
    f_W <- ifelse(Prec >= f_W_1_lower_p & Prec <= f_W_1_upper_p, 1,
                  ifelse(Prec < f_W_0.5_p, 0.5, 0.7))
    
    f_S <- 1 - vegetative_phase_stress
    
    delta_B <- LUE_onion_p * PAR * (1 - exp(-k * LAI_veg_p)) * f_T * f_W * f_S
    return(delta_B)
  }
  
  calculate_biomass_daily_bulbing <- function(PAR, LAI_bulbing_p, Tavg, Prec, bulbing_phase_stress, LUE_onion_p) {
    
    
    # Temperature effect
    f_T <- ifelse(Tavg >= f_T_1_lower_p & Tavg <= f_T_1_upper_p, 1,
                  ifelse(Tavg < f_T_0_lower_p | Tavg > f_T_0_upper_p, 0, 0.5))
    
    # Water availability
    f_W <- ifelse(Prec >= f_W_1_lower_p & Prec <= f_W_1_upper_p, 1,
                  ifelse(Prec < f_W_0.5_p, 0.5, 0.7))
    
    f_S <- 1 - bulbing_phase_stress
    
    delta_B <- LUE_onion_p * PAR * (1 - exp(-k * LAI_bulbing_p)) * f_T * f_W * f_S
    return(delta_B)
  }
  
  calculate_biomass_daily_maturation <- function(PAR, LAI_maturation_p, Tavg, Prec, maturation_phase_stress, LUE_onion_p) {
    
    
    # Temperature effect
    f_T <- ifelse(Tavg >= f_T_1_lower_p & Tavg <= f_T_1_upper_p, 1,
                  ifelse(Tavg < f_T_0_lower_p | Tavg > f_T_0_upper_p, 0, 0.5))
    
    # Water availability
    f_W <- ifelse(Prec >= f_W_1_lower_p & Prec <= f_W_1_upper_p, 1,
                  ifelse(Prec < f_W_0.5_p, 0.5, 0.7))
    
    f_S <- 1 - maturation_phase_stress
    
    delta_B <- LUE_onion_p * PAR * (1 - exp(-k * LAI_maturation_p)) * f_T * f_W * f_S
    return(delta_B)
  }
  
  #### Apply functions for all phases
  
  biomass_all_scenarios <- lapply(weather_scenario_list, function(df) {
    
    # Emergence phase
    df_emergence <- df[df$in_emergence_phase, ]
    biomass_emergence <- if (nrow(df_emergence) > 0) {
      sum(calculate_biomass_daily_emergence(
        PAR = df_emergence$PAR,
        LAI_emergence_p = LAI_emergence_p,
        Tavg = df_emergence$Tavg,
        Prec = df_emergence$Prec,
        emergence_phase_stress = df$emergence_phase_stress[1],
        LUE_onion_p = LUE_onion_p
      ), na.rm = TRUE)
    } else 0
    
    # Vegetative phase
    df_veg <- df[df$in_vegetative_phase, ]
    biomass_veg <- if (nrow(df_veg) > 0) {
      sum(calculate_biomass_daily_vegetative(
        PAR = df_veg$PAR,
        LAI_veg_p = LAI_veg_p,
        Tavg = df_veg$Tavg,
        Prec = df_veg$Prec,
        vegetative_phase_stress = df$vegetative_phase_stress[1],
        LUE_onion_p = LUE_onion_p
      ), na.rm = TRUE)
    } else 0
    
    # Bulbing phase
    df_bulbing <- df[df$in_bulbing_phase, ]
    biomass_bulbing <- if (nrow(df_bulbing) > 0) {
      sum(calculate_biomass_daily_bulbing(
        PAR = df_bulbing$PAR,
        LAI_bulbing_p = LAI_bulbing_p,
        Tavg = df_bulbing$Tavg,
        Prec = df_bulbing$Prec,
        bulbing_phase_stress = df$bulbing_phase_stress[1],
        LUE_onion_p = LUE_onion_p
      ), na.rm = TRUE)
    } else 0
    
    # Maturation phase
    df_maturation <- df[df$in_maturation_phase, ]
    biomass_maturation <- if (nrow(df_maturation) > 0) {
      sum(calculate_biomass_daily_maturation(
        PAR = df_maturation$PAR,
        LAI_maturation_p = LAI_maturation_p,
        Tavg = df_maturation$Tavg,
        Prec = df_maturation$Prec,
        maturation_phase_stress = df$maturation_phase_stress[1],
        LUE_onion_p = LUE_onion_p
      ), na.rm = TRUE)
    } else 0
    
    # Output as list
    # Corrected: Calculate total_biomass_current_scenario before using it for total_yield
    total_biomass_current_scenario <- biomass_emergence + biomass_veg + biomass_bulbing + biomass_maturation
    
    list(
      emergence   = biomass_emergence,
      vegetative  = biomass_veg,
      bulbing     = biomass_bulbing,
      maturation  = biomass_maturation,
      total       = total_biomass_current_scenario, # Using the calculated total
      total_yield_per_ha = (total_biomass_current_scenario * HI_onions_p * onions_per_ha_p * dry_onion_weight_p) / 1000000 # Using the calculated total
    )
  })
  
  return(biomass_all_scenarios) # This returns the list of biomass results for each scenario
  
}


# Run the Monte Carlo simulation using the model function
model_mc_simulation <- mcSimulation(estimate = as.estimate(input_variables),
                                    model_function = onion_climate_impact,
                                    numberOfModelRuns = num_simulations_c,
                                    functionSyntax = "plainNames")

saveRDS(model_mc_simulation, "MC_results/mc_onions.RDS")
#write.csv(model_mc_simulation, "MC_results/mc_onions.csv")

onions<-readRDS("MC_results/mc_onions.RDS")
#onions<-read.csv("MC_results/mc_onions.csv")

hist_all <- onions$y[, grepl("historical", names(onions$y))]
ssp1_all <- onions$y[, grepl("ssp1", names(onions$y))]
ssp2_all <- onions$y[, grepl("ssp2", names(onions$y))]
ssp3_all <- onions$y[, grepl("ssp3", names(onions$y))]
ssp5_all <- onions$y[, grepl("ssp5", names(onions$y))]

hist_all_long <- hist_all %>%
  pivot_longer(cols = everything(), names_to = "Name", values_to = "Value")
hist_all_long$Scenario <- rep("historical", length.out = nrow(hist_all_long))

ssp1_all_long <- ssp1_all %>%
  pivot_longer(cols = everything(), names_to = "Name", values_to = "Value")
ssp1_all_long$Scenario <- rep("ssp1", length.out = nrow(ssp1_all_long))

ssp2_all_long <- ssp2_all %>%
  pivot_longer(cols = everything(), names_to = "Name", values_to = "Value")
ssp2_all_long$Scenario <- rep("spp2", length.out = nrow(ssp2_all_long))

ssp3_all_long <- ssp3_all %>%
  pivot_longer(cols = everything(), names_to = "Name", values_to = "Value")
ssp3_all_long$Scenario <- rep("ssp3", length.out = nrow(ssp3_all_long))

ssp5_all_long <- ssp5_all %>%
  pivot_longer(cols = everything(), names_to = "Name", values_to = "Value")
ssp5_all_long$Scenario <- rep("ssp5", length.out = nrow(ssp5_all_long))

scenarios_long<-rbind(hist_all_long,ssp1_all_long,ssp2_all_long,ssp3_all_long,ssp5_all_long)
scenarios_long$Name <- gsub("weather_historical.", "", scenarios_long$Name, fixed = T)
scenarios_long$Name <- gsub("weather_ssp126.", "", scenarios_long$Name, fixed = T)
scenarios_long$Name <- gsub("weather_ssp245.", "", scenarios_long$Name, fixed = T)
scenarios_long$Name <- gsub("weather_ssp370.", "", scenarios_long$Name, fixed = T)
#scenarios_long$Name, fixed = T)
scenarios_long$Name <- gsub("weather_ssp585.", "", scenarios_long$Name, fixed = T)
summary(scenarios_long$Name)


ggplot(scenarios_long, aes(x=Name, y=Value))+
  geom_boxplot()+
  ylab("g/plant &  t/ha")+
  xlab("Biomass of Growth Phase")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 25))+
  facet_wrap(~ Scenario, ncol = 5)