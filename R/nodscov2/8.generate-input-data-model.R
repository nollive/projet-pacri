################################################################################
##                  Generate input data for the simulation of 
##                            epidemics in ICU wards
################################################################################

# Librairies
library(data.table)
library(tidyverse)
library(Rcpp)
rm(list=ls())

# Conditions
networks = c("herriot", "poincare")
thresholds = c(30*2, 60*2, 90*2)

conditions = expand.grid(networks, thresholds)


for (r in 1:nrow(conditions)) {
  
  network = conditions[r,1]
  threshold = conditions[r,2]
  
  ## Load data--------------------------------------------------------------------
  loc_path <- file.path("data","data-synthetic-graphs","loc", paste0(network, "-simulated-reconstructed-locations-", threshold,".rda"))
  load(loc_path)
  
  ## Modification of the global objects-------------------------------------------
  # Create admission_sim dataframe 
  admission_sim = admission %>%
    mutate(info = ifelse(status == 'PE', 1, 0)) %>%
    select(id,info) %>%
    left_join(., rooms, by = 'id') %>%
    mutate(room = ifelse(info == 1, as.integer(rooms$id_room[rooms$id == 'NS']), as.integer(room)),
           info = as.integer(info)) %>%
    select(id,info,room)
  
  # Create global status dataframe that will be used to store the epidemic dynamics
  global_status = data.frame(
    id = admission$id, 
    t_inf = as.integer(-1), 
    t_incub = as.integer(-1), 
    t_recover = as.integer(-1), 
    inf_by = '', 
    inf_room = as.integer(-1)
  )
  
  # Create global_environment list to store aerosol quantities at each time step
  env_t = rooms %>%
    group_by(room, id_room) %>%
    summarise(volume = max(volume), .groups = "drop") %>%
    mutate(env = 0, id_room = as.integer(id_room)) %>%
    select(room, id_room, env, volume)
  
  global_environment <- lapply(1:n_subdivisions, function(t){
    return (env_t)
  })
  rm(env_t)
  
  ## Truncation of the interaction objects----------------------------------------
  ## This chunk has to be removed after implementing the synthetic network algorithm
  end_date - begin_date ##end of last interaction
  t_begin = 1
  t_end = length(global_environment)
  
  ## TRUNCATED INFO EXTRACTED FROM TOTAL STUDY (KEEP 24HOURS OF THE 36HOURS)
  global_interaction = map(global_interaction, function(df) {
    df = df %>% 
      select(from, to, time)# ,location %>%
      #mutate(location = as.integer(location))
    return(df)
  }) 
  
  # truncated_data <- lapply(1:length(truncated_localization), function(t){
  #   df <- data.frame(id = truncated_localization[[t]]$id,
  #                    info = admission_sim$info,
  #                    #room = admission_sim$room,
  #                    #interacting = ifelse(id %in% truncated_interaction$from | id %in% truncated_interaction$to, T, F)
  #                    localization_ti = truncated_localization[[t]]$localization,
  #                    lambda_c = 0,
  #                    lambda_e = 0
  #   )
  #   df <- df %>% filter(localization_ti != -1)
  #   return(df)
  # })
  
  if (identical(colnames(paths), admission$id)) {
    global_data <- lapply(1:length(global_environment), function(t){
      df <- data.frame(id = colnames(paths),
                       info = admission_sim$info,
                       #room = admission_sim$room,
                       #interacting = ifelse(id %in% truncated_interaction$from | id %in% truncated_interaction$to, T, F)
                       location_ti = as.integer(paths[t, ]),
                       lambda_c = 0,
                       lambda_e = 0
      )
      df <- df %>% filter(location_ti > 0)
      return(df)
    }) 
  } else {
    warning("Ids are not ordered correctly")
  }
  
  ## List of model parameters-----------------------------------------------------
  # Fixed parameters 
  B <- 0.023 * 60 * 24 # Breathing rate (0.023 m3/min)
  # mu_air <- 1 # Mechanical ventilation in ISO 8 rooms (>10 ACH) but a bit high considering that there is no ventilation in Raymond Poincaré
  # Mechanical ventilation is not perfect either 
  #mu_ventilation <-  4 * 24 # Mechanical ventilation, Quanta removal (4-6 Air changes/h)
  mu_inac <- (log(2)/1.1) * 24 # Quanta inactivation (computed with viral half life -> 0.63 quanta inactivated/h) 
  mu <- mu_inac #+ (mu_ventilation - mu_air)
  dt <- 30 # Time step
  tau <- 60 * 60 *24 # Seconds in 1 day
  deltat <- dt/tau
  env_model = "linear"
  nu <- 1.8e7 * 2 * 24 # (1.8e7 RNA copies per 30-min) Lai et al., Clinical Infectious Diseases (2024)
  
  # Parameters explored in the simulation study (tested to give the same 20% SAR)
  beta_c <- 1/2 # Transmission rate for close contact interaction /h
  beta_e <- 1e-8
  
  
  ## Save the final object used as input of the epidemic simulation algorithm-----
  # Pay attention: it might take a while to save the 
  # the global objects when the number of days is high
  n_days = floor(as.numeric(difftime(end_date, begin_date, units = "day")))
  
  save(global_interaction,
       global_data,
       global_status,
       global_environment,
       admission,
       beta_c,
       beta_e,
       nu,
       mu,
       B,
       env_model,
       dt,
       n_subdivisions,
       n_days,
       t_begin,
       t_end,
       begin_date,
       end_date,
       file = file.path("out", paste0("parameters-synthetic-", network, "-", threshold, ".rda"))
  )
}



# ## FILTER INDIVIDUALS NOT INTERACTING
# id_interacting <- #do.call(rbind,truncated_interaction) %>%
#   do.call("rbind", global_interaction) %>%
#   pivot_longer(cols = c(from, to), names_to = 'direction', values_to = 'individual') %>%
#   distinct(individual) %>%
#   pull()
# 
# ## FILTER THE DFS TO REMOVE USELESS ROWS (INDIVIDUALS NOT INTERACTING)
# filter_interacting <- function(list_df, id_interacting) {
#   lapply(list_df, function(df) {
#     df %>% 
#       filter(id %in% id_interacting)
#   })
# }
# 
# # truncated_localization <- filter_interacting(truncated_localization, id_interacting )
# # truncated_status <- global_status %>% filter(id %in% id_interacting)
# global_location <- filter_interacting(global_location, id_interacting )
# global_status <- global_status %>% filter(id %in% id_interacting)
# admission <- admission %>% filter(id %in% id_interacting)
# admission_sim <- admission_sim  %>% filter(id %in% id_interacting)
# 
