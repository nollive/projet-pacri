################################################################################
##          Parallelized code to reconstruct individual locations
##
################################################################################

# Librairies
library(tidyverse)
library(foreach)
library(doParallel)
library(Rcpp)
rm(list = ls())

# Source helper functions
source("R/nodscov2/helper-functions.R")
source("R/nodscov2/helper-functions-cpp.R")

# Paths
nodscov2_synthetic_path <- file.path("data", "data-synthetic-graphs") 
nodscov2_path <- file.path("data", "data-nodscov2") 

# Network to study
networks = c("herriot-simulated", "poincare-simulated", "herriot-observed", "poincare-observed")

registerDoParallel(4)
verif = foreach (network=networks, .combine=c) %dopar% {
  ## Process verification---------------------------------------------------------
  checkpoints = paste0(networks, "\n\n")
  
  ## Load data--------------------------------------------------------------------
  net = gsub("-.*", "", network)
  db = gsub(".*-", "", network)
  
  if (db == "simulated") {
    # Admission data
    admission = read.csv2(file.path(nodscov2_synthetic_path, "input", paste0("admission_", net, ".csv"))) %>%
      mutate(firstDate = as_date(firstDate), lastDate = as_date(lastDate))
    
    # HCW schedule 
    agenda = read.csv2(file.path(nodscov2_synthetic_path, "input", paste0("agenda_", net,".csv"))) %>%
      mutate(firstDate = as_datetime(firstDate), lastDate = as_datetime(lastDate))
    
    # Interaction data 
    interactions = read.csv2(file.path(nodscov2_synthetic_path, "full", net, "matContactBuiltSimulatedCtcNetworks1_oneday.csv")) %>%
      mutate(date_posix = as_datetime(date_posix), from_status = gsub("-.*", "", from), to_status = gsub("-.*", "", to)) 
  }
  
  if (db == "observed") {
    # Admission data
    admission = read.csv2(file.path(nodscov2_path, "clean", paste0("admission_cleaned_", net, ".csv"))) %>%
      mutate(firstDate = as_date(firstDate), lastDate = as_date(lastDate))
    
    # HCW schedule 
    agenda = read.csv2(file.path(nodscov2_path, "clean", paste0("agenda_cleaned_", net,".csv"))) %>%
      mutate(firstDate = as_datetime(firstDate), lastDate = as_datetime(lastDate))
    
    # Interaction data 
    interactions = read.csv2(file.path(nodscov2_path, "clean", paste0("interaction_cleaned_", net, ".csv"))) %>%
      mutate(date_posix = as_datetime(date_posix), from_status = gsub("-.*", "", from), to_status = gsub("-.*", "", to))
  }
  
  print(network)
  
  # Patient rooms 
  patient_rooms = read.csv2(file.path(nodscov2_synthetic_path, "loc", paste0("patient_rooms_", net, ".csv")))
  
  ## Get category ids-------------------------------------------------------------
  # Ids of participants by category
  id_paramedical <- admission %>% filter(cat == "Paramedical") %>% pull(id)
  id_medical <- admission %>% filter(cat == "Medical") %>% pull(id)
  id_hcw <- admission %>% filter(status == "PE") %>% pull(id)
  id_patient <- admission %>% filter(status == "PA") %>% pull(id)
  
  ## Trim interactions to match schedule------------------------------------------
  # Cut interactions
  checkpoints = paste0(checkpoints, "Number of rows before cutting interactions: ", nrow(interactions), "\n")
  interactions = interactions %>%
    mutate(n = 1:n()) %>%
    nest(.by = n) %>%
    mutate(data = map(data, trim_interactions_agenda, admission, agenda)) %>%
    unnest(data) %>%
    filter(length > 0, before_schedule == F) %>%
    select(from, to, date_posix, length, from_status, to_status)
  checkpoints = paste0(checkpoints, "Number of rows after cutting interactions: ", nrow(interactions), "\n")
  
  # Verify that interactions are within the schedule of healthcare workers
  not_in_hcw_schedule = interactions %>%
    filter(from %in% id_hcw | to %in% id_hcw) %>%
    mutate(n = 1:n()) %>%
    nest(.by = n) %>%
    mutate(data = map(data, hcw_coherent_interaction, agenda)) %>%
    select(data) %>%
    unnest(data) %>%
    summarise(sum(!data))
  checkpoints = paste0(checkpoints, "Number of interactions not in HCW schedule: ", not_in_hcw_schedule, "\n")
  
  # Get number of interactions that are not within the hospitalization stays of patients
  not_in_patient_schedule = interactions %>%
    filter(from %in% id_patient | to %in% id_patient) %>%
    mutate(n = 1:n()) %>%
    nest(.by = n) %>%
    mutate(data = map(data, patient_coherent_interaction, admission)) %>%
    select(data) %>%
    unnest(data) %>%
    summarise(sum(!data))
  checkpoints = paste0(checkpoints, "Number of interactions not in patient schedule: ", not_in_patient_schedule, "\n")
  
  # Verify that all individuals interact
  id_interacting <- interactions %>%
    pivot_longer(cols = c(from, to), names_to = "direction", values_to = "individual") %>%
    distinct(individual) %>% 
    arrange(individual) %>%
    pull()
  id_total <- admission %>% distinct(id) %>% arrange(id) %>% pull()
  checkpoints = paste0(checkpoints, "Do all individuals have at least one interaction? ", identical(id_total, id_interacting))
  
  ## Create room dataframe--------------------------------------------------------
  # Ids of patient rooms
  id_room_patient = as.character(unique(patient_rooms$room))
  
  # Id of double room
  if (file.exists(file.path(nodscov2_path, "clean", paste0("double_rooms_", net, ".rda")))) {
    load(file.path(nodscov2_path, "clean", paste0("double_rooms_", net, ".rda")))
    double_rooms_id = sort(unique(patient_rooms$room[patient_rooms$id %in% unlist(double_rooms)]))
  } else {
    double_rooms_id = ""
  }
  
  # Dataframe with all rooms
  rooms = patient_rooms %>%
    mutate(
      volume = ifelse(room %in% double_rooms_id, 30*2.2, 20*2.2), #18 square meters * 2.2 height except for double room
      id_room = as.character(room),
      room = as.character(room)
    ) %>% 
    bind_rows(.,
              data.frame(room = "Medical Staff Room", id = "M-R", volume = 30 * 2.2, id_room = "50"),
              data.frame(room = "Paramedical Staff Room", id = "PM-R", volume =  40 * 2.2, id_room = "51"),
              data.frame(room = "Nursing station", id = "NS", volume = 20 * 2.2, id_room = "52"),
              data.frame(room = "Office", id = "M-O", volume = 20 * 2.2, id_room = "53"),
              data.frame(room = "Corridor", id = "C", volume = 200 * 2.2, id_room = "54")
    ) 
  
  
  ## Create interaction list------------------------------------------------------
  # Time information
  begin_date <- floor_date(min(interactions$date_posix), "minutes")
  end_date <- ceiling_date(max(interactions$date_posix + interactions$length), "minutes")
  time_spent <- difftime(end_date, begin_date, units = "secs")  
  
  n_subdivisions <- as.numeric(time_spent) / 30 ## Time spent * hours * minutes * 2 (number of time subdivisions)
  ## floor because of the manner we extract the interactions (=30sec, i.e. at  the end, interactions are <30sec -> we dont take these)
  
  # Global interaction object
  global_interaction <- lapply(1:n_subdivisions, get_global_interaction, 
                               interactions = interactions %>% mutate(date_posix = as.numeric(difftime(date_posix, begin_date, "sec"))) )
  
  # List of interaction clusters
  clusters <- lapply(global_interaction, get_clusters)
  
  
  ## Assign locations-------------------------------------------------------------
  # Threshold of the time spent outside ward while wearing a sensor
  thresholds = c(30*2, 60*2, 90*2)
  
  for (threshold in thresholds) {
    ## Assign locations when individuals are interacting
    paths = get_global_location(clusters, admission, n_subdivisions, rooms)
    
    ## Assign location when not interacting
    # Add when individuals are not present in the ward and add their room 
    # otherwise (Patients)
    identical(colnames(paths[, id_patient]), id_patient)
    for (i in id_patient) {
      paths[,i] = patient_locations(
        i, 
        paths[, i], 
        admission = admission, 
        begin_date = begin_date, 
        end_date = end_date, 
        patient_rooms = patient_rooms)
    }
    
    # paths[,id_patient] = as_tibble(mapply(
    #   patient_locations,
    #   patient_id = id_patient,
    #   patient_loc = paths[, id_patient],
    #   MoreArgs = list(admission = admission, begin_date = begin_date, end_date = end_date, patient_rooms = patient_rooms)
    # ))
    
    # Add when individuals are not present in the ward and smooth
    # their individual trajectory (HCW)
    identical(colnames(paths[, id_hcw]), id_hcw)
    patient_rooms_ids = as.character(sort(unique(patient_rooms$room)))
    
    for (i in id_hcw) {
      paths[,i] = hcw_locations(
        i, 
        paths[, i], 
        threshold = threshold, 
        agenda = agenda, 
        admission = admission, 
        begin_date = begin_date, 
        end_date = end_date, 
        rooms = rooms, 
        id_medical = id_medical)
    }
    
    # Herriot
    # PE-01381 has location assigned when not present in the ward
    # PE-01384 has location assigned when not present in the ward
    
    # locations[,id_hcw] = as_tibble(mapply(
    #   hcw_locations, hcw_id = id_hcw, hcw_loc = locations[,id_hcw], 
    #   MoreArgs= list(threshold = threshold, agenda = agenda, admission = admission, begin_date = begin_date, end_date = end_date, patient_rooms_ids = patient_rooms_ids, id_medical = id_medical)
    # ))
    
    # Convert locations into a list
    # write.csv2(locations, file.path(nodscov2_synthetic_path, "loc", paste0("locations_long_format_", threshold, ".csv")), quote = F, row.names = F)
    # locations = locations %>%
    #   mutate(time = 1:n()) %>%
    #   pivot_longer(-time, names_to = "id", values_to = "location")
    # global_location = split(locations, f = locations$time)
    # rm(locations)
    
    ## Save data for simulations----------------------------------------------------
    save(clusters,
         paths,
         global_interaction,
         admission,
         rooms,
         net, 
         threshold,
         id_hcw,
         id_patient,
         begin_date,
         end_date,
         n_subdivisions,
         file = file.path(nodscov2_synthetic_path, "loc", paste0(network, "-reconstructed-locations-", threshold, ".rda")))
    
  }
  
    return(checkpoints)
}

cat(paste0(verif, collapse = "\n\n"),file="out/loc-reconstruction/process-verification.txt")

