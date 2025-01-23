################################################################################
##          Support functions, dictionaries and variables
##                   to analyze interaction data
################################################################################

# Dictionaries------------------------------------------------------------------
dict_cat = c("aux nurse" = "Paramedical", 
             "nurse" = "Paramedical",
             "student nurse" = "Paramedical",
             "reeducation staff" = "Paramedical",
             "ext physician" = "Medical",
             "physician" = "Medical",
             "patient" = "Patient")

dict_cat_initial = c("administration" = "Other",
                     "aux nurse" = "HCW", 
                     "ext physician" = "HCW",
                     "investigation" = "Other",
                     "logistic" = "Other",
                     "nurse" = "HCW",
                     "other" = "Other",
                     "physician" = "HCW",
                     "reeducation staff" = "HCW",
                     "student nurse" = "HCW",
                     "visitor" = "Visitor"
                    )

dict_cat_num = c(
  "1" = "nurse",
  "2" = "aux nurse",
  "3" = "reeducation staff",
  "4" = "physician",
  "5" = "ext physician",
  "6" = "administration",
  "7" = "logistic",
  "8" = "investigation",
  "9" = "patient",
  "10" = "visitor",
  "11" = "other",
  "12" = "student nurse"
)

dict_ward = c(
  "AURÉLIEN DINH" = "Maladies infectieuses",
  "DENIS MALVY" = "Maladies infectieuses",
  "DJILLALI ANNANE" = "Reanimation",
  "ELSA KERMORVANT" = "Reanimation pediatrique", 
  "GÉRALDINE MARTIN-GAUJARD" = "Geriatrie",
  "GÉRARD CHERON" = "Urgence pediatrique",
  "JULIE TOUBIANA" = "Pediatrie",
  "KARIM TAZAROURTE" = "Urgences",
  "LAURENT ARGAUD" = "Reanimation",
  "MAGALI GUICHARDON" = "Geriatrie",
  "MARION DOUPLAT" = "Urgences", 
  "OLIVIER BOUCHAUD" = "Maladies infectieuses", 
  "OLIVIER LAMBOTTE" = "Medecine interne", 
  "OLIVIER SITBON" = "Pneumologie", 
  "SÉBASTIEN BEAUNE" = "Urgences", 
  "THOMAS RIMMELE" = "Reanimation chir" 
)

dict_hosp = c(
  "AURÉLIEN DINH" = "APHP - RAYMOND POINCARÉ",
  "DENIS MALVY" = "CHU DE BORDEAUX - PELLEGRIN",
  "DJILLALI ANNANE" = "APHP - RAYMOND POINCARÉ",
  "ELSA KERMORVANT" = "APHP - NECKER ENFANTS MALADES", 
  "GÉRALDINE MARTIN-GAUJARD" = "HC DE LYON - EDOUARD HERRIOT",
  "GÉRARD CHERON" = "APHP - NECKER ENFANTS MALADES",
  "JULIE TOUBIANA" = "APHP - NECKER ENFANTS MALADES",
  "KARIM TAZAROURTE" = "HC DE LYON - EDOUARD HERRIOT",
  "LAURENT ARGAUD" = "HC DE LYON - EDOUARD HERRIOT",
  "MAGALI GUICHARDON" = "APHP - PAUL BROUSSE",
  "MARION DOUPLAT" = "HC DE LYON - SUD", 
  "OLIVIER BOUCHAUD" = "APHP - AVICENNE", 
  "OLIVIER LAMBOTTE" = "APHP - BICÊTRE", 
  "OLIVIER SITBON" = "APHP - BICÊTRE", 
  "SÉBASTIEN BEAUNE" = "APHP - AMBROISE PARÉ", 
  "THOMAS RIMMELE" = "HC DE LYON - EDOUARD HERRIOT" 
)

dict_scenarios = c("sim_1-4_20" = "Scenario 1",
                   "sim_1-2_16" = "Scenario 2",
                   "sim_3-4_12" = "Scenario 3",
                   "sim_1_8" = "Scenario 4",
                   "sim_3-2_5" = "Scenario 5"
)

dict_rooms = c(
  "1" = "Patient Room 1",
  "2" = "Patient Room 2",
  "3" = "Patient Room 3",
  "4" = "Patient Room 4",
  "5" = "Patient Room 5",
  "6" = "Patient Room 6",
  "8" = "Patient Room 8",
  "9" = "Patient Room 9",
  "10" = "Patient Room 10",
  "11" = "Patient Room 11",
  "12" = "Patient Room 12",
  "13" = "Patient Room 13",
  "14" = "Patient Room 14",
  "15" = "Patient Room 15",
  "16" = "Patient Room 16",
  "17" = "Patient Room 17",
  "18" = "Medical Staff Room",
  "19" = "Paramedical Staff Room",
  "20" = "Nursing Station",
  "21" = "Office",
  "22" = "Corridor"
)

# Plots-------------------------------------------------------------------------
pal = c('Medical' = "#5CD6D6", 'Paramedical' = "#A9B9E8", 'Patient' = "#FFA766", 'Room' = "#666699")
env_pal = c('Contact' = 'darkorange', 'Environment' =  'orchid')
algo_syn_pal = c("Observed" = "darkorchid", "Reconstructed"= "darkorange")

# Variables---------------------------------------------------------------------
noon_day1 = as_datetime("2020-05-06 12:00:00")
midnight_day1 = as_datetime("2020-05-06 00:00:00")
noon_day2 = as_datetime("2020-05-07 12:00:00")
midnight_day2 = as_datetime("2020-05-07 00:00:00")
noon_last_day = noon_day1 + 90*3600*24
midnight_last_day = midnight_day1 + 90*3600*24

# Functions to analyze data-----------------------------------------------------
# Get vector of hours between two POSIXct dates 
unroll_dates = function(df) {
  if ("DATEREMISE" %in% colnames(df)) {
    out = data.frame(
      date_list=seq(df$DATEREMISE, df$DATEREC, 3600)
    )
  }
  
  if ("firstDate" %in% colnames(df)) {
    out = seq(floor_date(as_datetime(df$firstDate), "hour"), floor_date(as_datetime(df$lastDate), "hour"), 3600)
  }
  
  return(out)
}

# Get vector of 30s steps between two dates
unroll_30s = function(df) {
  out = seq(df$firstDate, df$lastDate-30, 30)
  return(out)
}

# Get vector of days between two dates
unroll_days = function(df) {
  out = seq.Date(df$firstDate, df$lastDate, 1)
  return(out)
}

# In a dataframe, concatenate rows that correspond to consecutive time periods
concatenate_schedules = function(df) {
  out = df %>%
    arrange(firstDate) %>%
    mutate(tomerge = ifelse(firstDate %in% lag(lastDate), 1, 0)) %>%
    mutate(lastDate = case_when(lastDate %in% lead(firstDate) ~ lead(lastDate), .default = lastDate)) %>%
    filter(tomerge == 0) %>%
    select(-tomerge)
  return(out)
}


# Expand.grid for dataframes (combine all rows of two distinct dataframes)
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

# Create a new variable that corresponds to the pair id
# i.e. the concatenation of the ids of the participants in interactions sorted by 
# alphanumerical order
get_pair_id = function(df) {
  if (nrow(df) == 1) {
    df$pair = paste0(sort(c(df$from, df$to)), collapse = "_")
    return(df) 
  }
}

# Detect for a set of interactions between two individuals, overlapping contacts
get_overlapping_contacts = function(df) {
  df$overlap = F
  
  if (nrow(df) > 1) {
    for (i in 1:nrow(df)) {
      overlap = df$date_posix[i] >= df$date_posix & df$date_posix[i] <= df$date_posix+df$length
      
      if (sum(overlap) > 1) {
        j = which(overlap)
        j = j[j!=i]
        df$overlap[i] = T
        df$overlap[j] = T
      }
    }
  }
  
  return(df)
}

# Detect and fusion overlapping contacts
fusion_overlapping_contacts = function(df) {
  
  if (nrow(df) > 1) {
    for (i in 1:nrow(df)) {
      
      overlap = df$date_posix[i] >= df$date_posix & df$date_posix[i] <= df$date_posix+df$length
      
      if (sum(overlap) > 1) {
        # Contacts with overlap
        j = which(overlap)
        j = j[j!=i]
        
        # Start of interaction
        new_date_posix = min(df$date_posix[c(i,j)])
        
        # End of interaction
        new_end = max(df$date_posix[c(i,j)] + df$length[c(i,j)])
        new_length = as.numeric(difftime(new_end, new_date_posix, units = "secs"))
        
        # Replace date_posix and length
        df$date_posix[c(i,j)] = new_date_posix
        df$length[c(i,j)] = new_length
      }
    }
  }
  
  df = distinct(df)
  
  return(df)
  
}

# Verify whether interactions involving a patient lay within the hospitalization 
# stay of the patients 
patient_coherent_interaction = function(df, admission) {
  ids = unname(unlist(df[, c("from", "to")]))
  ids = ids[grepl("^PA-", ids)]
  firstDates = as_datetime(paste(admission$firstDate[admission$id %in% ids], "00:00:00"))
  lastDates = as_datetime(paste(admission$lastDate[admission$id %in% ids], "23:59:30"))
  
  out = all(firstDates <= df$date_posix & lastDates >= df$date_posix+df$length) 
  return(out)
}

# Verify whether interactions involving a HCW lay within their presence in the ward 
hcw_coherent_interaction = function(df, agenda) {
  ids = unname(unlist(df[, c("from", "to")]))
  ids = ids[grepl("^PE-", ids)]
  presence_dates = agenda %>%
    filter(id %in% ids, floor_date(firstDate, "hour") <= floor_date(df$date_posix, "hour"), floor_date(lastDate, "hour") >= floor_date(df$date_posix, "hour"))
  
  if (length(ids) != nrow(presence_dates)) {
    out = F
  } else {
    out = all(presence_dates$firstDate <= df$date_posix & presence_dates$lastDate >= df$date_posix+df$length) 
  }
  return(out)
}

#
outside_schedule = function(df, sensor_nodscov2) {
  
  one_slot = 0
  d = df$date_posix
  e = df$date_posix + df$length
  
  for (f in unname(unlist(df[,c("from", "to")]))) {
    df_in = sensor_nodscov2 %>% 
      filter(id == f) %>%
      select(id, firstDate_sensor, lastDate_sensor) %>%
      mutate(
        firstDate_in = as.numeric(firstDate_sensor <= d & d<= lastDate_sensor),
        lastDate_in = as.numeric(e <= lastDate_sensor & e >= firstDate_sensor) 
      ) %>%
      mutate(date_in = firstDate_in+lastDate_in)
    
    if (any(df_in$date_in == 2)) {
      df$in_schedule = T
      
    } else if (any(df_in$date_in == 1)) {
      one_slot = one_slot + 1
      if (df_in$firstDate_in[df_in$date_in == 1] == 1) {
        new_firstDate = df$date_posix
        new_length = as.numeric(difftime(df_in$lastDate_sensor[df_in$date_in == 1], new_firstDate, units = "secs"))
      } else {
        new_firstDate = df_in$firstDate_sensor[df_in$date_in==1]
        new_length = as.numeric(difftime(df$date_posix+df$length, new_firstDate, units = "secs"))
      }
      
      df$date_posix = new_firstDate
      df$length = new_length
      df$in_schedule = T
      
    } else {
      one_slot = one_slot + 1
      df$in_schedule = F
    }
  }
  
  if (one_slot > 1) message(paste("Following interaction with both participants not present in the ward:", df$from, df$to, d))
  return(df)
}

# Get location during schedule 
get_location_during_schedule = function(df, start_cut, end_cut) {
  day_to_select = start_cut:end_cut
  locations = lapply(seq_along(day_to_select), function(x) {
    out = global_data[day_to_select][[x]]
    out$time = day_to_select[x]
    return(out)
  })
  locations = do.call("rbind", locations) %>%
    filter(id == df$id) %>%
    select(id, location_ti, time) %>%
    rename(room = location_ti)
  return(locations)
}

# Functions to analyze networks-------------------------------------------------
# Function to return a set of network metrics (from Leclerc et al., 2024)
# https://gitlab.pasteur.fr/qleclerc/network_algorithm/-/tree/main?ref_type=heads
get_net_metrics = function(graph_data, adm_data, iter = 0, db_type = "Observed", network = "Hospital"){
  
  days = unique(graph_data$date_posix)
  
  data = data.frame(degrees = rep(0, length(days)),
                    densities = rep(0, length(days)),
                    transitivities = rep(0, length(days)),
                    assortativities = rep(0, length(days)),
                    efficiencies = rep(0, length(days)),
                    temp_corr = rep(0, length(days)),
                    iter = iter,
                    data = db_type,
                    network = network,
                    day = days)
  
  
  for(i in 1:length(days)){
    
    # full network
    data_d = graph_data %>%
      filter(date_posix == days[i])
    
    graph_d = graph_from_data_frame(data_d, directed = F)
    graph_d = simplify(graph_d)
    vertex_atts = data.frame(id = vertex_attr(graph_d, "name")) %>%
      left_join(adm_data, "id")
    graph_d = graph_d %>%
      set_vertex_attr("cat", value = vertex_atts$cat)
    
    data$degrees[i] = mean(degree(graph_d))
    data$densities[i] = edge_density(graph_d)
    data$transitivities[i] = transitivity(graph_d)
    data$assortativities[i] = assortativity_degree(graph_d, directed = F)
    data$efficiencies[i] = global_efficiency(graph_d, directed = F)
    
  }
  
  if (length(days) > 1) {
    data$temp_corr = c(0, temporal_correlation(graph_data)) 
  }
  
  return(data)
  
}

# Get the assortativity by degree
get_assortativity_degree = function(graph_data, adm_data, iter = 0, network = "Observed"){
  
  days = unique(graph_data$date_posix)
  
  data = data.frame(assortativities = rep(0, length(days)),
                    iter = iter,
                    network = network,
                    day = days)
  
  for(i in 1:length(days)){
    
    # full network
    data_d = graph_data %>%
      filter(date_posix == days[i])
    
    graph_d = graph_from_data_frame(data_d, directed = F)
    graph_d = simplify(graph_d)
    vertex_atts = data.frame(id = vertex_attr(graph_d, "name")) %>%
      left_join(adm_data, "id")
    graph_d = graph_d %>%
      set_vertex_attr("cat", value = vertex_atts$cat) 
    
    data$assortativities[i] = assortativity_degree(graph_d, directed = F)
    
  }
  return(data)
}

# Function to calculate the temporal correlation of a network (from Leclerc et al., 2024)
# https://gitlab.pasteur.fr/qleclerc/network_algorithm/-/tree/main?ref_type=heads
temporal_correlation = function(graph_data){
  
  temp_cor = c()
  
  for(i in 1:(length(unique(graph_data$date_posix))-1)){
    
    # full network
    data_t1 = graph_data %>%
      filter(date_posix == unique(graph_data$date_posix)[i])
    graph_t1 = graph_from_data_frame(data_t1, directed = F)
    graph_t1 = simplify(graph_t1)
    
    data_t2 = graph_data %>%
      filter(date_posix == unique(graph_data$date_posix)[i+1])
    graph_t2 = graph_from_data_frame(data_t2, directed = F)
    graph_t2 = simplify(graph_t2)
    
    temp_cor_indiv = c()
    
    for(indiv in names(V(graph_t1))){
      if(indiv %in% names(V(graph_t2))){
        
        neighbours_t1 = names(neighbors(graph_t1, indiv))
        neighbours_t2 = names(neighbors(graph_t2, indiv))
        temp_cor_indiv = c(temp_cor_indiv,
                           length(intersect(neighbours_t1, neighbours_t2))/sqrt(length(neighbours_t1)*length(neighbours_t2)))
        
      } else temp_cor_indiv = c(temp_cor_indiv,0)
      
    }
    
    temp_cor = c(temp_cor, mean(temp_cor_indiv))
  }
  
  return(temp_cor)
}

# Function to get the number of contacts by hour for the different pairs of participants
contact_numbers = function(db, db_type, network, analysis_type = "Final") {
  out = data.frame()
  
  if (analysis_type == "Final") {
    dict_pair = c("PA-PA" = "Patient-patient", "PA-PE" = "Patient-staff", "PE-PE" = "Staff-staff") 
    
    # Input data
    db = db %>% 
      select(from, to, date_posix, length) %>%
      mutate(
        type = case_when(
          substr(from, 1, 2) != substr(to, 1, 2) ~ "PA-PE",
          .default = paste0(substr(from, 1, 2), "-", substr(to, 1, 2))
        )
      )
  } else {
    dict_pair = c("PA-PA" = "Patient-Patient", "PA-PE" = "Patient-HCW", "PE-PE" = "HCW-HCW", 
                  "PA-OT" = "Patient-Other", "PE-OT" = "HCW-Other", "OT-OT" = "Other-Other",
                  "PA-VI" = "Patient-Visitor", "PE-VI" = "HCW-Visitor", "OT-VI" = "Other-Visitor",
                  "VI-VI" = "Visitor-Visitor") 
    
    # Input data
    db = db %>%
      mutate(
        type = case_when(
          from_cat == to_cat & from_cat == "Patient" ~ 'PA-PA',
          from_cat == to_cat & from_cat == "HCW" ~ 'PE-PE',
          from_cat == to_cat & from_cat == "Other" ~ 'OT-OT',
          from_cat == to_cat & from_cat == "Visitor" ~ 'VI-VI',
          from_cat != to_cat & from_cat %in% c("Patient", "HCW") & to_cat %in% c("Patient", "HCW") ~ "PA-PE",
          from_cat != to_cat & from_cat %in% c("Patient", "Other") & to_cat %in% c("Patient", "Other") ~ "PA-OT",
          from_cat != to_cat & from_cat %in% c("Patient", "Visitor") & to_cat %in% c("Patient", "Visitor") ~ "PA-VI",
          from_cat != to_cat & from_cat %in% c("HCW", "Visitor") & to_cat %in% c("HCW", "Visitor") ~ "PE-VI",
          from_cat != to_cat & from_cat %in% c("Other", "Visitor") & to_cat %in% c("Other", "Visitor") ~ "OT-VI",
          from_cat != to_cat & from_cat %in% c("Other", "HCW") & to_cat %in% c("Other", "HCW") ~ "PE-OT",
        )
      ) %>%
      select(from, to, date_posix, length, type)
  }
  
  for (p in names(dict_pair)) {
    out_temp = db %>%
      filter(type == p) %>%
      mutate(date_posix = floor_date(as_datetime(date_posix), "hour")) %>%
      select(-length) %>%
      distinct() %>%
      count(date_posix) %>%
      mutate(day = wday(date_posix, week_start = 1)) %>%
      mutate(day = as.character(day)) %>%
      mutate(day = replace(day, day %in% c("6", "7"), "Weekend")) %>%
      mutate(day = replace(day, day != "Weekend", "Weekday")) %>%
      mutate(date_posix = hour(date_posix)) %>%
      group_by(date_posix, day) %>%
      summarise(med = median(n),
                q25 = quantile(n, 0.25),
                q75 = quantile(n, 0.75), .groups = "drop") %>%
      mutate(type = dict_pair[p])
    
    out = bind_rows(out, out_temp)
  }
  
  out = out %>%
    mutate(date_posix = paste0(date_posix, ":00")) %>%
    mutate(
      data = db_type, 
      network = network
      )
  
  return(out)
}


# Functions to reconstruct individual locations---------------------------------
# Function to trim synthetic networks when contacts occur when an individual is not 
# present in the ward
trim_interactions_agenda = function(df, admission, agenda) {
  pa = c(df$from, df$to)[grepl("^PA-", c(df$from, df$to))]
  pe = c(df$from, df$to)[grepl("^PE-", c(df$from, df$to))]
  allLastDates = c()
  allFirstDates = c()
  
  df$before_schedule = FALSE
  
  if (length(pa)>0) {
    allLastDates = c(as_datetime(paste(admission$lastDate[admission$id %in% pa], "23:59:30")))
    allFirstDates = c(as_datetime(paste(admission$firstDate[admission$id %in% pa], "00:00:00")))
  } 
  
  if (length(pe)>0) {
    pe_firstDates = agenda %>%
      filter(id %in% pe, floor_date(firstDate, "hour") <= floor_date(df$date_posix, "hour"), floor_date(df$date_posix, "day") <= floor_date(lastDate, "day")) %>% 
      pull(firstDate)
    
    pe_lastDates = agenda %>% 
      filter(id %in% pe, floor_date(firstDate, "hour") <= floor_date(df$date_posix, "hour"), floor_date(df$date_posix, "day") <= floor_date(lastDate, "day")) %>% 
      pull(lastDate)
    
    allLastDates = c(allLastDates, pe_lastDates)
    allFirstDates = c(allFirstDates, pe_firstDates)
  }
  
  allLastDates = as_datetime(allLastDates)
  allFirstDates = as_datetime(allFirstDates)
  
  if (any(allLastDates < df$date_posix + df$length)) df$length = as.numeric(difftime(min(allLastDates), df$date_posix, units = "secs"))
  if (any(allFirstDates > df$date_posix)) {
    if (any(allFirstDates > df$date_posix & allFirstDates > df$date_posix+df$length)) {
      df$before_schedule = T
    } else {
      newStart = max(allFirstDates)
      df$length = as.numeric(difftime(df$date_posix+df$length, newStart, units = "secs"))
      df$date_posix = newStart
    }
  }
  
  return(df)
}

# Function to reconstruct patient location when not interacting
patient_locations = function(patient_id, patient_loc, 
                             admission = admission, 
                             begin_date = begin_date,
                             end_date = end_date,
                             patient_rooms = patient_rooms
                             ) {
  # Assign "NOT HERE" when not present in the ward
  firstDate = as_datetime(paste(admission$firstDate[admission$id==patient_id], "00:00:00"))
  firstDate = case_when(begin_date>firstDate ~ begin_date, .default = firstDate)
  lastDate = as_datetime(paste(admission$lastDate[admission$id==patient_id], "23:59:30"))
  lastDate = case_when(end_date<lastDate ~ end_date, .default = lastDate)
  allDates = seq(begin_date, end_date-30, 30)
  present = seq(firstDate, lastDate-30, 30)
  not_present = which(!allDates %in% present)
  if (length(not_present) != sum(is.na(patient_loc[not_present]))) stop(paste(patient_id, "has assigned location when individual not in the ward"))
  patient_loc[not_present] = -1
  
  # Assign patient room when present in the ward and not interacting
  patient_loc[is.na(patient_loc)] = patient_rooms$room[patient_rooms$id==patient_id]
  return(patient_loc)
}

# Function to reconstruct HCW location when not interacting
hcw_locations = function(hcw_id, hcw_loc,
                         threshold = 90*2,
                         agenda,
                         admission, 
                         begin_date,
                         end_date,
                         rooms,
                         id_medical
                         ) {
  
  # Assign "NOT HERE" when not present in the ward
  allDates = seq(begin_date, end_date-30, 30)
  present = agenda %>% 
    filter(id == hcw_id) %>%
    select(firstDate, lastDate) %>%
    mutate(n = 1:n()) %>%
    nest(.by = n) %>%
    mutate(data=map(data, unroll_30s)) %>%
    unnest(data) %>%
    pull(data)
  not_present = which(!allDates %in% present)
  if (!all(is.na(hcw_loc[not_present]))) message(paste(hcw_id, "has location assigned when not present in the ward"))
  hcw_loc[not_present] = -1

  # Split vector into elements of consecutive identical values
  loc_split = split(hcw_loc, data.table::rleid(hcw_loc))

  # Rooms
  patient_rooms_ids = unique(rooms$id_room[rooms$room == rooms$id_room])
  
  # Assign locations to elements
  for (k in seq_along(loc_split)) {

    current = loc_split[[k]]
    if (!is.na(unique(current))) next()

    before = ifelse(k == 1, "-1", unique(loc_split[[k-1]]))
    after = ifelse(k == length(loc_split), "-1", unique(loc_split[[k+1]]))

    # If between patient room and less than 5 mins
    if (before == after & length(current) <= 10 & before %in% patient_rooms_ids) {
      loc_split[[k]] = rep(before, length(current))

    } else if (length(current) == 1) {
      # If between two different rooms and only 30 seconds --> corridor
      loc_split[[k]] = c(rooms$id_room[rooms$room == "Corridor"])

    } else if (length(current) > threshold) {
      # If more than a threshold (default = 90 min)
      loc_split[[k]] = rep("-2", length(current))

      if (length(unique(before)) == 1 & all(before %in% patient_rooms_ids)) loc_split[[k]][1] = rooms$id_room[rooms$room == "Corridor"]
      if (length(unique(after)) == 1 & all(after %in% patient_rooms_ids)) loc_split[[k]][length(current)] = rooms$id_room[rooms$room == "Corridor"]

    } else {
      # If between two different rooms and/or more than 5 mins
      resting_room = ifelse(hcw_id %in% id_medical, rooms$id_room[rooms$room == "Medical Staff Room"], rooms$id_room[rooms$room == "Paramedical Staff Room"])
      loc_split[[k]] = c(rooms$id_room[rooms$room == "Corridor"], rep(resting_room, length(current)-2), rooms$id_room[rooms$room == "Corridor"])
    }
  }

  # Return
  out = unlist(loc_split)
  return(out)
}

# Functions for plots-----------------------------------------------------------
# Get summary statistics of a simulated epidemic
get_ss = function(df, n1, n2) {
  
  epidemic_curve = rep(0, 90*24*3600/30)
  for (i in df$id[df$t_inf>0]) {
    t_start = df$t_incub[df$id == i] 
    t_end = df$t_recover[df$id == i]
    epidemic_curve[t_start:t_end] = epidemic_curve[t_start:t_end] + 1
  }
  peak_time = which(epidemic_curve == max(epidemic_curve))
  peak_time = ifelse(length(peak_time)>0, min(peak_time), peak_time) 
  
  out = data.frame(
    couple = n1, 
    sim_id = n2,
    nind = length(unique(df$id)),
    epidemic_duration = max(df$t_inf)/3600,
    peak_time = peak_time / (60*2*24),
    ninf_patients_c = sum(grepl("CONTACT", df$inf_by[df$id %in% id_patient])), 
    ninf_patients_e = sum(grepl("ENVIRONMENT", df$inf_by[df$id %in% id_patient])),
    ninf_para_c = sum(grepl("CONTACT", df$inf_by[df$id %in% id_paramedical])),
    ninf_para_e = sum(grepl("ENVIRONMENT", df$inf_by[df$id %in% id_paramedical])),
    ninf_med_c = sum(grepl("CONTACT", df$inf_by[df$id %in% id_medical])),
    ninf_med_e = sum(grepl("ENVIRONMENT", df$inf_by[df$id %in% id_medical])),
    ninf_c = sum(grepl("CONTACT", df$inf_by)),
    ninf_e = sum(grepl("ENVIRONMENT", df$inf_by))
  )
  
  return(out)
}
# Functions for animated gif----------------------------------------------------
distribute_points <- function(n) {
  theta <- seq(0, 2 * pi, length.out = n + 1)[-1]
  radius <- 0.2  # Radius of the circle where the points will be plotted
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  return(data.frame(offset_x = x, offset_y = y))
}

apply_offsets <- function(group_by_data) {
  n <- nrow(group_by_data)
  offsets <- distribute_points(n)
  group_by_data <- group_by_data %>%
    mutate(offset_x = offsets$offset_x, offset_y = offsets$offset_y)
  return(group_by_data)
}