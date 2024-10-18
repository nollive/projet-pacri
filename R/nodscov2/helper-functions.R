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
             "physician" = "Medical")

# Plots-------------------------------------------------------------------------
pal = c("Observed" = "darkorchid", "Reconstructed"= "darkorange")

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

# Get vector of days between two dates
unroll_days = function(df) {
  out = seq.Date(df$firstDate, df$lastDate, 1)
  return(out)
}

# In a dataframe, cocnatenate rows that correspond to consecutive time periods
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
      
      if ("date_posix_first" %in% colnames(df)) {
        overlap = df$date_posix_first[i] >= df$date_posix_first & df$date_posix_first[i] <= df$date_posix_first+df$length
      }      
      if ("date_posix" %in% colnames(df)) {
        overlap = df$date_posix[i] >= df$date_posix & df$date_posix[i] <= df$date_posix+df$length
      }
      
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
      
      overlap = df$date_posix_first[i] >= df$date_posix_first & df$date_posix_first[i] <= df$date_posix_first+df$length
      
      if (sum(overlap) > 1) {
        # Contacts with overlap
        j = which(overlap)
        j = j[j!=i]
        
        # Start of interaction
        new_date_posix = min(df$date_posix_first[c(i,j)])
        
        # End of interaction
        new_end = max(df$date_posix_first[c(i,j)] + df$length[c(i,j)])
        new_length = as.numeric(difftime(new_end, new_date_posix, units = "sec"))
        
        # Replace date_posix and length
        df$date_posix_first[c(i,j)] = new_date_posix
        df$length[c(i,j)] = new_length
      }
    }
  }
  
  df = distinct(df)
  
  return(df)
  
}

# Functions to analyze networks-------------------------------------------------
# Function to return a set of network metrics (from Leclerc et al., 2024)
# https://gitlab.pasteur.fr/qleclerc/network_algorithm/-/tree/main?ref_type=heads
get_net_metrics = function(graph_data, adm_data, iter = 0, network = "Observed"){
  
  days = unique(graph_data$date_posix)
  
  data = data.frame(degrees = rep(0, length(days)),
                    densities = rep(0, length(days)),
                    transitivities = rep(0, length(days)),
                    assortativities = rep(0, length(days)),
                    efficiencies = rep(0, length(days)),
                    temp_corr = rep(0, length(days)),
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
      set_vertex_attr("cat", value = vertex_atts$cat) %>%
      set_vertex_attr("cat_ag", value = vertex_atts$cat_ag) %>%
      set_vertex_attr("staff", value = vertex_atts$staff) %>%
      set_vertex_attr("ward", value = vertex_atts$ward)
    
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
get_assortativity_degree = function(graph_data, adm_data, iter = 0, network = "Observed", 
                                    time_window){
  
  days = unique(graph_data$date_posix)
  
  data = data.frame(assortativities = rep(0, length(days)),
                    iter = iter,
                    network = network,
                    day = days,
                    time_window = time_window)
  
  
  for(i in 1:length(days)){
    
    # full network
    data_d = graph_data %>%
      filter(date_posix == days[i])
    
    graph_d = graph_from_data_frame(data_d, directed = F)
    graph_d = simplify(graph_d)
    vertex_atts = data.frame(id = vertex_attr(graph_d, "name")) %>%
      left_join(adm_data, "id")
    graph_d = graph_d %>%
      set_vertex_attr("cat", value = vertex_atts$cat) %>%
      set_vertex_attr("cat_ag", value = vertex_atts$cat_ag) %>%
      set_vertex_attr("staff", value = vertex_atts$staff) %>%
      set_vertex_attr("ward", value = vertex_atts$ward)
    
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



