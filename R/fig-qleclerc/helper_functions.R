
# function to return a set of network metrics
get_net_metrics = function(graph_data, iter = 0, network = "Observed"){
  
  days = unique(graph_data$date_posix)
  
  data = data.frame(degrees = rep(0, length(days)),
                    densities = rep(0, length(days)),
                    transitivities = rep(0, length(days)),
                    assortativities = rep(0, length(days)),
                    assortativities_ward = rep(0, length(days)),
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
    vertex_atts = data.frame(id = get.vertex.attribute(graph_d, "name")) %>%
      left_join(adm_data, "id") %>%
      mutate(ward = replace(ward, ward == "Menard 1", "Neurologic (1)"),
             ward = replace(ward, ward == "Menard 2", "Neurologic (2)"),
             ward = replace(ward, ward == "Sorrel 0", "Nutrition"),
             ward = replace(ward, ward == "Sorrel 1", "Neurologic (3)"),
             ward = replace(ward, ward == "Sorrel 2", "Geriatric"),
             ward = replace(ward, ward == "Other", "Mobile"))
    graph_d = graph_d %>%
      set_vertex_attr("cat", value = vertex_atts$cat) %>%
      set_vertex_attr("cat_ag", value = vertex_atts$cat_ag) %>%
      set_vertex_attr("staff", value = vertex_atts$staff) %>%
      set_vertex_attr("ward", value = vertex_atts$ward)
    
    data$degrees[i] = mean(degree(graph_d))
    data$densities[i] = edge_density(graph_d)
    data$transitivities[i] = transitivity(graph_d)
    data$assortativities[i] = assortativity.degree(graph_d, directed = F)
    data$assortativities_ward[i] = assortativity.nominal(graph_d, as.factor(V(graph_d)$ward), directed = F)
    data$efficiencies[i] = global_efficiency(graph_d, directed = F)
    
  }
  
  data$temp_corr = c(0, temporal_correlation(graph_data))
  
  return(data)
  
}

# function to calculate the temporal correlation of a network
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
