
#####################
## SAR COMPUTATION ##
#####################

##################### SAR FOR ONE SIMULATION
SAR_compute <- function(global_status) {
  ## WE NEED INDEX TYPE FOR THE SAR BY TYPE COMPUTATION 
  index_id <- global_status %>% filter(inf_by == 'INDEX') %>% pull(id)
  index_type <- case_when(
    index_id %in% id_paramedical ~ 'PARAMEDICAL',
    index_id %in% id_medical ~ 'MEDICAL',
    index_id %in% id_patient ~ 'PATIENT'
  )
  ## SUSCEPTIBLE INDIVIDUALS DURING THE STUDY
  n_individual <- length(global_status %>% distinct(id) %>% pull()) -1
  n_medical <- ifelse(index_type == 'MEDICAL', length(id_medical) - 1, length(id_medical))
  n_paramedical <- ifelse(index_type == 'PARAMEDICAL', length(id_paramedical) - 1, length(id_paramedical))
  n_patient <- ifelse(index_type == 'PATIENT', length(id_patient) - 1, length(id_patient))
  n_hcw <- ifelse(index_type == 'PATIENT', length(id_hcw) - 1, length(id_hcw) - 1)
  
  ## SAR BY TYPE
  SAR_global <- (length(global_status %>% filter(t_inf != -1) %>% pull(id)) -1 )/(n_individual - 1)
  SAR_patient <- length(global_status %>% filter(id %in% id_patient & t_inf != -1 & inf_by != 'INDEX') %>% pull(id)) / (n_patient)
  SAR_medical <- length(global_status %>% filter(id %in% id_medical & t_inf != -1 & inf_by != 'INDEX') %>% pull(id)) / n_medical
  SAR_paramedical <- length(global_status %>% filter(id %in% id_paramedical & t_inf != -1 & inf_by != 'INDEX') %>% pull(id)) / n_paramedical
  SAR_hcw <- length(global_status %>% filter(id %in% id_hcw & t_inf != -1 & inf_by != 'INDEX') %>% pull(id)) / n_hcw
  SAR_environment <- length(global_status %>% filter(id %in% id_hcw & t_inf != -1 & inf_by != 'INDEX') %>% pull(id)) / n_hcw
  ##DF
  SAR_df <- data.frame(
    Global = SAR_global,
    HCW = SAR_hcw,
    Patient = SAR_patient,
    Paramedical = SAR_paramedical,
    Medical = SAR_medical)
  
  return(SAR_df)
}

##################### GET ALL SAR FOR ONE COUPLE (N SIMULATIONS)
get_SAR_couple <- function(couple, list_sim){
  lapply(list_sim[[couple]], function(global_status) {
    SAR_compute(global_status = global_status)
  })
}

##################### GET SAR FOR ALL BETA/NU COUPLES
get_all_SAR <- function(list_sim, id_patient, id_hcw, id_paramedical, id_medical) {
  list_SAR <- list()
  for(couple in names(list_sim)){
    list_SAR[[couple]] <- get_SAR_couple(couple, list_sim)
  }
  return(list_SAR)
}


#################
## SAR METRICS ##
#################

################# SAR METRICS FOR ONE COUPLE
get_SAR_metrics <- function(couple, list_SAR){
  couple_parts <- strsplit(couple, "_")[[1]]
  SAR_df <- rbindlist(list_SAR[[couple]])
  SAR_df_metrics <- data.frame(couple = couple,
                               beta = eval(parse(text = gsub(pattern = '-', replacement = '/', x = couple_parts[2]))),
                               nu = eval(parse(text = gsub(pattern = '-', replacement = '/', x = couple_parts[3]))),
                               SAR_min = min(SAR_df$Global),
                               SAR_max = max(SAR_df$Global),
                               SAR_median = median(SAR_df$Global),
                               SAR_mean = mean(SAR_df$Global),
                               SAR_sd = sd(SAR_df$Global))
  return(SAR_df_metrics)
}

################# SAR METRICS FOR EVERY COUPLES
get_all_SAR_metrics <- function(list_SAR){
  all_SAR_metrics <- data.frame()
  for(couple in names(list_SAR)){
    all_SAR_metrics <- rbind(all_SAR_metrics, get_SAR_metrics(couple, list_SAR) )
  }
  return(all_SAR_metrics)
  
}