################################################################################
##          Support functions, dictionaries and variables
##                   to analyze interaction data
################################################################################

# Dictionaries------------------------------------------------------------------
dict_scenarios = c("sim_1-4_20" = "Scenario 1",
                   "sim_1-2_16" = "Scenario 2",
                   "sim_3-4_12" = "Scenario 3",
                   "sim_1_8" = "Scenario 4",
                   "sim_3-2_5" = "Scenario 5"
)
pal = c('Medical' = "#5CD6D6", 'Paramedical' = "#A9B9E8", 'Patient' = "#FFA766", 'Room' = "#666699")
env_pal = c('Contact' = 'darkorange', 'Environment' =  'orchid')

# Variables---------------------------------------------------------------------
cat_paramedical <- c("nurse", "student nurse", "aux nurse")
cat_medical <- c("physician", "ext physician", "reeducation staff", "other")

# Functions---------------------------------------------------------------------
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

