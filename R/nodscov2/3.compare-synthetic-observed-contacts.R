################################################################################
##
##                Compare synthetic and real contact data
##
################################################################################

## Working environment----------------------------------------------------------
# Libraries
library(tidyverse)
library(viridis)
library(igraph)
library(ggnetwork)
library(ggtext)
library(ggpubr)
rm(list = ls())

# Helper functions
source("helper-functions.R")

# Working directories and paths
synthetic_path = file.path("..", "..", "data", "data-synthetic-graphs")
fig_path = file.path("..", "..", "fig", "network_comparison")
if (!dir.exists(fig_path)) dir.create(fig_path)

## Load real and synthetic contact data-----------------------------------------
networks = c("poincare", "herriot")
full_network = F
admission = list()
schedule = list()
interaction_real = list()
interaction_synthetic = list()

for (net in networks) {
  # List of individuals 
  admission[[net]] = read.csv2(file.path(synthetic_path, "input", paste0("admission_", net, ".csv"))) %>%
    mutate(
      firstDate = as_date(firstDate),
      lastDate = as_date(lastDate)
    )
  
  # Schedule of healthcare workers
  schedule[[net]] = read.csv2(file.path(synthetic_path, "input", paste0("agenda_", net, ".csv"))) %>%
    mutate(
      firstDate = as_datetime(firstDate),
      lastDate = as_datetime(lastDate)
    )
  
  # Real interaction data
  interaction_real[[net]] = read.csv2(file.path(synthetic_path, "input", paste0("interactions_", net, ".csv"))) %>%
    mutate(date_posix = as_datetime(date_posix))
  
  # Synthetic data
  interaction_synthetic[[net]] = data.frame()
  reconstructed_path = file.path(synthetic_path, "biased", net)
  if (full_network) reconstructed_path = file.path(synthetic_path, "full", net)
  
  for (f in list.files(reconstructed_path, pattern = ".*csv", full.names = T)) {
    interaction_synthetic[[net]] = bind_rows(
      interaction_synthetic[[net]],
      read.csv2(f) %>% 
        mutate(
          nSim = gsub("^.*Networks|_oneday.*$", "", f),
          firstDate = as_datetime(date_posix),
          lastDate = as_datetime(date_posix) + length
        )
    )
  } 
}

## Verify study periods for full networks---------------------------------------
if (full_network) {
  study_period = 90 # days
  
  # Herriot
  c(min(admission$herriot$firstDate), max(admission$herriot$lastDate))
  c(floor_date(min(as_datetime(interaction_synthetic$herriot$date_posix)), "day"), floor_date(max(as_datetime(interaction_synthetic$herriot$date_posix)+interaction_synthetic$herriot$length), "day")) 
  
  # Poincaré
  c(min(admission$poincare$firstDate), max(admission$poincare$lastDate))
  c(floor_date(min(as_datetime(interaction_synthetic$poincare$date_posix), na.rm = T), "day"), floor_date(max(as_datetime(interaction_synthetic$poincare$date_posix), na.rm = T), "day"))  
}

## Compare network characteristics----------------------------------------------
# Analyses are based on the codes from Leclerc et al. (2024) that are available
# at : https://gitlab.pasteur.fr/qleclerc/network_algorithm/-/tree/main/analysis?ref_type=heads

### Number of unique contacts and contact duration
# Extract number of unique contacts and contact duration from both networks (synthetic + observed)
durations = data.frame()
numbers = data.frame()

for (net in networks) {
  
  # Prepare observed data 
  durations = bind_rows(
    durations, 
    interaction_real[[net]] %>%
      select(from, to, date_posix, length) %>%
      mutate(date_posix = floor_date(date_posix, "hour")) %>%
      group_by(from, to, date_posix) %>%
      summarise(length = sum(length)/60, .groups = "drop") %>% # Convert in minutes
      select(length) %>%
      mutate(data = "Observed", network = net)
  )
  
  numbers = bind_rows(
    numbers,
    contact_numbers(interaction_real[[net]], "Observed", net)
  )
  
  # Prepare reconstructed data 
  reconstructed_path = file.path(synthetic_path, "biased", net)
  if (full_network) reconstructed_path = file.path(synthetic_path, "full", net)
  
  for (f in list.files(reconstructed_path, ".*.csv", full.names = T)) {
    
    data = read.csv2(f)
    
    # Contact durations
    db_type = gsub(".*Networks|_oneday.*", "", f)
    db_type = paste("Synthetic", db_type)
    
    durations = bind_rows(
      durations, 
      data %>%
        mutate(date_posix = floor_date(as_datetime(date_posix), "hour")) %>%
        group_by(from, to, date_posix) %>%
        summarise(length = sum(length)/60, .groups = "drop") %>%
        mutate(data = db_type, network = net) %>%
        select(length, network, data)
    )
    
    numbers = bind_rows(
      numbers,
      contact_numbers(data, db_type, net)
    )
  }
}

# Change level order
numbers = numbers %>%
  mutate(
    network = ifelse(network == "poincare", "Poincaré", "Herriot"),
    data = factor(data, c("Observed", paste("Synthetic", 1:10))),
    type = factor(type, levels = unique(type)),
    date_posix = factor(date_posix, levels = c("0:00","1:00","2:00","3:00","4:00",
                                               "5:00","6:00","7:00","8:00","9:00",
                                               "10:00","11:00","12:00","13:00","14:00",
                                               "15:00","16:00","17:00","18:00","19:00",
                                               "20:00","21:00","22:00","23:00"))
  )

# Plot number of unique contacts
pa = numbers %>%
  ggplot(., aes(x = date_posix, y = med, ymin = q25, ymax = q75, fill = data, col = data, group = interaction(data, type, network))) +
  geom_ribbon(alpha = 0.3) +
  geom_point() +
  geom_line() +
  facet_grid(cols = vars(type), rows = vars(network)) +
  scale_colour_discrete(type = c("orange", colorRampPalette(c("dodgerblue4", "lightskyblue"))(10))) +
  scale_fill_discrete(type = c("orange", colorRampPalette(c("dodgerblue4", "lightskyblue"))(10))) +
  scale_y_continuous(breaks = seq(0,300,75)) +
  scale_x_discrete(breaks = c("0:00", "4:00", "8:00", "12:00",
                              "16:00","20:00")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))  +
  guides(fill = FALSE, color = guide_legend(override.aes = list(fill = NA))) +
  labs(x = "Hour", y = "Number of unique contacts", colour = "")

# Plot contact duration 
pb = durations %>%
  mutate(network = ifelse(network == "poincare", "Poincaré", "Herriot"),
         data = factor(data, c("Observed", paste("Synthetic", 1:10)))) %>%
  ggplot(., aes(y=length, x = network, colour = data)) +
  geom_boxplot(outlier.shape = NA) +
  scale_colour_discrete(type = c("orange", colorRampPalette(c("dodgerblue4", "lightskyblue"))(10))) +
  coord_cartesian(ylim = c(0,35)) +
  # guides(colour="none") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  labs(y = "Contact duration (minutes)", x = "Network", col = "") 

# Combine plots
p = ggarrange(pa, pb, nrow = 2, labels = c("A", "B"), 
              heights = c(1, 0.8), common.legend = T, legend= "right")
p

if (full_network) ggsave(file.path(fig_path, "unique_contacts_full.png"), p, height = 8, width = 11)  
if (!full_network) ggsave(file.path(fig_path, "unique_contacts_biased.png"), p, height = 8, width = 11)  


### Network metrics (degree, global efficiency, transitivity, assortativity, temporal correlation, density)
# Get summary metrics on observed and reconstructed data as well as two ICU networks
summary_data = data.frame()

for (net in networks) {
  
  # Observed data 
  adm_data = admission[[net]] %>%
    mutate(cat = ifelse(cat == "Patient", cat, "Staff")) %>%
    select(id, cat, ward) %>%
    distinct()
  
  graph_data = interaction_real[[net]] %>%
    select(from, to, date_posix, length) %>%
    rename(date_posix = date_posix) %>%
    mutate(date_posix = floor_date(date_posix, "day")) %>%
    select(from, to, date_posix) %>%
    distinct() %>%
    arrange(date_posix)
  summary_data = bind_rows(
    summary_data, 
    get_net_metrics(graph_data, adm_data, 1, "Observed", net)
  )
  
  included_days = unique(graph_data$date_posix)
  
  # Synthetic networks
  reconstructed_path = file.path(synthetic_path, "biased", net)
  if (full_network) reconstructed_path = file.path(synthetic_path, "full", net)
  simu_files = list.files(path = reconstructed_path, pattern = ".*.csv", full.names = T)
  simu_data = data.frame()
  
  for(f in simu_files){
    iter = as.numeric(gsub(".*Networks|_oneday.*", "", f))
    graph_data = read.csv2(f) %>%
      mutate(date_posix = floor_date(as_datetime(date_posix), "day")) %>%
      filter(date_posix %in% included_days) %>%
      select(from, to, date_posix) %>%
      distinct %>%
      arrange(date_posix)
    
    summary_data = bind_rows(
      summary_data,
      get_net_metrics(graph_data, adm_data, iter, "Synthetic", net)
    )
  }
}

# Get median
summary_data = summary_data %>%
  filter(iter == 1) %>%
  mutate(
    data = factor(data, levels = c("Observed", "Synthetic")),
    network = ifelse(network == "herriot", "Herriot", "Poincaré")
  ) %>%
  group_by(day, data, network) %>%
  summarise(across(everything(), median), .groups = "drop") 

# Plots
dict_metric = c("assortativities" = "Assortativities", 
                "degrees" = "Degrees",
                "densities" = "Densities",
                "efficiencies" = "Efficiencies",
                "temp_corr" = "Temporal correlation",
                "transitivities" = "Transitivities"
)

p = summary_data %>%
  pivot_longer(-c(iter, data, network, day), names_to = "metric", values_to = "Value") %>%
  mutate(metric = recode(metric, !!!dict_metric)) %>%
  ggplot(., aes(x = network, y = Value, col = data)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5)) + 
  facet_wrap(facets = vars(metric), ncol = 3, scales = "free_y") +
  theme_bw() +
  labs(col = "Data", x = "")

if (full_network) ggsave(file.path(fig_path, "metrics_full.png"), p, height = 5, width = 10)
if (!full_network) ggsave(file.path(fig_path, "metrics_biased.png"), p, height = 5, width = 10)


### Troubleshooting assortativity by degree-------------------------------------
# adm_data = admission %>%
#   select(id, cat, ward) %>%
#   distinct()
# 
# ## Reconstructed full #####
# 
# simu_files = list.files(path = file.path(synthetic_path, "full_0-5agenda"), pattern = "matContact.*", full.names = T)
# simu_data = data.frame()
# iter=1
# 
# for(f in simu_files){
#   data = read.csv2(f)
#   
#   graph_data = data %>%
#     mutate(date_posix = floor_date(as_datetime(date_posix), "day")) %>%
#     select(from, to, date_posix) %>%
#     distinct %>%
#     arrange(date_posix)
#     
#   simu_data = rbind(simu_data,
#                     get_assortativity_degree(graph_data, adm_data, iter, "Reconstructed"))
#   
#   iter=iter+1
# } 
# 
# ## Observed #####
# graph_data = interaction_real %>%
#   rename(date_posix = date_posix) %>%
#   mutate(date_posix = floor_date(date_posix, "day")) %>%
#   select(from, to, date_posix) %>%
#   distinct() %>%
#   arrange(date_posix) 
# 
# observed_data = get_assortativity_degree(graph_data, adm_data, 1, "Observed")
# 
# ## Plot simulated assortativity
# ggplot(simu_data, aes(x = day, y = assortativities, col = as.factor(iter), group = as.factor(iter))) +
#   geom_line() +
#   geom_point(x = observed_data$day[1], y = observed_data$assortativities[1], col = "red") +
#   geom_point(x = observed_data$day[2], y = observed_data$assortativities[2], col = "red") +
#   scale_color_manual(values = colorRampPalette(c("darkblue", "lightblue"))(10)) +
#   theme_bw() +
#   expand_limits(y=c(min(simu_data$assortativities), max(observed_data$assortativities))) +
#   labs(col = "Simulation", y = "Assortativity degree (by degree)", x = "")


## Verifications on the synthetic networks
### Verify that all individuals in the admission file have at least one interaction
# data_real = interaction_real %>%
#   rename(date_posix = date_posix)
# ids_real = sort(unique(c(data_real$from, data_real$to)))
# 
# data_syn = read.csv2(file.path(synthetic_path, "matContactBuiltSimulatedCtcNetworks1_oneday.csv")) %>%
#   mutate(date_posix = floor_date(as_datetime(date_posix), "hour"))
# ids_syn = sort(unique(c(data_syn$from, data_syn$to)))
# 
# all(ids_syn %in% ids_real)
# sort(ids_syn[!ids_syn %in% ids_real & !grepl("-[0-9]$", ids_syn)])
# sort(admission$id[!admission$id %in% ids_real & !grepl("-[0-9]+$", admission$id)])

### Verify whether individuals in the synthetic network have at least one interaction per hour 
# # Do all individuals have at least one contact per hour in the synthetic data?
# at_least_one_contact = data.frame()
# all_contacts = data.frame()
# 
# for (i in admission$id) {
#   
#   # Admission date    
#   admissionDate = admission$firstDate[admission$id == i]
#   proceed = F
#   
#   if (admission$status[admission$id == i] == "PE") {
#     # Vector of hours of presence in the ward
#     proceed = T
#     # lastHour = as_datetime("2020-08-03 23:00:00")
#     hoursOfPresence = schedule_real %>%
#       filter(id == i) %>%
#       select(firstDate, lastDate) %>%
#       mutate(r = 1:n()) %>%
#       group_by(r) %>%
#       nest() %>%
#       mutate(seqHours = map(data, unroll_dates)) %>%
#       unnest(seqHours) %>%
#       ungroup() %>%
#       distinct(seqHours) %>%
#       pull(seqHours)
#   } 
#   
#   if (admission$status[admission$id == i] == "PA") {
#     proceed = T
#     # Patient discharge date
#     dischargeDate = admission$lastDate[admission$id == i]
# 
#     # Hours of presence in the ward
#     firstHour = as_datetime(paste0(admissionDate, "08:00:00"))
#     lastHour = as_datetime(paste0(dischargeDate, "20:00:00"))
#     hoursOfPresence = seq(firstHour, lastHour, 3600)
#   }
#   
#   if (proceed) {
#     # Hours with contacts in the ward 
#     all_contacts_temp = interaction_synthetic %>%
#       filter(from == i | to == i) %>%
#       mutate(date_posix = floor_date(as_datetime(date_posix), "hour")) %>%
#       group_by(date_posix, nSim) %>%
#       summarise(n = n(), .groups = "drop")
#     
#     hoursInteracting = all_contacts_temp %>%
#       group_by(nSim) %>%
#       summarise(atLeastOne = all(hoursOfPresence %in% date_posix), .groups = "drop") %>% 
#       mutate(id = i)
#     
#     all_contacts = bind_rows(all_contacts,
#                              expand.grid(nSim = as.character(1:10), date_posix = hoursOfPresence, id = i) %>%
#                                left_join(., all_contacts_temp %>% select(date_posix, nSim) %>% mutate(contact = 1), by = c("nSim", "date_posix")) %>%
#                                mutate(contact = ifelse(is.na(contact), 0, contact), moment = ifelse(hour(date_posix) < 12, "morning", "afternoon")) %>%
#                                group_by(id, nSim, moment) %>%
#                                summarise(contact = sum(contact==0), n_tot = n(), .groups = "drop")
#                                )
#     
#     at_least_one_contact = bind_rows(
#       at_least_one_contact,
#       hoursInteracting
#     ) 
#   }
# }
# 
# # Check 
# at_least_one_contact %>%
#   mutate(status = gsub("-.*$", "", id)) %>%
#   count(status, atLeastOne)
# 
# # Individuals without interactions
# inds = at_least_one_contact %>%
#   group_by(id) %>%
#   summarise(n=sum(atLeastOne), .groups = "drop") %>%
#   filter(n == 10) %>%
#   pull(id)
# schedule_real %>% filter(id %in% inds)
# 
# inds_to_test = schedule_real %>% 
#   filter(lastDate <= as_datetime("2020-05-07 12:00:00")) %>% 
#   group_by(id) %>% 
#   nest() %>% 
#   mutate(firstDayInteraction = map(data, testFun)) %>% 
#   select(-data) %>% 
#   unnest(firstDayInteraction) %>% 
#   filter(firstDayInteraction == F)
#   
# # Individuals with at least one contact 
# ids = at_least_one_contact %>%
#   filter(atLeastOne == T) %>%
#   distinct(id) %>%
#   pull(id)
# admission %>% filter(id %in% ids)
# schedule_real %>% filter(id %in% ids) %>% filter(lastDate <= "2020-05-09 00:00:00") %>% arrange(id, firstDate)
# interaction_real %>% filter(from %in% ids | to %in% ids, date_posix <= as_datetime("2020-05-07 12:00:00")) %>%  
#   mutate(lastDate = date_posix + length) 
# 
# # Individuals that are present only during the first afternoon 
# testFun = function(df) {
#   presentOnlyFirstDay = all(df$firstDate <= as_datetime("2020-05-07 00:00:00") & df$lastDate <= as_datetime("2020-05-07 00:00:00"))
#   return(presentOnlyFirstDay)
# }
# 
# presentFirstDayOnly = schedule_real %>% 
#   filter(lastDate <= as_datetime('2020-05-07 12:00:00')) %>%
#   group_by(id) %>%
#   nest() %>%
#   mutate(out = map(data, testFun)) %>% 
#   select(-data) %>% 
#   unnest(out) %>% 
#   ungroup() 
# count(presentFirstDayOnly, out)
# 
# presentFirstDayOnly$id[presentFirstDayOnly$out==T]
# 
# # Missing contacts
# all_contacts %>% 
#   ggplot(., aes(x = nSim, y = contact/n_tot, col = moment)) +
#   geom_point(position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   labs(x = "Simulation", y = "Proportion of hours with contact", col = "")

### Verification of the coherence between schedule and admission databases for healthcare workers
# # Schedule of healthcare workers 
# ids_pe = admission$id[admission$status == "PE"]
# interactions_in_schedule = data.frame()
# 
# for (i in ids_pe) {
#   # Schedule of the HCW
#   schedule_id = schedule_real %>%
#     filter(id == i, lastDate <= as_datetime("2020-05-07 12:00:00")) %>%
#     select(firstDate, lastDate)
#   
#   # Get all interactions
#   if (sum(interaction_real$to == i | interaction_real$from == i) > 0) {
#     interactions_id = interaction_real %>%
#       filter(from == i| to == i) %>%
#       select(date_posix, length) %>%
#       mutate(n = 1:n())
#   
#     # Are all interactions within the schedule ?
#     out_vec = expand.grid.df(
#       schedule_id, 
#       interactions_id
#     ) %>%
#       group_by(n) %>%
#       summarise(out = any(date_posix >= firstDate & date_posix+length <= lastDate), .groups = "drop") %>%
#       pull(out)
#     
#     interactions_in_schedule = bind_rows(
#       interactions_in_schedule, 
#       data.frame(
#         id = i, 
#         coherent_interactions = all(out_vec)
#       )
#     )
#   
#   } else {
#     interactions_in_schedule = bind_rows(
#       interactions_in_schedule, 
#       data.frame(
#         id = i, 
#         coherent_interactions = NA
#       )
#     )
#   }
# }
# 
# summary(interactions_in_schedule)


## Compute recurring contact probability 
# # Daily individual probability of recurring contacts and mean individual
# # probability
# simu_files = list.files(path = file.path(synthetic_path, "full_0"), pattern = "matContact.*", full.names = T)
# p_df = data.frame()
# 
# for (f in simu_files) {
#   pind = rep(NA, nrow(admission))
#   
#   for (i in seq_along(admission$id)) {
#     ii = admission$id[i]
#     
#     # Get n days spent in the ward with a sensor 
#     if (grepl("^PA-", ii)) {
#       all_d = admission %>%
#         filter(id == ii) %>%
#         nest() %>%
#         mutate(out = map(data, unroll_days)) %>%
#         unnest(out) %>%
#         .$out
#     }
#     
#     if (grepl("^PE-", ii)) {
#       all_d = schedule_real %>%
#         filter(id == ii) %>%
#         mutate(firstDate = floor_date(firstDate, "day"), lastDate = floor_date(lastDate, "day"), k = 1:n()) %>%
#         arrange(firstDate) %>%
#         group_by(k) %>%
#         nest() %>%
#         mutate(out = map(data, unroll_dates)) %>%
#         unnest(out) %>%
#         select(out) %>%
#         distinct(out) %>%
#         .$out
#     }
#       
#     if (length(all_d) > 1) {
#       interaction_d = vector("list", length(all_d))
#       proba_d = rep(NA, length(all_d)-1)
#       
#       # Get lists of contacts per hour
#       for (d in seq_along(all_d)) {
#         dd = all_d[d]
#         contacts = read.csv2(f) %>%
#           filter(from == ii | to == ii, date_posix >= dd, date_posix < dd+3600*24) %>%
#           select(from, to) %>%
#           unlist()
#         contacts = unique(contacts[contacts != ii])
#         names(contacts) = NULL
#         interaction_d[[d]] = contacts 
#       }
#       
#       # Get hourly probability of recurring contacts
#       if (length(all_d)>1) {
#         for (h in 2:length(all_d)) {
#           if (length(interaction_d[[d]]) > 0) {
#             proba_d[d-1] = sum(interaction_d[[d]] %in% unlist(interaction_d[1:(d-1)])) / length(interaction_d[[d]])
#           } else {
#             proba_d[d-1] = 0
#           }
#         }
#       }
#       
#       # Mean individual hourly probability of recurring contacts
#       pind[i] = sum(proba_d) / length(proba_d)
#       
#     } else {
#       pind[i] = 0
#     }
#   }
#   
#   p_df = bind_rows(p_df, data.frame(ids = admission$id, nSim = gsub("^.*Networks|_oneday.*$", "", f), pind = pind))
# }
# 
# 
# # Plot mean individual probability of recurring contacts
# admission_nodscov2 %>%
#   mutate(pind = pind) %>%
#   filter(!is.na(pind)) %>%
#   ggplot(., aes(x = cat, y = pind)) +
#   geom_boxplot() +
#   geom_jitter() +
#   theme_bw() +
#   labs(x = "", y = "Probability of daily recurring contact")
# 
# # Mean probability of recurring contacts by individual category
# admission_nodscov2 %>%
#   mutate(pind = pind) %>%
#   group_by(cat) %>%
#   summarise(m = mean(pind, na.rm = T), .groups = "drop")
# 
# # Comparison by category of individual
# admission_nodscov2 %>%
#   mutate(pind = pind) %>%
#   filter(!is.na(pind)) %>%
#   rstatix::wilcox_test(pind ~ cat, ref.group = "Paramedical")
# 
# # Mean probability when aggregating paramedical and medical staff
# admission_nodscov2 %>%
#   mutate(pind = pind, cat = ifelse(cat == "Patient", cat, "HCW")) %>%
#   group_by(cat) %>%
#   summarise(m = mean(pind, na.rm = T), .groups = "drop")
