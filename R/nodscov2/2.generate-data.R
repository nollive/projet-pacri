################################################################################
##
##              GENERATE DATA FOR SYNTHETIC CONTACT GENERATION 
##
################################################################################

## Working environment----------------------------------------------------------
# Libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(tidyr)
library(viridis)
library(purrr)
library(igraph)
library(ggtext)
library(ggnetwork)
rm(list = ls())

# Source functions, dictionaries and helpful variables
source("R/nodscov2/helper-functions.R")

# Hospital of interest 
net = "herriot"
if (net == "poincare") {
  hosp = "APHP - RAYMOND POINCARÉ"
  hosp_comp = "raymond_poincare"
  w = "Medical ICU #1"
}
if (net == "herriot") {
  hosp = "HC DE LYON - EDOUARD HERRIOT"
  hosp_comp = "edouard_herriot"
  w = "Medical ICU #2"
}

## Load Nods-CoV-2 data (not shared)--------------------------------------------
# Admission dates
admission = read.csv2("data/data-nodscov2/raw/cleanweb/NODS-COV2_BASE_FINALE_20201020_TDT_MOD20201105_utf8.csv") %>%
  select(M1F1_REFPARTICIPANT, M1F1_NOMINVESTIG, M1F1_TYPEPARTICIPANT) %>%
  rename(id = M1F1_REFPARTICIPANT, cat = M1F1_TYPEPARTICIPANT) %>%
  mutate(
    hospital = recode(M1F1_NOMINVESTIG, !!!dict_hosp),
    status = case_when(cat == 9 ~ "PA", cat %in% c(1,2,3,4,5,12) ~ "PE", .default = "OT"),
    ward = recode(M1F1_NOMINVESTIG, !!!dict_ward),
    sex = NA, 
    hospitalization = ifelse(NA),
    cat = recode(cat, !!!dict_cat_num)
    ) %>%
  select(-M1F1_NOMINVESTIG)

# Interactions 
list_ward = read.csv(paste0("data/data-nodscov2/raw/final_sym/", hosp_comp, "_rea/", hosp_comp, "_symetrise_finale_63-75.csv")) %>%
  rename(from = Nom1, to = Nom2, date_posix = temps.début, date_posix_end = temps.fin, length = durée) %>%
  # Round dates
  mutate(
    date_posix_round = floor_date(as_datetime(date_posix), "10 seconds"),
    date_posix_end_round = floor_date(as_datetime(date_posix_end), "10 seconds")
  ) %>%
  mutate(length = as.numeric(difftime(date_posix_end_round, date_posix_round, units = "secs"))) %>%
  select(from, to, date_posix_round, length) %>%
  rename(date_posix = date_posix_round)

# Sensor distribution 
sensor = read.csv(paste0("data/data-nodscov2/raw/final_sym/", hosp_comp, "_rea/calendar-", hosp_comp, "-Reanimation.csv"), 
                  header = F)
colnames(sensor) = c("id_sensor", "id", "x1", "DATEREMISE", "DATEREC", "length", "x2")
sensor = sensor %>%
  select(id, DATEREMISE, DATEREC) %>%
  mutate(DATEREMISE = as_datetime(DATEREMISE), DATEREC = as_datetime(DATEREC))

## Change individual ids to create "anonymized" data----------------------------
# Create index correspondance table 
id_correspondance = admission %>%
  arrange(hospital, ward, status, id) %>%
  group_by(status) %>%
  mutate(new_id = paste0(status, "-", sprintf("%05d", 1:n()))) %>%
  ungroup() %>%
  select(id, new_id)
identical(nrow(id_correspondance), length(unique(id_correspondance$new_id)))
write.csv2(id_correspondance, "data/data-nodscov2/raw/final_sym/id_correspondance.csv", quote = F, row.names = F)

# Change id in the 3 dataframes 
# Admission data
admission = admission %>% left_join(., id_correspondance, by = "id") %>% relocate(new_id, .before = id) %>% select(-id) %>% rename(id = new_id)

# Sensor wearing data
sensor = sensor %>% left_join(., id_correspondance, by = "id") %>% relocate(new_id, .before = id) %>% select(-id) %>% rename(id = new_id)

# Interaction data
list_ward = list_ward %>% left_join(., id_correspondance, by = c("from"="id")) %>% select(-from) %>% rename(from = new_id) %>%
  left_join(., id_correspondance, by = c("to"="id")) %>% select(-to) %>% rename(to = new_id)

# Raw categories
raw_cat = admission %>% select(id, cat) 

## Initial composition of adult ICU---------------------------------------------
# Number of participants
# Investigation staff are excluded as they were part of the staff deploying the sensors
# and do not belong to the medical wards per se
n1 = admission %>%
  filter(hospital == hosp, ward == "Reanimation", !cat %in% "investigation") %>%
  distinct(id) %>%
  nrow()
print(n1)

# Number of participants by individual category
n2 = admission %>%
  filter(hospital == hosp, ward == "Reanimation", !cat %in% "investigation") %>%
  group_by(status, cat) %>%
  summarize(n = n(), .groups = "drop")
print(n2)

## Create input dataframes for the synthetic network algorithm------------------
# Create dataframe with admission dates for the synthetic network algorithm
admission_nodscov2 = admission %>%
  filter(
    hospital == hosp, 
    ward == "Reanimation",
    cat %in% names(dict_cat),
    !id %in% id_correspondance$new_id[id_correspondance$id %in% c("005-0114-L-G", "005-0079-P-P")] # Sensors losts
  ) %>%
  mutate(cat = recode(cat, !!!dict_cat)) %>%
  select(-hospital)
nrow(admission_nodscov2) # Number of individuals in the database after filtering

# Sensor distribution data 
# Keep individuals from the ICU of interest
sensor_nodscov2 = sensor %>% filter(id %in% admission_nodscov2$id)

# Individuals with no sensor during the study period to remove
length(unique(admission_nodscov2$id)) - length(unique(sensor_nodscov2$id)) 
admission_nodscov2 %>%
  filter(!id %in% sensor_nodscov2$id)
admission_nodscov2 = admission_nodscov2 %>% filter(id %in% sensor_nodscov2$id)

# Number of individuals in the database after filtering individuals without sensor
nrow(admission_nodscov2)
admission_nodscov2 %>% count(cat) 
raw_cat %>% filter(id %in% admission_nodscov2$id) %>% count(cat)

## Explore interaction data-----------------------------------------------------
# Get number of pairs that are not stored with the correct ids 
list_ward %>%
  filter(from %in% admission_nodscov2$id, to %in% admission_nodscov2$id) %>%
  mutate(id = 1:n()) %>%
  nest(.by = id) %>%
  mutate(data = map(data, get_pair_id)) %>%
  unnest(data) %>%
  distinct(from, to, pair) %>%
  group_by(pair) %>%
  mutate(n = n()) %>%
  filter(n>1) %>%
  nrow()

# Get number of overlapping contacts (duplicates and starting at the same time are not included)
list_ward %>%
  filter(from %in% admission_nodscov2$id, to %in% admission_nodscov2$id) %>%
  nest(.by = c(from, to)) %>%
  mutate(data = map(data, get_overlapping_contacts)) %>%
  unnest(data) %>%
  filter(overlap == T) %>%
  nrow()

## Create cleaned admission, schedule and interaction data----------------------
# Interaction data from adult ICU
interaction_nodscov2 = list_ward %>%
  filter(from %in% admission_nodscov2$id, to %in% admission_nodscov2$id) %>%
  select(from, to, date_posix, length)

# Verify whether some interactions are out of the schedule
interaction_nodscov2 %>%
  mutate(n = 1:n()) %>%
  nest(.by = n) %>%
  mutate(data = map(data, outside_schedule, sensor_nodscov2 %>% rename(firstDate_sensor = DATEREMISE, lastDate_sensor = DATEREC))) %>%
  unnest(data) %>%
  filter(in_schedule == F) %>%
  nrow()

# Use sensor distribution data for the first and last date in admission
admission_nodscov2 = admission_nodscov2 %>%
  left_join(., sensor_nodscov2 %>% 
               rename(firstDate = DATEREMISE, lastDate = DATEREC) %>%
               group_by(id) %>%
               summarise(firstDate = as_date(floor_date(min(firstDate), "day")), 
                         lastDate = as_date(floor_date(max(lastDate), "day")), 
                         .groups = "drop"), 
             by = "id") %>%
  select(id, status, firstDate, lastDate, ward, sex, hospitalization, cat)

# Verify that we have the same participants in all three databases
identical(sort(unique(admission_nodscov2$id)), sort(unique(sensor_nodscov2$id)))
all(c(interaction_nodscov2$from, interaction_nodscov2$to) %in% admission_nodscov2$id)

# Create agenda of healthcare workers
agenda_nodscov2 = sensor_nodscov2 %>%
  filter(grepl("^PE-", id)) %>%
  arrange(DATEREMISE, DATEREC)

# Individuals with no interaction when wearing the sensor 
admission_nodscov2 %>% 
  filter(!id %in% c(interaction_nodscov2$from, interaction_nodscov2$to)) %>%
  group_by(status, cat) %>%
  summarise(n = n(), .groups = "drop")

admission_nodscov2 %>% filter(!id %in% c(interaction_nodscov2$from, interaction_nodscov2$to))
sensor_nodscov2 %>% filter(!id %in% c(interaction_nodscov2$from, interaction_nodscov2$to)) 

# ## Explore contact patterns of HCW at night-------------------------------------
# # First morning
# present_night1 = agenda_nodscov2 %>%
#   filter(status == "PE", firstDate < as_datetime("2020-06-09 08:00:00"))
# 
# with_interaction_night1 = interaction_nodscov2 %>%
#   filter(from %in% present_night1$id | to %in% present_night1$id, date_posix < as_datetime("2020-06-09 08:00:00")) %>%
#   select(from, to) %>%
#   pivot_longer(c(from, to), values_to = "ids") %>%
#   distinct(ids) %>%
#   filter(grepl("PE", ids))
# 
# nrow(with_interaction_night1) / present_night1 %>% distinct(id) %>% nrow()
# 
# graph_data = interaction_nodscov2 %>%
#   filter(from %in% present_night1$id | to %in% present_night1$id, date_posix < as_datetime("2020-06-09 08:00:00")) %>%
#   mutate(date_posix = floor_date(date_posix, "day")) %>%
#   group_by(from, to, date_posix) %>%
#   summarise(length = sum(length), .groups = "drop") %>%
#   arrange(date_posix)
# 
# graph_example = simplify(graph_from_data_frame(graph_data, directed = F))
# vertex_atts = data.frame(id = vertex_attr(graph_example, "name")) %>%
#   left_join(raw_cat, "id")
# 
# graph_example = graph_example %>%
#   set_vertex_attr("cat", value = vertex_atts$cat)
# 
# ggplot(ggnetwork(graph_example, layout = igraph::layout_with_kk(graph_example)),
#        aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_edges(alpha = 0.5) +
#   geom_nodes(aes(colour = cat), size = 4) +
#   theme_blank() +
#   scale_colour_discrete(type = pal) +
#   labs(colour = "Category:") +
#   theme(legend.text = element_text(size=12))
# 
# # Second night
# present_night2 = agenda_nodscov2 %>%
#   filter(status == "PE", 
#          (firstDate > as_datetime("2020-06-09 18:00:00") & lastDate < as_datetime("2020-06-10 08:00:00")) |
#            (firstDate < as_datetime("2020-06-10 08:00:00") & lastDate > as_datetime("2020-06-10 10:00:00")))
# 
# with_interaction_night2 = interaction_nodscov2 %>%
#   filter(from %in% present_night2$id | to %in% present_night2$id, date_posix > as_datetime("2020-06-09 18:00:00"), date_posix < as_datetime("2020-06-10 08:00:00")) %>%
#   select(from, to) %>%
#   pivot_longer(c(from, to), values_to = "ids") %>%
#   distinct(ids) %>%
#   filter(grepl("PE", ids))
# 
# nrow(with_interaction_night2) / present_night2 %>% distinct(id) %>% nrow()
# 
# graph_data = interaction_nodscov2 %>%
#   filter(from %in% present_night2$id | to %in% present_night2$id, date_posix > as_datetime("2020-06-09 18:00:00"), date_posix < as_datetime("2020-06-10 08:00:00")) %>%
#   mutate(date_posix = floor_date(date_posix, "day")) %>%
#   group_by(from, to, date_posix) %>%
#   summarise(length = sum(length), .groups = "drop") %>%
#   arrange(date_posix)
# 
# graph_example = simplify(graph_from_data_frame(graph_data, directed = F))
# vertex_atts = data.frame(id = vertex_attr(graph_example, "name")) %>%
#   left_join(admission_nodscov2, "id")
# 
# graph_example = graph_example %>%
#   set_vertex_attr("cat", value = vertex_atts$cat)
# 
# ggplot(ggnetwork(graph_example, layout = igraph::layout_with_kk(graph_example)),
#        aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_edges(alpha = 0.5) +
#   geom_nodes(aes(colour = cat), size = 4) +
#   theme_blank() +
#   scale_colour_discrete(type = pal) +
#   labs(colour = "Category:") +
#   theme(legend.text = element_text(size=12))
# 
# # Individuals with very few interactions with patients at night
# with_patient_night1 = interaction_nodscov2 %>%
#   filter(from %in% present_night1$id | to %in% present_night1$id, 
#          date_posix < as_datetime("2020-05-06 08:00:00"),
#          grepl("PA-", from) | grepl("PA-", to)) %>%
#   mutate(pe = ifelse(grepl("PE-", from), from, to)) %>%
#   group_by(pe) %>%
#   summarise(n = n())
# 
# nrow(with_patient_night1) / present_night1 %>% distinct(id) %>% nrow()
# 
# 
# 
# interaction_nodscov2 %>%
#   filter(date_posix < as_datetime("2020-06-09 08:00:00") |
#            (date_posix > as_datetime("2020-06-09 18:00:00") & date_posix < as_datetime("2020-06-10 08:00:00"))) %>%
#   mutate(date_posix = floor_date(date_posix, "hour"), 
#          pair = case_when(
#            grepl("PE", from) & grepl("PE", to) ~ "PE-PE",
#            grepl("PA", from) & grepl("PA", to) ~ "PA-PA",
#            .default = "PA-PE"
#          )) %>%
#   group_by(pair, date_posix) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   ggplot(., aes(x = date_posix, y = n)) +
#   facet_grid(cols = vars(pair)) +
#   geom_line() +
#   theme_bw()

## Save cleaned and selected real data------------------------------------------
if (!dir.exists("data/data-nodscov2/clean/")) dir.create("data/data-nodscov2/clean/")

# Admission data 
write.csv2(admission_nodscov2, paste0("data/data-nodscov2/clean/admission_cleaned_", net, ".csv"), row.names = F, quote = F)

# HCW schedule
write.csv2(agenda_nodscov2, paste0("data/data-nodscov2/clean/agenda_cleaned_", net, ".csv"), row.names = F, quote = F)

# Schedule for all individuals
write.csv2(sensor_nodscov2, paste0("data/data-nodscov2/clean/sensor_cleaned_", net, ".csv"), row.names = F, quote = F)

# Interaction data
write.csv2(interaction_nodscov2, paste0("data/data-nodscov2/clean/interaction_cleaned_", net, ".csv"), row.names = F, quote = F)

## Frequency of contacts of more than one hour----------------------------------
# Frequency of contacts that last more than 1hour
interaction_nodscov2 %>%
  summarise(
    n = n(), 
    n_more_1h = sum(length/3600 >= 1),
    percentage = sum(length/3600 >= 1) / n() * 100
  )

# Distribution of contact duration for contacts lasting more than 1h
interaction_nodscov2 %>%
  filter(length >= 3600) %>%
  mutate(length = length/3600) %>%
  select(length) %>%
  summary()

# Interactions of more than 1hour
interaction_nodscov2 %>%
  filter(length >= 3600) %>%
  mutate(length = length/3600)

## Contact recurrence probability per hour--------------------------------------
# Hourly individual probability of recurring contacts and mean individual
# probability
pind = rep(NA, nrow(admission_nodscov2))
for (i in seq_along(sensor_nodscov2 %>% distinct(id) %>% pull(id))) { #seq_along(admission_nodscov2$id)) {
  ii = unique(sensor_nodscov2$id)[i]#admission_nodscov2$id[i]
  
  # Get n hours spent in the ward with a sensor 
  all_h = sensor_nodscov2 %>%
    filter(id == ii) %>%
    mutate(DATEREMISE = floor_date(DATEREMISE, "hour"), DATEREC = floor_date(DATEREC, "hour"), k = 1:n()) %>%
    arrange(DATEREMISE) %>%
    group_by(k) %>%
    nest() %>%
    mutate(out = map(data, unroll_dates)) %>%
    unnest(out) %>%
    .$date_list
  
  if (length(all_h) > 1) {
    interaction_h = vector("list", length(all_h))
    proba_h = rep(NA, length(all_h)-1)
    
    # Get lists of contacts per hour
    for (h in seq_along(all_h)) {
      hh = all_h[h]
      contacts = interaction_nodscov2 %>%
        filter(from == ii | to == ii, date_posix >= hh, date_posix < hh+3600) %>%
        select(from, to) %>%
        unlist()
      contacts = unique(contacts[contacts != ii])
      names(contacts) = NULL
      interaction_h[[h]] = contacts 
    }
    
    # Get hourly probability of recurring contacts
    if (length(all_h)>1) {
      for (h in 2:length(all_h)) {
        if (length(interaction_h[[h]]) > 0) {
          proba_h[h-1] = sum(interaction_h[[h]] %in% unlist(interaction_h[1:(h-1)])) / length(interaction_h[[h]])
        } else {
          proba_h[h-1] = 0
        }
      }
    }
    
    # Mean individual hourly probability of recurring contacts
    pind[i] = sum(proba_h) / length(proba_h)
    
  } else {
    pind[i] = 0
  }
  
}

# Plot mean individual probability of recurring contacts
admission_nodscov2 %>%
  mutate(pind = pind) %>%
  ggplot(., aes(x = cat, y = pind)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  labs(x = "", y = "Probability of hourly recurring contact", title = net)

# Comparison by category of individual
admission_nodscov2 %>%
  mutate(pind = pind) %>%
  rstatix::wilcox_test(pind ~ cat, ref.group = "Paramedical")

# Mean probability of recurring contacts by individual category
mean_recurrent_hourly = bind_rows(
  admission_nodscov2 %>%
    mutate(pind = pind) %>%
    group_by(cat) %>%
    summarise(m = mean(pind), .groups = "drop"),
  admission_nodscov2 %>%
    mutate(pind = pind, cat = ifelse(cat == "Patient", cat, "HCW")) %>%
    group_by(cat) %>%
    summarise(m = mean(pind), .groups = "drop")
) %>%
  distinct() %>%
  arrange(cat)
mean_recurrent_hourly
write.csv2(mean_recurrent_hourly, paste0("data/data-nodscov2/recurring_contacts/mean_recurrent_hourly_", net, ".csv"), row.names = F, quote = F)

## Contact recurrence probability per day---------------------------------------
# Daily individual probability of recurring contacts and mean individual
# probability
pind = rep(NA, nrow(admission_nodscov2))
for (i in seq_along(admission_nodscov2$id)) {
  ii = admission_nodscov2$id[i]

  # Get n days spent in the ward with a sensor
  all_d = sensor_nodscov2 %>%
    filter(id == ii) %>%
    mutate(DATEREMISE = floor_date(DATEREMISE, "day"), DATEREC = floor_date(DATEREC, "day"), k = 1:n()) %>%
    arrange(DATEREMISE) %>%
    group_by(k) %>%
    nest() %>%
    mutate(out = map(data, unroll_dates)) %>%
    unnest(out) %>%
    .$date_list

  if (length(all_d) > 1) {
    interaction_d = vector("list", length(all_d))
    proba_d = rep(NA, length(all_d)-1)

    # Get lists of contacts per hour
    for (d in seq_along(all_d)) {
      dd = all_d[d]
      contacts = interaction_nodscov2 %>%
        filter(from == ii | to == ii, date_posix >= dd, date_posix < dd+3600*24) %>%
        select(from, to) %>%
        unlist()
      contacts = unique(contacts[contacts != ii])
      names(contacts) = NULL
      interaction_d[[d]] = contacts
    }

    # Get hourly probability of recurring contacts
    if (length(all_d)>1) {
      for (hi in 2:length(all_d)) {
        if (length(interaction_d[[d]]) > 0) {
          proba_d[d-1] = sum(interaction_d[[d]] %in% unlist(interaction_d[1:(d-1)])) / length(interaction_d[[d]])
        } else {
          proba_d[d-1] = 0
        }
      }
    }

    # Mean individual hourly probability of recurring contacts
    pind[i] = sum(proba_d) / length(proba_d)

  } else {
    pind[i] = 0
  }

}

# Plot mean individual probability of recurring contacts
admission_nodscov2 %>%
  mutate(pind = pind) %>%
  filter(!is.na(pind)) %>%
  ggplot(., aes(x = cat, y = pind)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  labs(x = "", y = "Probability of daily recurring contact")

# Comparison by category of individual
admission_nodscov2 %>%
  mutate(pind = pind) %>%
  filter(!is.na(pind)) %>%
  rstatix::wilcox_test(pind ~ cat, ref.group = "Paramedical")

# Mean probability of recurring contacts by individual category
mean_recurrent_daily = bind_rows(
  admission_nodscov2 %>%
    mutate(pind = pind) %>%
    group_by(cat) %>%
    summarise(m = mean(pind, na.rm = T), .groups = "drop"),
  admission_nodscov2 %>%
    mutate(pind = pind, cat = ifelse(cat == "Patient", cat, "HCW")) %>%
    group_by(cat) %>%
    summarise(m = mean(pind, na.rm = T), .groups = "drop")
) %>%
  distinct() %>%
  arrange(cat)
mean_recurrent_daily
write.csv2(paste0("data/data-nodscov2/recurring_contacts/mean_recurrent_daily_", net, ".csv"), row.names = F, quote = F)

## Patient augmentation---------------------------------------------------------
# Double rooms - Cumulative time
patients_ids = unique(admission_nodscov2$id[admission_nodscov2$status == "PA"])
long_pa_pa_contacts = interaction_nodscov2 %>% 
  filter(from %in% patients_ids, to %in% patients_ids) %>%
  mutate(pair = paste0(from, "_", to)) %>%
  group_by(pair) %>%
  summarise(length = sum(length)/3600, .groups = "drop") %>%
  filter(length > 2)

# List of ids in double rooms
if (nrow(long_pa_pa_contacts) > 0) {
  double_rooms = interaction_nodscov2 %>%
    filter(from %in% patients_ids, to %in% patients_ids) %>%
    distinct(from, to)
  double_rooms = as.list(as.data.frame(t(double_rooms)))
  save(double_rooms, file = paste("data/data-nodscov2/clean/double_rooms_", net, ".rda")) 
}

# Patient rooms
patient_rooms = data.frame(
  id = admission_nodscov2$id[admission_nodscov2$status == "PA"],
  room = 1:sum(admission_nodscov2$status == "PA")
)

if (exists("double_rooms")) {
  for (r in double_rooms) {
    patient_rooms$room[patient_rooms$id %in% r] = min(patient_rooms$room[patient_rooms$id %in% r])
  }
}

# Patient augmentation
set.seed(20240924)
first_day = min(admission_nodscov2$firstDate) 
last_day = first_day + 90 - 1

for (id in patients_ids) {
  
  # Stay length of original individual
  if (admission_nodscov2$firstDate[admission_nodscov2$id == id] != admission_nodscov2$lastDate[admission_nodscov2$id == id]) {
    stay_length = sample(2:11, 1) # At least 2 days - max 11 days 
    new_lastDate = first_day + stay_length - 1
    admission_nodscov2$lastDate[admission_nodscov2$id == id] = new_lastDate
  } else {
    new_lastDate = admission_nodscov2$lastDate[admission_nodscov2$id == id]
  }
  
  # Augment patients by keeping incrementing the original id
  i = 1
  new_room = data.frame()
  while(new_lastDate < last_day) {
    stay_length = sample(6:13, 1) # between 6 and 13 days
    to_rbind = data.frame(
      id = paste0(id, "-", i),
      status = "PA",
      firstDate = new_lastDate+1,
      lastDate = min(new_lastDate+1+stay_length-1, last_day),
      ward = "Reanimation",
      sex = NA,
      hospitalization = "patient",
      cat = "Patient"
    )  
    
    patient_rooms = rbind(
      patient_rooms, 
      data.frame(
        id = paste0(id, "-", i),
        room = patient_rooms$room[patient_rooms$id == id]
      )
    )
    
    admission_nodscov2 = rbind(admission_nodscov2, to_rbind)
    new_lastDate = min(new_lastDate+1+stay_length-1, last_day)
    i = i+1
  }
}

# Final set of individuals in the database
admission_nodscov2 %>%
  group_by(status, cat) %>%
  summarise(n = n(), .groups = "drop") 


## Repeat schedule of HCWs over the 90-day period-------------------------------
# Repeat schedule over 90 days for individuals present only one day out 
# of the two days of recording
# HCWs that have day+night shifts 
every_other_day = sensor_nodscov2 %>%
  filter(grepl("PE-", id), difftime(DATEREC, DATEREMISE, units = "secs") >= 20*3600) %>%
  pull(id)

# Repeat for HCW that do not have day+night shifts
new_agenda = agenda_nodscov2 %>% rename(firstDate = DATEREMISE, lastDate = DATEREC)
for (d in 1:44) {
  new_agenda = rbind(
    new_agenda, 
    agenda_nodscov2 %>% 
      rename(firstDate = DATEREMISE, lastDate = DATEREC) %>%
      filter(!id %in% every_other_day) %>%
      mutate(firstDate = firstDate+3600*48*d, lastDate = lastDate+3600*48*d))
}

# Add HCW with long night shifts
for (d in 1:29) {
  new_agenda = rbind(
    new_agenda, 
    agenda_nodscov2 %>% 
      rename(firstDate = DATEREMISE, lastDate = DATEREC) %>%
      filter(id %in% every_other_day) %>%
      mutate(firstDate = firstDate+3600*72*d, lastDate = lastDate+3600*72*d))
}

# Add additionnal variables
new_agenda = new_agenda %>%
  left_join(., admission_nodscov2 %>% select(id, status, cat, ward), by = "id") %>%
  select(id, status, cat, ward, firstDate, lastDate)

# Get new admission dates for HCWs
new_admission_dates = new_agenda %>%
  group_by(id) %>%
  summarise(firstDate_new = as_date(floor_date(min(firstDate), "day")), 
            lastDate_new = as_date(floor_date(max(lastDate), "day")), .groups = "drop")

# Verify coherence and change admission dates for HCWs
admission_nodscov2 = admission_nodscov2 %>%
  left_join(., new_admission_dates, by = "id") %>%
  mutate(
    firstDate = case_when(status == "PE" ~ firstDate_new, .default = firstDate),
    lastDate = case_when(status == "PE" ~ lastDate_new, .default = lastDate)
  ) %>%
  select(-c(firstDate_new, lastDate_new))

## Verify the concordance between contacts and the schedule of healthcare workers--
# Schedule of healthcare workers 
ids_pe = admission_nodscov2$id[admission_nodscov2$status == "PE"]
interactions_in_schedule = data.frame()

for (i in ids_pe) {
  # Schedule of the HCW
  schedule_id = agenda_nodscov2 %>%
    rename(firstDate = DATEREMISE, lastDate = DATEREC) %>%
    filter(id == i) %>%
    select(firstDate, lastDate)
  
  # Get all interactions
  if (sum(interaction_nodscov2$to == i | interaction_nodscov2$from == i) > 0) {
    interactions_id = interaction_nodscov2 %>%
      filter(from == i| to == i) %>%
      select(date_posix, length) %>%
      mutate(n = 1:n())
    
    # Are all interactions within the schedule ?
    out_vec = expand.grid.df(
      schedule_id, 
      interactions_id
    ) %>%
      group_by(n) %>%
      summarise(out = any(date_posix >= firstDate & date_posix+length <= lastDate), .groups = "drop") %>%
      pull(out)
    
    interactions_in_schedule = bind_rows(
      interactions_in_schedule, 
      data.frame(
        id = i, 
        coherent_interactions = all(out_vec)
      )
    )
    
  } else {
    interactions_in_schedule = bind_rows(
      interactions_in_schedule, 
      data.frame(
        id = i, 
        coherent_interactions = NA
      )
    )
  }
}

# Number of HCW
admission_nodscov2 %>% filter(status == "PE") %>% nrow()

# Individuals with all interactions within schedule
sum(interactions_in_schedule$coherent_interactions, na.rm = T)

# Individuals with missing interactions
sum(is.na(interactions_in_schedule$coherent_interactions))

## Save data for synthetic algorithm--------------------------------------------
# Interaction data
write.csv2(interaction_nodscov2, paste0("data/data-synthetic-graphs/input/interactions_", net, ".csv"), quote = F, row.names = F)

# Admission data
write.csv2(admission_nodscov2, paste0("data/data-synthetic-graphs/input/admission_", net, ".csv"), quote = F, row.names = F)

# HCW schedule
write.csv2(new_agenda, paste0("data/data-synthetic-graphs/input/agenda_", net, ".csv"),  quote = F, row.names = F)

# Patient rooms
write.csv2(patient_rooms, paste0("data/data-synthetic-graphs/loc/patient_rooms_", net, ".csv"), quote = F, row.names = F)

## Plot raw data for supplementary materials------------------------------------
### Participant description
# Create directory to store figures
if (!dir.exists("fig/nodscov2-data")) dir.create("fig/nodscov2-data")

# Number of individuals by category
p1 = admission_nodscov2 %>%
  group_by(cat) %>%
  summarise(n =n()) %>%
  ggplot(., aes(x = cat, y = n, fill = cat)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = n, y = n+5)) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "", y = "Number of individuals")

# Number of interactions
p2 = interaction_nodscov2 %>%
  left_join(., admission_nodscov2 %>% select(id, cat) %>% rename(cat_from = cat), by = c("from" = "id")) %>%
  left_join(., admission_nodscov2 %>% select(id, cat) %>% rename(cat_to = cat), by = c("to" = "id")) %>%
  mutate(cat_pair = case_when(
    cat_from == "Patient" & cat_to == "Patient" ~ "Patient-Patient",
    (cat_from == "Patient" & cat_to != "Patient") | (cat_from != "Patient" & cat_to == "Patient") ~ "Patient-HCW",
    cat_from != "Patient" & cat_to != "Patient" ~ "HCW-HCW"
  ),
  id_pair = 1:n()) %>%
  filter(cat_pair != "Patient-Patient") %>%
  group_by(id_pair, from, to, cat_pair, date_posix) %>%
  nest() %>%
  mutate(out = map2(data, date_posix, function(.x, .y) data.frame(interacting = seq.POSIXt(.y, .y+.x$length, 10)))) %>%
  ungroup() %>%
  select(-data) %>%
  unnest(cols = out) %>%
  group_by(cat_pair, interacting) %>%
  summarise(n =n(), .groups = "drop") %>%
  ggplot(., aes(x= interacting, y = n, col = cat_pair)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("Patient-HCW" = "orange", "HCW-HCW" = "orchid")) +
  labs(x = "Time", y="Number of interactions", col = "Type of interaction")

# Combined plot
p = ggarrange(p1, p2, ncol = 2, align = "hv", widths = c(1,4), labels = c("A", "B"))
ggsave(paste0("fig/nodscov2-data/interactions_combined_", net, ".png"), p, height = 4, width = 10)

## Plot interaction network-----------------------------------------------------
# Create network objects
graph_data = interaction_nodscov2 %>%
  mutate(date_posix = floor_date(date_posix, "day")) %>%
  group_by(from, to, date_posix) %>%
  summarise(length = sum(length), .groups = "drop") %>%
  arrange(date_posix)

graph_data_PA_PA = graph_data %>%
  filter(grepl("PA-", from) & grepl("PA-", to))

graph_data_PE_PE = graph_data %>%
  filter(grepl("PE-", from) & grepl("PE-", to))

graph_data_PA_PE = graph_data %>%
  filter((grepl("PE-", from) & grepl("PA-", to)) | (grepl("PA-", from) & grepl("PE-", to)))

# Get all metrics
all_metrics = get_net_metrics(graph_data, adm_data = admission_nodscov2, network = "Full")
PA_PA_metrics = get_net_metrics(graph_data_PA_PA, adm_data = admission_nodscov2, network = "PA-PA")
PE_PE_metrics = get_net_metrics(graph_data_PE_PE, adm_data = admission_nodscov2, network = "PE-PE")
PA_PE_metrics = get_net_metrics(graph_data_PA_PE, adm_data = admission_nodscov2, network = "PA-PE")

# Summary statistics of full grpah and subgraphs
cols = colnames(all_metrics)[1:7]
summary_tab = rbind(all_metrics, PA_PA_metrics, PE_PE_metrics, PA_PE_metrics) %>%
  mutate(temp_corr = replace(temp_corr, temp_corr==0, NA)) %>%
  group_by(network) %>%
  summarise(across(all_of(cols), list(mean=mean, sd=sd), na.rm=T))

summary_tab[,-1] = round(summary_tab[,-1], 2)
View(summary_tab)

# Get all degrees
all_degrees = c()
for(d in unique(graph_data$date_posix)){
  data_d = graph_data %>%
    filter(date_posix == d)
  
  graph_d = graph_from_data_frame(data_d, directed = F)
  graph_d = simplify(graph_d)
  all_degrees = c(all_degrees, degree(graph_d))
}

# Distribution of daily number of degrees (dashed red line = mean degree)
# CV: coefficient of variation (standard deviation/mean)
pe = ggplot() +
  geom_histogram(aes(all_degrees, after_stat(density)), binwidth = 1, colour = "grey") +
  geom_vline(xintercept = mean(all_degrees), linetype = "dashed", linewidth = 1, colour = "red3") +
  geom_richtext(aes(x = 35, y = 0.05,
                    label = paste0("CV<sup>2</sup> = ",
                                   round((sd(all_degrees)/mean(all_degrees))^2, 3))),
                colour = "red3", size = 6) +
  theme_bw() +
  labs(x = "Node degree", y = "Frequency") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_x_continuous(breaks = seq(0,40,10))


# Full graph
example_date = min(graph_data$date_posix)
example_data = graph_data %>%
  filter(date_posix == example_date)
graph_example = graph_from_data_frame(example_data, directed = F)
graph_example = simplify(graph_example)

vertex_atts = data.frame(id = vertex_attr(graph_example, "name")) %>%
  left_join(admission_nodscov2, "id")

graph_example = graph_example %>%
  set_vertex_attr("cat", value = vertex_atts$cat)

pa = ggplot(ggnetwork(graph_example, layout = igraph::layout_with_kk(graph_example)),
            aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(colour = cat), size = 4) +
  theme_blank() +
  scale_colour_discrete(type = pal) +
  labs(colour = "Category:") +
  theme(legend.text = element_text(size=12))

# PA-PA subgraph
graph_example_PA_PA = graph_from_data_frame(graph_data_PA_PA %>% filter(date_posix == example_date),
                                            directed = F)
graph_example_PA_PA = simplify(graph_example_PA_PA)

vertex_atts = data.frame(id = vertex_attr(graph_example_PA_PA, "name")) %>%
  left_join(admission_nodscov2, "id")

graph_example_PA_PA = graph_example_PA_PA %>%
  set_vertex_attr("cat", value = vertex_atts$cat)

pb = ggplot(ggnetwork(graph_example_PA_PA, layout = igraph::layout_with_kk(graph_example_PA_PA)),
            aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(colour = cat), size = 3) +
  theme_blank() +
  scale_colour_discrete(type = pal["Patient"]) +
  labs(colour = "Category:") +
  theme(legend.text = element_text(size=12)) +
  guides(colour = "none")

# PE-PE subgraph
graph_example_PE_PE = graph_from_data_frame(graph_data_PE_PE %>% filter(date_posix == example_date),
                                            directed = F)
graph_example_PE_PE = simplify(graph_example_PE_PE)

vertex_atts = data.frame(id = vertex_attr(graph_example_PE_PE, "name")) %>%
  left_join(admission_nodscov2, "id")

graph_example_PE_PE = graph_example_PE_PE %>%
  set_vertex_attr("cat", value = vertex_atts$cat)

pc = ggplot(ggnetwork(graph_example_PE_PE, layout = igraph::layout_with_kk(graph_example_PE_PE)),
            aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(colour = cat), size = 3) +
  theme_blank() +
  scale_colour_discrete(type = pal[c("Medical", "Paramedical")]) +
  labs(colour = "Category:") +
  theme(legend.text = element_text(size=12)) +
  guides(colour = "none")

# PA-PE subgraph
graph_example_PA_PE = graph_from_data_frame(graph_data_PA_PE %>% filter(date_posix == example_date),
                                            directed = F)
graph_example_PA_PE = simplify(graph_example_PA_PE)

vertex_atts = data.frame(id = vertex_attr(graph_example_PA_PE, "name")) %>%
  left_join(admission_nodscov2, "id")

graph_example_PA_PE = graph_example_PA_PE %>%
  set_vertex_attr("cat", value = vertex_atts$cat)

pd = ggplot(ggnetwork(graph_example_PA_PE, layout = igraph::layout_with_kk(graph_example_PA_PE)),
            aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(colour = cat), size = 3) +
  theme_blank() +
  scale_colour_discrete(type = pal) +
  labs(colour = "Category:") +
  theme(legend.text = element_text(size=12)) +
  guides(colour = "none")

# Final plot
p = ggarrange(pa,
              ggarrange(pb,pc,pd, nrow=1, labels = c("B", "C", "D")),
              pe,
              heights = c(1,0.7,0.5),
              ncol = 1,
              labels = c("A", "", "E"))

ggsave(paste0("fig/nodscov2-data/network_subnetworks_", net, ".png"), p, height = 12, width = 10)

