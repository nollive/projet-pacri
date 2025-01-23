################################################################################
##
##                        COMPARE SYMMETRIZED SENSOR DATA 
##
################################################################################

## Working environment----------------------------------------------------------
rm(list = ls())

# Libraries
library(tidyverse)
library(ggpubr)
source("R/nodscov2/helper-functions.R")

# Working directories and paths
nodscov2_path = file.path("data", "data-nodscov2") 
fig_path = file.path("fig")

# Hospital of interest 
net = "herriot"

if (net == "poincare") {
  hosp = "APHP - RAYMOND POINCARÉ"
  hosp_min = "raymond_poincare"
  w = "Medical ICU #1"
}
if (net == "herriot") {
  hosp = "HC DE LYON - EDOUARD HERRIOT"
  hosp_min = "edouard_herriot"
  w = "Medical ICU #2"
}

## Load Nods-CoV-2 data (not shared)--------------------------------------------
# Admission dates with all participants (no filter based on individual category)
load(file.path(nodscov2_path, "raw", "initial_sym", "Post-processed", "admission_ctc_nodscov2.RData"))
Encoding(admission$hospital) <- 'latin1'
admission = admission %>%
  filter(
    hospital == hosp, 
    ward == "Reanimation"
  ) %>%
  mutate(cat_recoded = recode(cat, !!!dict_cat_initial)) %>%
  mutate(cat_recoded = ifelse(is.na(cat), "Patient", cat_recoded)) %>%
  select(-hospital)

# Interactions - George
  # Filters used :
  # 1. max RSSI > -75 during the contact
  # 2. true contacts are separated by 63 seconds without contact 
list_ward = read.csv(file.path(nodscov2_path, "raw", "initial_sym", "Post-processed", "list_ward_complete.csv")) %>%
  rename(date_posix = date_posix_first) %>%
  filter(from %in% admission$id, to %in% admission$id) %>%
  select(-c(ward_id, wardType, newID)) %>%
  mutate(date_posix = as_datetime(date_posix))

# Sensor distribution - George
sensor = read.csv(file.path(nodscov2_path, "raw", "initial_sym", "Post-processed", "met_timespent_mins_complete.csv")) %>%
  filter(id %in% admission$id)

# Load Maximilien's symmetrized data (same filters as George's)
list_ward_m = read.csv(paste0("data/data-nodscov2/raw/final_sym/", hosp_min, "_rea/", hosp_min, "_symetrise_finale_63-75.csv")) %>%
  select(Nom1, Nom2, temps.début, temps.fin, durée) %>%
  rename(from = Nom1, to = Nom2, date_posix = temps.début, date_posix_end = temps.fin, length_m = durée) %>%
  mutate(date_posix = as_datetime(date_posix),
         date_posix_end = as_datetime(date_posix_end)) %>%
  filter(!from %in% c("005-0079-P-P", "005-0114-L-G"), !to %in% c("005-0079-P-P", "005-0114-L-G"))

# Load Antoine's data
list_ward_a = read.csv(paste0("data/data-nodscov2/raw/initial_sym/Symmetrized/finalized-or-", hosp_min, ".csv"), header = F)
colnames(list_ward_a) = c("sensor1", "sensor2", "from", "to", "date_posix", "date_posix_end", "length", "rssi1", "rssi2", "rssi3")

list_ward_a = list_ward_a %>%
  select(from, to, date_posix, date_posix_end, length) %>%
  mutate(date_posix = as_datetime(date_posix), date_posix_end = as_datetime(date_posix_end),
         from = gsub(" ", "", from), to = gsub(" ", "", to))

# Load symmetrisation calendar
sensor_m = read.csv(paste0("data/data-nodscov2/raw/final_sym/", hosp_min, "_rea/calendar-", hosp_min, "-Reanimation.csv"), 
                    header = F)
colnames(sensor_m) = c("id_sensor", "id", "x1", "DATEREMISE", "DATEREC", "length", "x2")
sensor_m = sensor_m %>%
  mutate(DATEREMISE = as_datetime(DATEREMISE), DATEREC = as_datetime(DATEREC)) %>%
  select(id, DATEREMISE, DATEREC) %>%
  filter(!id %in% c("005-0079-P-P", "005-0114-L-G"))

# Data symmetrized by Maximilien with same assumptions as Antoine
if (net == "herriot") {
  list_ward_mbis = read.csv("data/data-nodscov2/raw/intial_sym_modified/edouard_herriot_symetrise_final-63-75_bilateral.csv") %>%
    select(Nom1, Nom2, temps.début, temps.fin, durée) %>%
    rename(from = Nom1, to = Nom2, date_posix = temps.début, date_posix_end = temps.fin, length_mbis = durée) %>%
    mutate(date_posix = as_datetime(date_posix), date_posix_end = as_datetime(date_posix_end)) #%>%
    #filter(!from %in% c("005-0079-P-P", "005-0114-L-G"), !to %in% c("005-0079-P-P", "005-0114-L-G"))
}

## Rearrange pair names---------------------------------------------------------
# George's data
list_ward %>% 
  distinct(from, to) %>%
  mutate(i = 1:n()) %>%
  nest(.by = i) %>%
  mutate(data = map(data, get_pair_id)) %>%
  unnest(data) %>%
  group_by(pair) %>%
  mutate(n = n()) %>%
  filter(n>1)

# Maximilien's data
list_ward_m %>% 
  distinct(from, to) %>%
  mutate(i = 1:n()) %>%
  nest(.by = i) %>%
  mutate(data = map(data, get_pair_id)) %>%
  unnest(data) %>%
  group_by(pair) %>%
  mutate(n = n()) %>%
  filter(n>1)

# Antoine's data
list_ward_a %>% 
  distinct(from, to) %>%
  mutate(i = 1:n()) %>%
  nest(.by = i) %>%
  mutate(data = map(data, get_pair_id)) %>%
  unnest(data) %>%
  group_by(pair) %>%
  mutate(n = n()) %>%
  filter(n>1)

# Maximilien's data with same hypotheses as Antoine
list_ward_mbis %>% 
  distinct(from, to) %>%
  mutate(i = 1:n()) %>%
  nest(.by = i) %>%
  mutate(data = map(data, get_pair_id)) %>%
  unnest(data) %>%
  group_by(pair) %>%
  mutate(n = n()) %>%
  filter(n>1)

## Modify data to make them comparable------------------------------------------
# Fusion overlapping contacts in George's data
list_ward_nooverlap = list_ward %>% 
  # Keep longest contacts
  group_by(from, to, date_posix) %>%
  summarise(length = max(length), .groups = "drop") %>%
  # Remove duplicates 
  distinct() %>%
  # Fusion overlapping contacts
  nest(.by = c(from, to)) %>%
  mutate(data = map(data, fusion_overlapping_contacts)) %>%
  unnest(data) %>%
  distinct()

# Round dates in Maximilien's data
list_ward_mbis_rounded = list_ward_mbis %>%
  mutate(date_posix_end = date_posix + length_mbis) %>%
  mutate(
    date_posix_round = floor_date(date_posix, "10 seconds"),
    date_posix_end_round = floor_date(date_posix_end, "10 seconds")
  ) %>%
  mutate(length_m = as.numeric(difftime(date_posix_end_round, date_posix_round, units = "secs"))) %>%
  select(from, to, date_posix_round, length_mbis) %>%
  rename(date_posix = date_posix_round)

## Compare antoine's and maximilien's data--------------------------------------
# Number of individuals
part_a = sort(unique(c(list_ward_a$from, list_ward_a$to)))
part_m = sort(unique(c(list_ward_m$from, list_ward_m$to)))
part_mbis = sort(unique(c(list_ward_mbis$from, list_ward_mbis$to)))
identical(part_a, part_m)
identical(part_a, part_mbis)

length(part_a)
length(part_mbis)
part_a[!part_a %in% part_mbis]

# Number of unique pairs 
pairs_a = list_ward_a %>% mutate(p = gsub(" ", "", paste0(from, "_", to))) %>% distinct(p) %>% arrange(p) %>% .$p
pairs_m = list_ward_m %>% mutate(p = paste0(from, "_", to)) %>% distinct(p) %>% arrange(p) %>% .$p
pairs_mbis = list_ward_mbis %>% mutate(p = paste0(from, "_", to)) %>% distinct(p) %>% arrange(p) %>% .$p
identical(pairs_a, pairs_m)
identical(pairs_a, pairs_mbis)
length(pairs_a)
length(pairs_m)
length(pairs_mbis)
head(pairs_m[!pairs_m %in% pairs_a])

list_ward_m %>% 
  filter(from == "001-0001-D-J", to == "001-0004-V-A") %>%
  mutate(date_posix = as.numeric(date_posix))

list_ward_a %>% 
  filter(to == "001-0001-D-J", from == "001-0004-V-A")

# Total number of contacts that are not shared
list_ward_a %>%
  select(-date_posix_end) %>%
  distinct() %>%
  anti_join(., list_ward_mbis %>% select(-date_posix_end) %>% rename(length = length_mbis), by = c("from", "to", "date_posix", "length")) %>%
  nrow()
list_ward_a %>% distinct() %>% nrow()
list_ward_mbis %>% distinct() %>% nrow()

# 
pairs_df_a = list_ward_a %>% distinct(from, to)
pairs_df_m = list_ward_m %>% distinct(from, to)

pairs_df_m %>%
  anti_join(., pairs_df_a, by = c("from", 'to')) %>%
  pivot_longer(c(from, to), names_to = "origin", values_to = "p") %>%
  count(p) %>%
  arrange(desc(n))
  
# Overlap in Antoine's data
list_ward_a %>%
  nest(.by = c(from, to)) %>%
  mutate(data = map(data, get_overlapping_contacts)) %>%
  unnest(data) %>%
  filter(overlap == T) %>%
  nrow()

# Overlap in Maximilien's data symmatrized with same hypotheses as Antoine
list_ward_mbis %>%
  rename(length = length_mbis) %>%
  nest(.by = c(from, to)) %>%
  mutate(data = map(data, get_overlapping_contacts)) %>%
  unnest(data) %>%
  filter(overlap == T) %>%
  nrow()

## Compare number of participants-----------------------------------------------
all_m = sort(unique(c(list_ward_m$from, list_ward_m$to)))
all(all_m %in% admission$id)

all_mbis = sort(unique(c(list_ward_mbis$from, list_ward_mbis$to)))
all(all_mbis %in% admission$id)

all(admission$id %in% all_m)
all(admission$id %in% all_mbis)
admission$id[!admission$id %in% all_mbis]

all(all_m %in% sensor_m$id)
sensor_m$id[!sensor_m$id %in% all_m] # sensors 005-0079-P-P and 005-0114-L-G

## Compare sensor wearing data--------------------------------------------------
# Clean George's sensor data
sensor = sensor %>%
  select(id, DATEREMISEJ1, DATERECJ1, DATEREMISEJ2, DATERECJ2) %>%
  pivot_longer(-id, names_to = c(".value", "set"), names_pattern = "(.*)J(.$)") %>%
  filter(!is.na(DATEREMISE), !is.na(DATEREC)) %>%
  mutate(DATEREMISE = as_datetime(gsub("Z", "", gsub("T", " ", DATEREMISE))),
         DATEREC = as_datetime(gsub("Z", "", gsub("T", " ", DATEREC)))
  ) %>%
  select(-set)

# Compare sensor data 
symdiff(sensor, sensor_m) %>%
  arrange(id, DATEREMISE) %>%
  left_join(., sensor %>% mutate(data = "George"), by = c("id", "DATEREMISE", "DATEREC")) %>%
  mutate(data = ifelse(is.na(data), "Maximilien", data)) %>%
  write.csv2(., paste0("data/data-nodscov2/raw/comparison_calendars", net, ".csv"), quote = F, row.names = F)

## Compare contact data---------------------------------------------------------
# Number of recorded contacts
nrow(list_ward)
nrow(list_ward_a)
nrow(list_ward_m)

nrow(list_ward_mbis)
list_ward_a %>%
  nest(.by = c(from, to)) %>%
  mutate(data = map(data, get_overlapping_contacts)) %>%
  unnest(data) %>%
  filter(overlap == F) %>%
  nrow()

# Number of unique pairs of contacts
list_ward %>% distinct(from, to) %>% nrow()
list_ward_m %>% distinct(from, to) %>% nrow()
list_ward_mbis %>% distinct(from, to) %>% nrow()

# Contacts in agreement
list_ward_m %>%
  left_join(., list_ward, by = c("from", "to", "date_posix")) %>%
  summary()

# Number of duplicated contacts ? overlapping contacts ? 
nrow(list_ward) - nrow(distinct(list_ward))
list_ward %>%
  nest(.by = c(from, to)) %>%
  mutate(data = map(data, get_overlapping_contacts)) %>%
  unnest(data) %>%
  filter(overlap == T) %>%
  nrow()

## Impact of date rounding in symmetrized Maximilien's data---------------------
# Are there duplicates in Maximilien's data ? overlapping contacts ? 
nrow(list_ward_m) == nrow(distinct(list_ward_m))
list_ward_m %>%
  rename(length = length_m) %>%
  nest(.by = c(from, to)) %>%
  mutate(data = map(data, get_overlapping_contacts)) %>%
  unnest(data) %>%
  filter(overlap == T)

# Date rounding 
list_ward_round_m = list_ward_m %>%
  mutate(date_posix_end = date_posix + length_m) %>%
  mutate(
    date_posix_round = floor_date(date_posix, "10 seconds"),
    date_posix_end_round = floor_date(date_posix_end, "10 seconds")
    ) %>%
  mutate(length_m = as.numeric(difftime(date_posix_end_round, date_posix_round, units = "secs"))) %>%
  select(from, to, date_posix_round, length_m) %>%
  rename(date_posix = date_posix_round)

# Does date rounding impact the number of duplicates ? overlapping contacts ?
nrow(list_ward_round_m) == nrow(distinct(list_ward_round_m))
list_ward_round_m %>%
  rename(length = length_m) %>%
  nest(.by = c(from, to)) %>%
  mutate(data = map(data, get_overlapping_contacts)) %>%
  unnest(data) %>%
  filter(overlap == T)

## Number of distinct contacts and total hours of interactions------------------
# Total number of hours in interaction
total_h = list_ward %>%
  pivot_longer(c(from, to), values_to = "id", names_to = "names") %>%
  left_join(., admission %>% select(id, cat_recoded), by = c("id" = "id")) %>%
  group_by(id, cat_recoded) %>%
  summarise(length = sum(length)/3600, .groups = "drop") %>%
  filter(cat_recoded %in% c("Patient", "Visitor", "HCW")) %>%
  mutate(data = "George")

total_hnooverlap = list_ward_nooverlap %>%
  pivot_longer(c(from, to), values_to = "id", names_to = "names") %>%
  left_join(., admission %>% select(id, cat_recoded), by = c("id" = "id")) %>%
  group_by(id, cat_recoded) %>%
  summarise(length = sum(length)/3600, .groups = "drop") %>%
  filter(cat_recoded %in% c("Patient", "Visitor", "HCW")) %>%
  mutate(data = "George without overlaps")

total_h_m = list_ward_mbis_rounded %>%
  rename(length = length_mbis) %>%
  pivot_longer(c(from, to), values_to = "id", names_to = "names") %>%
  left_join(., admission %>% select(id, cat_recoded), by = c("id" = "id")) %>%
  group_by(id, cat_recoded) %>%
  summarise(length = sum(length)/3600, .groups = "drop") %>%
  filter(cat_recoded %in% c("Patient", "Visitor", "HCW")) %>%
  mutate(data = "Maximilien")

pa = bind_rows(total_h_m, total_hnooverlap, total_h) %>%
  ggplot(., aes(x = length, y = data)) +
  facet_grid(rows = vars(cat_recoded)) +
  geom_violin() +
  geom_point(data = bind_rows(total_h_m, total_hnooverlap, total_h) %>% group_by(data, cat_recoded) %>% summarise(length = median(length), .groups = "drop"), aes(x = length, y = data), col = "red") +
  geom_point(data = bind_rows(total_h_m, total_hnooverlap, total_h) %>% group_by(data, cat_recoded) %>% summarise(length = mean(length), .groups = "drop"), aes(x = length, y = data), col = "red", shape = 4) +
  theme_bw() +
  labs(x = "Total hours in contact over study period", y = "", title = ifelse(net == "poincare", "Poincaré", "Herriot"))
  
# Number of distinct contacts
total_c = list_ward %>%
  distinct(from, to) %>%
  pivot_longer(c(from, to), values_to = "id", names_to = "names") %>%
  left_join(., admission %>% select(id, cat_recoded), by = c("id" = "id")) %>%
  group_by(id, cat_recoded) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(cat_recoded %in% c("Patient", "Visitor", "HCW")) %>%
  mutate(data = "George")

total_cnooverlap = list_ward_nooverlap %>%
  distinct(from, to) %>%
  pivot_longer(c(from, to), values_to = "id", names_to = "names") %>%
  left_join(., admission %>% select(id, cat_recoded), by = c("id" = "id")) %>%
  group_by(id, cat_recoded) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(cat_recoded %in% c("Patient", "Visitor", "HCW")) %>%
  mutate(data = "George without overlaps")

total_c_m = list_ward_mbis_rounded %>%
  distinct(from, to) %>%
  pivot_longer(c(from, to), values_to = "id", names_to = "names") %>%
  left_join(., admission %>% select(id, cat_recoded), by = c("id" = "id")) %>%
  group_by(id, cat_recoded) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(cat_recoded %in% c("Patient", "Visitor", "HCW")) %>%
  mutate(data = "Maximilien")

pb = bind_rows(total_c_m, total_cnooverlap, total_c) %>%
  ggplot(., aes(x = n, y = data)) +
  facet_grid(rows = vars(cat_recoded)) +
  geom_violin() +
  geom_point(data = bind_rows(total_c_m, total_cnooverlap, total_c) %>% group_by(data, cat_recoded) %>% summarise(n = median(n), .groups = "drop"), aes(x = n, y = data), col = "red") +
  geom_point(data = bind_rows(total_c_m, total_cnooverlap, total_c) %>% group_by(data, cat_recoded) %>% summarise(n = mean(n), .groups = "drop"), aes(x = n, y = data), col = "red", shape = 4) +
  theme_bw() +
  labs(x = "Total unique contacts over study period", y = "", title = ifelse(net == "poincare", "Poincaré", "Herriot"))
ggsave(paste0("fig/symmetrization/total_hours_contacts_", hosp_min, ".png"), ggarrange(pa, pb), height = 6, width = 10)

## Compare number of interactions per hour--------------------------------------
# Extract number of unique contacts and contact duration from both networks (synthetic + observed)
durations = data.frame()
numbers = data.frame()
sources = c("George", "George without overlaps", "Maximilien")

for (s in sources) {
  if (s == "George") contacts = list_ward
  if (s == "George without overlaps") contacts = list_ward_nooverlap 
  if (s == "Maximilien") contacts = list_ward_mbis_rounded %>% rename(length = length_mbis)
  
  # Prepare observed data 
  durations = bind_rows(
    durations, 
    contacts %>%
      select(from, to, date_posix, length) %>%
      mutate(date_posix = floor_date(date_posix, "hour")) %>%
      group_by(from, to, date_posix) %>%
      summarise(length = sum(length)/60, .groups = "drop") %>% # Convert in minutes
      select(length) %>%
      mutate(data = s)
  )
  
  numbers = bind_rows(
    numbers,
    contact_numbers(contacts %>%
                      left_join(., admission[, c("id", "cat_recoded")], by = c("from" = "id")) %>%
                      rename(from_cat = cat_recoded) %>%
                      left_join(., admission[, c("id", "cat_recoded")], by = c("to" = "id")) %>%
                      rename(to_cat = cat_recoded)
                      , s, NA, "Initial")
  )
  
}

# Change level order
numbers = numbers %>%
  mutate(date_posix = factor(date_posix, levels = c("0:00","1:00","2:00","3:00","4:00",
                                                    "5:00","6:00","7:00","8:00","9:00",
                                                    "10:00","11:00","12:00","13:00","14:00",
                                                    "15:00","16:00","17:00","18:00","19:00",
                                                    "20:00","21:00","22:00","23:00"))
         )

# Plot number of unique contacts
pa = ggplot(numbers, aes(x = date_posix, y = med, ymin = q25, ymax = q75, fill = data, col = data, group = interaction(data, type))) +
  geom_ribbon(alpha = 0.3) +
  geom_point() +
  geom_line() +
  facet_wrap(facets = vars(type), ncol = 3, scales = "free_y") +
  scale_colour_discrete(type = c("dodgerblue4", "orange")) +
  scale_fill_discrete(type = c( "dodgerblue4", "orange")) +
  scale_x_discrete(breaks = c("0:00", "4:00", "8:00", "12:00",
                              "16:00","20:00")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  guides(fill = FALSE, color = guide_legend(override.aes = list(fill = NA))) +
  labs(x = "Hour", y = "Number of unique contacts", colour = "",
       title = ifelse(net == "poincare", "Poincaré", "Herriot"))
ggsave(paste0("fig/symmetrization/n_unique_contacts_", hosp_min, ".png"), pa, height = 6, width = 10)

# Plot contact duration 
pb = ggplot(durations, aes(y=length, x = data, colour = data)) +
  geom_boxplot(outlier.shape = NA) +
  scale_colour_discrete(type = c("dodgerblue4", "orange")) +
  coord_cartesian(ylim = c(0,35)) +
  theme_bw() +
  labs(y = "Contact duration (minutes)", x = "", col = "", 
       title= ifelse(net == "poincare", "Poincaré", "Herriot")) 
ggsave(paste0("fig/symmetrization/contact_duration_", hosp_min, ".png"), pb, height = 3, width = 4)

