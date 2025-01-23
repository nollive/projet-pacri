################################################################################
##
##              FIGURES OF LOCATION RECONSTRUCTION PROCEDURE
##
################################################################################

## Working environment----------------------------------------------------------
# Libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(igraph)
library(tidyr)
library(viridis)
library(purrr)
library(stringr)
library(extrafont)
library(extrafontdb)
library(cowplot)
rm(list=ls())

# Helper functions
source("helper-functions.R")

# Figure paths
fig_loc_path <- file.path("..","..", "fig", "loc-reconstruction")
ind_paths_path <- file.path("..","..", "fig", "individual-paths")
if (!dir.exists(fig_loc_path)) dir.create(fig_loc_path, recursive = TRUE)
if (!dir.exists(ind_paths_path)) dir.create(ind_paths_path, recursive = TRUE)

# Input data paths
threshold = 90*2
network = "herriot"

## Load data--------------------------------------------------------------------
all_data = list()
data_path = file.path("..", "..", "data", "data-synthetic-graphs") 
## Observed data 
load(file.path(data_path, "loc", paste0(network, "-observed-reconstructed-locations-", threshold, ".rda")))
all_data[["Observed network"]] = list(
  clusters = clusters,
  paths = paths,
  global_interaction = global_interaction,
  admission = admission,
  rooms = rooms,
  threshold = threshold,
  id_hcw = id_hcw,
  id_patient = id_patient,
  begin_date = begin_date,
  end_date = end_date,
  n_subdivisions = n_subdivisions
)
rm(list = c("clusters", "paths", "global_interaction", "admission", "rooms",
            "id_hcw", "id_patient", "begin_date", "end_date", "n_subdivisions"))

## Synthetic network
load(file.path(data_path, "loc", paste0(network, "-simulated-reconstructed-locations-", threshold, ".rda")))
all_data[["Synthetic network"]] = list(
  clusters = clusters,
  paths = paths,
  global_interaction = global_interaction,
  admission = admission,
  rooms = rooms,
  threshold = threshold,
  id_hcw = id_hcw,
  id_patient = id_patient,
  begin_date = begin_date,
  end_date = end_date,
  n_subdivisions = n_subdivisions
)
rm(list = c("clusters", "paths", "global_interaction", "admission", "rooms",
            "id_hcw", "id_patient", "begin_date", "end_date", "n_subdivisions"))


## Plot distribution of time spent in resting rooms/corridor--------------------
## Room dictionary
dict_rooms_simplified = all_data$`Observed network`$rooms %>%
  select(room, id_room) %>%
  mutate(room = ifelse(room == id_room, "Patient Room", room)) %>%
  distinct()
dict_rooms_simplified = setNames(as.character(dict_rooms_simplified$room), dict_rooms_simplified$id_room)
dict_rooms_simplified = c(dict_rooms_simplified, "-2" = "Outside ward")

# Load data and get cumulative time in each room per day
cum_time_data = data.frame()

for (d in names(all_data)) {
  # Get data
  paths = all_data[[d]]$paths
  admission = all_data[[d]]$admission
  begin_date = all_data[[d]]$begin_date
  end_date = all_data[[d]]$end_date
  
  # Dictionaries 
  dict_times = setNames(as.character(seq(begin_date, end_date-30, 30)),
                        as.character(1:nrow(paths)))
  
  dict_id_cat = setNames(admission %>% select(id, cat) %>% pull(cat),
                         admission %>% select(id, cat) %>% pull(id))
  
  ## Boxplot
  cumulative_time = do.call("rbind", mapply(
    function(x,y) {
      out = data.frame(
        time = floor_date(seq(begin_date, end_date-30, 30), "day"),
        loc = x
      ) %>%
        filter(loc != "-1") %>%
        mutate(loc = recode(loc, !!!dict_rooms_simplified)) %>%
        group_by(time, loc) %>%
        # Convert time spent in each room into hours 
        summarise(time_spent = n()*30/3600, .groups = "drop") 
      
      if (grepl("^PE-", y)) {
        time_per_day_in_ward = data.frame(
          time = floor_date(seq(begin_date, end_date-30, 30), "day"),
          loc = x
        ) %>%
          filter(loc != "-1") %>%
          group_by(time) %>%
          # Convert time spent in each room into hours 
          summarise(time_in_ward = n()*30/3600, .groups = "drop")
        
        out = out %>%
          left_join(., time_per_day_in_ward, by = "time") %>%
          mutate(time_spent = time_spent / time_in_ward) %>%
          group_by(loc) %>%
          summarise(time_spent = mean(time_spent), .groups = "drop")
        
      } else {
        out = out %>%
          group_by(time) %>%
          mutate(time_spent = time_spent/sum(time_spent)) %>%
          ungroup() %>%
          group_by(loc) %>%
          summarise(time_spent = mean(time_spent), .groups = "drop")
      }
      out$id = y
      return(out)
    },
    lapply(seq_len(ncol(paths)), function(i) paths[,i]),
    colnames(paths),
    SIMPLIFY = F
  ))
  
  # Same data but displayed in a different way (all locations for each individual category)
  cum_time_data = bind_rows(
    cum_time_data, 
    cumulative_time %>%
      mutate(
        id = gsub("\\.", "-", id),
        loc = gsub("Restroom", "resting room", loc),
        db = d
      ) %>%
      left_join(., admission %>% select(id, cat), by = "id")
  )
}

# Plots
n = ifelse(network == "herriot", "Herriot", "Poincaré")
th = setNames(c("30 min", "1h", "1.5h"), c(30*2, 60*2, 90*2))

# Plot for medical staff
pa = cum_time_data %>%
  filter(cat == "Medical") %>%
  complete(nesting(db, id, cat), loc, fill = list(time_spent = 0)) %>%
  ggplot(., aes(x = loc, y = time_spent, fill = db)) +
  geom_boxplot(position = position_dodge(width = 0.8), outliers = F) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 0.1) +
  expand_limits(y = c(0,1)) +
  facet_grid(cols = vars(cat)) +
  scale_fill_manual(values = c("Observed network" = "orange", "Synthetic network" = "darkorchid")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  labs(y = "Average proportion of daily time spent in room")

# Plot for paramedical staff
pb = cum_time_data %>%
  filter(cat == "Paramedical") %>%
  complete(nesting(db, id, cat), loc, fill = list(time_spent = 0)) %>%
  ggplot(., aes(x = loc, y = time_spent, fill = db)) +
  geom_boxplot(position = position_dodge(width = 0.8), outliers = F) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 0.1) +
  expand_limits(y = c(0,1)) +
  facet_grid(cols = vars(cat)) +
  scale_fill_manual(values = c("Observed network" = "orange", "Synthetic network" = "darkorchid")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(title = paste0(n, " - ", th[as.character(threshold)]))

# Plot for patients 
pc = cum_time_data %>%
  filter(cat == "Patient") %>%
  complete(nesting(db, id, cat), loc, fill = list(time_spent = 0)) %>%
  ggplot(., aes(x = loc, y = time_spent, fill = db)) +
  geom_boxplot(position = position_dodge(width = 0.8), outliers = F) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 0.1) +
  expand_limits(y = c(0,1)) +
  facet_grid(cols = vars(cat)) +
  scale_fill_manual(values = c("Observed network" = "orange", "Synthetic network" = "darkorchid")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Combined plot
p = ggarrange(pa, pb, pc, ncol = 3, align = "h", widths = c(1,0.8,0.4), common.legend = T, legend = "bottom")
ggsave(file.path(fig_loc_path, paste0("cumulative_time_", network, "_", threshold, ".png")), p, height = 4, width = 11)

## Plot individual trajectories for randomly selected individuals and days------
# Example of individual trajectories 
set.seed(123)
if (data_to_load == "observed") agenda = read.csv2(file.path(data_path, "clean", "agenda_cleaned.csv"))
if (data_to_load == "simulated") agenda = read.csv2(file.path(data_path, "input", "agenda.csv"))
individuals_to_show = c(
  sample(admission$id[admission$cat == "Paramedical"], 5),
  sample(admission$id[admission$cat == "Medical"], 5),
  sample(id_patient, 5)
)

for (y in individuals_to_show) {
  if (grepl("^PA-", y)) {
    days_of_presence = seq(admission %>% filter(id == y) %>% pull(firstDate), 
                           admission %>% filter(id == y) %>% pull(lastDate), 
                           1)
  } else {
    days_of_presence = agenda %>% 
      filter(id == y) %>% 
      mutate(k = 1:n(), firstDate = as_date(floor_date(as_datetime(firstDate), "day")), lastDate = as_date(floor_date(as_datetime(lastDate), "day"))) %>%
      nest(.by = k) %>%
      mutate(data = map(data, unroll_days)) %>%
      unnest(data) %>%
      pull(data)
  }
  
  random_day = sample(unique(days_of_presence), 1)
  individual_path = data.frame(
    time = seq(begin_date, end_date-30, 30),
    loc = paths[[gsub("-", ".", y)]],
    cat = admission$cat[admission$id == y] 
  ) %>%
    filter(floor_date(time, "day") == random_day) %>%
    mutate(loc = recode(loc, !!!c(dict_rooms, "-1" = "Absent")))
  individual_cat = unique(individual_path$cat)
  
  p = ggplot(individual_path, aes(x = time, y = loc, group=cat)) +
    geom_path(col = pal[individual_cat]) +
    theme_bw() +
    labs(x = "", y = paste("Location of", y))
  ggsave(file.path(ind_paths_path, data_to_load, paste0(gsub("\\.", "-", y), ".png")), p, height = 4, width = 10)
}

# Plot one week for the same individual
y = individuals_to_show[1]
if (grepl("^PA-", y)) {
  days_of_presence = seq(admission %>% filter(id == y) %>% pull(firstDate), 
                         admission %>% filter(id == y) %>% pull(lastDate), 
                         1)
} else {
  days_of_presence = agenda %>% 
    filter(id == y) %>% 
    mutate(k = 1:n(), firstDate = as_date(floor_date(as_datetime(firstDate), "day")), lastDate = as_date(floor_date(as_datetime(lastDate), "day"))) %>%
    nest(.by = k) %>%
    mutate(data = map(data, unroll_days)) %>%
    unnest(data) %>%
    pull(data)
}

random_day = sample(unique(days_of_presence), 1)
individual_path = data.frame(
  time = seq(begin_date, end_date-30, 30),
  loc = paths[[gsub("-", ".", y)]],
  cat = admission$cat[admission$id == y] 
) %>%
  filter(floor_date(time, "day") %in% seq(random_day, random_day+7, 1)) %>%
  mutate(loc = recode(loc, !!!c(dict_rooms, "-1" = "Absent")))
individual_cat = unique(individual_path$cat)

p = ggplot(individual_path, aes(x = time, y = loc, group=cat)) +
  geom_path(col = pal[individual_cat]) +
  theme_bw() +
  labs(x = "", y = paste("Location of", y))
# ggsave(file.path(ind_paths_path, paste0(gsub("\\.", "-", y), ".png")), p, height = 4, width = 10)

## Plot contact-related information---------------------------------------------
# Patient - Patient contact rate per hour 
interactions_synthetic = read.csv2(file.path(data_path, "input", "interactions_trimmed.csv")) %>%
  mutate(
    interaction_type = case_when(
      from_status == "PA" & to_status == "PA" ~ "PA-PA",
      from_status != to_status ~ "PA-PE",
      from_status == "PE" & to_status == "PE" ~ "PE-PE"
    ),
    date_hour = floor_date(as_datetime(date_posix), "hour")
  ) %>%
  group_by(interaction_type, date_hour) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(data = "synthetic")

interactions_observed = read.csv2(file.path("..", "..", "data", "data-nodscov2", "clean", "interaction_cleaned.csv")) %>%
  mutate(
    interaction_type = case_when(
      grepl("^PA-", from) & grepl("^PA-", to) ~ "PA-PA",
      grepl("^PA-", from) != grepl("^PA-", to) ~ "PA-PE",
      grepl("^PE-", to) & grepl("^PE-", to) ~ "PE-PE"
    ),
    date_hour = floor_date(as_datetime(date_posix), "hour")
  ) %>%
  group_by(interaction_type, date_hour) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(data = "observed")

# Plot
bind_rows(interactions_observed, interactions_synthetic) %>%
  ggplot(., aes(x = data, y = n)) +
  ggh4x::facet_grid2(cols = vars(interaction_type), scales = "free_y", independent  ="y") +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 0.5) +
  theme_bw() +
  labs(x = "Contact dataset", y = "Distribution of the number of interactions per hour")

# Load cluster data 
id_paramedical = admission$id[admission$cat =="Paramedical"]
id_medical = admission$id[admission$cat =="Medical"]

cluster_composition = lapply(seq_along(clusters), function(x) {
  out = sapply(clusters[[x]], function(y) {
    out2 = data.frame(patient = sum(y %in% id_patient), 
                      medical = sum(y %in% id_medical), 
                      paramedical = sum(y %in% id_paramedical),
                      n = length(y))
    return(out2)
  })
  out=data.frame(t(out))
  out$time = x
  return(out)
}
)

# Plot custer composition
p = do.call("rbind", cluster_composition) %>%
  mutate(cluster_type = case_when(
    paramedical > 0 & medical == 0 & patient == 0 ~ "Paramedical",
    paramedical == 0 & medical > 0 & patient == 0 ~ "Medical",
    paramedical == 0 & medical == 0 & patient > 0 ~ "Patient",
    paramedical > 0 & medical == 0 & patient > 0 ~ "Paramed-Patient",
    paramedical == 0 & medical > 0 & patient > 0 ~ "Medical-Patient",
    .default = "All"
  ), 
  time = floor_date(as_datetime(recode(time, !!!dict_times)), "hour")
  ) %>%
  group_by(time, cluster_type) %>%
  summarise(n=n(), .groups = "drop") %>%
  ggplot(., aes(x = time, y = n, col = cluster_type)) +
  geom_line() +
  theme_bw()
