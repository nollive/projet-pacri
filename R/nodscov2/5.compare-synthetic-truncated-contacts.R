################################################################################
##          Parallelized code to reconstruct individual locations
##
################################################################################

# Librairies
library(tidyverse)
library(viridis)
library(igraph)
library(ggnetwork)
library(ggtext)
library(ggpubr)
library(foreach)
library(doParallel)
rm(list = ls())

# Source helper functions
source("R/nodscov2/helper-functions.R")

# Paths
synthetic_path <- file.path("data", "data-synthetic-graphs") 
fig_path = file.path("..", "..", "fig", "network_comparison")

# Network to study
networks = c("herriot", "poincare")
nsim = 1:10
conditions = expand.grid("nsim" = nsim, "network" = networks)

# Truncate data-----------------------------------------------------------------
# registerDoParallel(4)
# foreach (r=1:nrow(conditions)) %dopar% {
# 
#   net = as.character(conditions[r,"network"])
#   i = conditions[r,"nsim"]
# 
#   ## Load data
#   # Admission data
#   admission = read.csv2(file.path(synthetic_path, "input", paste0("admission_", net, ".csv"))) %>%
#     mutate(firstDate = as_date(firstDate), lastDate = as_date(lastDate))
# 
#   # HCW schedule
#   agenda = read.csv2(file.path(synthetic_path, "input", paste0("agenda_", net,".csv"))) %>%
#     mutate(firstDate = as_datetime(firstDate), lastDate = as_datetime(lastDate))
# 
#   # Interaction data
#   interactions = read.csv2(file.path(synthetic_path, "full", net, paste0("matContactBuiltSimulatedCtcNetworks", i, "_oneday.csv"))) %>%
#     mutate(date_posix = as_datetime(date_posix), from_status = gsub("-.*", "", from), to_status = gsub("-.*", "", to))
# 
#   ## Trim interactions to match schedule
#   interactions = interactions %>%
#     mutate(n = 1:n()) %>%
#     nest(.by = n) %>%
#     mutate(data = map(data, trim_interactions_agenda, admission, agenda)) %>%
#     unnest(data) %>%
#     filter(length > 0, before_schedule == F) %>%
#     select(from, to, date_posix, length, from_status, to_status)
#   write.csv2(interactions, file.path(synthetic_path, "full", net, paste0("truncatedContacts", i, ".csv")),
#              row.names = F, quote = F)
# }

# Compare networks--------------------------------------------------------------
## Load real and synthetic contact data
networks = c("poincare", "herriot")

admission = list()
schedule = list()
interaction_real = list()
interaction_synthetic = list()
interaction_truncated = list()

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
  reconstructed_path = file.path(synthetic_path, "full", net)
  interaction_synthetic[[net]] = data.frame()
  
  for (f in list.files(reconstructed_path, pattern = "matContact*", full.names = T)) {
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
  
  # Synthetic truncated data
  interaction_truncated[[net]] = data.frame()
  
  for (f in list.files(reconstructed_path, pattern = "truncatedContacts*", full.names = T)) {
    interaction_truncated[[net]] = bind_rows(
      interaction_truncated[[net]],
      read.csv2(f) %>% 
        mutate(
          nSim = gsub("^.*truncatedContacts|.csv$", "", f),
          firstDate = as_datetime(date_posix),
          lastDate = as_datetime(date_posix) + length
        )
    )
  } 
}


## Compare number of unique contacts and contact duration 
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
  reconstructed_path = file.path(synthetic_path, "full", net)
  for (f in list.files(reconstructed_path, ".*.csv", full.names = T)) {
    
    data = read.csv2(f)
    
    # Contact durations
    if (grepl("_oneday", f)) db_type = "Synthetic"
    if (grepl("truncated", f)) db_type = "Truncated"
    db_type = paste(db_type, gsub(".*Networks|_oneday.*$|.*truncatedContacts|.csv$", "", f))
    
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
    data = factor(data, c("Observed", paste("Synthetic", 1:10), paste("Truncated", 1:10))),
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
  scale_colour_discrete(type = c("orange", colorRampPalette(c("dodgerblue4", "lightskyblue"))(10), colorRampPalette(c("darkmagenta", "mistyrose"))(10))) +
  scale_fill_discrete(type = c("orange", colorRampPalette(c("dodgerblue4", "lightskyblue"))(10), colorRampPalette(c("darkmagenta", "mistyrose"))(10))) +
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
         data = factor(data, c("Observed", paste("Synthetic", 1:10), paste("Truncated", 1:10)))) %>%
  ggplot(., aes(y=length, x = network, colour = data)) +
  geom_boxplot(outlier.shape = NA) +
  scale_colour_discrete(type = c("orange", colorRampPalette(c("dodgerblue4", "lightskyblue"))(10), colorRampPalette(c("darkmagenta", "mistyrose"))(10))) +
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
ggsave(file.path(fig_path, "unique_contacts.png"), p, height = 8, width = 11)  


## Get summary metrics on observed and reconstructed data as well as two ICU networks
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
  reconstructed_path = file.path(synthetic_path, "full", net)
  simu_files = list.files(path = reconstructed_path, pattern = ".*.csv", full.names = T)
  simu_data = data.frame()
  
  for(f in simu_files){
    if (grepl("_oneday", f)) db_type = "Synthetic"
    if (grepl("truncated", f)) db_type = "Truncated"
    iter = as.numeric(gsub(".*Networks|_oneday.*$|.*truncatedContacts|.csv$", "", f))
    graph_data = read.csv2(f) %>%
      mutate(date_posix = floor_date(as_datetime(date_posix), "day")) %>%
      filter(date_posix %in% included_days) %>%
      select(from, to, date_posix) %>%
      distinct %>%
      arrange(date_posix)
    
    summary_data = bind_rows(
      summary_data,
      get_net_metrics(graph_data, adm_data, iter, db_type, net)
    )
  }
}

# Get median
summary_data = summary_data %>%
  filter(iter == 1) %>%
  mutate(
    data = factor(data, levels = c("Observed", "Synthetic", "Truncated")),
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

ggsave(file.path(fig_path, "metrics_full.png"), p, height = 5, width = 10)
