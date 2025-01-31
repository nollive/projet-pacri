################################################################################
##
##                  Analysis of the simulated scenarios
##
################################################################################
  
## Working environment----------------------------------------------------------
#### LIBRARIES
library(tidyverse)
library(cowplot)
library(viridis)
library(ggpubr)
library(ggborderline)
library(rstatix)
library(DescTools)
library(RColorBrewer)
library(ggpubr)
rm(list=ls())

## Suppporting functions, dictionaries and variables
source("helper-functions-simulations.R")
source('helper-functions.R')

## Paths
out_path <- file.path('..','..', 'out')
simulations_path <- file.path('..','..', 'out', 'final-simulations')
fig_path <- file.path('..','..', 'fig', 'simulations')

## Load data--------------------------------------------------------------------
# LOAD DATA
load(file.path(out_path, "parameters-synthetic-data.rda"))

## IDS BY CATEGORY
id_patient <- admission$id[admission$cat == "Patient"]
id_hcw <- admission$id[admission$cat != "Patient"]
id_paramedical <- admission$id[admission$cat == "Paramedical"]
id_medical <- admission$id[admission$cat == "Medical"]

## Selected parameter pairs
# Conditions to load
conditions = c("sim_1-4_1-12000", "sim_3-4_1-18000", "sim_9-10_1-20000", "sim_5-4_1-40000", "sim_6-4_1-60000")
dict_scenarios = paste("Scenario", 1:length(conditions))
names(dict_scenarios) = gsub("sim_", "", conditions)

# list_dirs = list.dirs(path = simulations_path, recursive = T, full.names = T)
# list_dirs = list_dirs[sapply(list_dirs, function(y) {any(sapply(paste0(conditions, "$"), grepl, x = y))})]
# 
# ## READ ALL SIMULATION FOR ALL AVAILABLE COUPLES
# list_sim <- list()
# for (dir in list_dirs){
#   dir_name = gsub(".*simulations/", "", dir)
#   dir_name = gsub("/", "_", dir_name)
#   list_sim[[dir_name]] <- list()
#   rdata_files <- list.files(path = dir, pattern = "^sim_.*\\.rda$", full.names = F)
#   invisible(lapply(rdata_files, function(file) load_rdata_to_list(file, list_sim, dir)))
# }
# save(list_sim, file = file.path(simulations_path, "list_sim.rda"))
load(file.path(simulations_path, "list_sim.rda"))

# # Get SAR
# list_SAR <- get_all_SAR(list_sim = list_sim, 
#                         id_patient = id_patient, 
#                         id_hcw = id_hcw, 
#                         id_paramedical = id_paramedical, 
#                         id_medical = id_medical)
# 
# all_SAR = do.call("rbind", lapply(list_SAR, function(x) do.call("rbind", x)))
# rownames(all_SAR) = NULL
# write.csv2(all_SAR, file.path(simulations_path, "all_SAR.csv"), quote = F, row.names = F)
all_SAR = read.csv2(file.path(simulations_path, "all_SAR.csv"))

# ## Get epidemic duration 
# list_epidemic_duration <- get_all_epidemic_duration(list_sim = list_sim)
# all_durations = do.call("rbind", lapply(list_epidemic_duration, function(x) do.call("rbind", x))) %>%
#   mutate(model = gsub("_.*$", "", couple), scenario = recode(gsub(".*sim_", "", couple), !!!dict_scenarios)) %>%
#   select(-couple)
# rownames(all_durations) = NULL
# write.csv2(all_durations, file.path(simulations_path, "all_durations.csv"), quote = F, row.names = F)
all_durations = read.csv2(file.path(simulations_path, "all_durations.csv"))

## Plot probability of infection for the short-range transmission route--------- 
# Info 
beta_c = c(1/4, 3/4, 9/10, 5/4, 6/4)
scaling_factor = 30/(24*60*60)
days = 1:10

dict_scenario_beta_c = paste("Scenario", 1:length(beta_c))
names(dict_scenario_beta_c) = beta_c

# Plot
p = expand.grid(b_c = beta_c, day = days) %>%
  mutate(
    p = 1-exp(-b_c*scaling_factor)^(day*24*60*2),
    scenario = recode(b_c, !!!dict_scenario_beta_c)
  ) %>%
  ggplot(., aes(x = factor(day), y = p, col = scenario, group = scenario)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(n = 5, name = "Set2")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(y = "Probability of infection by contact after x days", x = "Exposure in days")
ggsave(file.path(fig_path, "transmission_proba_contact.png"), p, height = 4, width = 5)

## Plots with all models--------------------------------------------------------
### Plot SAR
all_SAR = read.csv2(file.path(simulations_path, "all_SAR.csv"))
dict_scenarios = c(
  "1-4_1-12000" = "Scenario 1",
  "3-4_1-18000" = "Scenario 2",
  "5-4_1-40000" = "Scenario 3",
  "7-4_1-100000" = "Scenario 4"
)

## Plot global SAR
SAR_test = all_SAR %>%
  mutate(scenario = recode(paste0(beta_c, "_", beta_e), !!!dict_scenarios)) %>%
  group_by(scenario) %>%
  wilcox_test(Global ~ model) %>% 
  add_xy_position(x = "supp")

all_SAR %>%
  mutate(
    scenario = recode(paste0(beta_c, "_", beta_e), !!!dict_scenarios)
  ) %>%
  ggplot(., aes(x = scenario, y = Global, col = model)) +
  geom_boxplot(outliers = F, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  labs(y = "Global secondary attack rate") #+
# stat_pvalue_manual(SAR_test, label = "p") +
# scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

## Plot SAR by individual category
all_SAR %>%
  mutate(scenario = recode(paste0(beta_c, "_", beta_e), !!!dict_scenarios)) %>%
  select(scenario, model, Patient, Paramedical, Medical) %>%
  pivot_longer(c(Patient, Paramedical, Medical), values_to = "value", names_to = "Category") %>%
  ggplot(., aes(x = Category, y = value, col = scenario)) +
  facet_grid(rows = vars(model)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  labs(y = "Secondary attack rate by individual category")

## Plot SAR by source of infection - Patients
all_SAR %>%
  mutate(scenario = recode(paste0(beta_c, "_", beta_e), !!!dict_scenarios)) %>%
  select(scenario, model, matches("Patient_")) %>%
  pivot_longer(matches("Patient_"), values_to = "value", names_to = "Source") %>%
  mutate(Source = gsub("Patient_", "", Source)) %>%
  ggplot(., aes(x = Source, y = value, col = scenario)) +
  facet_grid(rows = vars(model)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  labs(y = "Secondary attack rate by individual category")

all_SAR %>%
  filter(Patient > 0) %>%
  mutate(
    scenario = recode(paste0(beta_c, "_", beta_e), !!!dict_scenarios),
    ratio = Patient_Contact / Patient_Environment
  ) %>%
  ggplot(., aes(x = scenario, y = ratio)) +
  facet_grid(rows = vars(model)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  labs(y = "Ratio of contact / environment transmission")

## Plot SAR by source of infection - Paramedical
all_SAR %>%
  mutate(scenario = recode(paste0(beta_c, "_", beta_e), !!!dict_scenarios)) %>%
  select(scenario, model, matches("Paramedical_")) %>%
  pivot_longer(matches("Paramedical_"), values_to = "value", names_to = "Source") %>%
  mutate(Source = gsub("Paramedical_", "", Source)) %>%
  ggplot(., aes(x = Source, y = value, col = scenario)) +
  facet_grid(rows = vars(model)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  labs(y = "Secondary attack rate among paramedical staff")


## Plot SAR by source of infection - Medical
all_SAR %>%
  mutate(scenario = recode(paste0(beta_c, "_", beta_e), !!!dict_scenarios)) %>%
  select(scenario, model, matches("^Medical_")) %>%
  pivot_longer(matches("^Medical_"), values_to = "value", names_to = "Source") %>%
  mutate(Source = gsub("^Medical_", "", Source)) %>%
  ggplot(., aes(x = Source, y = value, col = scenario)) +
  facet_grid(rows = vars(model)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  labs(y = "Secondary attack rate among medical staff")

### Plots for linear model------------------------------------------------------
#### SAR 
comparisons = c("Scenario 1_Scenario 2", "Scenario 2_Scenario 3", "Scenario 3_Scenario 4", "Scenario 4_Scenario 5")

# Extract linear model results
linear_SAR = read.csv2(file.path(simulations_path, "all_SAR.csv")) %>%
  filter(model == "linear") %>%
  mutate(scenario = recode(paste0(beta_c, "_", beta_e), !!!dict_scenarios)) %>%
  filter(scenario %in% dict_scenarios) %>%
  arrange(scenario)

# Create composite figure 
# Global SAR
global_test = linear_SAR %>% 
  wilcox_test(Global ~ scenario) %>%
  add_xy_position(x = "scenario") %>%
  filter(p.adj <= 0.05)

p1 = ggboxplot(linear_SAR, x = "scenario", y = "Global", fill = "scenario") +
  stat_pvalue_manual(global_test, label = "p.adj") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Set2")) +
  geom_jitter(size = 0.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  labs(y = "Global SAR")

# SAR stratified by source of infection
linear_SAR_source = linear_SAR %>%
  select(scenario, Environment, Contact) %>%
  pivot_longer(cols = c(Environment, Contact), names_to = "Source", values_to = "SAR") %>%
  arrange(Source)

source_test = linear_SAR_source %>%
  group_by(scenario) %>%
  wilcox_test(SAR ~ Source) %>%
  add_xy_position(x = "Source") %>%
  add_significance() %>%
  filter(p <= 0.05)

p2 = ggboxplot(linear_SAR_source, x = "Source", y = "SAR", fill = "Source", facet.by = "scenario", ncol = 5) +
  stat_pvalue_manual(source_test, label = "p.signif") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.6))) +
  scale_fill_manual(values = env_pal) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(hjust = 1, angle = 30),
        legend.position = "none") +
  labs(y = "SAR by source")

p = ggarrange(p1, p2, nrow = 2)
ggsave(file.path(fig_path, "sar.png"), p, height = 5, width = 6)

## SAR stratified by individual category and source-----------------------------
# SAR for patients
# patient_test = linear_SAR %>% 
#   wilcox_test(Patient ~ scenario) %>%
#   add_xy_position(x = "scenario") %>%
#   filter(p.adj <= 0.05, paste0(group1, "_", group2) %in% comparisons)
p_patient = data.frame(
  p = JonckheereTerpstraTest(split(linear_SAR$Patient, linear_SAR$scenario), alternative = "increasing", exact = F)$p.value,
  scenario = "Scenario 3",
  Patient = 0.4
) %>%
  mutate(p = paste0("increasing trend p-value = ", signif(p, digits = 2)))

p1 = #ggboxplot(linear_SAR, x = "scenario", y = "Patient", fill = "scenario") +
  #stat_pvalue_manual(patient_test, label = "p.adj.signif") +
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ggplot(linear_SAR, aes(x = scenario, y = Patient, fill = scenario)) +
  geom_boxplot() +
  geom_text(data = p_patient, aes(label = p)) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Set2")) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(y = "SAR among patients")

p2 = linear_SAR %>%
  select(scenario, matches("^Patient_")) %>%
  pivot_longer(matches("^Patient_"), names_to = "Source", values_to = "value") %>%
  mutate(Source = gsub("Patient_", "", Source)) %>%
  ggplot(., aes(x = Source, y = value, col = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1)) +
  scale_color_manual(values = brewer.pal(n = 5, name = "Set2")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.title = element_blank(), legend.position = "bottom") +
  labs(y = "SAR among patients")

# SAR for paramedical staff
# paramedical_test = linear_SAR %>% 
#   wilcox_test(Paramedical ~ scenario) %>%
#   add_xy_position(x = "scenario") %>%
#   filter(p.adj <= 0.05, paste0(group1, "_", group2) %in% comparisons)
p_paramedical = data.frame(
  p = JonckheereTerpstraTest(split(linear_SAR$Paramedical, linear_SAR$scenario), alternative = "decreasing", exact = F)$p.value,
  scenario = "Scenario 3",
  Paramedical = 0.9
) %>%
  mutate(p = paste0("decreasing trend p-value = ", signif(p, digits = 2)))

p3 = #ggboxplot(linear_SAR, x = "scenario", y = "Paramedical", fill = "scenario") +
  # stat_pvalue_manual(paramedical_test, label = "p.adj.signif") +
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ggplot(linear_SAR, aes(x = scenario, y = Paramedical, fill = scenario)) +
  geom_boxplot() +
  geom_text(data = p_paramedical, aes(label = p)) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Set2")) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(y = "SAR among paramedical staff")

p4 = linear_SAR %>%
  select(scenario, matches("^Paramedical_")) %>%
  pivot_longer(matches("^Paramedical_"), names_to = "Source", values_to = "value") %>%
  mutate(Source = gsub("Paramedical_", "", Source)) %>%
  ggplot(., aes(x = Source, y = value, col = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1)) +
  scale_color_manual(values = brewer.pal(n = 5, name = "Set2")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  labs(y = "SAR among paramedical staff")

# SAR for medical staff
# Medical_test = linear_SAR %>% 
#   wilcox_test(Medical ~ scenario) %>%
#   add_xy_position(x = "scenario") %>%
#   filter(p.adj <= 0.05, paste0(group1, "_", group2) %in% comparisons)
p_medical = data.frame(
  p = JonckheereTerpstraTest(split(linear_SAR$Medical, linear_SAR$scenario), alternative = "decreasing", exact = F)$p.value,
  scenario = "Scenario 3",
  Medical = 0.9
) %>%
  mutate(p = paste0("decreasing trend p-value = ", signif(p, digits = 2)))

p5 = # ggboxplot(linear_SAR, x = "scenario", y = "Medical", fill = "scenario") +
  # stat_pvalue_manual(Medical_test, label = "p.adj.signif") +
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ggplot(linear_SAR, aes(x = scenario, y = Medical, fill = scenario)) +
  geom_boxplot() +
  geom_text(data = p_medical, aes(label = p)) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Set2")) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(y = "SAR among medical staff")

p6 = linear_SAR %>%
  select(scenario, matches("^Medical_")) %>%
  pivot_longer(matches("^Medical_"), names_to = "Source", values_to = "value") %>%
  mutate(Source = gsub("^Medical_", "", Source)) %>%
  ggplot(., aes(x = Source, y = value, col = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1)) +
  scale_color_manual(values = brewer.pal(n = 5, name = "Set2")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  labs(y = "SAR among medical staff")

# Combine all plots
p = ggarrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2, legend = F)
p = ggarrange(p, get_legend(p2), nrow = 2, heights = c(1,0.1))
ggsave(file.path(fig_path, "sar_stratified.png"), height = 8, width = 8)

## Epidemic duration, epidemic curves, and time to peak-------------------------
# Plot of epidemic curves
all_curves = get_all_SEIR(list_sim, t_end)
all_curves = do.call("rbind", lapply(all_curves, function(x) do.call("rbind", x))) %>%
  mutate(model = gsub("_.*", "", couple), scenario = recode(gsub(".*sim_", "", couple), !!!dict_scenarios)) %>%
  filter(model == "linear") 

dict_time = seq(new_begin_date, new_end_date, by = 3000)
names(dict_time) = seq(1, t_end, 100)

p1 = all_curves %>%
  filter(status == "Infectious") %>%
  mutate(time = recode(time, !!!dict_time)) %>%
  ggplot(., aes(x = time, y = count, group = id_sim)) +
  geom_line(col = "grey80", linewidth = 0.1) +
  facet_grid(cols = vars(scenario)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(y = "Number of infectious")

# Wilcoxon test to compare the distribution across scenarios (p-values adjusted for multiple testing)
p_duration = all_durations %>%
  filter(model == "linear") %>%
  wilcox_test(duration ~ scenario) %>%
  add_xy_position(x = "scenario") %>%
  filter(p.adj <= 0.05)

# Plot of epidemic duration
p2 = all_durations %>%
  filter(model == "linear") %>%
  arrange(scenario) %>%
  ggboxplot(., x = "scenario", y = "duration", fill = "scenario") +
  stat_pvalue_manual(p_duration, label = "p.adj.signif") +
  geom_jitter(size = 0.8) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Set2")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  labs(y = "Epidemic duration (in days)")

# Plot of time to peak
time_to_peak = all_curves %>%
  nest(.by = c(scenario, id_sim)) %>%
  mutate(time_to_peak = map(data, get_time_to_peak, dict_time)) %>%
  unnest(time_to_peak) 

p_time_to_peak = data.frame(
  p = format(signif(JonckheereTerpstraTest(split(time_to_peak$time_to_peak, time_to_peak$scenario), alternative = "decreasing", exact = F)$p.value, 2), scientific = T),
  scenario = "Scenario 3",
  time_to_peak = 90
)

p3 = time_to_peak %>%
  ggplot(., aes(x = scenario, y = time_to_peak, fill = scenario)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 0.8) +
  geom_text(data = p_time_to_peak, aes(label = paste0("Decreasing trend p = ", p))) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Set2")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  labs(y = "Time to peak (in days)")

# Combine plots
p = ggarrange(p1, p2, p3, nrow = 3)
ggsave(file.path(fig_path, "epidemic_timing.png"), p, height = 6, width = 6)

## Example of aerosol concentration dynamics------------------------------------
# Simulations to load
load("../../out/final-simulations/linear/sim_1-4_1-12000/sim_1-4_1-12000_1.rda")
all_simulations = do.call("rbind", `sim_1-4_1-12000_1`[["global_environment"]]) %>%
  mutate(time = rep(1:259263, each = 21))

# First day
start_cut = 1
end_cut = start_cut + 24*60*2 - 1
dict_time = seq(new_begin_date, new_begin_date+24*3600-30, by = 30)
names(dict_time) = start_cut:end_cut
all_dates_locations = expand.grid(room = unname(dict_rooms), time = unname(dict_time))

first_day_infectious = `sim_1-4_1-12000_1`[["global_status"]] %>%
  filter((start_cut >= t_incub & start_cut <= t_recover) | 
           (end_cut >= t_incub & end_cut <= t_recover) ) %>%
  mutate(k = 1:n()) %>%
  nest(.by = k) %>%
  mutate(data = map(data, get_location_during_schedule, start_cut, end_cut)) %>%
  unnest(data) %>%
  mutate(
    room = recode(room, !!!dict_rooms),
    time = recode(time, !!!dict_time)
  ) %>%
  group_by(time, room) %>%
  summarise(n = n(), .groups = "drop") %>%
  full_join(., all_dates_locations, by = c("room", "time")) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(room %in% c("Corridor", "Nursing Station", "Office", "Medical Staff Room", "Paramedical Staff Room", "Patient Room 1", "Patient Room 6"))


p1 = all_simulations %>% 
  filter(time <= end_cut, time >= start_cut, id_room %in% c(1, 6, 18:22)) %>%
  mutate(
    room = recode(id_room, !!!dict_rooms),
    time = recode(time, !!!dict_time)
  ) %>%
  mutate(
    room = factor(room, c("Corridor", "Nursing Station", "Office", "Medical Staff Room", "Paramedical Staff Room", paste("Patient Room", c(1:6,8:17))))
  ) %>%
  ggplot(., aes(x = time, y = env/volume)) +
  facet_grid(cols = vars(room)) +
  scale_x_datetime(labels = scales::date_format("%H:%M")) +
  geom_line(data = first_day_infectious, aes(y = n*1e5), col = "red") +
  geom_line() +
  theme_bw() +
  scale_y_continuous(
    "Aerosol concentration (per m3)", 
    sec.axis = sec_axis(~ . /1e5, name = "Number of infectious")
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "First day")

# One week later
set.seed("123")
random_day = sample(seq(new_begin_date+24*3600*7, new_end_date, 24*3600), 1)
start_cut = as.integer(difftime(random_day, new_begin_date, units = "secs")/30)
end_cut = start_cut + 24*60*2 - 1
dict_time = seq(random_day, random_day+24*3600-30, by = 30)
names(dict_time) = start_cut:end_cut
all_dates_locations = expand.grid(room = unname(dict_rooms), time = unname(dict_time))

random_day_infectious = `sim_1-4_1-12000_1`[["global_status"]] %>%
  filter((start_cut >= t_incub & start_cut <= t_recover) | 
           (end_cut >= t_incub & end_cut <= t_recover)) %>%
  mutate(k = 1:n()) %>%
  nest(.by = k) %>%
  mutate(data = map(data, get_location_during_schedule, start_cut, end_cut)) %>%
  unnest(data) %>%
  mutate(
    room = recode(room, !!!dict_rooms),
    time = recode(time, !!!dict_time)
  ) %>%
  group_by(time, room) %>%
  summarise(n = n(), .groups = "drop") %>%
  full_join(., all_dates_locations, by = c("room", "time")) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(room %in% c("Corridor", "Nursing Station", "Office", "Medical Staff Room", "Paramedical Staff Room", "Patient Room 1", "Patient Room 6"))

p2 = all_simulations %>% 
  filter(time <= end_cut, time >= start_cut, id_room %in% c(1, 6, 18:22)) %>%
  mutate(
    room = recode(id_room, !!!dict_rooms),
    time = recode(time, !!!dict_time)
  ) %>%
  mutate(
    room = factor(room, c("Corridor", "Nursing Station", "Office", "Medical Staff Room", "Paramedical Staff Room", paste("Patient Room", c(1:6,8:17))))
  ) %>%
  ggplot(., aes(x = time, y = env/volume)) +
  facet_grid(cols = vars(room)) +
  scale_x_datetime(labels = scales::date_format("%H:%M")) +
  geom_line(data = random_day_infectious, aes(y = n*2.5e5), col = "red") +
  geom_line() +
  scale_y_continuous(
    "Aerosol concentration (per m3)", 
    labels = scales::label_scientific(),
    sec.axis = sec_axis(~ . /2.5e5, name = "Number of infectious")
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "Random day one week after the start of the epidemic")

# Combine plots
p = ggarrange(p1, p2, nrow = 2)
ggsave(file.path(fig_path, "dynamics_aerosols.png"), p, height = 6, width = 12)
