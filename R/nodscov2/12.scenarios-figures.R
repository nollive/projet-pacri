################################################################################
##
##                  Figures of scenarios comparison
##
################################################################################

## Working environment----------------------------------------------------------
## Libraries
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(lubridate)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)
library(ggborderline)
rm(list = ls())

## Suppporting functions, dictionaries and variables
source("scenarios-figures-supp.R")
source('scenarios-analysis_fun.R')

##PATHS
loc_path <- file.path('..','..', 'out','loc-nodscov2')
scenarios_path <- file.path('..','..', 'out', 'scenarios-analysis')

# LOAD DATA
load(file.path(scenarios_path, "parameters-admission-nodscov2.RData"))

## IDS BY CATEGORY
id_patient <- admission_sim %>% filter(info == 0) %>% distinct(id) %>% pull()
id_hcw <- admission_sim %>% filter(info == 1) %>% distinct(id) %>% pull()

## FOR MEDICAL/PARAMEDICAL WE USE ADMISSION (INDIVIDUALS BEFORE ADDING NEW PATIENTS)
id_paramedical <- admission %>% filter(cat %in% cat_paramedical) %>% distinct(id) %>% pull()
id_medical <- admission %>% filter(cat %in% cat_medical) %>% distinct(id) %>% pull()

## ADJUST THE BEGIN/END DATES OF SIMULATIONS
# new_begin_date <- begin_date + t_begin*30 ## OFFSET TO START THE DAY AT 12AM
# new_end_date <- new_begin_date + new_n_subdivisions*30
# print(new_end_date - new_begin_date)


#### WHEN COUPLES ARE SELECTED--------------------------------------------------
all_metrics <- read.csv2(file = file.path(scenarios_path, 'scenarios-all_metrics.csv'))
couple_same_SAR <- all_metrics %>% filter(dplyr::between(Global_median, 0.20, 0.30)) %>% pull(couple)

##SELECT COUPLES FOR SAR +/- 20%
# all_SAR_g_metrics %>% filter(dplyr::between(SAR_median, 0.10, 0.30))
# all_SAR_g_metrics %>% filter(dplyr::between(SAR_mean, 0.10, 0.30)) %>% select(beta, nu , SAR_mean, SAR_median, SAR_sd)


### READ ONLY RDATA FOR SELECTED BETA/NU COUPLES
list_sim <- list()
for (couple in couple_same_SAR){
  list_sim[[couple]] <- list()
  sim_path <- file.path(scenarios_path, 'scenarios-simulations', couple, 'results')
  rdata_files <- list.files(sim_path, pattern = "^sim_.*\\.RData$", full.names = F)
  invisible(lapply(rdata_files, function(file) load_rdata_to_list(file, list_sim, couple)))
}

## COMPUTE ALL THE SAR FOR EACH COUPE/SIMULATION
list_SAR <- get_all_SAR(list_sim = list_sim, 
                        id_patient = id_patient, 
                        id_hcw = id_hcw, 
                        id_paramedical = id_paramedical, 
                        id_medical = id_medical)

list_epidemic_duration <- get_all_epidemic_duration(list_sim = list_sim)
all_epidemic_duration_metrics <- get_all_epidemic_duration_metrics(list_epidemic_duration = list_epidemic_duration)

## CONCATENATE THE NSIMULATION DATAFRAMES  FOR 1 COUPLE TO 1 DATAFRAME FOR 1 COUPLE
list_df <- list()
for(couple in couple_same_SAR){
  list_df[[couple]] <- rbindlist(list_SAR[[couple]], idcol = 'id_sim') %>% mutate(couple = couple)
}

df_SAR <- rbindlist(list_df) 
n_ext <- df_SAR %>% group_by(couple) %>% filter(Global == 0) %>% summarise(n = n())


df_SAR_long <- df_SAR %>%
  pivot_longer(cols = c("Global", "HCW", "Patient", "Paramedical", "Medical", "Environment", "Contact" , 
                        "Patient_Contact", "Patient_Environment", "Paramedical_Contact", "Paramedical_Environment",
                        "Medical_Contact", "Medical_Environment"), names_to = "Type_SAR", values_to = "SAR")
df_SAR_long$Type_SAR <- factor(df_SAR_long$Type_SAR, levels = c("Global", "HCW", "Patient", "Paramedical", "Medical", "Environment", "Contact",  "Patient_Contact", "Patient_Environment", "Paramedical_Contact", "Paramedical_Environment", "Medical_Contact", "Medical_Environment"))

selected_couples <- c(paste(rep('sim', 5), c('1-4_20', '1-2_16', '3-4_12', '1_8', '3-2_5'), sep = '_'))
selected_df <-  df_SAR_long %>% filter(couple %in% selected_couples) 
selected_df$couple <- factor(selected_df$couple, levels = selected_couples)


## t TEST/ Wilcoxon test (significative difference between SAR contact/env)-----
for (id_couple in selected_couples){
  print(id_couple)
  
  print(t.test(x = df_SAR %>% filter(couple == id_couple) %>% pull(Contact),
               y = df_SAR %>% filter(couple == id_couple) %>% pull(Environment),
               paired = T))
  
  print(wilcox.test(x = df_SAR %>% filter(couple == id_couple) %>% pull(Contact),
                    y = df_SAR %>% filter(couple == id_couple) %>% pull(Environment),
                    paired = T))
  
}

## PLOTS------------------------------------------------------------------------
p_all_c <- df_SAR_long %>%
  ggplot( aes(x=couple, y=SAR, fill=couple)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggh4x::facet_grid2(cols = vars(Type_SAR), scales = "free_y", independent  ="y") +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p_selected_c <- selected_df %>%
  ggplot( aes(x=couple, y=SAR, fill=couple)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggh4x::facet_grid2(cols = vars(Type_SAR), scales = "free_y", independent  ="y") +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

print(p_all_c)
print(p_selected_c)

## FINAL PLOTS - CONTRIBUTION---------------------------------------------------
selected_df <- selected_df %>%
  mutate(scenario = recode(couple, !!!dict_scenarios))
selected_df$scenario <- factor(selected_df$scenario, levels = c(paste0(rep('Scenario ',5), 1:5)))

p1 = selected_df %>%
  filter(Type_SAR %in% c("Contact", "Environment")) %>%
  ggplot( aes(x=scenario, y=SAR, col = Type_SAR)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.5) +
  geom_jitter(size=0.4, alpha=0.9, position = position_jitterdodge()) +
  scale_color_manual(values = env_pal) +
  theme_bw() +
  theme(
    #axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(size=11)
  ) + 
  labs(x = "", col = "")


p2 = selected_df %>%
  filter(Type_SAR %in% c("Patient", "Paramedical", "Medical")) %>%
  mutate(Type_SAR = factor(Type_SAR, c("Paramedical", "Medical", "Patient"))) %>%
  ggplot( aes(x=scenario, y=SAR, col = Type_SAR)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.5) +
  geom_jitter(size=0.4, alpha=0.9, position = position_jitterdodge()) +
  scale_color_manual(values = pal[1:3]) +
  theme_bw() +
  theme(
    #axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(size=11)
    
  ) + 
  labs(x = "", col = "") 


p3 = selected_df %>%
  filter(Type_SAR %in% c("Patient", "Paramedical", "Medical")) %>%
  pivot_wider(names_from = Type_SAR, values_from = SAR) %>%
  mutate(
    Patient = Patient * length(id_patient),
    Paramedical = Paramedical * length(id_paramedical), 
    Medical = Medical * length(id_medical)
  ) %>%
  pivot_longer(cols =  c(Patient, Paramedical, Medical), names_to = "Type_SAR", values_to = 'n') %>%
  group_by(scenario, Type_SAR) %>%
  summarise(n = mean(n), .groups = "drop") %>%
  group_by(scenario) %>%
  mutate(n = n/sum(n)) %>%
  ungroup() %>%
  ggplot(., aes(x = scenario, y = n, fill = Type_SAR)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  #theme(axis.text.x = element_text(hjust = 1, angle = 30)) +
  scale_fill_manual(values = pal[1:3]) +
  labs(x  = "", y = "Proportion of infected\nindividuals by type", fill = "")


p4 = selected_df %>%
  filter(grepl("_", Type_SAR)) %>% # Pour avoir uniquement les SAR contact + environment pour chaque type d'individu
  mutate(
    n = case_when(
      grepl("Patient", Type_SAR) ~ SAR * length(id_patient),
      grepl("Paramedical", Type_SAR) ~ SAR * length(id_paramedical),
      grepl("Medical", Type_SAR) ~ SAR * length(id_medical)
    ),
    ind_type = gsub("_.*$", "", Type_SAR),
    source_inf = gsub("^.*_", "", Type_SAR)
  ) %>%
  group_by(scenario, Type_SAR, ind_type, source_inf) %>%
  summarise(n = mean(n), .groups = "drop") %>%
  group_by(scenario, ind_type) %>%
  mutate(n = n/sum(n)) %>%
  ungroup() %>%
  ggplot(., aes(x = scenario, y = n, fill = source_inf)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = env_pal) +
  facet_grid(cols = vars(ind_type)) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 30)) +
  labs(x  = "", y = "Contribution of modes of transmission", fill = "")


p_contrib <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), nrow = 2, ncol = 2, align = "hv")
ggsave(filename = file.path(scenarios_path, 'fig', 'scenario-contribution.png'), plot = p_contrib, height = 7, width = 14)


## Plots of epidemic duration---------------------------------------------------
# Plot
epi_durations = do.call(
  "rbind", 
  lapply(list_epidemic_duration[selected_couples], function(x) do.call("rbind", x))
) %>%
  mutate(scenario = recode(couple, !!!dict_scenarios))

ggplot(epi_durations, aes(x = scenario, y = duration)) +
  geom_boxplot() +
  geom_jitter() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey") +
  geom_hline(yintercept = 90, linetype = "dashed", col = "grey") +
  theme_bw() +
  labs(x = "", y = "Epidemic duration (in days)")

# ANOVA to compare epidemic duration 
epi_durations_aov = aov(duration ~ scenario, data = epi_durations)
summary(epi_durations_aov) # No significant difference

# Get contributions
selected_scenarios_output = list_sim[selected_couples]
all_info = data.frame()
for (i in seq_along(selected_scenarios_output)) {
  for (j in seq_along(selected_scenarios_output[[i]])) {
    all_info = bind_rows(all_info, 
                         get_ss(selected_scenarios_output[[i]][[j]], 
                                names(selected_scenarios_output)[[i]], 
                                j))
  }
}

# Plot contributions
relative_contribution = all_info %>%
  mutate(
    scenario = recode(couple, !!!dict_scenarios),
    r_c = ifelse(ninf_c+ninf_e == 0, 0, ninf_c/(ninf_c + ninf_e)*100) 
  ) 

ggplot(relative_contribution, aes(x = scenario, y = r_c)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  labs(x = "", y = "Relative contribution of contacts over environment")

# Wilcoxon test 
mean(relative_contribution$r_c[relative_contribution$scenario == "Scenario 1"])
mean(relative_contribution$r_c[relative_contribution$scenario == "Scenario 5"])
wilcox.test(relative_contribution$r_c[relative_contribution$scenario == "Scenario 1"],
            relative_contribution$r_c[relative_contribution$scenario == "Scenario 5"])


# Patient SAR
all_info$ninf_patients = all_info$ninf_patients_c + all_info$ninf_patients_e
wilcox.test(all_info$ninf_patients[all_info$couple == "sim_1-4_20"],
            all_info$ninf_patients[all_info$couple == "sim_3-2_5"])
c(mean(all_info$ninf_patients[all_info$couple == "sim_1-4_20"]),
  mean(all_info$ninf_patients[all_info$couple == "sim_3-2_5"])
)

all_info$ninf_para = all_info$ninf_para_c + all_info$ninf_para_e
wilcox.test(all_info$ninf_para[all_info$couple == "sim_1-4_20"],
            all_info$ninf_para[all_info$couple == "sim_3-2_5"])
c(mean(all_info$ninf_para[all_info$couple == "sim_1-4_20"]),
  mean(all_info$ninf_para[all_info$couple == "sim_3-2_5"]))


all_info$ninf_med = all_info$ninf_med_c + all_info$ninf_med_e
wilcox.test(all_info$ninf_med[all_info$couple == "sim_1-4_20"],
            all_info$ninf_med[all_info$couple == "sim_3-2_5"])

# Peak time
all_info %>%
  mutate(scenario = recode(couple, !!!dict_scenarios)) %>%
  ggplot(., aes(x = scenario, y = peak_time)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  labs(x = "", y = "Peak time")
peak_time_aov = aov(peak_time ~ couple, data = all_info)
summary(peak_time_aov) # No significant difference

# Peak time
all_info %>%
  mutate(scenario = recode(couple, !!!dict_scenarios)) %>%
  ggplot(., aes(x = scenario, y = epidemic_duration)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  labs(x = "", y = "Epidemic duration")
epidemic_duration_aov = aov(epidemic_duration ~ couple, data = all_info)
summary(epidemic_duration_aov) # No significant difference


all_SAR_metrics <- compute_all_SAR_metrics(list_SAR = list_SAR)
list_gs_df <- list()
for(couple in selected_couples){
  list_gs_df[[couple]] <- rbindlist(list_sim[[couple]], idcol = 'id_sim') %>% mutate(couple = couple)
}

df_status <- rbindlist(list_gs_df) 
test = df_status %>%
  filter(t_inf > 0) %>%
  mutate(
    source = gsub("-.*$", "", inf_by),
    type_ind = case_when(
      id %in% id_patient ~ "Patient",
      id %in% id_paramedical ~ "Paramedical",
      id %in% id_medical ~ "Medical"
    ),
    # start_exposed = as.Date(new_begin_date + t_inf * 30), 
    # end_exposed = as.Date(new_begin_date + (t_inf + t_incub) * 30)#,
    start_infected = as.Date(new_begin_date + seconds((t_inf + t_incub) * 30)),
    end_infected = as.Date(new_begin_date + seconds((t_inf + t_incub + t_recover) * 30))
  ) %>%
  dplyr::select(couple, id_sim, id, type_ind, matches("start_|end_")) %>%
  group_by(couple, id_sim, id, type_ind) %>% 
  reframe(date = seq(start_infected, end_infected, by="1 day")) %>%
  group_by(couple, id_sim, date, type_ind) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(couple, id_sim, date, type_ind, fill = list(n = 0))

test <- test%>%
  mutate(scenario = case_when(
    couple == 'sim_1-4_20'~ 'Scenario 1',
    couple == 'sim_1-2_16'~ 'Scenario 2',
    couple == 'sim_3-4_12'~ 'Scenario 3',
    couple == 'sim_1_8'~ 'Scenario 4',
    couple == 'sim_3-2_5'~ 'Scenario 5'))
test$scenario <- factor(test$scenario, levels = c(paste0(rep('Scenario ',5), 1:5)))


test_mean = test %>% 
  group_by(scenario, date, type_ind) %>%
  summarise(n = mean(n), .groups = "drop") %>%
  mutate(id_sim = NA)

p_seir <- test %>%
  ggplot(., aes(x = date, y = n, group = id_sim)) +
  geom_line(linewidth = 0.5, col = "gray80") +
  geom_line(data = test_mean, col = "red") +
  facet_grid(cols = vars(scenario), rows = vars(type_ind)) +
  scale_x_date(date_labels = "%b/%d") +
  labs(
    x = 'Time',
    y = 'Infectious individuals') +
  theme_bw()

print(p_seir)
ggsave(filename = file.path(scenarios_path, 'fig', 'seir-cat.png'), plot = p_seir, width = 14, height = 7)

## without cat facetting

test_global = df_status %>%
  filter(t_inf > 0) %>%
  mutate(
    source = gsub("-.*$", "", inf_by),
    type_ind = case_when(
      id %in% id_patient ~ "Patient",
      id %in% id_paramedical ~ "Paramedical",
      id %in% id_medical ~ "Medical"
    ),
    # start_exposed = as.Date(new_begin_date + t_inf * 30), 
    # end_exposed = as.Date(new_begin_date + (t_inf + t_incub) * 30)#,
    start_infected = as.Date(new_begin_date + seconds((t_inf + t_incub) * 30)),
    end_infected = as.Date(new_begin_date + seconds((t_inf + t_incub + t_recover) * 30))
  ) %>%
  dplyr::select(couple, id_sim, id, type_ind, matches("start_|end_")) %>%
  group_by(couple, id_sim, id) %>% 
  reframe(date = seq(start_infected, end_infected, by="1 day")) %>%
  group_by(couple, id_sim, date) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(couple, id_sim, date, fill = list(n = 0))
test_global <- test_global%>%
  mutate(scenario = case_when(
    couple == 'sim_1-4_20'~ 'Scenario 1',
    couple == 'sim_1-2_16'~ 'Scenario 2',
    couple == 'sim_3-4_12'~ 'Scenario 3',
    couple == 'sim_1_8'~ 'Scenario 4',
    couple == 'sim_3-2_5'~ 'Scenario 5'))
test_global$scenario <- factor(test_global$scenario, levels = c(paste0(rep('Scenario ',5), 1:5)))


test_mean_global = test_global %>% 
  group_by(scenario, date) %>%
  summarise(n = mean(n), .groups = "drop") %>%
  mutate(id_sim = NA)

p_seir_global <- test_global %>%
  ggplot(., aes(x = date, y = n, group = id_sim)) +
  geom_line(linewidth = 0.5, col = "gray80") +
  geom_line(data = test_mean_global, col = "red") +
  facet_grid(cols = vars(scenario)) +
  scale_x_date(date_labels = "%b/%d") +
  labs(
    x = 'Time',
    y = 'Infectious individuals') +
  theme_bw()

print(p_seir_global)
ggsave(filename = file.path(scenarios_path, 'fig', 'seir-global-scenario.png'), plot = p_seir_global, width = 9, height = 3)


### x = beta < 4
ggplot(all_metrics %>% filter(beta < 4), aes(x = beta, y = Global_median, color = as.factor(nu), group = nu)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "darkgrey") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.3, alpha = 0.3, fill = "lightgrey") +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "darkgrey") +
  labs(title = 'SAR Median vs Beta for Different Nu Values',
       x = 'Beta',
       y = 'SAR Median',
       color = 'Nu') +
  theme_minimal()


### x = nu
ggplot(all_metrics, aes(x = nu, y = Global_median, color = as.factor(beta), group = beta)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "darkgrey") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.3, alpha = 0.3, fill = "lightgrey") +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "darkgrey") +
  labs(title = 'SAR Median vs Nu for Different Beta Values',
       x = 'nu',
       y = 'SAR Median',
       color = 'Beta') +
  theme_minimal()

### x = beta < 4 AND GREY ZONE FOR SAR
beta <- all_metrics$beta
p_grid <- ggplot(all_metrics , aes(x = beta, y = Global_median, color = as.factor(nu), group = nu)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "darkgrey") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.3, alpha = 0.3, fill = "lightgrey") +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "darkgrey") +
  scale_x_continuous("Beta", labels = as.character(all_metrics$beta), breaks = all_metrics$beta) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6)) +
  labs(x = 'Beta',
       y = 'SAR Median',
       # title = 'SAR Median vs Beta for Different Nu Values',
       color = 'Nu') +
  theme_bw() + 
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))

print(p_grid)
ggsave(filename = file.path(scenarios_path, 'fig', 'grid-search.png'), plot = p_grid, width = 21, height = 7)

##### HEATMAP
ggplot(all_metrics, aes(x = beta, y = nu, fill = Global_median)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("red", "yellow", "green", "blue")) +
  labs(title = "2D Surface Plot of Grid Search Results",
       x = "Beta",
       y = "Nu",
       fill = "SAR Median") +
  theme_minimal()

## SEIR PLOTS 
couple = 'sim_3-4_12'
list_SEIR <- get_all_SEIR(list_sim, n_subdivisions = new_n_subdivisions)

all_SEIR_metrics <- get_all_SEIR_metrics(list_SEIR = list_SEIR)

#Temp
begin_date <- new_begin_date
plot_SEIR_n(couple = couple, list_SEIR = list_SEIR, all_SEIR_metrics = all_SEIR_metrics)

## VARIABILITY BETWEEN SIMULATION - SEIR, ED, SAR FOR BASELINE
## Baseline scenario SEIR with epidemic duration and SAR
scenario <- 'sim_3-4_12'
sar_scenario <- df_SAR_long %>% filter(couple == scenario) %>% filter(Type_SAR == 'Global')
epidemic_duration <- list_epidemic_duration[[scenario]] %>% rbindlist()
list_sim_scenario <- list()
list_sim_scenario[[scenario]] <-  list_sim[[scenario]]
list_SEIR_scenario <- get_all_SEIR(list_sim_scenario, new_n_subdivisions)
seir_metrics_scenario <- get_all_SEIR_metrics(list_SEIR_scenario)


## VARIABILITY BETWEEN SIMULATION - PLOTS (2 jittered boxplots)
SEIR_colors <- c("Susceptible" = "green", "Exposed" = "pink", "Infectious" = "red", "Recovered" = "blue")
n_sim <- length(list_SEIR_scenario[[scenario]])
title <- paste0('SEIR - ', couple)
counts_list <- list_SEIR_scenario[[scenario]]
couple_id <- couple
metrics_df <- seir_metrics_scenario %>% filter(couple == couple_id)
## SEIR: NUMBER OF INDIVIDUAL
SEIR_n <- ggplot() + 
  labs(x = "Time",
       y = "Number of individuals",
       # title = title,
       color = "Status") +
  scale_color_manual(values = SEIR_colors) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18))

for (i in 1:n_sim){
  SEIR_n <- SEIR_n + 
    geom_line(data = counts_list[[i]],
              aes(x = time * 30 + new_begin_date, y = count, color = status),
              alpha = 0.5,
              linewidth = 0.3,
              linetype = "solid")
}

# SEIR
p1 <- SEIR_n +
  geom_borderline(data = metrics_df,
                  aes(x = time * 30 + new_begin_date, y = count_median , color = status),
                  linewidth = 2,
                  linetype = "solid",
                  bordercolour = "white") +
  geom_line(data = metrics_df,
            aes(x = time * 30 + new_begin_date, y = count_min, color = status),
            linewidth = 0.5,
            linetype = "solid") +
  geom_line(data = metrics_df,
            aes(x = time * 30 + new_begin_date, y = count_max, color = status),
            linewidth = 0.5,
            linetype = "solid")

# SAR
p2 <- ggplot(sar_scenario, aes(x = "", y = SAR, fill = "SAR")) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Boxplot
  geom_jitter(color = "black", size = 0.4, alpha = 0.7) +  # Jittered points
  geom_hline(yintercept = median(sar_scenario$SAR), linetype = "dashed", linewidth = 1, color = "red3") +
  scale_fill_manual(values = c("SAR" = "#1f77b4")) +  
  theme_bw() +
  labs(x = "", y = "SAR (Secondary Attack Rate)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# Epidemic duration
p3 <- ggplot(epidemic_duration, aes(x = "", y = duration, fill = "Epidemic Duration")) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Boxplot
  geom_jitter(color = "black", size = 0.4, alpha = 0.7) +  # Jittered points
  geom_hline(yintercept = median(epidemic_duration$duration), linetype = "dashed", linewidth = 1, color = "red3") +
  scale_fill_manual(values = c("Epidemic Duration" = "#ff7f0e")) +  
  theme_bw() +
  labs(x = "", y = "Epidemic Duration (days)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# p_scenario <- plot_grid(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3, rel_widths = c(3, 0.5, 0.5))
p_scenario <- plot_grid(p2, p3, p1, labels = c("A", "B", "C"), ncol = 3, rel_widths = c(0.5, 0.5, 3))
ggsave(filename = file.path(scenarios_path, 'fig', 'seir-sar-ed-scenario-3.png'), plot = p_scenario, height = 7, width = 21)


## VARIABILITY BETWEEN SIMULATION - PLOTS (2 hstograms)
# SAR plot
p2bis <- ggplot(sar_scenario, aes(x = SAR)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.01, color = "grey", fill = "lightblue") +
  geom_vline(xintercept = median(sar_scenario$SAR), linetype = "dashed", linewidth = 1, color = "red3") +
  theme_bw() +
  labs(x = "SAR (Secondary Attack Rate)", y = "Density") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))

# Epidemic Duration plot
p3bis <- ggplot(epidemic_duration, aes(x = duration)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, color = "grey", fill = "lightblue") +
  geom_vline(xintercept = median(epidemic_duration$duration), linetype = "dashed", linewidth = 1, color = "red3") +
  theme_bw() +
  labs(x = "Epidemic Duration (days)", y = "Density") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))

final_plot <- plot_grid(p2bis, p3bis, p1, labels = c("A", "B", "C"), ncol = 3, rel_widths = c(0.5, 0.5, 3))
ggsave(filename = file.path(scenarios_path, 'fig', 'final_plot.png'), plot = final_plot, height = 7, width = 21)

print(final_plot)
