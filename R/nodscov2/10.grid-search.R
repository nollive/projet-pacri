################################################################################
##
##                  Grid-search procedure
##
################################################################################

## Working environment----------------------------------------------------------
## DATA MANAGEMENT
rm(list=ls())
library(tidyverse)

##PATHS
scenarios_path <- file.path('..','..', 'out', 'grid-search')

## LOAD FUNCTIONS
source('helper-functions-simulations.R')

## Load data--------------------------------------------------------------------
#load(file.path(loc_path, 'dev-localization-nodscov2.RData'))
## Load synthetic contact data
load(file.path("..", "..", "out", "parameters-synthetic-data.rda"))

## Participant ids by category
id_patient <- admission$id[grepl("^PA-", admission$id)]
id_hcw <- admission$id[admission$cat != "Patient"]
id_paramedical <- admission$id[admission$cat == "Paramedical"]
id_medical <- admission$id[admission$cat == "Medical"]

## BETA/NU COUPLES -> READ SIMULATION RESULTS
## ALL COUPLES INSIDE A LIST, EVERY SIMULATION INSIDE SUB-LIST NAMED sim<beta>_<nu> --> WILL BE NAMED A COUPLE
list_dirs = list.dirs(path = scenarios_path, recursive = T, full.names = T)
list_dirs = list_dirs[grepl("sim_", list_dirs)]

## READ ALL SIMULATION FOR ALL AVAILABLE COUPLES
list_sim <- list()
for (dir in list_dirs){
  dir_name = gsub(".*grid-search/", "", dir)
  dir_name = gsub("/", "_", dir_name)
  list_sim[[dir_name]] <- list()
  rdata_files <- list.files(path = dir, pattern = "^sim_.*\\.rda$", full.names = F)
  invisible(lapply(rdata_files, function(file) load_rdata_to_list(file, list_sim, dir)))
}

## FOR EACH SIMULATION -> COMPUTE SAR-------------------------------------------
## GET SARs FOR ALL COUPLES
list_SAR <- get_all_SAR(list_sim = list_sim, 
                        id_patient = id_patient, 
                        id_hcw = id_hcw, 
                        id_paramedical = id_paramedical, 
                        id_medical = id_medical)
all_SAR = do.call("rbind", lapply(list_SAR, function(x) do.call("rbind", x)))

# Grid search values for all models
dict_beta_c = c(
  "0" = "0", 
  "1-4" = "0.25", 
  "1-2" = "0.5", 
  "3-4" = "0.75", 
  "1" = "1", 
  "5-4" = "1.25", 
  "6-4" = "1.5",
  "7-4" = "1.75", 
  "2" = "2", 
  "9-4" = "2.25", 
  "10-4" = "2.5", 
  "11-4" = "2.75", 
  "3" = "3"
)

p = all_SAR %>%
  group_by(model, beta_e, beta_c) %>%
  summarise(median = median(Global), mean = mean(Global), .groups = "drop") %>%
  pivot_longer(c(median, mean), values_to = "SAR", names_to = "metric") %>%
  mutate(
    beta_e = factor(gsub("-", "/", beta_e), c("0", "1/120000", "1/110000", "1/100000", "1/90000", "1/80000", "1/70000", "1/60000", "1/50000", "1/40000", "1/30000", "1/20000", "1/18000", "1/16000", "1/14000", "1/12000", "1/10000", "1/8000", "1/6000", "1/4000", "1/2000")),
    beta_c = recode(beta_c, !!!dict_beta_c),
  ) %>%
  ggplot(., aes(x = beta_e, y = SAR, group = beta_c, col = beta_c)) +
  geom_rect(ymin = 0.2, ymax = 0.3, xmax = Inf, xmin = -Inf, fill = "grey", col = "grey") + 
  facet_grid(cols = vars(model), rows = vars(metric)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  labs(x = expression(beta[e]), y = "Secondary attack rate", col = expression(beta[c]))
ggsave(file.path("..", "..", "fig", "simulations", "grid-search.png"), p, height = 7, width = 14)

# Grid search values for linear model only
p = all_SAR %>%
  filter(model == "linear") %>%
  mutate(
    beta_e = factor(gsub("-", "/", beta_e), c("0", "1/120000", "1/110000", "1/100000", "1/90000", "1/80000", "1/70000", "1/60000", "1/50000", "1/40000", "1/30000", "1/20000", "1/18000", "1/16000", "1/14000", "1/12000", "1/10000", "1/8000", "1/6000", "1/4000", "1/2000")),
    beta_c = recode(beta_c, !!!dict_beta_c),
  ) %>%
  mutate(scenario = paste0(beta_c, "_", beta_e)) %>%
  select(scenario, Environment, Contact) %>%
  pivot_longer(c(Environment, Contact), names_to = "Source", values_to = "SAR") %>%
  ggplot(., aes(x = Source, y = SAR)) +
  facet_wrap(facets = vars(scenario), ncol = 10) +
  geom_rect(ymin = 0.2, ymax = 0.3, xmax = Inf, xmin = -Inf, fill = "grey", col = "grey") + 
  geom_boxplot() +
  geom_point() +
  theme_bw() #+
# theme(
#   axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#   axis.text.y = element_text(size = 14),
#   strip.text = element_text(size = 14),
#   axis.title = element_text(size = 16),
#   legend.title = element_text(size = 16),
#   legend.text = element_text(size = 14)
# ) +
# labs(x = expression(beta[e]), y = "Secondary attack rate", col = expression(beta[c]))
ggsave(file.path("..", "..", "fig", "simulations", "grid-search_env.png"), p, height = 30, width = 20)

## SAR METRICS------------------------------------------------------------------
# all_SAR_g_metrics <- get_all_SAR_g_metrics(list_SAR = list_SAR) ## If only interested in global SAR
all_SAR_metrics <- compute_all_SAR_metrics(list_SAR = list_SAR)

## EPIDEMIC DURATION
list_epidemic_duration <- get_all_epidemic_duration(list_sim = list_sim)
all_epidemic_duration_metrics <- get_all_epidemic_duration_metrics(list_epidemic_duration = list_epidemic_duration)

## ALL METRICS IN ONE DF
all_metrics <- left_join(all_SAR_metrics, all_epidemic_duration_metrics, by="couple", suffix=c("",".y")) %>%
  select(-ends_with(".y"))

write.csv2(all_metrics, file = file.path(scenarios_path, 'scenarios-all_metrics.csv'), row.names = F)

