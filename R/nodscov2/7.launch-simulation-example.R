################################################################################
##
##                  Generate input data for the simulation of 
##                            epidemics in ICU wards
##
################################################################################
  
## Working environment-----------------------------------------------------------
# Libraries
library(Rcpp)
library(tidyverse)
library(ggpubr)
rm(list = ls())

# Load functions and model
source("helper-functions-simulations.R")
sourceCpp(file.path("..", "..", "cpp", "dev-model-nodscov2_fun.cpp"))
sourceCpp(file.path("..", "..", "cpp", "dev-sensibility-analysis.cpp"))

## Prepare simulation-----------------------------------------------------------
# Load data 
load(file.path("..", "..", "out", "parameters-synthetic-data.rda"))

## Change parameter values  
beta_c <- 3/4
beta_e <- 1e-4
nu <- 1.8e7 * 2 * 24
env_model = "linear" # linear exponential log-linear

## Draw random index case
id_index <- sample(x = admission$id[admission$firstDate == "2020-05-06" & admission$firstDate != admission$lastDate], size = 1)

## Get paramed and med ids
id_paramedical = admission$id[admission$cat == "Paramedical"]
id_medical = admission$id[admission$cat == "Medical"]

## Update global status of the index case (they are infectious at the start of the simulation)
global_status <- global_status %>%
  mutate(t_inf = ifelse(id == id_index,
                        as.integer(1),
                        t_inf),
         t_incub = ifelse(id == id_index,
                          as.integer(1),
                          t_recover),
         t_recover = ifelse(id == id_index,
                            as.integer((t_incub)  + runif(1, min = 2880*3, max = 2880*7)),
                            t_recover),
         inf_by = ifelse(id == id_index,
                         "INDEX",
                         inf_by))

sum(global_status$t_inf>0)

## Conversion factor 
cf = (24*60*2)

## Cumulative risk - Contact
# 1 day
1-exp(-beta_c/cf)^cf
# 5 days
1-exp(-beta_c/cf)^(5*cf)

## Cumulative risk - Environment - Room of 60 m3
# 1 day
1-exp(-(beta_e/cf)*(B/cf)/60*(nu/cf))^cf
# 5 days
1-exp(-(beta_e/cf)*(B/cf)/60*(nu/cf))^(5*cf)


## SIMULATION USING RCPP--------------------------------------------------------
result <- simulation(global_interaction = global_interaction,
                     global_environment = global_environment,
                     global_data = global_data,
                     global_status = global_status,
                     beta_c = beta_c,
                     beta_e = beta_e,
                     B = B,
                     nu = nu,
                     mu = mu,
                     env_model = env_model,
                     dt = dt)
print(sum(result$t_inf>1))

i = 1
while(sum(result$t_inf > 1) < 10) {
  print(i)
  result <- simulation(global_interaction = global_interaction,
                       global_environment = global_environment,
                       global_data = global_data,
                       global_status = global_status,
                       beta_c = beta_c,
                       beta_e = beta_e,
                       B = B,
                       nu = nu,
                       mu = mu,
                       env_model = env_model,
                       dt = dt)
  print(paste0(i, " - ", sum(result$t_inf>1)))
  i = i+1
}

## Plot output------------------------------------------------------------------
summary(result)
sum(result$t_inf>1) / (nrow(result)-1)

p1 = compute_SEIR(result, 90*24*60*2) %>%
  mutate(time = time*30/(3600*24)) %>% # In days
  ggplot(., aes(x = time, y = count, col = status)) +
  geom_line() +
  theme_bw() +
  labs(x = "Time (in days)", y = "Counts", col =  "")

p2 = result %>% 
  summarise(
    n_contact = sum(grepl("CONTACT", inf_by)), 
    n_env = sum(grepl("ENVIRONMENT", inf_by)),
    p_contact = sum(grepl("CONTACT", inf_by)) / (n()-1),
    p_env = sum(grepl("ENVIRONMENT", inf_by)) / (n()-1)
  ) %>%
  pivot_longer(cols = everything(), names_to = c("metric", "inf_source"), values_to = "value", names_pattern = "(.*)_(.*)") %>%
  mutate(metric = ifelse(metric == "n", "Count", "Proportion")) %>%
  ggplot(., aes(x = inf_source, y = value)) +
  ggh4x::facet_grid2(cols = vars(metric), scales = "free_y", independent = "y") +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  theme_bw() +
  labs(x = "Source of infection", y = "")


p3 = data.frame(
  inf_source = c("Paramedical staff", "Medical staff", "Patient"),
  value = c(sum(gsub("CONTACT-", "", result$inf_by) %in% id_paramedical), sum(gsub("CONTACT-", "", result$inf_by) %in% id_medical), 
            sum(grepl("CONTACT-PA", result$inf_by)))
) %>%
  ggplot(., aes(x = inf_source, y = value)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(y = "Counts")


p4 = data.frame(
  inf_source = c("Corridor", "Patient room", "Nursing station", "Office", "Paramedical resting room", "Medical resting room"),
  value = c(sum(result$inf_room==22 & grepl("ENVIRONMENT", result$inf_by)), sum(result$inf_room<=17 & grepl("ENVIRONMENT", result$inf_by)),
            sum(result$inf_room==20 & grepl("ENVIRONMENT", result$inf_by)), sum(result$inf_room==21 & grepl("ENVIRONMENT", result$inf_by)),
            sum(result$inf_room==19 & grepl("ENVIRONMENT", result$inf_by)), sum(result$inf_room==18 & grepl("ENVIRONMENT", result$inf_by))
  )
) %>%
  ggplot(., aes(x = inf_source, y = value)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(y = "Counts")


ggarrange(p1, 
          ggarrange(p2, p3, p4, ncol = 3, labels = c("B", "D", "E")), 
          nrow = 2, labels = c("A", ""))


###########
## PATHS ##
###########
beta <- 3
beta <- 'mean'
beta <- 'median'

## sensibility analysis
loc_path <- file.path('..','..', 'out','loc-nodscov2','dev-localization-nodscov2.RData')
sensibility_path <- file.path('..','..', 'out', 'sensibility-analysis')
RData_path <- file.path(sensibility_path, paste0('sim-beta-', beta), paste0('parameters-model-beta-', beta, '.RData'))
cpp_path <- file.path('..', '..', 'cpp', 'nodscov2', 'dev-sensibility-analysis.cpp')

## scenarios analysis
scenarios_path <- file.path("..","..","out", 'scenarios-analysis')
RData_path <-  file.path(scenarios_path, "parameters-admission-nodscov2.RData")


####################
## RDATA AND RCPP ##
####################
load(RData_path)
sourceCpp(cpp_path)

################################
## CREATE ENDLESS DAY OBJECTS ##
################################
test_interaction <- rep(truncated_interaction, n_days)
test_environment <- rep(truncated_environment, n_days)
test_data <- rep(truncated_data, n_days)
test_status <- truncated_status
new_n_subdivisions <- (t_end-t_begin +1)*n_days


#######################
## RANDOM INDEX CASE ##
#######################
id_index <- sample( x = admission_sim %>% filter(id %in% admission$id) %>% distinct(id) %>% pull(), size = 1)
# UPDATE STATUS
# global_status = test_status
global_status <- global_status %>%
  mutate(t_inf = ifelse(id == id_index,
                        as.integer(1),
                        t_inf),
         t_incub = ifelse(id == id_index, 
                          as.integer(t_inf + runif(1, min = 2880*1, max = 2880*3)),
                          t_recover),
         t_recover = ifelse(id == id_index,
                            as.integer(t_incub + runif(1, min = 2880*3, max = 2880*7)),
                            t_recover),
         inf_by = ifelse(id == id_index,
                         'INDEX',
                         inf_by))

################################# SIMULATION (LOCAL) #################################  -->
## SINGLE SIMULATION (dev-model)
cpp_path <- file.path('..', '..', 'cpp', 'nodscov2')
sourceCpp(file.path(cpp_path, 'dev-model-nodscov2.cpp'))

## Sensibility analysis (only 17 pqtients)
# data_sim <- simulation(global_interaction = test_interaction,
#                      global_environment = test_environment,
#                      global_data = test_data,
#                      global_status = test_status,
#                      beta_c = beta_c,
#                      beta_e = beta_e,
#                      B = B,
#                      nu = nu,
#                      mu = mu,
#                      env_threshold = env_threshold,
#                      dt = dt)

## Scenarios analysis (no replicated truncated data, new patients)
beta = 3/4
nu = 12
data_sim <- simulation(global_interaction = global_interaction,
                       global_environment = global_environment,
                       global_data = global_data,
                       global_status = global_status,
                       beta_c = beta_c,
                       beta_e = beta_e,
                       B = B,
                       nu = nu,
                       mu = mu,
                       env_model = "linear",
                       dt = dt)

args = commandArgs(trailingOnly=TRUE)
id_sim <- ifelse(length(args)==0, "dev-test-endless-day", args[1])
id_sim_path <- file.path(wd, '..', '..', 'out', 'sim-nodscov2' , id_sim)

if (!dir.exists(id_sim_path)) {
  dir.create(id_sim_path, recursive = TRUE)
}

## SAVE FOR SCENARIOS ANALYSIS SINGLE SIMULATION
save(admission,
     admission_sim,
     beta,
     nu,
     mu,
     B,
     dt,
     new_n_subdivisions,
     n_days,
     t_begin,
     t_end,
     data_sim,
     file = file.path(scenarios_path, "scenario-3-unique-simulation-nodscov2.RData")
)

save(truncated_interaction,
     truncated_localization,
     truncated_environment,
     admission,
     admission_sim,
     beta,
     nu,
     mu,
     B,
     dt,
     new_n_subdivisions,
     n_days,
     new_begin_date,
     new_end_date,
     t_begin,
     t_end,
     data_sim,
     file = file.path(id_sim_path, paste0(id_sim, "-", "unique-simulation-nodscov2.RData"))
)

## MULTIPLE SIMULATIONS (C/E/C+E)
n_sim <- 1
sourceCpp(file.path(cpp_path, 'dev-sensibility-analysis.cpp'))

sim_C <- replicate(n_sim, {
  simulation(
    global_interaction = test_interaction,
    global_environment = test_environment,
    global_data = test_data,
    global_status = test_status,
    beta = beta,
    B = 0,
    nu = nu,
    mu = mu,
    dt = dt)
}, simplify = FALSE)


sim_E <- replicate(n_sim, {
  simulation(
    global_interaction = test_interaction,
    global_environment = test_environment,
    global_data = test_data,
    global_status = test_status,
    beta = 0,
    B = B,
    nu = nu,
    mu = mu,
    dt = dt)
}, simplify = FALSE)


sim_C_E <- replicate(n_sim, {
  simulation(
    global_interaction = test_interaction,
    global_environment = test_environment,
    global_data = test_data,
    global_status = test_status,
    beta = beta,
    B = B,
    nu = nu,
    mu = mu,
    dt = dt)
}, simplify = FALSE)


##SAVE .RDATA
save_path <- file.path(id_sim_path, paste0(id_sim, "-", "multiple-simulation-nodscov2.RData"))
save(sim_C,
     sim_E,
     sim_C_E,
     truncated_interaction,
     truncated_localization,
     truncated_data,
     admission,
     admission_sim,
     n_sim,
     rooms,
     beta,
     nu,
     mu,
     B,
     dt,
     new_n_subdivisions,
     n_days,
     new_begin_date,
     new_end_date,
     file = save_path
)


## BENCHMARK--------------------------------------------------------------------
library(microbenchmark)
a <- microbenchmark(simulation(global_interaction = test_interaction,
                               global_environment = test_environment,
                               global_data = test_data,
                               global_status = test_status,
                               beta = beta,
                               B = B,
                               nu = nu,
                               mu = mu,
                               dt = dt), times = 10)
  