## LIBRARIES
library(dplyr)
library(Rcpp)
library(tidyr)

wd <- getwd()
##PATHS
loc_path <- file.path(wd,"out","loc-nodscov2","dev-localization-nodscov2.RData")
##########################
## LOAD DATA & CPP CODE ##
##########################
load(loc_path)

##################################
## GLOBAL OBJECTS MODIFICATIONS ## 
##################################
admission_sim <- admission %>%
  mutate(info = ifelse(status == "PE", 1,0)) %>%
  select(id,info) %>%
  mutate(id_ind = id) %>% 
  left_join(rooms, by = c("id_ind" = "id")) %>%
  mutate(room = ifelse(info == 1, rooms[rooms$id == "M-PM", "id_room"], room)) %>%
  select(id,info,room)

admission_sim$info <- as.integer(admission_sim$info)
admission_sim$room <- as.integer(admission_sim$room)

env_t <- rooms %>%
  distinct(room) %>%
  mutate(env = 0) %>%
  mutate(id_room = rooms %>% distinct(id_room) %>% pull() ) %>%
  left_join(rooms %>% select(id_room, volume), by = "id_room") %>%
  select(room, id_room, env, volume)
env_t$id_room <- as.integer(env_t$id_room)

global_environment <- lapply(1:n_subdivisions, function(t){
  return (env_t)
})
rm(env_t)

global_status = data.frame(id = admission$id, t_inf = as.integer(-1), t_incub = as.integer(-1), t_recover = as.integer(-1), inf_by = "", inf_room = as.integer(-1))

## TRUNCATE TO KEEP ONLY 24H (24/05/06 12:00 TO 24/05/07 12:00)
end_date - begin_date ##end of last interaction
new_begin_date <- begin_date + (5*60*60) ## OFFSET TO START THE DAY AT 12AM
new_end_date <- new_begin_date + (24*60*60)
print(new_end_date - new_begin_date)
t_begin <- 5*60*2 ## OFFSET t
t_end <- (t_begin + 24*60*2) - 1


## TRUNCATED INFO
truncated_interaction <- global_interaction[(u <- seq_along(global_interaction)) %in% t_begin:(t_end)]
truncated_localization <- global_localization[(u <- seq_along(global_localization)) %in% t_begin:(t_end)]
truncated_environment <- global_environment[(u <- seq_along(global_environment)) %in% t_begin:(t_end)]
truncated_status <- global_status


## FILTER INDIVIDUALS NOT INTERACTING
id_interacting <- do.call(rbind,truncated_interaction) %>%
  pivot_longer(cols = c(from, to), names_to = "direction", values_to = "individual") %>%
  distinct(individual) %>%
  pull()
id_total <- admission %>% distinct(id) %>% pull()
id_not_interacting <- setdiff(id_total, id_interacting)

filter_interacting <- function(list_df, id_interacting) {
  lapply(list_df, function(df) {
    df %>% filter(id %in% id_interacting)
  })
}

admission_sim <- admission_sim %>% filter(id %in% id_interacting)
truncated_localization <- filter_interacting(truncated_localization, id_interacting )
truncated_status <- truncated_status %>% filter(id %in% id_interacting)
truncated_data <- lapply(1:length(truncated_localization), function(t){
  df <- data.frame(id = truncated_localization[[t]]$id,
                   info = admission_sim$info,
                   #room = admission_sim$room,
                   #interacting = ifelse(id %in% truncated_interaction$from | id %in% truncated_interaction$to, T, F)
                   localization_ti = truncated_localization[[t]]$localization,
                   lambda_c = 0,
                   lambda_e = 0
  )
  df <- df %>% filter(localization_ti != -1)
  return(df)
})



## ENDLESS DAY OBJECTS
n_days <- 7*10

##PARAMETERS
B <- 0.48 *24 # Breathing rate ( 0.48 m3/h)
mu_air <- 4 * 24 # Quanta removal (4-6 Air changes/h)
mu_inac <- (log(2)/1.1) * 24 # Quanta inactivation (computed with viral half life -> 0.63 quanta inactivated/h) 
mu <- mu_air + mu_inac
nu <- 10 * 24 # (10 quanta/h)
dt <- 30 # Time step
tau <- 60 * 60 *24 # Seconds in 1 day
deltat <- dt/tau
env_threshold <- 0 # Quanta threshold above which the environment is infectious

# BETA
p_PA_PA <- 2.3e-5
p_PA_PE <- 1.19e-4
p_PE_PA <- 7.89e-4
p_PE_PE <- 1.66e-4
#p <- median(c(p_PA_PA,p_PA_PE,p_PE_PA,p_PE_PE))
p <- mean(c(p_PA_PA,p_PA_PE,p_PE_PA,p_PE_PE))
beta <- (-log(1-p))*(86400/30)


## SAVE ALL BETA CONFIGURATIONS
for (beta in seq(from = 0.5, to = 10, by = 0.5)){
  save(n_days,
       truncated_interaction,
       truncated_environment,
       truncated_data,
       truncated_status,
       admission_sim,
       beta,
       B,
       nu,
       mu,
       env_threshold,
       dt,
       t_begin,
       t_end,
       file = file.path(paste0("parameters-model-beta-", beta, ".RData"))
  )
}
## MEDIAN BETA
p_median <- median(c(p_PA_PA,p_PA_PE,p_PE_PA,p_PE_PE))
beta<- (-log(1-p_median))*(86400/30) ##SEC IN A DAY/30 SEC (30sec time-step)

save(n_days,
     truncated_interaction,
     truncated_environment,
     truncated_data,
     truncated_status,
     admission_sim,
     beta,
     B,
     nu,
     mu,
     env_threshold,
     dt,
     t_begin,
     t_end,
     file = file.path("parameters-model-beta-median.RData")
)
## MEAN
p_mean <- mean(c(p_PA_PA,p_PA_PE,p_PE_PA,p_PE_PE))
beta <- (-log(1-p_mean))*(86400/30) ##SEC IN A DAY/30 SEC (30sec time-step)
save(n_days,
     truncated_interaction,
     truncated_environment,
     truncated_data,
     truncated_status,
     admission_sim,
     beta,
     B,
     nu,
     mu,
     env_threshold,
     dt,
     t_begin,
     t_end,
     file = file.path("parameters-model-beta-mean.RData")
)


###### launch_sensibility.R (CLUSTER)
load(file.path("parameters-model-beta-10.0.RData"))
cpp_path <- file.path(wd,'cpp', 'nodscov2', 'dev-sensibility-analysis.cpp')
sourceCpp(cpp_path)

## CREATE ENDLESS DAY OBJECTS
test_interaction <- rep(truncated_interaction, n_days)
test_localization <- rep(truncated_localization, n_days)
test_lambda <- rep(truncated_lambda, n_days)
test_environment <- rep(truncated_environment, n_days)
test_status <- truncated_status

new_n_subdivisions <- (t_end-t_begin +1)*n_days


## RANDOM INDEX CASE
id_index <- sample( x = admission_sim %>% distinct(id) %>% pull(), size = 1)
# UPDATE STATUs
test_status <- test_status %>%
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
                         "INDEX",
                         inf_by))

##SIMULATION
n_sim <- 50

sim_C <- replicate(n_sim, {
  simulation(
    global_interaction = test_interaction,
    global_localization = test_localization,
    global_environment = test_environment,
    global_lambda = test_lambda,
    global_status = test_status,
    admission = admission_sim,
    beta = beta,
    B = 0,
    nu = nu,
    mu = mu,
    env_threshold = env_threshold,
    dt = dt
  )
}, simplify = FALSE)

sim_E <- replicate(n_sim, {
  simulation(
    global_interaction = test_interaction,
    global_localization = test_localization,
    global_environment = test_environment,
    global_lambda = test_lambda,
    global_status = test_status,
    admission = admission_sim,
    beta = 0,
    B = B,
    nu = nu,
    mu = mu,
    env_threshold = env_threshold,
    dt = dt
  )
}, simplify = FALSE)

sim_C_E <- replicate(n_sim, {
  simulation(
    global_interaction = test_interaction,
    global_localization = test_localization,
    global_environment = test_environment,
    global_lambda = test_lambda,
    global_status = test_status,
    admission = admission_sim,
    beta = beta,
    B = B,
    nu = nu,
    mu = mu,
    env_threshold = env_threshold,
    dt = dt
  )
}, simplify = FALSE)



##SAVE .RDATA
save_path <- file.path(wd,"dev-sensibility-analysis.RData")
save(truncated_interaction,
     truncated_localization,
     admission,
     admission_sim,
     n_sim,
     rooms,
     beta,
     nu,
     mu,
     B,
     env_threshold,
     dt,
     new_n_subdivisions,
     n_days,
     new_begin_date,
     new_end_date,
     file = save_path
)

