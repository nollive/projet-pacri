## LIBRARIES
library(dplyr)
library(Rcpp)

wd <- getwd()
##PATHS
loc_path <- file.path(wd,"out","loc-nodscov2","dev-localization-nodscov2.RData")
cpp_path <- file.path(wd,'cpp', 'nodscov2', 'dev-sensibility-analysis.cpp')
##LOAD DATA & CPP CODE
load(loc_path) ##overwrite admission
sourceCpp(cpp_path)

## GLOBAL OBJECTS MODIFICATIONS
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

lambda_t <- admission %>%
  distinct(id) %>%
  mutate(lambda_c = 0, lambda_e = 0)
global_lambda <- lapply(1:n_subdivisions, function(t){
  return(lambda_t)
})
rm(lambda_t)


##PARAMETERS
p_PA_PA <- 2.3e-5
p_PA_PE <- 1.19e-4
p_PE_PA <- 7.89e-4
p_PE_PE <- 1.66e-4

#p <- median(c(p_PA_PA,p_PA_PE,p_PE_PA,p_PE_PE))
p <- mean(c(p_PA_PA,p_PA_PE,p_PE_PA,p_PE_PE))
beta <- (-log(1-p))*(86400/30)
B <- 0.48 *24 # Breating rate ( 0.48 m3/h)
mu_air <- 4 * 24 # Quanta removal (4-6 Air changes/h)
mu_inac <- (log(2)/1.1) * 24 # Quanta inactivation (computed with viral half life -> 0.63quanta inactivated/h) 
mu <- mu_air + mu_inac
nu <- 10 * 24 #( 10 quanta/h)
dt <- 30 # Time step
tau <- 60 * 60 *24 # Seconds in 1 day
deltat <- dt/tau
env_threshold <- 0 # Quanta threshold above which the environment is infectious


## INDEX
id_index <- "001-0038-B-S"
status <- global_status %>%
  mutate(t_inf = ifelse(id == id_index,
                        as.integer(1),
                        t_inf),
         t_incub = ifelse(id == id_index, 
                          as.integer(t_inf + runif(1, min = 2880*1, max = 2880*5)),
                          t_recover),
         t_recover = ifelse(id == id_index,
                            as.integer(t_incub + runif(1, min = 2880*3, max = 2880*7)),
                            t_recover),
         inf_by = ifelse(id == id_index,
                         "INDEX",
                         inf_by))



## endless day
end_date - begin_date ##end of last interaction
new_begin_date <- begin_date + (5*60*60) ## OFFSET TO START THE DAY AT 12AM
new_end_date <- new_begin_date + (24*60*60)
print(new_end_date - new_begin_date)
t_begin <- 5*60*2 ## OFFSET t
t_end <- (t_begin + 24*60*2) - 1


## TRUNCATED INFO
truncated_interaction <- global_interaction[(u <- seq_along(global_interaction)) %in% t_begin:(t_end)]
truncated_localization <- global_localization[(u <- seq_along(global_localization)) %in% t_begin:(t_end)]
truncated_lambda <- global_lambda[(u <- seq_along(global_lambda)) %in% t_begin:(t_end)]
truncated_environment <- global_environment[(u <- seq_along(global_environment)) %in% t_begin:(t_end)]
truncated_status <- status

## ENDLESS DAY OBJECTS
n_days <- 7*7
test_interaction <- rep(truncated_interaction, n_days)
test_localization <- rep(truncated_localization, n_days)
test_lambda <- rep(truncated_lambda, n_days)
test_environment <- rep(truncated_environment, n_days)
test_status <- truncated_status

new_n_subdivisions <- (t_end-t_begin +1)*n_days

##SIMULATION
n_sim <- 50

save(n_sim,
     n_days,
     test_interaction,
     test_localization,
     test_environment,
     test_lambda,
     test_status,
     admission_sim,
     id_index,
     beta,
     B,
     nu,
     mu,
     env_threshold,
     dt,
     file = file.path("parameters-model.RData")
     )


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

