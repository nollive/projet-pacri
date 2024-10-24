library(dplyr)
library(Rcpp)

wd <- getwd()
##correct<d
wd <- file.path(wd,"R", "nodscov2")
cpp_path <-  file.path(wd,"..", "..", "cpp")
cpp_model_path <- file.path(cpp_path,"nodscov2")
sim_nodscov2_path <- file.path(wd, "..", "..", "out", "sim-nodscov2")
loc_nodscov2_path <- file.path(wd, "..", "..", "out", "loc-nodscov2")


## NodsCov2 data (not shared)
args = commandArgs(trailingOnly=TRUE)
#load("Y:/ogaufres/pacri-non-sync/data/data-nodscov2/admission_ctc_nodscov2.RData")
#load("Y:/ogaufres/pacri-non-sync/data/data-nodscov2/list_ward.RData")
load(file.path(loc_nodscov2_path,"localization-nodscov2.RData")) ##overwrite admission


## ADMISSION
admission_sim <- admission %>%
  mutate(info = ifelse(status == "PE", 1,0)) %>%
  select(id,info) %>%
  mutate(id_ind = id) %>% 
  left_join(rooms, by = c("id_ind" = "id")) %>%
  mutate(room = ifelse(info == 1, rooms[rooms$id == "PE", "id_room"], room)) %>%
  select(id,info,room)
admission_sim$info <- as.integer(admission_sim$info)
admission_sim$room <- as.integer(admission_sim$room)


## ROOMS
#head(rooms)


## INTERACTIONS
#head(global_interaction)


## LOCALIZATIONS
global_localization <- lapply( 1:n_subdivisions, function(t){
  loc_t <- global_localization[[t]]
  loc_t <- loc_t %>% mutate(localization = ifelse(localization == "NOT HERE", -1, localization))
  loc_t$localization <- as.integer(loc_t$localization)
  return(loc_t)
})


## ENVIRONMENT
env_t <- rooms %>%
  distinct(room) %>%
  mutate(env = 0) %>%
  mutate(id_room = rooms %>% distinct(id_room) %>% pull() ) %>%
  select(room, id_room, env)
env_t$id_room <- as.integer(env_t$id_room)

global_environment <- lapply(1:n_subdivisions, function(t){
  return (env_t)
})

rm(env_t)


## STATUS
global_status = data.frame(id = admission$id, t_inf = as.integer(-1), t_recover = as.integer(-1), inf_by = "", inf_room = as.integer(-1))


## LAMBDA
lambda_t <- admission %>%
  distinct(id) %>%
  mutate(lambda_c = 0, lambda_e = 0)
global_lambda <- lapply(1:n_subdivisions, function(t){
  return(lambda_t)
})
rm(lambda_t)



### SIMULATION
## Rcpp
#load("~/CNAM-PASTEUR-FIXE/projet-pacri/data-model-nodscov.RData")
sourceCpp(file.path(cpp_model_path,"model-nodscov2.cpp"))


## PARAMETERS
beta <- 200
epsilon <- 50
mu <- 5
nu <- 5
dt <- 30
tau <- 60 * 60 *24
deltat <- dt/tau
env_threshold <- 0


## INDEX
id_index <- "001-0038-B-S"
status <- global_status %>% mutate(t_inf = ifelse(id == id_index, as.integer(1), t_inf ),
                                   t_recover = ifelse(id == id_index, as.integer(1000), t_recover ),
                                   inf_room = ifelse(id ==id_index, rooms[rooms$id == id_index, "id_room"], inf_room),
                                   inf_by = ifelse(id == id_index, "INDEX", inf_by))


## Simulation
data_sim <- simulation(global_interaction = global_interaction,
                       global_localization = global_localization,
                       global_environment = global_environment,
                       global_lambda = global_lambda,
                       global_status = status,
                       admission = admission_sim,
                       beta = beta,
                       epsilon = epsilon,
                       nu = nu,
                       mu = mu,
                       dt = dt
)


#SAVE TO .RData FILES
id_sim <- ifelse(length(args)==0, "01", args[1])
save(global_interaction,
     global_localization,
     admission,
     admission_sim,
     rooms,
     beta,
     nu,
     mu,
     dt,
     n_subdivisions,
     begin_date,
     end_date,
     data_sim,
     clusters,
     file = file.path(sim_nodscov2_path, paste0(id_sim, "-", "simulation-nodscov2.RData"))
)




