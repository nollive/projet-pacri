library(Rcpp)
library(dplyr)

wd <- getwd()

args <- commandArgs(trailingOnly = TRUE)

load(file.path("/pasteur/appa/homes/ogaufres/scenarios-analysis/RData/", 'parameters-admission-nodscov2.RData'))
sourceCpp("/pasteur/appa/homes/ogaufres/scenarios-analysis/cpp/dev-sensibility-analysis.cpp")

## EDIT BETA AND NU
sim_id <- as.character(args[1])
beta_type <- gsub(pattern = '/', replacement = '-', x = as.character(args[2]))
nu_type <- gsub(pattern = '/', replacement = '-', x = as.character(args[3]))

beta <- eval(parse(text = args[2]))
nu <- eval(parse(text = args[3]))

## RANDOM INDEX CASE
id_index <- sample(x= admission_sim %>% filter(id %in% admission$id) %>% distinct(id) %>% pull(), size = 1)
# Update global status (index is infectious when the simulation starts)
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

## SIMULATION USING RCPP
result <- simulation(global_interaction = global_interaction,
                     global_environment = global_environment,
                     global_data = global_data,
                     global_status = global_status,
                     beta = beta,
                     B = B,
                     nu = nu,
                     mu = mu,
                     env_threshold = env_threshold,
                     dt = dt)


## assign file name
assign(paste0("sim_", beta_type, '_', nu_type, "_",  sim_id), result)
save_path <- file.path(wd, "out", 'scenarios-analysis', 'out', paste0('sim_', beta_type, '_', nu_type ), 'results')
dir.create(save_path, recursive = T, showWarnings = F)

# Save simulation results
save(list = paste0("sim_", beta_type, '_', nu_type, "_",  sim_id), file = file.path(save_path, paste0('sim_', beta_type, '_', nu_type, '_', sim_id, ".RData")))


