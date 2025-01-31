### SIMULATION MODEL AVEC NEXTFLOW
library(Rcpp)
library(dplyr)
library(lubridate)
library(MASS)

# On recupere la ligne avec le set de parametre
args <- commandArgs(trailingOnly = TRUE)
cat(args, file="test.txt", append = TRUE)
cat("\n", file="test.txt", append = TRUE)


# Fichier avec le code du modele
rlib=args[1] 
sourceCpp(rlib)

# Fichier avec les codes pour l'analyse de la simulation
rlib2 = args[2]
source(rlib2)

# Other arguments
sim_id = as.character(args[4])
beta_c <- as.numeric(args[5])
beta_e <- as.numeric(args[6])
threshold = as.numeric(args[7])
network <- args[8]
model = args[9]

# Load synthetic data
rinput = args[3] # dossier avec les fichiers input
data = paste0(rinput, "/", network, "-simulated-reconstructed-locations-", threshold, ".rda")
load(data)

## RANDOM INDEX CASE
firstDay = as_date(floor_date(begin_date, "day"))
id_index <- sample(x = admission$id[admission$firstDate == firstDay & admission$firstDate != admission$lastDate], size = 1)
# id_index <- sample(x= admission %>% filter(id %in% admission$id) %>% distinct(id) %>% pull(), size = 1)

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
                     beta_c = beta_c,
                     beta_e = beta_e,
                     B = B,
                     nu = nu,
                     mu = mu,
                     env_model = model,
                     dt = dt)


## assign file name
beta_c_type = fractions(beta_c)
beta_e_type = fractions(beta_e)

assign(paste0("sim_", sim_id), result)
save_path <- file.path('grid-search', model, threshold, network, paste0('sim_', beta_c_type, '_', beta_e_type ))
# if (!dir.exists(save_path)) dir.create(save_path, recursive = T, showWarnings = F)

# Save simulation results
save(list = paste0("sim_", sim_id), 
     file = file.path(save_path, paste0('sim_', sim_id, ".rda")))


# Save simulation summary statistics
stats = data.frame(
  network = network,
  model = model,
  threshold = threshold,
  beta_c = beta_c,
  beta_e = beta_e,
  nSim = sim_id
  ) %>%
write.csv2(., paste0("summary_stat_", model, "_", threshold, "_", network, "_", beta_c, "_", beta_e, "_", sim_id, ".csv"), row.names = F)





