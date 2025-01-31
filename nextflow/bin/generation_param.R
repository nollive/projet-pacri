### GENERATION DES PARAMETRES ###

args <- commandArgs(trailingOnly = TRUE)
nb_simu = seq(1,as.numeric(args[1]))
beta_c = 1 #seq(0, 3, 1/4)
beta_e = 1/12000 #seq(0, 1/120000, 1/10000)
networks = "herriot" #c("herriot", "poincare")
threshold = c(60,120,180)
models = c("linear", "log-linear", "exponential")


params = expand.grid(
  nb_simu = nb_simu, 
  beta_c = beta_c, 
  beta_e = beta_e, 
  threshold = threshold,
  network = networks,
  model = models
  ) 

write.table(params, file="param_grid.txt")
