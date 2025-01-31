### Script of a nextflow workflow: stats on all simulation indicators
library(tidyverse)

## Collect all indicators in each argument of args and put them together in a dataframe
args <- commandArgs(trailingOnly = TRUE)

simu_count = args[1]
params = data.frame(read.table(args[2], header = TRUE))
data = data.frame(read.table(args[3], header = TRUE, sep = "\t"))


if(length(args>=4)){
    for (i in 4:length(args)){
        data = rbind(data, data.frame(read.table(args[i], header = TRUE, sep = "\t")))
    }
}

#data test
#data = read.table("resu_simu_all.txt", header = TRUE, sep="\t")
#params = read.table("param_grid.txt", header = TRUE)

## Do statistics on the indicators of all simulation

# mean and sd
stats_mean_sd = function(data, params){
  
  data = data %>%
    left_join(params, by="ID_SIMU") 
  
  results_analyse = data %>%
    group_by(group, Indicator, Individual, Population) %>%
    summarise(Mean_value = mean(Value),
              Sd_value = sd(Value))

  return(results_analyse)
}
results_analyse = stats_mean_sd(data, params)

write.table(results_analyse, file = "resu_stat_mean_sd.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# convergence of indicators
# on va regarder mean et sd de chaque indicator et de chaque set de parametre en le calculant pour 1 simulation, 2 simulations, etcc et voir si ça converge au bout d'un certains nombres de simulations ou non
# permet de savoir si on a réalisé assez de simulations de chaque set de parametres

evaluate_convergence = function(data, params, simu_count){
  # evalue la convergence de chaque groupe de simulation et pour chaque indicateur
  # ARGUMENTS
    # data : dataframe avec en colonne : group, Indicator, Individual, Population, Value
    # params : dataframe avec en colonne : nb_simu, input, ... (tous les parametres) ... , init, seed, ID_SIMU, group
    # simu_count : nombre de simulations effectuees avec le meme jeu de parametre
  

  data = data %>%
    left_join(params, by="ID_SIMU") %>% 
    mutate(Indicator_complet = paste(Indicator, Individual, Population, sep=" "))  %>% 
    select(-Indicator, -Individual, -Population)
  
  n_indicators = length(unique(data$Indicator))
  n_groups = length(unique(data$group))


  # calcul de la moyenne pour chaque groupe de même parametres pour i simulation
  res_conv = data.frame()
  for (i in 1:simu_count){
    results_convergence = data %>%
      group_by(group, Indicator_complet) %>%
      slice_head(n = i) %>%
      summarise(Mean_value = mean(Value)) %>% 
      mutate(nb_simu = i) %>% 
      ungroup()
    
    res_conv = rbind(res_conv, results_convergence) 
  }
    
  # calcul des ecarts de moyenne 2 a 2
  res_conv2 = res_conv %>% 
    arrange(Indicator_complet, group, nb_simu) %>% 
    group_by(Indicator_complet, group) %>% 
    mutate(Ecart_Mean_value = Mean_value - lag(Mean_value))
    
  # evaluation de la convergence, critere : ecart des moyennes entre j et j-1 <= 0.1 * moyenne j
  # le n° de la simulation ou il y a convergence correspond au dernier où il n'y a pas convergence +1
  res_conv3 = res_conv2 %>% 
    mutate(threshold = 0.1*Mean_value) %>% 
    mutate(Condition = ifelse(Ecart_Mean_value > threshold, TRUE, FALSE)) %>% 
    summarise(
      nb_simu_cv = if (any(Condition, na.rm = TRUE)) {max(nb_simu[Condition], na.rm = TRUE)+1 } 
      else { min(nb_simu, na.rm = TRUE)},
      .groups = "drop") %>% 
    mutate(convergence = nb_simu_cv <= simu_count )
    
  # evaluation du nombre d'indicateurs qui convergent pour chaque groupe
  res_conv4 = res_conv3 %>% 
    group_by(group) %>% 
    filter(convergence) %>% 
    summarise(nb_indicator_cv = n()/n_indicators*100,
              max_nb_simu_cv = max(nb_simu_cv))
    
  # plot pour chaque indicateurs le nombre de groupes qui convergent
  res_conv5 = res_conv3 %>% 
    group_by(Indicator_complet) %>% 
    filter(convergence) %>% 
    summarise(nb_indicator_cv = n()/n_groups*100,
              max_nb_simu_cv = max(nb_simu_cv))
  
  
  return(list(res_per_group=res_conv4, res_per_indicator=res_conv5, res_conv = res_conv))
  
}



  
res = evaluate_convergence(data, params, simu_count)

write.table(res$res_per_group, file = "resu_stat_convergence_per_group.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(res$res_per_indicator, file = "resu_stat_convergence_per_indicator.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(res$res_conv, file = "resu_stat_convergence.txt", sep = "\t", row.names = FALSE, col.names = TRUE)




