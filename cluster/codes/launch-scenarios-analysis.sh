#!/bin/bash

#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=maylis.layan@pasteur.fr
#SBATCH --output=/pasteur/appa/homes/maylayan/Pacri/Hospitals_modes_transmission/github_codes/cluster/logs_loop/slurm-%j.out
#SBATCH --error=/pasteur/appa/homes/maylayan/Pacri/Hospitals_modes_transmission/github_codes/cluster/logs_loop/slurm-%j.err

module purge
module load R/4.4.0


n_sim=$1
#list_beta_c=("0" "1/6" "1/4" "1/3" "1/2" "1" "3/2")
list_beta_c=("3/2")
#list_beta_e=("0" "5" "10" "15" "20" "25" "30" "35")
list_beta_e=("1/100")
list_models=("linear", "log-linear", "exponential")


for model in "${list_models[@]}"; do
  for beta_c in "${list_beta_c[@]}"; do
    for beta_e in "${list_beta_e[@]}"; do
      for i in $(seq 1 $n_sim); do
        echo "Launching simulation $i using model = ${model}, beta_c = $beta_c and beta_e = $beta_e"
        sbatch --qos=fast --job-name="sim_${model}_${beta_c}_${beta_e}_${i}" --export=ALL,beta_c=$beta_c,beta_e=$beta_e,sim_id=$i,model=$model run_simulation_scenario.sh
      done
    done
  done
done
