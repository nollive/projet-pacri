#!/bin/bash

#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=olivier.gaufres@pasteur.fr
#SBATCH --output=/pasteur/appa/homes/ogaufres/logs/scenarios-analysis/slurm-%j.out
#SBATCH --error=/pasteur/appa/homes/ogaufres/logs/scenarios-analysis/slurm-%j.err

module purge
module load R/4.4.0


n_sim=$1
list_beta=("0" "1/6" "1/4" "1/3" "1/2" "1" "3/2")
list_nu=("0" "5" "10" "15" "20" "25" "30" "35")

for beta in "${list_beta[@]}"; do
  for nu in "${list_nu[@]}"; do
    for i in $(seq 1 $n_sim); do
      echo "Launching simulation $i using beta = $beta and nu = $nu"
      sbatch --qos=fast --job-name="sim_${beta}_${nu}_${i}" --export=ALL,beta=$beta,nu=$nu,sim_id=$i run_simulation_scenario.sh
    done
  done
done
