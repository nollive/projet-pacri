#!/bin/bash

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=olivier.gaufres@pasteur.fr
#SBATCH --output=/pasteur/appa/homes/ogaufres/logs/scenarios-analysis/simulation-%j.out
#SBATCH --error=/pasteur/appa/homes/ogaufres/logs/scenarios-analysis/simulation-%j.err

module purge
module load R/4.4.0

# vars from launch-scenarios-analysis.sh
beta=${beta:-"ERROR-NO-BETA"}
nu=${nu:-"ERROR-NO-NU"}
sim_id=${sim_id:-"ERROR-NO-ID"}

# pass vars to R script
Rscript $HOME/scenarios-analysis/run_simulation_scenarios.R $sim_id $beta $nu