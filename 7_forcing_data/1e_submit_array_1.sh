#!/bin/bash
#SBATCH --job-name=derived_vars_1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=1GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=//gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/slurm-%A_%a.out
#SBATCH --array=0-999

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 geo-stack/2022a
source /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/camels-spat-env/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data
python 1e_add_derived_variables.py $SLURM_ARRAY_TASK_ID


