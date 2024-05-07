#!/bin/bash
#SBATCH --job-name=to_lst_1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=2GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=//gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/time_to_lst/slurm-%A_%a.out
#SBATCH --array=0-999 # NOTE: needs to run from 1-1698 because of 1-based indexing in line below

module load StdEnv/2020 libspatialindex/1.8.5 gcc/9.3.0  openmpi/4.0.3 geo-stack/2022a
source /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/camels-spat-env/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data
python 4_update_time_index.py $SLURM_ARRAY_TASK_ID
