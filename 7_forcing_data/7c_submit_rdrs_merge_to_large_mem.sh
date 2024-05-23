#!/bin/bash
#SBATCH --job-name=merge_rdrs
#SBATCH --ntasks=1
#SBATCH --time=0:45:00
#SBATCH --mem=50GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_mem_clark
#SBATCH --output=//gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/rdrs/7b_merge/slurm-%A_%a.out
#SBATCH --array=64-467

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gdal/3.5.1 libspatialindex/1.8.5 geo-stack/2022a
source /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/camels-spat-env/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data
python 7b_merge_rdrs_to_month_and_drop_vars.py $SLURM_ARRAY_TASK_ID