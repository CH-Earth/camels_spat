#!/bin/bash
#SBATCH --job-name=winddir_rr
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=//gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/wind_dir/slurm-rerun-%A_%a.out
#SBATCH --array=0-918

module load StdEnv/2020 libspatialindex/1.8.5 gcc/9.3.0  openmpi/4.0.3 geo-stack/2022a
source /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/camels-spat-env/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data
python 3b_add_derived_wind_direction_RERUNS.py $SLURM_ARRAY_TASK_ID


