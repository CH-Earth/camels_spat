#!/bin/bash
#SBATCH --job-name=cleanup_rdrs
#SBATCH --exclude=cnicgiwscpu009 # seems broken right now
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=4GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=//gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/rdrs/7d_cleanup/slurm-%A_%a.out
#SBATCH --array=1-1697 # start at 1 because we did basin 0 during tests

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gdal/3.5.1 libspatialindex/1.8.5 geo-stack/2022a
source /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/camels-spat-env/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data
python 7d_tar_and_remove_day_files.py $SLURM_ARRAY_TASK_ID