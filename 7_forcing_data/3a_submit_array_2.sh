#!/bin/bash
#SBATCH --job-name=easymore_2
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=2GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=//gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/slurm-%A_%a.out
#SBATCH --array=1000-1697

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gdal/3.5.1 libspatialindex/1.8.5 geo-stack/2022a
source /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/camels-spat-env/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data
python 3a_average_forcing_to_hrus.py $SLURM_ARRAY_TASK_ID > ./easymore_logs/easymore_2_$SLURM_ARRAY_TASK_ID.out


