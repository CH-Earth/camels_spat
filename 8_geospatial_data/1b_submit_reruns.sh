#!/bin/bash
#SBATCH --job-name=soilgrids
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=16GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=/gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/8_geospatial_data/slurm_logs/slurm-%A_%a.out
#SBATCH --array=34,39,44,49,54,59

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gdal/3.5.1 libspatialindex/1.8.5 geo-stack/2022a
source /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/camels-spat-env/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/8_geospatial_data
python 1b_soil.py $SLURM_ARRAY_TASK_ID > ./logs/1b_RERUNS_{$SLURM_ARRAY_TASK_ID}.out


