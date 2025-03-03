#!/bin/bash
#SBATCH --job-name=ss_shp2
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=64GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=/gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/8_geospatial_data/slurm_logs/slurm-2b-%A_%a.out
#SBATCH --array=1005

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gdal/3.5.1 libspatialindex/1.8.5 geo-stack/2022a
source /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/camels-spat-env/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/8_geospatial_data
python 2b_subset_continental_shps_to_basins.py $SLURM_ARRAY_TASK_ID > ./logs/2b_{$SLURM_ARRAY_TASK_ID}.out

