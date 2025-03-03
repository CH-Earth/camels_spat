#!/bin/bash
#SBATCH --job-name=restructure
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=10GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=/gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/10_analysis/logs/restructure_3/slurm-%A_%a.out
#SBATCH --array=0-1425 

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gdal/3.5.1 libspatialindex/1.8.5 geo-stack/2022a
source /globalhome/wmk934/HPC/camels_spat/camels-spat-env-jupyter/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/10_analysis
python 8b_restructure_basin.py $SLURM_ARRAY_TASK_ID