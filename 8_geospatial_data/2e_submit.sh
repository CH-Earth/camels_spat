#!/bin/bash
#SBATCH --job-name=ss_cfvo_unc
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=128GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=/gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/8_geospatial_data/slurm_logs/slurm-2e-%A_%a.out

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gdal/3.5.1 libspatialindex/1.8.5 geo-stack/2022a scipy-stack/2023b
source /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/camels-spat-env/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/8_geospatial_data

# Loop from 0 to 1697 basins
for ((i=0; i<=1697; i++)); do
    BASIN_IX=$i
    python 2e_subset_continental_tifs_to_basins_CFVO_UNCERTAINTY_RERUNS.py $BASIN_IX > ./logs/2e_{$BASIN_IX}.out
done




