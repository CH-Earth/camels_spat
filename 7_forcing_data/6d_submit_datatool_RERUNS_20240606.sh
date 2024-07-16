#!/bin/bash
#SBATCH --job-name=subset_daymet
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=20GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=//gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/daymet/6d_subset/slurm-%A_%a.out
#SBATCH --array=19,28,33,70,173,187,193,197,203,204,207,216,218-227,229-234,236,240,242,244,251,253,255,258,260,265,267,271,276,286-296,298-329,331-334,336-374

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gdal/3.5.1 libspatialindex/1.8.5 geo-stack/2022a
source /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/camels-spat-env/bin/activate
cd /gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data
python 6d_daymet_to_basins.py $SLURM_ARRAY_TASK_ID $SLURM_TMPDIR