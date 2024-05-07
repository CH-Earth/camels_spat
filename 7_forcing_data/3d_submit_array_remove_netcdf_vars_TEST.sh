#!/bin/bash
#SBATCH --job-name=delvars_TEST
#SBATCH --ntasks=1
#SBATCH --time=0:10:00
#SBATCH --mem=1GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=//gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/del_vars/slurm-%A_%a.out
#SBATCH --array=1-2 # NOTE: needs to run from 1-1698 because of 1-based indexing in line below

module load nco

# This contains a list of basin forcing folders; i.e. /full/path/to/CAN_01AD002/forcing
FOLDER_LIST="/globalhome/wmk934/HPC/camels_spat/7_forcing_data/3d_basin_forcing_folder_list.txt"

# Select the correct folder at the current job array ID
FOLDER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $FOLDER_LIST) 

# Per forcing folder, find the ERA5 files in sub-folders 'raw', 'lumped' and 'dist'
for SUBFOLDER in /raw /lumped /dist; do
 FILES=$(find $FOLDER$SUBFOLDER -name 'ERA5_*.nc' | sort)

 # Per ERA5 file, remove the two variables we don't want, while overwriting the original files
 for FILE in $FILES; do
  ncks -x -v reflected_sw,net_radiation -O $FILE $FILE
 done
done

echo "Completed processing: $FOLDER"

# BACKUP
# This contains a full loop (basins, aggregation levels), but this is too slow to be feasible
#
## Loop over all camels-spat forcing folders and remove net_radiation and reflected_sw variables from netcdf files
#
## Top folder with basin folders in it
#MAIN="/project/gwf/gwf_cmt/wknoben/camels_spat/camels-spat-data/basin_data/"
#
## Find folders named 'forcing' at any depth inside $MAIN
#FORCING=$(find $MAIN -type d -name forcing | sort)
#
## Loop over the forcing folders
#for FOLDER in $FORCING; do
# echo "Processing: $FOLDER"
#
# # Per forcing folder, find the ERA5 files in sub-folders 'raw', 'lumped' and 'dist'
# for SUBFOLDER in /raw /lumped /dist; do
#  FILES=$(find $FOLDER$SUBFOLDER -name 'ERA5_*.nc' | sort)
#
#  # Per ERA5 file, remove the two variables we don't want, while overwriting the original files
#  for FILE in $FILES; do
#   ncks -x -v reflected_sw,net_radiation -O $FILE $FILE
#  done
# done
#done