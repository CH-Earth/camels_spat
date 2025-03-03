#!/bin/bash
#SBATCH --job-name=delvars_rr
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=1GB
#SBATCH --mail-user=wouter.knoben@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=hpc_c_giws_clark
#SBATCH --output=//gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/del_vars/slurm-rerun-%A_%a.out
#SBATCH --array=1-498 # Length of the file

module load nco

# Define a text file where we'll save the links to the forcing files that still need to be processed
FILE_LIST="/globalhome/wmk934/HPC/camels_spat/7_forcing_data/3d_basin_forcing_file_rerun_list.txt"

# Get the right line in the file
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $FILE_LIST)

# Process
ncks -x -v reflected_sw,net_radiation -O $FILE $FILE


## I ran the code below only once, to get a full list of files we still need to process. This is easier to parallelize
# This contains a list of basin forcing folders; i.e. /full/path/to/CAN_01AD002/forcing
#FOLDER_LIST="/globalhome/wmk934/HPC/camels_spat/7_forcing_data/3d_basin_forcing_folder_list.txt"

# Select the correct folder at the current job array ID
# We know basins 761 and 987 did not complete (3d_check_what_to_rerun.sh)
# NOTE: we also know these are really large basins that we probably don't want to use anyway
#FOLDER1=$(sed -n "761p" < $FOLDER_LIST)
#FOLDER2=$(sed -n "987p" < $FOLDER_LIST)



# Create a list of specific files that we still need to process
# BASIN 761: manual checks indicate only part of /distributed is left to complete
# ---------------------------------------------------------------------------------
#SUBFOLDER="/distributed"
#FILES=$(find $FOLDER1$SUBFOLDER -name 'ERA5_*.nc' | sort)
#LAST_COMPLETE="/ERA5_dist_remapped_2015-11-01-00-00-00.nc"

# Per ERA5 file, remove the two variables we don't want, while overwriting the original files
#for FILE in $FILES; do
# if [[ $FILE > $FOLDER1$SUBFOLDER$LAST_COMPLETE ]]; then
#  echo "$FILE" >> $FILE_LIST
# fi
#done

# BASIN 987: manual checks indicate only part of /distributed is left to complete
# ---------------------------------------------------------------------------------
#SUBFOLDER="/distributed"
#FILES=$(find $FOLDER2$SUBFOLDER -name 'ERA5_*.nc' | sort)
#LAST_COMPLETE="/ERA5_dist_remapped_1988-09-01-00-00-00.nc"

# Per ERA5 file, remove the two variables we don't want, while overwriting the original files
#for FILE in $FILES; do
# if [[ $FILE > $FOLDER2$SUBFOLDER$LAST_COMPLETE ]]; then
#  echo "$FILE" >> $FILE_LIST
# fi
#done

# Manually check how many lines this file has
# wc -l $FILE_LIST
# $498

# We thus need to run an array of length 1-498 

# Create some testing code
#for i in {1..498}; do
# FILE=$(sed -n "${i}p" < $FILE_LIST)
# echo $FILE
#done