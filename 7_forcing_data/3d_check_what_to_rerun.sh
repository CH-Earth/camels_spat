#!/bin/bash
# Loops over some SLURM outputs and looks for the word 'completed', which should only be in there if my script actually finished and printed that to the terminal

FOLDER="/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/del_vars"
for FILE in $FOLDER/*.out; do
 CONTENTS=$(<$FILE)
 if [[ $CONTENTS != *"Completed"* ]]; then
  echo "DNF: $FILE"
 fi
done