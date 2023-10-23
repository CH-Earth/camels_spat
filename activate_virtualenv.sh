#!/bin/bash

# Activates a virtual environment for CAMELS-spat computations.
# Virtualenv is created by script `1_make_venv.sh` in folder `1_Python_setup`.

# To run Bash files under Windows, use a dedicated terminal such as:
# - Windows Subsystem for Linux: https://docs.microsoft.com/en-us/windows/wsl/install
# - Git Bash (part of Git for Windows): https://gitforwindows.org/
# - Cygwin: https://www.cygwin.com/

# Define the location of the config file
# --------------------------------------
config="./0_config/config.txt"

# Figure out what OS we are on
# ----------------------------
# Source: https://stackoverflow.com/a/3466183
# Note: just adding a few common ones. More options in table here: https://en.wikipedia.org/wiki/Uname
# Implicit assumption: follow-up commands are identical for Cygwin and Mingw. Tested only on Mingw.
# Implicit assumption: follow-up commands are identical for Linux and MacOS. Not tested.
uname_out="$(uname -s)"
case "${uname_out}" in
    Linux*)     machine=Unix;;
    Darwin*)    machine=Unix;;
    CYGWIN*)    machine=Windows;; 
    MINGW*)     machine=Windows;;
    *)          machine="UNKNOWN:${uname_out}" 
esac

# Check if we're ready for the OS we're on
# ----------------------------------------
if [[ "UNKNOWN" == *${machine}* ]]; then
    echo "OS ${machine}"
    exit 1 # Catch-all error code: https://tldp.org/LDP/abs/html/exitcodes.html
fi

# Get all the relevant info from the config file
# ----------------------------------------------
# Config file always has the setting we want as the 2nd item in a row, separated by \ or #
echo "Attempting to read from ${config} using ${machine} command."
if [[ "Unix" == *${machine}* ]]; then
    venv_path=$(sed -n 's/.*|\(.*\)|.*/\1/p' <<< $(grep -m 1 "^venv_path" $config) | tr -d ' ')
    venv_name=$(sed -n 's/.*|\(.*\)|.*/\1/p' <<< $(grep -m 1 "^venv_name" $config) | tr -d ' ')
elif [[ "Windows" == *${machine}* ]]; then
    while IFS='|#' read -ra LINE; do venv_path=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^venv_path" $config) # Where the venv should go
    while IFS='|#' read -ra LINE; do venv_name=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^venv_name" $config) # What the venv will be called
fi

# Activate the virtualenv
# -----------------------
# Source: https://docs.python.org/3/library/venv.html
echo "Attempting to activate ${venv_name} using ${machine} command."
if [[ "Unix" == *${machine}* ]]; then
    source "${venv_path}/${venv_name}/bin/activate"
elif [[ "Windows" == *${machine}* ]]; then
    # Source: https://medium.com/@presh_onyee/activating-virtualenv-on-windows-using-git-bash-python-3-7-1-6b4b21640368
    chmod +x "${venv_path}/${venv_name}/Scripts/activate.bat"
    . "${venv_path}/${venv_name}/Scripts/activate"
fi

# Launch Jupyter notebook
# -----------------------
jupyter notebook