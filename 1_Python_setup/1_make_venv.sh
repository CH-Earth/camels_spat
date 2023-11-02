#!/bin/bash

# Creates a virtual environment for CAMELS-spat computations.
# Uses pip and the `0_requirements.txt` file in this folder.

# To run Bash files under Windows, use a dedicated terminal such as:
# - Windows Subsystem for Linux: https://docs.microsoft.com/en-us/windows/wsl/install
# - Git Bash (part of Git for Windows): https://gitforwindows.org/
# - Cygwin: https://www.cygwin.com/

# Note: using pip over conda on Windows also required downloading Microsoft Visual Studio C++ Build Tools.
# VIsual Studio needed to be further updated with <..>

# Define the location of the config file
# --------------------------------------
config="../0_config/config.txt"

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
    echo "Unknown OS ${machine}. Aborting."
	exit 1 # Catch-all error code: https://tldp.org/LDP/abs/html/exitcodes.html
fi

# Get all the relevant info from the config file
# ----------------------------------------------
# Config file always has the setting we want as the 2nd item in a row, separated by \ or #
echo "Attempting to read from ${config} using ${machine} command."
if [[ "Unix" == *${machine}* ]]; then
    venv_path=$(sed -n 's/.*|\(.*\)|.*/\1/p' <<< $(grep -m 1 "^venv_path" $config) | tr -d ' ')
    venv_name=$(sed -n 's/.*|\(.*\)|.*/\1/p' <<< $(grep -m 1 "^venv_name" $config) | tr -d ' ')
    code_path=$(sed -n 's/.*|\(.*\)|.*/\1/p' <<< $(grep -m 1 "^code_path" $config) | tr -d ' ')
    reqs_path=$(sed -n 's/.*|\(.*\)|.*/\1/p' <<< $(grep -m 1 "^reqs_path" $config) | tr -d ' ')
    reqs_file=$(sed -n 's/.*|\(.*\)|.*/\1/p' <<< $(grep -m 1 "^reqs_file" $config) | tr -d ' ')
elif [[ "Windows" == *${machine}* ]]; then
	while IFS='|#' read -ra LINE; do venv_path=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^venv_path" $config) # Where the venv should go
	while IFS='|#' read -ra LINE; do venv_name=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^venv_name" $config) # What the venv will be called
	while IFS='|#' read -ra LINE; do code_path=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^code_path" $config) # Root path to code
	while IFS='|#' read -ra LINE; do reqs_path=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^reqs_path" $config) # Requirements file sub-folder(s)
	while IFS='|#' read -ra LINE; do reqs_file=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^reqs_file" $config) # Requirements filename
fi

# Ensure a clean install
# ----------------------
echo "Attempting to install virtualenv ${venv_path}/${venv_name}. Will remove venv folder if it exists to ensure a clean install."
mkdir -p "${venv_path}/${venv_name}" # Ensure the destination directory exists
rm -rf "${venv_path}/${venv_name}" # Remove the virtualenv if it exists

# Make the virtualenv
# -------------------
# Source: https://docs.python.org/3/library/venv.html
# Note: command should work on both Windows and Unix(-like) OS
# Note: assumes Python 3 is available on the system and added to the PATH
python3 -m venv "${venv_path}/${venv_name}"

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

# Check if we activated the virtualenv
# ------------------------------------
# Source: https://stackoverflow.com/a/3063887
if [[ -z "$VIRTUAL_ENV" ]]; then
    echo "Virtual environment not successfully activated. Aborting."
	exit 1
else
    echo "Virtual environment successfully activated in: ${VIRTUAL_ENV}"
	echo # empty line for clarity on the terminal - ugly, but no idea how to do better
fi

# Update the basics
# -----------------
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools

# Make sure we have the required libraries (GDAL)
# -----------------------------------------------
echo "Attempting to activate/install libraries using ${machine} command."
if [[ "Unix" == *${machine}* ]]; then
    if gdalinfo --version &> /dev/null; then
    	gdal_version=$(gdalinfo --version) # E.g. 'GDAL 3.7.2, released 2023/09/05'
    	gdal_version=$(echo $gdal_version | awk '{print $2}' | tr -d ',')
    	echo "GDAL ${gdal_version} install found on system. Continuing."
    else
    	echo "Attempting to install GDAL using homebrew."
    	brew install gdal
    fi
elif [[ "Windows" == *${machine}* ]]; then
    echo "Using pre-built wheels to install GDAL and Fiona for Windows x64."
	# Install GDAL 3.4.3 with a pre-built wheel for Windows
	# Wheel source: https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal
	python -m pip install GDAL-3.4.3-cp38-cp38-win_amd64.whl # GDAL 3.4.3, built for Python 3.8.x, on Windows 64-bit
	export GDAL_DATA="${venv_path}/${venv_name}/Lib/site-packages/osgeo/data/gdal"
	export GDAL_VERSION="3.4.3" # Needed so Fiona knows what GDAL it's working with
	
	# Install Fiona with a pre-built wheel for Windows
	# Wheel source: https://www.lfd.uci.edu/~gohlke/pythonlibs/#fiona
	python3 -m pip install Fiona-1.8.21-cp38-cp38-win_amd64.whl
fi

# Install the required packages
# -----------------------------
python3 -m pip install -r "${code_path}/${reqs_path}/${reqs_file}"

# Force reinstall of GDAL, to ensure numpy connection is correct
# See: https://gis.stackexchange.com/questions/83138/cannot-import-gdal-array
# Note: specifying numpy as a dependency for gdal in requirements.txt (gdal [numpy]) does not seem to work
# Note: specifying the version string without purging cache isn't enough
# Note: Without purging the cache it still fails
python3 -m pip cache purge
gdal_str="python3 -m pip install -v --force-reinstall gdal[numpy]==${gdal_version}"
eval $gdal_str

# Ensure the virtualenv is available as a notebook kernel
# -------------------------------------------------------
ipython kernel install --user --name=${venv_name}

# Run test to check if we have everything we need
# -----------------------------------------------
# Source: https://stackoverflow.com/a/45474387
echo # empty line for clarity on the terminal
echo "Checking packages installed in: ${venv_name}"
cd $code_path/python_tests # needed because using relative paths to a parent directory doesn't play nice with next line
python3 -m unittest test_requirements

# Clean-up
# --------
# Note: virtualenv will be automatically deactivated upon script exit, 
# because executing a script opens a new virtual terminal. Everything
# happens in there and hence in the current terminal (from which we
# activate this script) no virtualenv is activated at all.
