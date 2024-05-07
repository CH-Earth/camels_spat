# Load the geostack so we have most of what we need
#module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 geo-stack/2022a # This works
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 gdal/3.5.1 libspatialindex/1.8.5 geo-stack/2022 geo-stack/2022a
# Define the location of the config file
# --------------------------------------
config="../0_config/config.txt"

# Get all the relevant info from the config file
# ----------------------------------------------
# Config file always has the setting we want as the 2nd item in a row, separated by \ or #
while IFS='|#' read -ra LINE; do venv_path=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^venv_path" $config) # Where the venv should go
while IFS='|#' read -ra LINE; do venv_name=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^venv_name" $config) # What the venv will be called
while IFS='|#' read -ra LINE; do code_path=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^code_path" $config) # Root path to code
while IFS='|#' read -ra LINE; do reqs_path=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^reqs_path" $config) # Requirements file sub-folder(s)
while IFS='|#' read -ra LINE; do reqs_file=$(echo ${LINE[1]}); done <<< $(grep -m 1 "^reqs_file" $config) # Requirements filename

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
python -m venv "${venv_path}/${venv_name}"

# Activate the virtualenv
# -----------------------
# Source: https://docs.python.org/3/library/venv.html
echo "Attempting to activate ${venv_name} using ${machine} command."
source "${venv_path}/${venv_name}/bin/activate"

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
python -m pip install --upgrade pip
python -m pip install --upgrade wheel
python -m pip install --upgrade setuptools
python -m pip install cdsapi baseflow rasterstats
python -m pip install --upgrade pandas xarray 
python -m pip install --upgrade baseflow rasterstats
python -m pip install --upgrade numpy numba
python -m pip install --upgrade geopandas