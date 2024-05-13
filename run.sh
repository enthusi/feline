#!/bin/bash

cube_link="martinwendt.de/cube.fits"
cube_name="cube.fits"
#cube_link="https://amused.univ-lyon1.fr/data/UDF/HUDF/download/DATACUBE_UDF-10.fits"
#cube_name="DATACUBE_UDF-10.fits"
# Initialize variables with default values
verbose=0
zlow="0"
zhigh="1.9"
max_match="20"
ignore_below="7"
help_option=""
extension='/dev/null'

# Define function to print usage
usage() {
    echo "Usage: $0 [-v] [--zlow <value>] [--zhigh <value>] [--max <value>] [--ign <value>] [-h | --help]"
    echo "Options:"
    echo "  -v, --verbose              Enable verbose mode"
    echo "      --zlow <value>         Set zlow to <value>"
    echo "      --zhigh <value>        Set zhigh to <value>"
    echo "      --max <value>          Set max_match to <value>"
    echo "      --ign <value>          Set ignore_below to <value>"
    echo "  -h, --help                 Display this help message"
}


# Parse command-line options
VALID_ARGS=$(getopt -o vh --long verbose,zlow:,zhigh:,max:,ign:,help -n "$0" -- "$@")
if [ $? -ne 0 ]; then
    exit 1;
fi

eval set -- "$VALID_ARGS"
while [ "$1" != "--" ]; do
  case "$1" in
    -v | --verbose)
      verbose=1
      shift
      ;;
    --zlow)
      zlow="$2"
      shift 2
      ;;
    --zhigh)
      zhigh="$2"
      shift 2
      ;;
    --max)
      max_match="$2"
      shift 2
      ;;
    --ign)
      ignore_below="$2"
      shift 2
      ;;
    -h | --help)
      help_option="1"
      shift
      ;;
    *)
      echo "Invalid option: $1" >&2
      exit 1
      ;;
  esac
done
shift $((OPTIND - 1))

# If help option is provided, print usage and exit
if [ -n "$help_option" ]; then
    usage
    exit 0
fi

# Your script logic goes here
if [ $verbose -eq 1 ]; then
  echo "Verbose mode enabled"
  extension='/dev/tty'
fi

# Function to download the cube file if not already downloaded
download_cube_file() {
    if [ ! -f $cube_name ]; then
        echo "Downloading cube file..."
        wget --no-check-certificate $cube_link
    else
        echo "Cube file already exists."
    fi
}

draw_progress_bar() {
    local width=50 # Width of the progress bar
    local progress=$1
    local completed=$((progress * width / 100))
    local remaining=$((width - completed))

    # Save cursor position
    tput sc

    # Move cursor to the last line of the terminal window
    tput cup $(tput lines) 0

    # Draw progress bar
    printf "["
    printf "%${completed}s" | tr ' ' '='
    printf ">"
    printf "%${remaining}s" | tr ' ' ' '
    printf "] %d%%" "$progress"

    # Restore cursor position
    tput rc
}
draw_progress_bar 0
# Install Requirements
echo "Creating Virtual Python Environment..."
python3 -m venv venv
. venv/bin/activate
echo "Installing Python dependencies..."
pip install -r requirements.txt >"$extension"
draw_progress_bar 10
export PYTHONWARNINGS="ignore"

# Check and download cube file
download_cube_file

# Set up cube file path
export CUBEFILE=$(realpath $cube_name)
export CUBEFILENAME=$cube_name
# Set number of cores
export CORES=$(nproc)
echo "Number of cores on the current System: $CORES"

echo "Starting preprocessing..."
echo ""
echo "Create med_filt.fits..."
python src/preprocessing/median-filter-cube.py $CUBEFILE --signalHDU=1 --varHDU=2 --num_cpu=$CORES --width=151 --output=data/processed/med_filt.fits >"$extension"
draw_progress_bar 15
echo "Create spatial_cc.fits..."
python src/preprocessing/lsd_cc_spatial.py --input=data/processed/med_filt.fits --SHDU=1 --NHDU=2 --threads=$CORES --gaussian --lambda0=7050 -pc 0.7 --classic --output=data/processed/spatial_cc.fits --overwrite >"$extension"
draw_progress_bar 17
echo "Create spectral_cc.fits..."
python src/preprocessing/lsd_cc_spectral.py --input=data/processed/spatial_cc.fits --threads=$CORES --FWHM=250 --SHDU=1 --NHDU=2 --classic --output=data/processed/spectral_cc.fits --overwrite >"$extension"
draw_progress_bar 22
echo "Create s2n_v250.fits..."
time python src/preprocessing/s2n-cube.py --input=data/processed/spectral_cc.fits --output=data/processed/s2n_v250.fits --clobber --NHDU=2 --SHDU=1 >"$extension" 2> /dev/null

rm data/processed/spatial_cc.fits
rm data/processed/spectral_cc.fits
cp $CUBEFILE data/raw/

cd src/preprocessing || exit
echo "Create Masking Plot and transpose Cube for better Cache Access..."
python combination.py $CUBEFILENAME s2n_v250.fits
draw_progress_bar 25
echo "Starting FELINE..."

# Compile and run FELINE
cd ../../ || exit
make >"$extension"
./feline.bin $zlow $zhigh $max_match $ignore_below
draw_progress_bar 70
echo "Starting Postprocessing and creating PDF..."
# Detect objects and plot results
cd src/postprocessing || exit
python detect_objects.py s2n_v250.fits > catalog.txt
draw_progress_bar 72
python create_final_plots.py $CUBEFILENAME s2n_v250.fits sorted_catalog.txt med_filt.fits J0014m0028 >"$extension" 2>&1
python create_pdf.py
draw_progress_bar 100
echo "Your results are saved to data/final_plots/"

cd ../../ || exit
# Clean Everything (Optional)
read -p "Do you want to clean everything for the next run? (y/n): " clean_choice
if [ "$clean_choice" = "y" ]; then
  find data/raw ! -name '.gitkeep' -type f -delete

	# Delete all files except .gitkeep in data/processed
	find data/processed ! -name '.gitkeep' -type f -delete

	# Remove all files with .bmp and .raw extension from the current directory
	rm -f *.bmp *.raw

	# Delete all files with .txt and .fits extension in src/postprocessing
	find src/postprocessing -type f \( -name '*.txt' -o -name '*.fits' -o -name '*.png' -o -name '*.log' -o ! -name "*.*" \) -delete
fi

exit
