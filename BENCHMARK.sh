#!/bin/bash

# Ergebnisdatei erstellen
RESULT_FILE="benchmark_results.txt"
echo "Benchmarking Results" > $RESULT_FILE
echo "====================" >> $RESULT_FILE
echo "" >> $RESULT_FILE

# Hardware Informationen sammeln
echo "Collecting hardware information..."

# CPU Informationen
CPU_NAME=$(lscpu | grep "Model name:" | awk -F: '{print $2}' | sed 's/^ *//')
CPU_CORES=$(lscpu | grep "^CPU(s):" | awk '{print $2}')
TOTAL_RAM=$(grep MemTotal /proc/meminfo | awk '{print $2 / 1024 " MB"}')

# GPU Informationen (NVIDIA)
if command -v nvidia-smi &> /dev/null
then
    GPU_MODEL=$(nvidia-smi --query-gpu=name --format=csv,noheader)
    GPU_VRAM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader)
else
    GPU_MODEL="N/A"
    GPU_VRAM="N/A"
fi

# Hardware Informationen in die Ergebnisdatei schreiben
echo "CPU Name: $CPU_NAME" >> $RESULT_FILE
echo "CPU Cores: $CPU_CORES" >> $RESULT_FILE
echo "Total RAM: $TOTAL_RAM" >> $RESULT_FILE
echo "GPU Model: $GPU_MODEL" >> $RESULT_FILE
echo "GPU VRAM: $GPU_VRAM MB" >> $RESULT_FILE
echo "" >> $RESULT_FILE

# CUBENAME aus der Makefile auslesen
CUBENAME="DATACUBE_UDF-10.fits"
ZLOW="0"
ZHIGH="1.9"
MAX_MATCH="20"
IGNORE_BELOW="7"

# Preprocessing der Cubefile und Erstellen der Binärdatei
#cp $CUBENAME data/raw/$CUBENAME
echo "Starting preprocessing of the cubefile..."
echo ""
echo "apply median filter to 'flat' out the data (emission lines remain)..."
#python3 src/preprocessing/median-filter-cube.py data/raw/$CUBENAME --signalHDU=1 --varHDU=2 --num_cpu=$CPU_CORES --width=151 --output=data/processed/med_filt.fits > /dev/null
echo ""
echo "filter the data cube with a 'typical line' template in spatial dimension first..."
#python3 src/preprocessing/lsd_cc_spatial.py --input=data/processed/med_filt.fits --SHDU=1 --NHDU=2 --threads=$CPU_CORES --gaussian --lambda0=7050 -pc 0.7 --classic --output=data/processed/spatial_cc.fits --overwrite > /dev/null
echo "filter the data cube with a 'typical line' template in spectral dimension..."
#python3 src/preprocessing/lsd_cc_spectral.py --input=data/processed/spatial_cc.fits --threads=$CPU_CORES --FWHM=250 --SHDU=1 --NHDU=2 --classic --output=data/processed/spectral_cc.fits --overwrite > /dev/null
echo "construct a signal-to-noise cube..."
#python3 src/preprocessing/s2n-cube.py --input=data/processed/spectral_cc.fits --output=data/processed/s2n_v250.fits --clobber --NHDU=2 --SHDU=1 > /dev/null

echo "deleting tmp files..."
#rm data/processed/spatial_cc.fits
#rm data/processed/spectral_cc.fits
#ln -s data/raw/$CUBENAME data/raw/

#echo "Create Masking Plot and transpose Cube for better Cache Access..."
#cd src/preprocessing
#python3 combination.py $CUBENAME s2n_v250.fits
#cd ../..

# Compile the CUDA code
#echo "Starting FELINE..."
#if [ -e feline.bin ]; then
#    rm -f feline.bin
#fi

#make

# Benchmarking feline.bin
echo "Benchmarking feline.bin..."
for i in $(seq 1 10)
do
    echo "Run #$i" >> $RESULT_FILE
    
    # Zeitmessung starten
    START_TIME=$(date +%s.%N)
    
    # Feline ausführen
    time ./feline.bin $ZLOW $ZHIGH $MAX_MATCH $IGNORE_BELOW
    
    # Zeitmessung beenden
    END_TIME=$(date +%s.%N)
    
    # Laufzeit berechnen
    RUN_TIME=$(echo "$END_TIME - $START_TIME" | awk '{print $1 - $2}')
    
    # Laufzeit in die Ergebnisdatei schreiben
    echo "Run time: $RUN_TIME seconds" >> $RESULT_FILE
    echo "" >> $RESULT_FILE
done
rm data/raw/$CUBENAME
echo "Benchmark completed. Results are saved in $RESULT_FILE."
