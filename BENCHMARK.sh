#!/bin/bash

# Ergebnisdatei erstellen
RESULT_FILE="benchmark_results.txt"
echo "Benchmarking Results" > $RESULT_FILE
echo "====================" >> $RESULT_FILE
echo "" >> $RESULT_FILE

# Hardware Informationen sammeln
echo "Collecting hardware information..."

CPU_NAME=$(lscpu | grep "Model name:" | awk -F: '{print $2}' | sed 's/^ *//')
CPU_CORES=$(lscpu | grep "^CPU(s):" | awk '{print $2}')
CPU_MAX_FREQ=$(lscpu | grep "CPU max MHz:" | awk -F: '{print $2}' | sed 's/^ *//')
TOTAL_RAM=$(grep MemTotal /proc/meminfo | awk '{print $2 / 1024 " MB"}')

# Hardware Informationen in die Ergebnisdatei schreiben
echo "CPU Name: $CPU_NAME" >> $RESULT_FILE
echo "CPU Cores: $CPU_CORES" >> $RESULT_FILE
echo "CPU Max Frequency: $CPU_MAX_FREQ MHz" >> $RESULT_FILE
echo "Total RAM: $TOTAL_RAM" >> $RESULT_FILE
echo "" >> $RESULT_FILE

# CUBENAME aus der Makefile auslesen
CUBENAME=$(grep -oP '^CUBENAME := "\K[^"]+' Makefile)

if [ -n "$CUBENAME" ]; then
    # Dateigröße ermitteln
    if [ -f "$CUBENAME" ]; then
        CUBE_SIZE=$(du -h "$CUBENAME" | awk '{print $1}')
        echo "Cube File: $CUBENAME" >> $RESULT_FILE
        echo "Cube File Size: $CUBE_SIZE" >> $RESULT_FILE
    else
        echo "Cube File: $CUBENAME not found" >> $RESULT_FILE
    fi
    echo "" >> $RESULT_FILE
else
    echo "CUBENAME not found in Makefile" >> $RESULT_FILE
    echo "" >> $RESULT_FILE
fi

# Benchmark durchführen
echo "Starting the benchmark..."
for i in {1..10}
do
    echo "Run #$i" >> $RESULT_FILE
    
    # Initialisieren der Variablen zur Berechnung der Durchschnittswerte
    FREQ_SUM=0
    FREQ_COUNT=0
    
    # Zeitmessung starten
    START_TIME=$(date +%s.%N)
    
    # make run ausführen und währenddessen Taktrate überwachen
    make run &
    PID=$!
    
    while kill -0 $PID 2> /dev/null; do
        # Alle CPU-Kern-Taktraten sammeln
        FREQ=$(cat /proc/cpuinfo | grep "MHz" | awk '{print $4}')
        
        # Durchschnitt der Taktraten berechnen
        AVG_CORE_FREQ=$(echo "$FREQ" | awk '{sum+=$1} END {print sum/NR}')
        
        # Taktrate summieren
        FREQ_SUM=$(echo "$FREQ_SUM + $AVG_CORE_FREQ" | bc)
        
        # Anzahl der Messungen erhöhen
        FREQ_COUNT=$((FREQ_COUNT + 1))
        
        sleep 1  # 1 Sekunde warten, bevor die Werte erneut überprüft werden
    done
    
    wait $PID
    
    # make clean ausführen
    make clean
    
    # Zeitmessung beenden
    END_TIME=$(date +%s.%N)
    
    # Laufzeit berechnen
    RUN_TIME=$(echo "$END_TIME - $START_TIME" | bc)
    
    # Durchschnittstaktrate berechnen
    if [ $FREQ_COUNT -gt 0 ]; then
        AVG_FREQ=$(echo "$FREQ_SUM / $FREQ_COUNT" | bc -l)
    else
        AVG_FREQ="N/A"
    fi
    
    # Laufzeit und Durchschnittstaktrate in die Ergebnisdatei schreiben
    echo "Run time: $RUN_TIME seconds" >> $RESULT_FILE
    echo "Average CPU Clock Speed: $AVG_FREQ MHz" >> $RESULT_FILE
    echo "" >> $RESULT_FILE
done
python3 benchmark_plot.py
echo "Benchmark completed. Results are saved in $RESULT_FILE."

