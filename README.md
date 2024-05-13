# Find Emmision Lines (FELINE) #
<img src="https://github.com/enthusi/feline/assets/115990442/3f7e12a4-6ca2-4a22-a5f0-b358296e80b2" width="390" height="300" />

## 1. Clone Repository and Install Requirements ##
### Clone Repository ###
```bash
git clone git@github.com:enthusi/feline.git
```
### Install Requirements ###
If not already done create a Python3 Virtual Environment
```bash
cd feline
```
```bash
python3 -m venv venv
```
```bash
source venv/bin/activate
```
install requirements:
```bash
pip install -r requirements.txt
```
## 2. Usage ##
### Preprocessing (LSDcat: BSD-3 License Christian Herenz) ###
set Path of Cubefile which should be used:
```bash
export CUBEFILE=$(realpath <cubefile>.fits)
```
set the Number of Cores that should be used for preprocessing (default: 4)
```bash
export CORES=<num_cores>
```
do all preprocessing
```bash
sh preprocessing.sh
```
### some tweaking after the preprocessing ###
```bash
cd src/postprocessing
```
create optional masking plot 
```bash
python create_masking_plot.py <cubefile>.fits
```
transpose the cube data so for Cache Optimization
```bash
python transpose_cube.py s2n_v250.fits none
```

### Compile and run FELINE ###
(go back to project root)
```bash
make
```
```bash
./feline.bin <zlow> <zhigh> <max_match> <ignore_below>
```
### detect objects and plot results ###
```bash
cd src/postprocessing
```
```bash
python detect_objects.py s2n_v250.fits > catalog.txt
```
```bash
sort -rn -k5 catalog.txt > sorted_catalog.txt
```
```bash
python create_final_plots.py <cubefile>.fits s2n_v250.fits sorted_catalog.txt med_filt.fits J0014m0028
```
<<Your results are saved to data/final_plots/*>>
### Clean Everything for next run (Optional) ###
```bash
sh cleanup.sh
```





