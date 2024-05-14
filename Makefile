CC = gcc
CFLAGS = -O3 -ffast-math -fopenmp -g -std=c99
LDFLAGS = -lm
TARGET = feline.bin
SOURCE = src/feline.c



# Check for SDL2 availability
SDL_CONFIG := sdl2-config
SDL2_CFLAGS := $(shell $(SDL_CONFIG) --cflags 2>/dev/null)
SDL2_LIBS := $(shell $(SDL_CONFIG) --libs 2>/dev/null)

ifeq ($(SDL2_CFLAGS)$(SDL2_LIBS),)
    SDLavailable = 0
else
    SDLavailable = 1
endif

# Add the SDL2 availability macro to the compilation flags
CFLAGS += $(SDL2_CFLAGS)
LDFLAGS += $(SDL2_LIBS)

# Define the SDLavailable macro for use in your source code
CFLAGS += -D SDLavailable=$(SDLavailable)

#CUBELINK := "martinwendt.de/cube.fits"
#CUBENAME := "cube.fits"
CUBE_LINK := "https://amused.univ-lyon1.fr/data/UDF/HUDF/download/DATACUBE_UDF-10.fits"
CUBE_NAME := "DATACUBE_UDF-10.fits"

ZLOW="0"
ZHIGH="1.9"
MAX_MATCH="20"
IGNORE_BELOW="7"
                                                                                      



all: clean $(TARGET)                                   

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS)


run:

	@echo "Downloading Cube File..."
	@if [ -f $(CUBE_NAME) ]; then \
		echo "File exists";\
	else\
		echo "File does not exist";\
		wget --no-check-certificate $(CUBE_LINK) ;\
	fi

	CUBEFILE:=$(shell realpath $(CUBENAME))
	CUBEFILENAME=$(CUBENAME)
	CORES=$(shell nproc)


	@echo "Setting up environment..."
	python3 -m venv venv
	@. venv/bin/activate; pip install -Ur requirements.txt > /dev/null
	export PYTHONWARNING="ignore"

	@echo "Starting preprocessing..."
	@echo ""
	@echo "Create med_filt.fits..."
	@. venv/bin/activate && \
	python src/preprocessing/median-filter-cube.py $(CUBEFILE) --signalHDU=1 --varHDU=2 --num_cpu=$(CORES) --width=151 --output=data/processed/med_filt.fits > /dev/null
	@echo ""
	@echo "Create spatial_cc.fits..."
	@. venv/bin/activate && \
	python src/preprocessing/lsd_cc_spatial.py --input=data/processed/med_filt.fits --SHDU=1 --NHDU=2 --threads=$(CORES) --gaussian --lambda0=7050 -pc 0.7 --classic --output=data/processed/spatial_cc.fits --overwrite > /dev/null
	@echo "Create spectral_cc.fits..."
	@. venv/bin/activate && \
	python src/preprocessing/lsd_cc_spectral.py --input=data/processed/spatial_cc.fits --threads=$(CORES) --FWHM=250 --SHDU=1 --NHDU=2 --classic --output=data/processed/spectral_cc.fits --overwrite > /dev/null
	@echo "Create s2n_v250.fits..."
	@. venv/bin/activate && \
	python src/preprocessing/s2n-cube.py --input=data/processed/spectral_cc.fits --output=data/processed/s2n_v250.fits --clobber --NHDU=2 --SHDU=1 > /dev/null

	@echo "deleting tmp files..."
	rm data/processed/spatial_cc.fits
	rm data/processed/spectral_cc.fits
	cp $(CUBEFILE) data/raw/

	@echo "Create Masking Plot and transpose Cube for better Cache Access..."
	@. venv/bin/activate ; \
	cd src/preprocessing ; \
	python combination.py $(CUBEFILENAME) s2n_v250.fits
	echo "Starting FELINE..."

	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS)
	./feline.bin $(ZLOW) $(ZHIGH) $(MAX_MATCH) $(IGNORE_BELOW)

	@echo "Starting Postprocessing and creating PDF..."
	# Detect objects and plot results
	@. venv/bin/activate ; \
	cd src/postprocessing || exit ; \
	python detect_objects.py s2n_v250.fits ; \
	python create_final_plots.py $(CUBEFILENAME) s2n_v250.fits sorted_catalog.txt med_filt.fits J0014m0028 ; \
	python create_pdf.py


clean:

	rm -f $(TARGET)
	find data/raw ! -name '.gitkeep' -type f -delete
	find data/processed ! -name '.gitkeep' -type f -delete
	rm -f *.bmp *.raw
	find src/postprocessing -type f \( -name '*.txt' -o -name '*.fits' -o -name '*.png' -o -name '*.log' -o ! -name "*.*" \) -delete

