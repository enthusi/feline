CC = clang
CFLAGS = -O3 -ffast-math -fopenmp=libomp -g -std=c99
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


CUBELINK := "https://amused.univ-lyon1.fr/data/UDF/HUDF/download/DATACUBE_UDF-10.fits"
#CUBENAME := "DATACUBE_UDF-10.fits"
CUBENAME := "cube.fits"


ZLOW="0"
ZHIGH="1.9"
MAX_MATCH="20"
IGNORE_BELOW="7"
CORES := $(shell uname -s| awk '{if ($$0 == "Darwin") print "sysctl -n hw.physicalcpu"; else print "nproc"}' | sh)


all: $(TARGET)

$(TARGET): $(SOURCE)

	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS)


run:

	@if [ -f data/raw/$(CUBENAME) ]; then \
		echo "File exists";\
	else\
		echo "Downloading Cube File...";\
		wget --no-check-certificate $(CUBELINK) ;\
		mv $(CUBENAME) data/raw/$(CUBENAME);\
	fi

	@if [ ! -d "venv" ]; then \
		echo "Setting up environment..."; \
		python3 -m venv venv; \
	fi

	@. venv/bin/activate; pip install --prefer-binary -r requirements.txt > /dev/null
	@export PYTHONWARNING="ignore"

	@echo "Starting preprocessing of the cubefile..."
	@echo ""
	@echo "apply median filter to 'flat' out the data (emission lines remain)..."
	@. venv/bin/activate && \
	python src/preprocessing/median-filter-cube.py data/raw/$(CUBENAME) --signalHDU=1 --varHDU=2 --num_cpu=$(CORES) --width=151 --output=data/processed/med_filt.fits > /dev/null
	@echo ""
	@echo "filter the data cube with a 'typical line' template in spatial dimension first..."
	@. venv/bin/activate && \
	python src/preprocessing/lsd_cc_spatial.py --input=data/processed/med_filt.fits --SHDU=1 --NHDU=2 --threads=$(CORES) --gaussian --lambda0=7050 -pc 0.7 --classic --output=data/processed/spatial_cc.fits --overwrite > /dev/null
	@echo "filter the data cube with a 'typical line' template in spectral dimension..."
	@. venv/bin/activate && \
	python src/preprocessing/lsd_cc_spectral.py --input=data/processed/spatial_cc.fits --threads=$(CORES) --FWHM=250 --SHDU=1 --NHDU=2 --classic --output=data/processed/spectral_cc.fits --overwrite > /dev/null
	@echo "construct a signal-to-noise cube..."
	@. venv/bin/activate && \
	python src/preprocessing/s2n-cube.py --input=data/processed/spectral_cc.fits --output=data/processed/s2n_v250.fits --clobber --NHDU=2 --SHDU=1 > /dev/null

	@echo "deleting tmp files..."
	@rm data/processed/spatial_cc.fits
	@rm data/processed/spectral_cc.fits

	@echo "Create Masking Plot and transpose Cube for better Cache Access..."
	@. venv/bin/activate ; \
	python -m src.preprocessing.combination $(CUBENAME) s2n_v250.fits
	@echo "Starting FELINE..."
	@if [ -e feline.bin ]; then \
		rm -f feline.bin; \
	else \
		echo ""; \
	fi
	@$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS)
	./feline.bin $(ZLOW) $(ZHIGH) $(MAX_MATCH) $(IGNORE_BELOW)

	@echo "Starting Postprocessing and creating PDF..."
	@. venv/bin/activate ; \
	python -m src.postprocessing.detect_objects s2n_v250.fits ; \
	python -m src.postprocessing.create_final_plots $(CUBENAME) s2n_v250.fits sorted_catalog.txt med_filt.fits J0014m0028 ; \
	python -m src.postprocessing.create_pdf

cuda:

	@if [ -f data/raw/$(CUBENAME) ]; then \
		echo "File exists";\
	else\
		echo "Downloading Cube File...";\
		wget --no-check-certificate $(CUBELINK) ;\
		mv $(CUBENAME) data/raw/$(CUBENAME) ;\
	fi

	@if [ ! -d "venv" ]; then \
		echo "Setting up environment..."; \
		python3 -m venv venv; \
	fi

	@. venv/bin/activate; pip install --prefer-binary -r requirements.txt > /dev/null
	@export PYTHONWARNING="ignore"

	@echo "Starting preprocessing of the cubefile..."
	@echo ""
	@echo "apply median filter to 'flat' out the data (emission lines remain)..."
	@. venv/bin/activate && \
	python src/preprocessing/median-filter-cube.py data/raw/$(CUBENAME) --signalHDU=1 --varHDU=2 --num_cpu=$(CORES) --width=151 --output=data/processed/med_filt.fits > /dev/null
	@echo ""
	@echo "filter the data cube with a 'typical line' template in spatial dimension first..."
	@. venv/bin/activate && \
	python src/preprocessing/lsd_cc_spatial.py --input=data/processed/med_filt.fits --SHDU=1 --NHDU=2 --threads=$(CORES) --gaussian --lambda0=7050 -pc 0.7 --classic --output=data/processed/spatial_cc.fits --overwrite > /dev/null
	@echo "filter the data cube with a 'typical line' template in spectral dimension..."
	@. venv/bin/activate && \
	python src/preprocessing/lsd_cc_spectral.py --input=data/processed/spatial_cc.fits --threads=$(CORES) --FWHM=250 --SHDU=1 --NHDU=2 --classic --output=data/processed/spectral_cc.fits --overwrite > /dev/null
	@echo "construct a signal-to-noise cube..."
	@. venv/bin/activate && \
	python src/preprocessing/s2n-cube.py --input=data/processed/spectral_cc.fits --output=data/processed/s2n_v250.fits --clobber --NHDU=2 --SHDU=1 > /dev/null

	@echo "deleting tmp files..."
	@rm data/processed/spatial_cc.fits
	@rm data/processed/spectral_cc.fits

	@echo "Create Masking Plot and transpose Cube for better Cache Access..."
	@. venv/bin/activate ; \
	python -m src.preprocessing.combination $(CUBENAME) s2n_v250.fits
	@echo "Starting FELINE..."
	@if [ -e feline.bin ]; then \
                rm -f feline.bin; \
        else \
                echo ""; \
        fi
	nvcc -O3 --use_fast_math -o feline.bin src/feline.cu
	./feline.bin $(ZLOW) $(ZHIGH) $(MAX_MATCH) $(IGNORE_BELOW)

	@echo "Starting Postprocessing and creating PDF..."
	@. venv/bin/activate ; \
	python -m src.postprocessing.detect_objects s2n_v250.fits ; \
	python -m src.postprocessing.create_final_plots $(CUBENAME) s2n_v250.fits sorted_catalog.txt med_filt.fits J0014m0028 ; \
	python -m src.postprocessing.create_pdf

clean:

	rm -f $(TARGET)
	find data/raw -type f ! -name '*.fits' ! -name '.gitkeep' -delete
	find data/processed ! -name '.gitkeep' -type f -delete
	find src/postprocessing -type f \( -name '*.txt' -o -name '*.fits' -o -name '*.png' -o -name '*.log' -o ! -name "*.*" \) -delete
	find data/runtime_files ! -name '.gitkeep' -type f -delete

debug:

