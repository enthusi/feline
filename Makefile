CC = gcc
CFLAGS = -O3 -ffast-math -fopenmp -g -std=c99
LDFLAGS = -lm
TARGET = feline.bin
SOURCE = src/feline.c
#CUBE_LINK = "martinwendt.de/cube.fits"
#CUBE_NAME = "cube.fits"
CUBE_LINK = "https://amused.univ-lyon1.fr/data/UDF/HUDF/download/DATACUBE_UDF-10.fits"
CUBE_NAME = "DATACUBE_UDF-10.fits"

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

all: clean $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS)



# whole runthrough
#

download_file:
	@if [ -f $(CUBE_NAME) ]; then \
		echo "File exists";\
	else\
		echo "File dows not exist";\
		wget --no-check-certificate $(CUBE_LINK);\
	fi

clean: 
	rm -f $(TARGET)

