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

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
