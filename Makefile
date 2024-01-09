CC = gcc
CFLAGS = -O3 -ffast-math -fopenmp -g $(shell sdl2-config --cflags) -std=c99
LDFLAGS = $(shell sdl2-config --libs) -lm
TARGET = feline.bin
SOURCE = src/feline.c

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
