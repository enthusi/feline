CC = gcc
CFLAGS = -O3 -ffast-math -fopenmp -g $(shell sdl-config --cflags) -std=c99
LDFLAGS = $(shell sdl-config --libs) -lm
TARGET = feline.bin
SOURCE = src/feline.c

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)

