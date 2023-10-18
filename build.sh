#!/bin/sh
gcc -O3 --fast-math -fopenmp -g `sdl-config --cflags` -lm -o feline07.bin felinev07.c `sdl-config --libs` -std=c99 

