#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#if SDLavailable
    #include <SDL.h>
#endif

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define DATA_PATH_PROCESSED "data/processed/"
#define DATA_PATH_RUNTIME_FILES "data/runtime_files/"
#define FELINE_VERSION "Feline v1.0.0"


/**
 * Generates a color value based on an input integer.
 *
 * Calculates a color represented in RGB format based on the value
 * of the input integer 'v'. Assigns different values to the red,
 * green, and blue components depending on the range in which 'v' falls:
 *
 * @param v: An integer value between 0 and 255 that determines the resulting RGB color.
 *
 * @return: An integer representing the RGB color, packed in the format 0xRRGGBB.
 */
int getcolor(int v) {
    int red, green, blue;
    red = 255;
    green = 255;
    blue = 255;
    if (v < 64) {
        red = 0;
        green = 1024 * (v) / 256;
    } else if (v < 128) {
        red = 0;
        blue = 256 + 1024 * (64 - v) / 256;
    } else if (v < 192) {
        red = 1024 * (v - 128) / 256;
        blue = 0;
    } else {
        green = 256 + 1024 * (192 - v) / 256;
        blue = 0;
    }
    return blue + green * 256 + red * 256 * 256;
}


/**
 * Counts the number of set bits (1s) in an integer.
 *
 * Takes an integer `v` as input and returns the number
 * of bits that are set to 1 in its binary representation.
 *
 * @param v The integer whose set bits are to be counted.
 * @return The number of set bits (1s) in the integer.
 */
inline int countbits(int v) {
    int c;
    for (c = 0; v; v >>= 1) { c += v & 1; }
    return c;
}


const float lines_first[11] = {6564.61, 4862.72, 4341.68, 4102.89, 3727.09, 4960.30, 6549.86, 6718.29, 3869.81, 1908.73,
                               1215.67};
const float lines_second[11] = {0, 0, 0, 0, 3729.88, 5008.24, 6585.27, 6732.67, 3968.53, 1906.68, 0};
const float scale = 10.0;
const float significance = 1.5;
float ignore_below = 1.5;
const int min_used = 2;
const float lmin = 4699.59;
const float lmax = 9350.84;
//how many frames do we store?
const int layer = 4;
const float redshift_step = 0.0002;
//const float max_match=20.0;
//for VERY weak ones
float max_match = 8.0;

//do we load a previously processed cube?

int previous = 0;
//this parameter de-values single super strong peaks in the S2N cube
//and thus favors multi-line-solutions


int main(int argc, char *argv[]) {
    setbuf(stdout, NULL);
    float *temp;

    float *res_i;
    float *res_z;
    float *res_t;
    float *res_u;

    float *prev_array;
#if SDLavailable
    SDL_Surface *screen;
    SDL_Window *window;
    SDL_Renderer *renderer;
    SDL_Event event;

    Uint8 *keys;
    Uint8 pixelr, pixelg, pixelb;
#endif
    int i, x, y;

    printf("%s\n", FELINE_VERSION);
    if (argc < 5) {
        printf("Syntax: %s zlow zhigh max_match ignore_below\n", argv[0]);
        return 1;
    }
    float zlow = atof(argv[1]);
    float zqso = atof(argv[2]);
    float max_match = atof(argv[3]);
    float ignore_below = atof(argv[4]);

    printf("Checking for emission line objects between z=%f and z=%f\n", zlow, zqso);

    FILE *f;
    FILE *image;
    temp = malloc(sizeof(float) * 4);
    if ((f = fopen(DATA_PATH_PROCESSED "raw_reordered_s2ncube.dat", "rb")) == NULL) {
        perror("No File found");
        free(temp);
        return -1;
    }

    int num = fread(temp, sizeof(float) * 4, 1, f);
    if (!num) {
        perror("Couldn't read from file");
        fclose(f);
        free(temp);
        return -1;
    }
    fclose(f);

    int size_header = 4;
    int dz = (int) (temp[0]);
    int dy = (int) (temp[1]);
    int dx = (int) (temp[2]);
    float lmin = (float) (temp[3]);
    float lmax = lmin + dz * 1.25;
    printf("%dx%dx%d (z,y,x) cube dimensions, start: %.2f end: %.2f\n", dz, dy, dx, lmin, lmax);

    if (argc==4){
        printf("Loading previously generated file '%s'\n" , argv[3]);
        FILE *prev_res_i_file;
        prev_array=malloc(sizeof(float)*dy*dx*layer);
        prev_res_i_file = fopen(argv[3], "rb");
        if (prev_res_i_file == NULL) {
            perror("Failed to open file");
            free(prev_res_i_file);
            free(temp);
            free(prev_array); // Ensure to free allocated memory to avoid memory leak
            return -1; // Or handle the error as appropriate
        }
        int r = fread(prev_array, (sizeof(float)*dy*dx*layer), 1, prev_res_i_file);
        if (!r) {
            perror("Failed to read from file");
            fclose(prev_res_i_file);
            free(prev_res_i_file);
            free(temp);
            free(prev_array); // Ensure to free allocated memory to avoid memory leak
            return -1; // Or handle the error as appropriate
        }
        fclose(prev_res_i_file);
        previous=1; //set proper flag
    }

    printf("Reading in full cube (%.1f MB)... ", (dx * dy * dz * sizeof(float) / 1048576.0));
    free(temp);
    temp = malloc(sizeof(float) * dz * dy * dx + sizeof(float) * size_header);
    if ((f = fopen(DATA_PATH_PROCESSED "raw_reordered_s2ncube.dat", "rb")) == NULL) {
        perror("No File found");
        free(temp);
        return -1;
    }

    num = fread(temp, (sizeof(float) * dz * dy * dx + sizeof(float) * size_header), 1, f);
    if (!num) {
        perror("Could not read from File!");
    }
    //printf("%d\n",sizeof(float));
    printf("DONE\n");

#if SDLavailable
    // Init SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fputs(SDL_GetError(), stderr);
        exit(1);
    }

    atexit(SDL_Quit);

    window = SDL_CreateWindow(FELINE_VERSION, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, dx * 4, dy + 1, SDL_WINDOW_SHOWN);
    if (window == NULL) {
        fputs(SDL_GetError(), stderr);
        exit(1);
    }

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (renderer == NULL) {
        fputs(SDL_GetError(), stderr);
        exit(1);
    }

// Access pixels directly from the window surface
    screen = SDL_CreateRGBSurface(0, dx * 4, dy + 1, 32, 0, 0, 0, 0);
    if (screen == NULL) {
        fputs(SDL_GetError(), stderr);
        exit(1);
    }
    SDL_Texture *texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, dx * 4, dy + 1);
    if (SDL_MUSTLOCK(screen)) {
        SDL_LockSurface(screen);
    }

    Uint32 *pixels = (Uint32 *)screen->pixels;
#endif


    res_i = malloc(sizeof(float) * dy * dx * layer);

    printf("Seeking for Emitters... ");
#if SDLavailable
#pragma omp parallel for schedule(static) shared(temp, dx, dy, dz, size_header, pixels, texture, res_i, prev_array, previous)  default(shared)
#else
#pragma omp parallel for schedule(static) shared(temp, dx, dy, dz, size_header, res_i, prev_array, previous)  default(shared)
#endif
    for (int i = 0; i < dy * dx; i++) {

        float v;

        int y = i / (dx);
        int x = i % (dx);
        int z;
        //printf("%d %d\n",y,x);
        float all_lines[20];
        int all_pix[20];
        //set to 0 as mask indicator by preprocessor
        //skip whole spaxel if first bin is 0, this might be because we masked it before
        //or its an area of NANs, i.e. in rotated cubes
        if (temp[0 + y * dx * dz + x * dz + size_header] == 0.0) {
            res_i[y * dx + x + dx * dy * 0] = 0.0;
            res_i[y * dx + x + dx * dy * 1] = 0.0;
            res_i[y * dx + x + dx * dy * 2] = 0;
            res_i[y * dx + x + dx * dy * 3] = 0;
            continue;
        }

        float sum = 0;
        //build up all templates
        //test all z-ranges
        //howto: generate array of n pixels
        //add up those n pixel values and store: template + z
        int template, toggle;
        int lines_ptr;
        float l1;
        float l2;
        float zmin;
        float zmax;
        float emission;
        float redshift;
        int pix;
        float best_redshift = 0;
        int best_template = 0;
        float best_sum = 0;
        int best_used = 0;
        best_template = 0;
        float avoid_z = 0;
        if (previous > 0) {
            if (!prev_array) {
                exit(EXIT_FAILURE); // Or handle the error as appropriate
            }
            // Now you can safely access prev_array
            avoid_z = prev_array[y * dx + x + dx * dy * 1]; // Accessing prev_array after ensuring it's initialized
        }
        for (template = 1; template < 512; template++) {
            zmin = zlow;
            zmax = zqso;
            lines_ptr = 0;
            toggle = template;

            for (int c = 0; c < 11; c++) {
                if (toggle & 1) {
                    l1 = lines_first[c];
                    l2 = lines_second[c];
                    //add doublets of available
                    if (l1 > 1.0) {
                        all_lines[lines_ptr] = l1;
                        lines_ptr++;
                    }
                    if (l2 > 1.0) {
                        all_lines[lines_ptr] = l2;
                        lines_ptr++;
                    }
                }
                toggle /= 2;
            }//optional lines loop

            if (lines_ptr < min_used) continue; //try next template

            for (int k = 0; k < lines_ptr; k++) {
                emission = all_lines[k];
                zmin = MAX(lmin / emission - 1, zmin);
                zmax = MIN(lmax / emission - 1, zmax);
            }
            //printf("\n%f %f\n",zmin,zmax);

            //now check full redshift range
            int samples = (int) ((zmax - zmin) / redshift_step);
            //printf("%f %f %d\n",zmin,zmax,samples);

            for (int s = 0; s < samples; s++) {
                redshift = zmin + redshift_step * s;
                if (fabs(redshift - avoid_z) < 0.05) {
                    continue;
                }
                //construct template
                sum = 0.0;

                for (int k = 0; k < lines_ptr; k++) {
                    emission = all_lines[k] * (redshift + 1);
                    z = (int) ((emission - lmin) / 1.25);
                    v = temp[z + y * dx * dz + x * dz + size_header];//*scale;
                    if (v < ignore_below) v = 0;

                    if (v > max_match) v = max_match + log(v); //crop intensity spikes
                    sum += v * scale;
                }//Emission loop

                sum = sum - lines_ptr * scale * significance;
                if (sum > best_sum) {
                    //printf("new best %d %d\n",best_sum,sum);
                    best_sum = sum;
                    best_redshift = redshift;
                    best_template = template;
                    best_used = lines_ptr;
                }
                //iter++;

            }//redshift loop

        }//templates loop


#if SDLavailable
        int width = dx * 4;


        pixels[y * width + x + dx * 0] = getcolor((((int) (best_sum)) >> 2));
        pixels[y * width + x + dx * 1] = getcolor(((best_used - min_used + 1) << 4));
        if ((((int) (best_sum)) >> 2) > 10)
            pixels[y * width + x + dx * 2] = getcolor((Uint8)(best_redshift / zqso * 256.0));
        else
            pixels[y * width + x + dx * 2] = 0xffffff;
        pixels[y * width + x + dx * 3] = getcolor((Uint8)(best_template));
#endif
        res_i[y * dx + x + dx * dy * 0] = best_sum;
        res_i[y * dx + x + dx * dy * 1] = best_redshift;
        res_i[y * dx + x + dx * dy * 2] = best_template;
        res_i[y * dx + x + dx * dy * 3] = best_used;

        int tid = omp_get_thread_num();
        if (tid == 0) {

            {
#if SDLavailable
                SDL_UpdateTexture(texture, NULL, pixels, dx * 4 * sizeof(Uint32));
                SDL_RenderClear(renderer);
                SDL_RenderCopy(renderer, texture, NULL, NULL);
                SDL_RenderPresent(renderer);
#endif
                if (i % 10 == 0 && (i * omp_get_num_threads()) < (dx*dy)) {
                    //printf("\rNumber of Threads:   %d", omp_get_num_threads());
                    //printf("\rIteration i:   %d", i);
                    printf("\rSeeking for Emitters... %.1f%%", ((omp_get_num_threads() * i * 100.0) / (dy * dx)));
                    fflush(stdout); // Ensure immediate output
                }
            }

        }



    }//xy-loop


    printf("\r                                 ");
    printf("\rSeeking for Emitters... 100%%\n");
    fflush(stdout);

    FILE *res_i_file;
    res_i_file = fopen(DATA_PATH_RUNTIME_FILES "float32_array_omp4.raw", "wb");
    fwrite(res_i, (sizeof(float) * dy * dx * layer), 1, res_i_file);


#if SDLavailable
    char *file = DATA_PATH_RUNTIME_FILES "map_omp4.bmp";
    SDL_RenderReadPixels(renderer, NULL, SDL_PIXELFORMAT_ARGB8888, screen->pixels, screen->pitch);
    if (SDL_SaveBMP(screen, file) != 0) {
        fprintf(stderr, "Could not write %s!\n", file);
    }

    SDL_RenderClear(renderer);
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);

    if (SDL_MUSTLOCK(screen)) {
        SDL_UnlockSurface(screen);
    }

#endif
    return 0;
}
