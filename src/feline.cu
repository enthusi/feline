#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define DATA_PATH "data/raw/"
#define DATA_PATH_PROCESSED "data/processed/"

const float lines_first[11] = {6564.61, 4862.72, 4341.68, 4102.89, 3727.09, 4960.30, 6549.86, 6718.29, 3869.81, 1908.73,
                               1215.67};
const float lines_second[11] = {0, 0, 0, 0, 3729.88, 5008.24, 6585.27, 6732.67, 3968.53, 1906.68, 0};
const float scale = 10.0;
const float significance = 1.5;
float ignore_below = 1.5;
const int min_used = 2;
float lmin = 4699.59;
float lmax = 9350.84;
//how many frames do we store?
const int layer = 4;
const float redshift_step = 0.0002;
//const float max_match=20.0;
//for VERY weak ones
float max_match = 8.0;

//do we load a previously processed cube?
int previous = 0;


__global__ void seekEmitters(float *temp, float *res_i, int dz, int dy, int dx, int size_header, float zlow, float zqso, float max_match, float ignore_below, float lmin, float lmax, float redshift_step, int min_used, const float* lines_first, const float* lines_second) {

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= dx*dy){return;} // Check if index is out of range
    float v;
    int y = index / dx;
    int x = index % dx;
    int z;

    float all_lines[20];

    if (temp[0 + y * dx * dz + x * dz + size_header] == 0.0) {
        res_i[y * dx + x + dx * dy * 0] = 0.0;
        res_i[y * dx + x + dx * dy * 1] = 0.0;
        res_i[y * dx + x + dx * dy * 2] = 0;
        res_i[y * dx + x + dx * dy * 3] = 0;
        return;
    }

    float sum = 0;
    int template_id, toggle;
    int lines_ptr;
    float l1;
    float l2;
    float zmin;
    float zmax;
    float emission;
    float redshift;
    float best_redshift = 0;
    int best_template = 0;
    float best_sum = 0;
    int best_used = 0;
    float avoid_z = 0;
    for (template_id = 1; template_id < 512; template_id++) {
        zmin = zlow;
        zmax = zqso;
        lines_ptr = 0;
        toggle = template_id;
        for (int c = 0; c < 11; c++) {
            if ((toggle & 0x1) == 1) {
                l1 = lines_first[c];
                l2 = lines_second[c];
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
        }
        if (lines_ptr < min_used) continue;
        for (int k = 0; k < lines_ptr; k++) {
            emission = all_lines[k];
            zmin = MAX(lmin / emission - 1, zmin);
            zmax = MIN(lmax / emission - 1, zmax);
        }
        int samples = (int) ((zmax - zmin) / redshift_step);

        for (int s = 0; s < samples; s++) {
            redshift = zmin + redshift_step * s;
            if (abs(redshift - avoid_z) < 0.05) {continue;}
            sum = 0.0;
            for (int k = 0; k < lines_ptr; k++) {
                emission = all_lines[k] * (redshift + 1);
                z = (int) ((emission - lmin) / 1.25);
                v = temp[z + y * dx * dz + x * dz + size_header];
                if (v < ignore_below) v = 0;
                if (v > max_match) v = max_match + logf(v);
                sum += v * scale;
            }
            sum = sum - lines_ptr * scale * significance;
            if (sum > best_sum) {
                best_sum = sum;
                best_redshift = redshift;
                best_template = template_id;
                best_used = lines_ptr;
            }
        }
    }
    res_i[y * dx + x + dx * dy * 0] = best_sum;
    res_i[y * dx + x + dx * dy * 1] = best_redshift;
    res_i[y * dx + x + dx * dy * 2] = best_template;
    res_i[y * dx + x + dx * dy * 3] = best_used;
}

int main(int argc, char *argv[]) {
    setbuf(stdout, NULL);
    float *temp;

    float *res_i;

    printf("Feline v1.0.0 - 2024/09/12\n");
    printf("Martin Wendt\n");
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
    temp = (float*) malloc(sizeof(float) * 4);
    if ((f = fopen(DATA_PATH_PROCESSED "raw_reordered_s2ncube.dat", "rb")) == NULL) {
        perror("No File found");
        return -1;
    }

    int num = fread(temp, sizeof(float) * 4, 1, f);
    if (!num) {
        perror("Couldn't read from file");
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
    printf("Reading in full cube (%.1f MB)... ", (dx * dy * dz * sizeof(float) / 1048576.0));
    temp = (float*) malloc(sizeof(float) * dz * dy * dx + sizeof(float) * size_header);
    if ((f = fopen(DATA_PATH_PROCESSED "raw_reordered_s2ncube.dat", "rb")) == NULL) {
        perror("No File found");
        return -1;
    }
    num = fread(temp, (sizeof(float) * dz * dy * dx + sizeof(float) * size_header), 1, f);
    if (!num) {
        perror("Could not read from File!");
    }
    printf("\n%d SIZE OF DX * DY\n", dx*dy);
    //printf("%d\n",sizeof(float));
    printf("DONE\n");
    res_i = (float*) malloc(sizeof(float) * dy * dx * layer);

    printf("Seeking for Emitters... ");

    /*
     * CUDA Part
     *
     *
     *
     */

    // Get GPU memory information
    cudaSetDevice(0);
    // CUDA Runtime Version
    int runtimeVersion;
    cudaRuntimeGetVersion(&runtimeVersion);
    printf("CUDA Runtime Version: %d.%d\n", runtimeVersion / 1000, (runtimeVersion % 100) / 10);

    // CUDA Driver Version
    int driverVersion;
    cudaDriverGetVersion(&driverVersion);
    printf("CUDA Driver Version: %d.%d\n", driverVersion / 1000, (driverVersion % 100) / 10);
    float *d_temp;
    float *d_res_i;
    cudaError_t c1 = cudaMalloc((void **)&d_temp, sizeof(float) * dz * dy * dx + sizeof(float) * size_header);
    cudaError_t c2 = cudaMemcpy(d_temp, temp, sizeof(float) * dz * dy * dx + sizeof(float) * size_header, cudaMemcpyHostToDevice);
    cudaError_t c3 = cudaMalloc((void **)&d_res_i, sizeof(float) * dy * dx * layer);
    cudaError_t c4 = cudaMemcpy(d_res_i, res_i, sizeof(float)*dy*dx*layer, cudaMemcpyHostToDevice);
    float *d_lines_first, *d_lines_second;

    // Allocate memory on the device
    cudaMalloc((void **)&d_lines_first, sizeof(lines_first));
    cudaMalloc((void **)&d_lines_second, sizeof(lines_second));

    // Copy data from host to device
    cudaMemcpy(d_lines_first, lines_first, sizeof(lines_first), cudaMemcpyHostToDevice);
    cudaMemcpy(d_lines_second, lines_second, sizeof(lines_second), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    printf("\n%s\n",cudaGetErrorString(c1));
    printf("\n%s\n",cudaGetErrorString(c2));
    printf("\n%s\n",cudaGetErrorString(c3));
    printf("\n%s\n", cudaGetErrorString(c4));

    // Define grid and block dimensions
    // right now for the cube with 6004 pixels manually 
    // we want blockdim and gridDim to be as close as possible to each other

    dim3 blockDim(64, 1, 1);
    dim3 gridDim((dx*dy + blockDim.x)/blockDim.x, 1, 1);
    printf("\n%u\n", gridDim.x);
    cudaError_t err1 = cudaGetLastError();
    if (err1 != cudaSuccess){
        printf("\nhere Error: %s\n", cudaGetErrorString(err1));
    }
    // Launch CUDA kernel
    seekEmitters<<<gridDim, blockDim>>>(d_temp, d_res_i, dz, dy, dx, size_header, zlow, zqso, max_match, ignore_below, lmin, lmax, redshift_step, min_used, d_lines_first, d_lines_second);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("here2 Error: %s\n", cudaGetErrorString(err));
    }
    free(res_i);
    res_i = (float*) malloc(sizeof(float)*dy*dx*layer);
    cudaError_t h2d = cudaMemcpy(res_i, d_res_i, sizeof(float) * dy * dx * layer, cudaMemcpyDeviceToHost);
    printf("\n%s\n", cudaGetErrorString(h2d));

    /*
     * End of CUDA Part
     *
     *
     */
    printf("%dx%dx%d (z,y,x) cube dimensions, start: %.2f end: %.2f\n", dz, dy, dx, lmin, lmax);

    FILE *res_i_file;
    res_i_file = fopen("float32_array_omp4.raw", "wb");
    fwrite(res_i, (sizeof(float) * dy * dx * layer), 1, res_i_file);

    free(temp);
    free(res_i);
    fclose(res_i_file);
    return 0;
}

