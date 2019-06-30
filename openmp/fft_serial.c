// COMPILE: gcc -o fft_serial fft_serial.c -lm
// RUN: ./fft_serial

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#define PI 3.1415926
#define Size 4//Problem Size

double real_sum[Size], img_sum[Size];
double sin_table[Size], cos_table[Size];
double complex e_table[Size];
double table[Size][3];

void run() {
    FILE * fin, * fout;
    char filename[100];
    sprintf(filename, "%d.txt", Size);
    fin = fopen(filename, "r");
    for(int i = 0; i < Size; i++) {
        table[i][0] = i;
        fscanf(fin, "%lf%lf", &table[i][1], &table[i][2]);
    }
    fclose(fin);
    clock_t start = clock();
    for (int k = 0; k < Size; k++) {
        //compute sin, cos table
        double tmp = 2*PI * k / Size;
        sin_table[k] = sin(tmp);
        cos_table[k] = cos(tmp);
        e_table[k] = cos_table[k] - (sin_table[k]*I);
    }
    for (int k = 0; k < Size / 2; k++) {
        double real_even = 0.0, img_even = 0.0, real_odd = 0.0, img_odd = 0.0;
        for(int i = 0; i < Size; i++) {
            double complex val1 = (table[i][1] + table[i][2] * I);
            int tmp = (i*k) % Size;
            val1 *= e_table[tmp];
            if (i & 1) {
                real_odd += creal(val1);
                img_odd += cimag(val1);
            }
            else {
                real_even += creal(val1);
                img_even += cimag(val1);
            }
        }

        real_sum[k] = real_even + real_odd;
        img_sum[k]  = img_even + img_odd;
        real_sum[k + Size/2] = real_even - real_odd;
        img_sum[k + Size/2] = img_even - img_odd;
    }
    clock_t end = clock();
    printf("Time: %ld ms\n", (end-start)*1000/CLOCKS_PER_SEC);
    fout = fopen("out.txt", "w");
    for (int k = 0; k < Size ; k++) {
        fprintf( fout, "%lf %lf\n", real_sum[k], img_sum[k]);
    }
    fclose(fout);
}

int main() {
    run();
    return 0;
}
