// COMPILE: mpicc -o fft_serial fft_serial.c -lm
// RUN: ./fft_serial

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#define PI 3.1415926
#define Size 8192 //Problem Size

void run() {
    clock_t start = clock();
    double complex even[(Size / 2)], odd[(Size / 2)];
    double real_sum[Size], img_sum[Size];
    double table[Size][3];
    for(int i = 0; i < Size; i++) {
        table[i][0] = i;
        for(int j = 1; j < 3; j++) {
            table[i][j] = 1.0 / j;
        }
    }
    for (int k = 0; k < Size / 2; k++) {
        double real_even = 0.0, img_even = 0.0, read_odd = 0.0, img_odd = 0.0;
        for(int i = 0; i < Size/2; i++) {
            for (int j = 0; j < 2; j++) {
                double complex val1 = (table[2*i+j][1] + table[2*i+j][2] * I);
                double tmp = ((2*PI) * ((2*i+j)*k)) / Size;
                double complex val2 = (cos(tmp) - (sin(tmp)*I));
                if (j)  odd[i] = (val1 * val2);
                else    even[i] = (val1 * val2);
            }
        }

        for(int i = 0; i < Size / 2; i++) {
                real_even += creal(even[i]);
                img_even += cimag(even[i]);
                read_odd += creal(odd[i]);
                img_odd += cimag(odd[i]);
            }
            real_sum[k] = real_even + read_odd;
            img_sum[k]  = img_even + img_odd;
            real_sum[k + Size/2] = real_even - read_odd;
            img_sum[k + Size/2] = img_even - img_odd;
    }
    clock_t end = clock();
    printf("Time: %d ms\n", (end-start)*1000/CLOCKS_PER_SEC);
}

int main() {
    run();
    return 0;
}
