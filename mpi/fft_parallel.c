// COMPILE: mpicc -o fft_parallel fft_parallel.c -lm
// RUN: mpirun -np 8 ./fft_parallel

#include <stdio.h>
#include <mpi.h>
#include <complex.h>
#include <math.h>

#define PI 3.1415926
#define Size 8192 //Problem Size

void run() {
    int my_rank, comm_size;
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    double start, end;
    if(my_rank == 0) {
        start = MPI_Wtime();
    }

    double complex even[(Size / comm_size / 2)], odd[(Size / comm_size / 2)];
    double complex even_leader[ (Size / comm_size / 2) * comm_size], odd_leader[ (Size / comm_size / 2) * comm_size];
    double real_sum[Size], img_sum[Size];
    double table[Size][3];
    double rev_table[(Size / comm_size)][3];
    for(int i = 0; i < Size; i++) {
        table[i][0] = i;
        for(int j = 1; j < 3; j++) {
            table[i][j] = 1.0 / j;
        }
    }
    int send_receive_count = (Size / comm_size) * 3;
    MPI_Scatter(table, send_receive_count, MPI_DOUBLE, rev_table,
                send_receive_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int k = 0; k < Size / 2; k++) {
        double real_even = 0.0, img_even = 0.0, read_odd = 0.0, img_odd = 0.0;
        for(int i = 0; i < (Size/comm_size)/2; i++) {
            for (int j = 0; j < 2; j++) {
                double tmp = 0.0;
                int val0 = my_rank * rev_table[2*i+j][0];
                double complex val1 = (rev_table[2*i+j][1] + rev_table[2*i+j][2] * I);
                if(my_rank == 0) {
                    tmp = ((2*PI) * ((2*i)*k)) / Size;
                } else {
                    tmp = ((2*PI) * ((val0)*k)) / Size;
                }
                double complex val2 = (cos(tmp) - (sin(tmp)*I));
                if (j)  odd[i] = (val1 * val2);
                else    even[i] = (val1 * val2);
            }
        }
        MPI_Gather(even, (Size / comm_size / 2), MPI_DOUBLE_COMPLEX, even_leader, (Size / comm_size / 2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
        MPI_Gather(odd, (Size / comm_size / 2), MPI_DOUBLE_COMPLEX, odd_leader, (Size / comm_size / 2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

        if(my_rank == 0) {
            for(int i = 0; i < (Size / comm_size / 2) * comm_size; i++) {
                real_even += creal(even_leader[i]);
                img_even += cimag(even_leader[i]);
                read_odd += creal(odd_leader[i]);
                img_odd += cimag(odd_leader[i]);
            }
            real_sum[k] = real_even + read_odd;
            img_sum[k]  = img_even + img_odd;
            real_sum[k + Size/2] = real_even - read_odd;
            img_sum[k + Size/2] = img_even - img_odd;
        }
    }
    if(my_rank == 0) {
        end = MPI_Wtime();
        printf("Time: %f ms\n", (end-start)*1000);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

int main() {
    run();
    return 0;
}
