// COMPILE: mpicc -o fft_parallel fft_parallel.c -lm
// RUN: mpirun -np 8 ./fft_parallel

#include <stdio.h>
#include <mpi.h>
#include <complex.h>
#include <math.h>

#define PI 3.1415926
#define Size 65536//Problem Size

double real_sum[Size], img_sum[Size];
double table[Size][3];
double sin_table[Size], cos_table[Size];
double complex e_table[Size];

void run() {

    int my_rank, comm_size;
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    double start, end;
    if(my_rank == 0) {
        FILE * fin;
        char filename[100];
        sprintf(filename, "%d.txt", Size);
        fin = fopen(filename, "r");
        for(int i = 0; i < Size; i++) {
            table[i][0] = i;
            fscanf(fin, "%lf%lf", &table[i][1], &table[i][2]);
        }
        fclose(fin);
        start = MPI_Wtime();
    }

    for (int k = 0; k < Size; k++) {
        //compute sin, cos table
        double tmp = 2*PI * k / Size;
        sin_table[k] = sin(tmp);
        cos_table[k] = cos(tmp);
        e_table[k] = cos_table[k] - (sin_table[k]*I);
    }

    //double complex even[(Size / comm_size / 2)], odd[(Size / comm_size / 2)];
    //double complex even_leader[ (Size / comm_size / 2) * comm_size], odd_leader[ (Size / comm_size / 2) * comm_size];
    double rev_table[(Size / comm_size)][3];

    int send_receive_count = (Size / comm_size) * 3;
    MPI_Scatter(table, send_receive_count, MPI_DOUBLE, rev_table,
            send_receive_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //printf("%d, %lf, %lf\n", my_rank, rev_table[0][0], rev_table[0][1]);
    for (int k = 0; k < Size / 2; k++) {
        double real_even_leader[comm_size], img_even_leader[comm_size], real_odd_leader[comm_size], img_odd_leader[comm_size];
        double real_even = 0.0, img_even = 0.0, real_odd = 0.0, img_odd = 0.0;
        for(int i = 0; i < (Size/comm_size); i++) {
            int idx = (int)rev_table[i][0];
            int tmp = (idx * k ) % Size;
            double complex val1 = (rev_table[i][1] + rev_table[i][2] * I);
            //printf("k:%d, rank: %d, input: %lf + %lfj\n", k, my_rank, creal(val1), cimag(val1));
            val1 *= e_table[tmp];
            //printf("k:%d, rank: %d, res: %lf + %lfj\n", k, my_rank, creal(val1), cimag(val1));
            if (idx & 1) {
                real_odd += creal(val1);
                img_odd += cimag(val1);
            }
            else {
                real_even += creal(val1);
                img_even += cimag(val1);
            }
        }
        //printf("k:%d, rank: %d, %lf, %lf\n", k, my_rank, real_even, real_odd);
        MPI_Gather(&real_even, 1, MPI_DOUBLE, real_even_leader, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&real_odd, 1, MPI_DOUBLE, real_odd_leader, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&img_even, 1, MPI_DOUBLE, img_even_leader, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&img_odd, 1, MPI_DOUBLE, img_odd_leader, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(my_rank == 0) {
            double real_even_sum = 0.0, img_even_sum = 0.0, real_odd_sum = 0.0, img_odd_sum = 0.0;
            for(int i = 0; i < comm_size; i++) {
                //printf("k: %d, rank: root, %lf, %lf\n", k, real_even_leader[i], real_odd_leader[i]);
                real_even_sum += real_even_leader[i];
                real_odd_sum += real_odd_leader[i];
                img_even_sum += img_even_leader[i];
                img_odd_sum += img_odd_leader[i];
            }
            real_sum[k] = real_even_sum + real_odd_sum;
            img_sum[k]  = img_even_sum + img_odd_sum;
            real_sum[k + Size/2] = real_even_sum - real_odd_sum;
            img_sum[k + Size/2] = img_even_sum - img_odd_sum;
        }
    }
    // Synchronize before timing
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0) {
        end = MPI_Wtime();
        printf("Time: %f ms\n", (end-start)*1000);
        FILE * fout;
        fout = fopen("out.txt", "w");
        for (int k = 0; k < Size ; k++) {
            fprintf( fout, "%lf %lf\n", real_sum[k], img_sum[k]);
        }
        fclose(fout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

int main() {
    run();
    return 0;
}
