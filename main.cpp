#include <complex>
#include <iostream>
#include <valarray>
#include <omp.h>
#include <chrono>
#include <time.h>
#include <math.h>
typedef std::chrono::high_resolution_clock Clock;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
const double PI = 3.1415926;


void fft(CArray& x) {
    const size_t N = x.size();
    if (N <= 1) return;

    CArray even;
    CArray  odd;

#pragma omp parallel for
    for(int i = 0; i <= 1; i++) {
        CArray temp = x[std::slice(i, N / 2, 2)];
        fft(temp);
        if (i) odd = temp;
        else even = temp;
    }

    size_t k;
    if (N > pow(2,4)) {
    //#pragma omp parallel for num_threads(8) schedule(dynamic)
        for (k = 0; k < N / 2; ++k) {
            Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
            x[k] = even[k] + t;
            x[k + N / 2] = even[k] - t;
        }
    }
    else {
        for (k = 0; k < N / 2; ++k) {
            Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
            x[k] = even[k] + t;
            x[k + N / 2] = even[k] - t;
        }
    }

}

void fac() {
    int sum = 1;
    for (int i = 0; i < 10000; ++i) {
        sum *= i;
    }
}

void test1() {
    int k, N;
    N=100000;
#pragma omp parallel
#pragma omp for
    for (k = 0; k < N / 2; ++k)
    {
        fac();
    }
}

void test2() {
    int k, N;
    N=100000;
    Complex data[N];
    for (long k = 0; k < N; ++k) {
        data[k].real(1);
        data[k].imag(1);
    }
    CArray arr(data, N);
    CArray even;
    CArray  odd;

    auto start = Clock::now();

    for(int i = 0; i <= 1; i++) {
        CArray temp = arr[std::slice(i, N / 2, 2)];
    }
    auto end = Clock::now();
    std::cout<<"T0: "<< std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " nanoseconds \n\n";
    start = Clock::now();
    start = Clock::now();

#pragma omp parallel for
    for(int i = 0; i <= 1; i++) {
        CArray temp = arr[std::slice(i, N / 2, 2)];
    }
    end = Clock::now();
    std::cout<<"T1: "<< std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " nanoseconds \n\n";
    start = Clock::now();

}
int main()
{
    long size = pow(2, 15);
    srand (time(NULL));

    Complex data[size];
    for (long k = 0; k < size; ++k) {
        data[k].real(1);
        data[k].imag(1);
    }
    CArray arr(data, size);
    auto start = Clock::now();

    // forward fft
    // fft(arr);
    // test1();
    test2();

    auto end = Clock::now();

    std::cout<<"cost: "<< std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
                         << " nanoseconds \n\n";
    for (int i = 0; i < 5; ++i) {
        std::cout << arr[i] << std::endl;
    }
    return 0;
}
