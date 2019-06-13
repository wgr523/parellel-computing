#include <iostream>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <ctime>
#include <cmath>
#define lint 8
#define LINT (1<<lint)
#define LIMIT 8800001

using namespace std;
typedef std::chrono::high_resolution_clock Clock;
const double PI = 3.1415926;

unsigned char rev[LINT];

double re[LIMIT], im[LIMIT];

void compute_rev() {
    for (int i=0; i<LINT; i++) {
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (lint - 1));
    }
}

unsigned int reverse(const unsigned int idx, const unsigned int n) {
    if (n <= 8) {
        return (unsigned int)rev[(idx & 0xFF)] >> (8 - n);
    } else if (n <= 16) {
        return (unsigned int)(rev[((idx >> 8) & 0xFF)] | (rev[(idx & 0xFF)] << 8) ) >> (16 - n);
    } else if (n <= 24) {
        return (unsigned int)(rev[((idx >> 16) & 0xFF)] | (rev[((idx >> 8) & 0xFF)] << 8) | (rev[(idx  & 0xFF)] << 16) ) >> (24 - n);
    } else { return 0;}
}

void reverse_input(double* re, double* im, const unsigned int n) {
    double tmp;
    unsigned int j;
#pragma omp parallel for num_threads(4)
    for (int i=0; i< (1<<n); i++) {
        j = reverse(i, n);
        if (i < j) {
            tmp = re[i];
            re[i] = re[j];
            re[j] = tmp;
            tmp = im[i];
            im[i] = im[j];
            im[j] = tmp;
        }
    }
}
        
int main() {
    compute_rev();

    double input_re, input_im;
    ifstream fin;
    fin.open("rand23.txt");
    unsigned int cnt = 0;
    while (fin>>input_re>>input_im) {
        re[cnt] = input_re;
        im[cnt] = input_im;
        cnt++;
    }
    fin.close();
    cout << cnt << endl;
    unsigned int n =  (unsigned int)ceil(log2((double)cnt));
    auto start = Clock::now();
    reverse_input(re,im,n);
    auto end = Clock::now();
    cout<<"time: "<< chrono::duration_cast<chrono::nanoseconds>(end - start).count() << " nanoseconds" << endl;

    for (int i=cnt-1;i<cnt;i++) {
        cout<<re[i]<<' '<<im[i]<<endl;
    }
    return 0;
}
