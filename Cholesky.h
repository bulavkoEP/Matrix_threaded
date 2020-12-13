#ifndef MATR
#include "functions.h"
#define MATR
#endif

#include <math.h>
#include <iostream>

using namespace std;

int get_inverse(double* mat, double* res, double* d, int n, int n_threads, int thread_i, int* res_code);