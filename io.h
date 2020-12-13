#ifndef MATR
#include "functions.h"
#define MATR
#endif
#include <iostream>
#include <fstream>

using namespace std;

double fun(int k, int n, int i, int j);

int read_matrix_from_file(double* matrix, int n, string filename);

void generate_matrix_from_formula(double* matrix, int n, int k);

void print(double* matrix, int n, int m, ostream& of);

int ind(int i, int j);