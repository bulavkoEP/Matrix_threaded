#include "io.h"

int read_matrix_from_file(double* matrix, int n, string filename) {
    ifstream in;
    in.open(filename);
    if (!in) {
        cout << "Error opening file" << endl;
        return -1;
    }
    double tmp;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (!(in >> tmp)) {
                cout << "Error reading file" << endl;
                return -1;
            }
        }
        for (int j = i; j < n; ++j) {
            if (!(in >> tmp)) {
                cout << "Error reading file" << endl;
                return -1;
            }
            matrix[ind(i, j)] = tmp;
        }
    }
    in.close();
    return 1;
}

void generate_matrix_from_formula(double* matrix, int n, int k) {
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            matrix[ind(i, j)] = fun(k, n, i, j);
        }
    }
}

void print(double* matrix, int n, int m, ostream& of) {
    for (int i = 0; i < min(n, m); ++i) {
        for (int j = 0; j < min(n, m); ++j) {
            of << matrix[ind(i, j)] << " ";
        }
        of << endl;
    }
}

double fun(int k, int n, int i, int j) {
    switch(k) {
        case 1: {
            return n - max(i, j) + 1;
        }
        case 2: {
            return max(i, j) + 1;
        }
        case 3: {
            return abs(i - j);
        }
        case 4: {
            return 1.0 / (i + j + 1);
        }
        default: {
            return 0;
        }
    }
} 