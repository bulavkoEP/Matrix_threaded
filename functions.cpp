#include "functions.h"
#include <math.h>


int ind(int i, int j) {
    if (i < j) return j * (j + 1) / 2 + i;
	else return i * (i + 1) / 2 + j;
}

double mult_err(double* mat1, double* mat2, int n) {
    double norm = 0;
    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			double s = 0;
			for (int k = 0; k < n; ++k)
				s += mat1[ind(i, k)] * mat2[ind(k, j)];
            
            if (i == j)
                 s -= 1;

			norm += s * s;
		}
	}
	return sqrt(norm);
}

double norm(double* mat, int n) {
    double max = 0, sum = 0;
    for (int j = 0; j < n; ++j) {
        sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += abs(mat[ind(i, j)]);
        }
        if (sum > max) max = sum;
    }
    return max;
}