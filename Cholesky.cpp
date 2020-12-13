#include "Cholesky.h"
#include "io.h"


#define EPS 1e-15

void synchronize(int total_threads) {
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;

	pthread_mutex_lock(&mutex);

	threads_in++;
	if (threads_in >= total_threads)
	{
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	} else
		while (threads_in < total_threads)
			pthread_cond_wait(&condvar_in,&mutex);

	threads_out++;
	if (threads_out >= total_threads)
	{
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	} else
		while (threads_out < total_threads)
			pthread_cond_wait(&condvar_out,&mutex);

	pthread_mutex_unlock(&mutex);
}

/*
int ind(int i, int j) {
    if (i < j) return j * (j + 1) / 2 + i;
	else return i * (i + 1) / 2 + j;
}
*/

int get_inverse(double* mat, double* res, double* d, int n, int n_threads, int thread_i, int* res_code) {
    double s = 0;
	double no = 0;
	int row_from = 0, row_to = 0;

	if (thread_i == 0){
		s = 0;
		no = norm(mat, n);
		row_from = 0; row_to = 0;
	}

	row_from = n * thread_i / n_threads;
	row_to = n * (thread_i +  1) / n_threads;

	for (int i = row_from; i < row_to; ++i) {
		for (int j = i; j < n; ++j) {
			res[j * (j + 1) / 2 + i] = mat[j * (j + 1) / 2 + i];
		}
	}

	synchronize(n_threads);

	for (int k = 0; k < n; ++k) {
		if (thread_i == 0) {
			
			if (thread_i == 0 & fabs(res[k * (k + 1) / 2 + k] / no) < EPS) {
				*res_code = -1;
			}

			if (res[k * (k + 1) / 2 + k] < 0) {
				d[k] = -1;
				res[k * (k + 1) / 2 + k] = -res[k * (k + 1) / 2 + k];
			} else {
				d[k] = 1;
			}
			res[k * (k + 1) / 2 + k] = sqrt(res[k * (k + 1) / 2 + k]);

			for (int j = k + 1; j < n; ++j) {
				res[j * (j + 1) / 2 + k] /= res[k * (k + 1) / 2 + k] * d[k];
			}
		}

		
		synchronize(n_threads);
		if (*res_code == -1) return -1;

		row_from = (n - k - 1) * thread_i / n_threads + k + 1;
    	row_to = (n - k - 1) * (thread_i + 1) / n_threads + k + 1;
		
		for (int i = row_from; i < row_to; ++i) {
			for (int j = i; j < n; ++j) {
				res[j * (j + 1) / 2 + i] -= res[i * (i + 1) / 2 + k] * res[j * (j + 1) / 2 + k] * d[k];
			}
		}
		synchronize(n_threads);		
	}

	synchronize(n_threads);
    
	row_from = n * thread_i / n_threads;
    row_to = n * (thread_i + 1) / n_threads;
    
    for (int i = row_from; i < row_to; ++i) {
		
		for (int j = i; j >= 0; --j) {
			s = (double)(i == j);
			for (int k = j + 1; k <= i; ++k)
				s -= mat[i * (i + 1) / 2 + k] * res[k * (k + 1) / 2 + j];
				
			mat[i * (i + 1) / 2 + j] = s / res[j * (j + 1) / 2 + j];
		}
	}


	synchronize(n_threads);

	for (int i = 0; i < n; ++i) {
		//synchronize(n_threads);
		row_from = (i + 1) * thread_i / n_threads;
    	row_to = (i + 1) * (thread_i + 1) / n_threads;

		for (int j = row_from; j <= row_to; ++j)
		{
			s = 0.0;
			for (int k = i; k < n; ++k)
				s += d[k] * mat[k * (k + 1) / 2 + j] * mat[k * (k + 1) / 2 + i];
			res[i * (i + 1) / 2 + j] = s;
		}

		synchronize(n_threads);
		
		row_from = (n - (i + 1)) * thread_i / n_threads + i + 1;
    	row_to = (n - (i + 1)) * (thread_i + 1) / n_threads + i + 1;
		for (int j = row_from; j < row_to; ++j)
		{
			s = 0.0;
			for (int k = j; k < n; ++k)
				s += d[k] * mat[k * (k + 1) / 2 + j] * mat[k * (k + 1) / 2 + i];
			res[j * (j + 1) / 2 + i] = s;
		}
	}
	synchronize(n_threads);
    return 1;
}

int cholesky_decomp(double* mat, double* res, double* d, int n, int n_threads, int thread_i) {
    double s = 0;
    double no = norm(mat, n);
    for (int i = 0; i < n; ++i)
		for (int j = i; j < n; ++j)
		{
			s = mat[j * (j + 1) / 2 + i];
			for (int k = 0; k < i; ++k)
				s -= res[i * (i + 1) / 2 + k] * res[j * (j + 1) / 2 + k] * d[k];


			if (i == j)
			{
				if (s < 0)
				{
					s = -s;
					d[i] = -1.0;
				}
				else
					d[i] = 1.0;

				if (fabs(s / no) < EPS)
					return -1;

				s = sqrt(s);
				res[i * (i + 1) / 2 + i] = s;
			}
			else
				res[j * (j + 1) / 2 + i] = s / (res[i * (i + 1) / 2 + i] * d[i]);
		}
    return 1;
}