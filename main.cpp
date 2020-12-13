#ifndef MATR
#include "functions.h"
#define MATR
#endif
#include "io.h"
#include <time.h>

#include "Cholesky.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <pthread.h>

using namespace std;

long int get_full_time(void)
{
	struct timeval buf;
	gettimeofday(&buf, 0);
	return buf.tv_sec * 100 + buf.tv_usec/10000;
}

typedef struct
{
	int n;
	double* mat;
	double* d;
    double* res;
	int n_threads;
    int thread_i;
    int* code;
} ARGS;

void* solve(void* arg) {
    ARGS* args = (ARGS*) arg;
    get_inverse(args->mat, args->res, args->d, args->n, args->n_threads, args->thread_i, args->code);
	return NULL;
}

int main(int argc, char * argv[]) {
    string filename = "";
    int n, m, k, n_threads;
    if (argc == 6) {
        n = stoi(argv[1]);
        m = stoi(argv[2]);
        k = stoi(argv[3]);
        n_threads = stoi(argv[4]);
        filename = argv[4];
    } else if (argc == 5) {
        n = stoi(argv[1]);
        m = stoi(argv[2]);
        k = stoi(argv[3]);
        n_threads = stoi(argv[4]);
        
    } else {
        cout << "Usage: ./main n m k filename or ./main n m k when k != 0";
        return 0;
    }

    double* matrix;
    double* d;
    double* res;
    pthread_t* threads;
    ARGS* args;
    try {
        res = new double[n * (n + 1) / 2];
        matrix = new double[n * (n + 1) / 2];
        d = new double[n];
        threads = (pthread_t*) malloc(n_threads * sizeof(pthread_t));
        args = (ARGS*) malloc(n_threads * sizeof(ARGS));
    } catch (const bad_alloc) {
        cout << "memory error" << endl;
        return 0;
    }

    if (k != 0) {
        generate_matrix_from_formula(matrix, n, k);
    } else {
        if (read_matrix_from_file(matrix, n, filename) != 1) {
            delete[] matrix; delete[] res; delete[] d;
            free(threads); free(args);
            return 0;
        }
    }

    cout << "MATRIX: " << endl;
    print(matrix, n, m, cout);
    cout << endl;

    int res_code = 0;
    for (int i = 0; i < n_threads; i++) {
		args[i].n = n;
		args[i].mat = matrix;
		args[i].d = d;
        args[i].res = res;
		args[i].thread_i = i;
		args[i].n_threads = n_threads;
        args[i].code = &res_code;
	}

    long int t_start = get_full_time();

    
    for (int i = 0; i < n_threads; i++)
		if (pthread_create(threads + i, 0, solve, args + i))
		{
			cout << "Error creating thread" << endl;
			delete[] matrix; delete[] res; delete[] d;
			if (threads) free(threads);
			if (args) free(args);
			return 0;
		}

    
	for (int i = 0; i < n_threads; i++)
		if (pthread_join(threads[i], 0))
		{
			cout << "Error joining thread" << endl;
			delete[] matrix; delete[] res; delete[] d; 
            if (threads) free(threads);
			if (args) free(args);
			return 0;
		}
    
   //get_inverse(matrix, res, d, n, n_threads, 0);

    cout << "code: " << res_code << endl;
    if (res_code == -1) {
        cout << "det is zero" << endl;
         delete[] matrix; delete[] res;
         delete[] d; 
         free(threads); free(args);
         return 0;
    }

    //get_inverse(matrix, res, d, n, n_threads, 0);
    cout << "THREADS: " << n_threads << endl;
    cout << "TIME: " << (double) (get_full_time() - t_start) << endl;

    cout << "A: " << endl;
    print(matrix, n, m, cout);
    cout << "RES: " << endl;
    print(res, n, m, cout);
    cout << endl;

    if (k != 0) {
        generate_matrix_from_formula(matrix, n, k);
    } else {
        if (read_matrix_from_file(matrix, n, filename) != 1) {
            delete[] matrix; delete[] res; delete[] d;
            free(threads); free(args);
            return 0;
        }
    }

    cout << "NORM: " << mult_err(matrix, res, n) << endl;
    delete[] matrix; delete[] res; delete[] d; 
    free(threads); free(args);
}