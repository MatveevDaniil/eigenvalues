#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <iostream>
#include "read_print.h"
#include "valuesalgorithm.h"

int main(int argc, char** argv) {
	int n, m, k;
	std::string file_name;
	int iterations;
	double *matrix, *x1, *x2, *values, EPS;
	_CRT_DOUBLE _EPS;

	if (argc < 4 || argc > 5) {
		std::cout << "wrong argument's number = " << argc << std::endl;
		std::cout << "usage: 'prog n m EPS 0 filename' or 'prog n m EPS k'" << std::endl;
		return -1;
	}
	else {
		n = atoi(argv[1]);
		m = atoi(argv[2]);
		int retval = _atodbl(&_EPS, argv[3]);
		if (retval != 0) {
			std::cout << "usage: 'prog n m EPS 0 filename' or 'prog n m EPS k'" << std::endl;
			std::cout << "where EPS must be double" << std::endl;
		}
		EPS = _EPS.x;
		k = atoi(argv[4]);
		if (k == 0) {
			if (argc == 5) {
				file_name = argv[5];
			}
			else {
				std::cout << "if k == 0, u need to give file with matrix: 'prog n m EPS 0 filename'" << std::endl;
				return -2;
			}
		}
		if (n <= 0 || m <= 0 || k < 0 || k > 6 || EPS < 1e-16) {
			std::cout << "usage: 'prog n m 0 filename' or 'prog n m k', where n > 0, m > 0, 4 >= k > 0, EPS > 0" << std::endl;
			return -2;
		}
	}
	
	matrix = new double[n * n];
	x1 = new double[n];
	x2 = new double[n];
	values = new double[n];
	if(matrix == NULL || x1 == NULL || x2 == NULL || values == NULL) {
		std::cout << "not enough memory" << std::endl;
		if(matrix) delete[] matrix;
		if(x1) delete[] x1;
		if(x2) delete[] x2;
		if (values) delete[] values;
		return -5;
	}

	if (k == 0) {
		int read_err = read_matrix(file_name.c_str(), matrix, n);
		switch (read_err) {
		case (-1):
			std::cout << "can not open file " << file_name.c_str() << std::endl;
			delete[] matrix, x1, x2, values;
			return -4;
		case (-2):
			std::cout << "can not read matrix from file " << file_name.c_str() << std::endl;
			delete[] matrix, x1, x2, values;
			return -5;
		case (-3):
			std::cout << "error with input array allocation" << std::endl;
			delete[] matrix, x1, x2, values;
			return -6;
		default:
			break;
		}
	}
	else
		fill_matrix(k, n, matrix);

	std::cout << "inputed matrix:" << std::endl;
	print_matrix(m, m, n, matrix);
	double trace(0), len(0);
	for (int i = 0; i < n; i++) {
		trace += matrix[i * n + i];
		for (int j = 0; j < n; j++)
			len += matrix[i * n + j] * matrix[j * n + i];
	}
	std::clock_t total_time = 0, start = 0;	
	values_search(n, matrix, values, x1, x2, EPS);
	total_time += std::clock() - start;
	std::cout << "total time: " << (double)total_time / CLOCKS_PER_SEC << " sec.\n" << std::endl;
	std::cout << "values: " << std::endl;
	for (int i = 0; i < n; i++) {
		trace -= values[i];
		len -= values[i] * values[i];
	}
	std::cout << std::scientific;
	std::cout << "trace residual: " << trace << std::endl;
	std::cout << "L2_len^2 residual: " << len << std::endl;
	print_matrix(m, 1, n, values);
	delete[] matrix, x1, x2, values;
	return 0;
}
