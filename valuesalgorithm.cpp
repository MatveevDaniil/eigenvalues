#include <cmath>
#include <cstdio>
#include "read_print.h"
#include "valuesalgorithm.h"

void semitriangle_rotation(int n, double* matrix) {
	double x1, x2;
	for (int i = 1; i < n; i++)
		for (int j = i + 1; j < n; j++) {
			x1 = matrix[i * n + i - 1], x2 = matrix[j * n + i - 1];
			if (fabs(x2) < 1.11e-16)
				continue;
			double r = sqrt(x1 * x1 + x2 * x2);
			if (r < 1.11e-16)
				continue;
			double cos_phi = x1 / r;
			double sin_phi = -x2 / r;
			matrix[i * n + i - 1] = r;
			matrix[j * n + i - 1] = 0;
			for (int k = i; k < n; k++) {
				x1 = matrix[i * n + k], x2 = matrix[j * n + k];
				matrix[i * n + k] = x1 * cos_phi - x2 * sin_phi;
				matrix[j * n + k] = x1 * sin_phi + x2 * cos_phi;
			}
			for (int k = 0; k < n; k++) {
				x1 = matrix[k * n + i], x2 = matrix[k * n + j];
				matrix[k * n + i] = x1 * cos_phi - x2 * sin_phi;
				matrix[k * n + j] = x1 * sin_phi + x2 * cos_phi;
			}
		}
}

void QR2RQ(int n, double* matrix, int loc_n, double* x1, double* x2) {
	for (int i = 0; i < loc_n - 1; i++)
		for (int j = 0; j < i + 2; j++) {
			double foo = (matrix[j * n + i] * x1[i] + matrix[j * n + i + 1] * x2[i]) * 2;
			matrix[j * n + i] -= foo * x1[i];
			matrix[j * n + i + 1] -= foo * x2[i];
		}
}

double inf_norm(int n, int m, double* matrix) {
	double norm = 0;
	for (int i = 0; i < n; i++) {
		double local_norm = 0;
		for (int j = 0; j < m; j++)
			local_norm += fabs(matrix[i * n + j]);
		norm = std::max(norm, local_norm);
	}
	return norm;
}

void QR_reflection(int n, double* matrix, int loc_n, double* x1, double* x2) {
	for (int i = 0; i < loc_n - 1; i++) {
		double foo = matrix[(i + 1) * n + i] * matrix[(i + 1) * n + i], tmp;
		if (foo < 1e-16) {
			tmp = fabs(matrix[i * n + i]);
			if (matrix[i * n + i] > 0) x1[i] = 1;
			else x1[i] = -1;
			x2[i] = 0;
		}
		else {
			tmp = sqrt(matrix[i * n + i] * matrix[i * n + i] + foo);
			matrix[i * n + i] -= tmp;
			foo = sqrt(matrix[i * n + i] * matrix[i * n + i] + foo);
			x1[i] = matrix[i * n + i] / foo;
			x2[i] = matrix[(i + 1) * n + i] / foo;
		}
		for (int j = i + 1; j < loc_n; j++) {
			foo = (x1[i] * matrix[i * n + j] + x2[i] * matrix[(i + 1) * n + j]) * 2;
			matrix[i * n + j] -= foo * x1[i];
			matrix[(i + 1) * n + j] -= foo * x2[i];
		}
		matrix[i * n + i] = tmp;
		matrix[(i + 1) * n + i] = 0;
	}
}

void values_search(int n, double* matrix, double* values, double* x1, double* x2, double EPS) {
	double s_k, stop_param = inf_norm(n, n, matrix) * EPS;

	semitriangle_rotation(n, matrix);
	//std::cout << std::endl;
	//print_matrix(9,9, n, matrix);
	int i;
	double check_param;
	for (int loc_n = n; loc_n > 2; loc_n--) {
		check_param = fabs(matrix[(loc_n - 1) * n + loc_n - 2]);
		while (check_param > stop_param) {
			std::cout << std::endl << check_param << std::endl;
			std::cout << std::endl;
			print_matrix(6, 6, n, matrix);
			s_k = matrix[(loc_n - 1) * n + loc_n - 1];
			for (i = 0; i < loc_n; i++) matrix[i * n + i] -= s_k;
			QR_reflection(n, matrix, loc_n, x1, x2);
			std::cout << "qr" << std::endl;
			print_matrix(6, 6, n, matrix);
			std::cout << std::endl;
			QR2RQ(n, matrix, loc_n, x1, x2);
			std::cout << "RQ" << std::endl;
			print_matrix(6, 6, n, matrix);
			std::cout << std::endl;
			for (i = 0; i < loc_n; i++) matrix[i * n + i] += s_k;
			check_param = fabs(matrix[(loc_n - 1) * n + loc_n - 2]);
		}
	}
	double trace = matrix[0] + matrix[n + 1];
	double det = matrix[0] * matrix[n + 1] - matrix[1] * matrix[n];
	double D = trace * trace - 4 * det;
	matrix[0] = (trace + sqrt(D)) / 2;
	matrix[n + 1] = (trace - sqrt(D)) / 2;

	for (int i = 0; i < n; i++)
		values[i] = matrix[i * n + i];
}
