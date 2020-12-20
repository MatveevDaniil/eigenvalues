void values_search(int n, double* matrix, double* values, double* x1, double* x2, double EPS = 1e-10);
void semitriangle_rotation(int n, double* matrix);
double inf_norm(int n, int m, double* matrix);
void QR_reflection(int n, double* matrix, int loc_n, double* x1, double* x2);
void QR2RQ(int n, double* matrix, int loc_n, double* x1, double* x2);