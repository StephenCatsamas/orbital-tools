

int expsys(double* z0, double* dz);
int rK4(double* z0, double* z, double h, int (*system)(double*, double*));