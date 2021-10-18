#include <vector>
#include <array>


int rK4(double* z0, double* z, double h, int (*system)(double*, double*));
int stepper(std::vector<std::array<double,SYSDIM>>& path, double t_stop, double* error, int (*system)(double*, double*));