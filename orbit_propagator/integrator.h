#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <array>

bool time_stop(const std::array<double,SYSDIM>& z);
int set_t_stop(const double t);
bool landed(const std::array<double,SYSDIM>& z);
int rK4(double* z0, double* z, double h, int (*system)(double*, double*));
int stepper(std::vector<std::array<double,SYSDIM>>& path, bool (*check_stop)(const std::array<double,SYSDIM>&), double* error, int (*system)(double*, double*));
double new_step_size(double hdyn, double d_norm, double e_norm);
double Q_fqrt(double);
double Q_qqrt(double);
float Q_rsqrt(float);

#endif