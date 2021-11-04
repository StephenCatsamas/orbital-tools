#ifndef INTERFACE_H
#define INTERFACE_H

#include "util.h"

std::array<double, SYSDIM> orbit_to_position(double apoapsis, double periapsis, double sea_altitude, double inclination, bool accending);

int load_inputs(craft& rocket, body& planet);
int write_path(std::vector<std::array<double,SYSDIM>>& path, std::vector<std::array<double,AUXDIM>>& aux_path);
int write_meta(void);

double get_s_maj_a(double apoapsis, double periapsis);
double get_sin_to_slr(double Lsq, double mu, double esqm1, double r);
double get_specific_ke(double E, double mu, double r);
double get_specific_energy(double apoapsis, double periapsis, double mu);
double get_sqr_specific_anglular_momentum(double apoapsis, double periapsis, double mu);
double get_eccentrity_sqrm1(double E, double Lsq, double mu);


#endif