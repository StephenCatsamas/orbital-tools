#ifndef INTERFACE_H
#define INTERFACE_H

#include "util.h"

std::array<double, SYSDIM> orbit_to_position(double apoapsis, double periapsis, double sea_altitude, double inclination, bool accending);

int write_path(std::vector<std::array<double,SYSDIM>>& path);
int write_meta(void);

#endif