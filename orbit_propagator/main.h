#ifndef MAIN_H
#define MAIN_H



int expsys(double* z0, double* dz);
int gravsys(double* z0, double* dz);
int solve_BVP(std::vector<std::array<double,SYSDIM>>& path);
int main(int argc, char **argv);



#endif