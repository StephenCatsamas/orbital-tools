#ifndef MAIN_H
#define MAIN_H



int expsys(double* z0, double* dz);
int gravsys(double* z0, double* dz);
int main(int argc, char **argv);
int write_path(std::vector<std::array<double,SYSDIM>>& path);
int write_meta(void);


#endif