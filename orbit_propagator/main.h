#ifndef MAIN_H
#define MAIN_H


int expsys(std::array<double,SYSDIM>& z0, std::array<double,SYSDIM>& dz, std::array<double,AUXDIM>& aux_r0);
int gravsys(std::array<double,SYSDIM>& z0, std::array<double,SYSDIM>& dz, std::array<double,AUXDIM>& aux_r0);
int solve_BVP(std::vector<std::array<double,SYSDIM>>& path, std::vector<std::array<double,AUXDIM>>& aux_path);
int main(int argc, char **argv);



#endif