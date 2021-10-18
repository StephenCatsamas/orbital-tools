#include <vector>
#include <array>
#include <stdio.h>
#include "util.h"
#include "integrator.h"

int expsys(double* z0, double* dz){
    dz[1] = -z0[2];
    dz[2] = z0[1];
    return 0;    
}

int main(int argc, char **argv){
	printf("Quaternions!\n");
    
    double_v3 q = {3,4,5};
    double_v3 p = {1,2,3};

    q = q+p;
    
    PRINTDV3(q);
    
    const int ts = 1000000;
    double* error;

    std::vector<std::array<double,SYSDIM>> path;
    std::array<double, SYSDIM> z0 = {0,1};
    path.push_back(z0);
    
    stepper(path, ts, error, expsys);

    FILE* fp = fopen("vis.dat", "w");
    // for(int i = 0; i < path.size(); i++){
    for(int j = 0; j < SYSDIM; j++){
        fprintf(fp, "%f,", path[0][j]);
    }
        fprintf(fp, "\n");
    // }

	return 0;
}






