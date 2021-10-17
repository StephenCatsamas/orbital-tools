#include <stdio.h>
#include "util.h"
#include "integrator.h"

int main(int argc, char **argv){
	printf("Quaternions!\n");
    
    double_v3 q = {3,4,5};
    double_v3 p = {1,2,3};

    q = q+p;
    
    PRINTDV3(q);
    
    const int ts = 100;
    double h = 0.1;
    double zs[ts][2];
    double z0[] = {0,1};
    veccpy(zs[0],z0);
    for(int i = 1; i < ts; i++){
        rK4(zs[i-1], zs[i], h, expsys);
    }
    
    FILE* fp = fopen("vis.dat", "w");
    for(int i = 0; i < ts; i++){
    for(int j = 0; j < SYSDIM; j++){
        fprintf(fp, "%f,", zs[i][j]);
    }
        fprintf(fp, "\n");
    }

	return 0;
}






