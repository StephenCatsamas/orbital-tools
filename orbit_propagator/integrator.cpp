#include "util.h"
#include "integrator.h"

int expsys(double* z0, double* dz){

    dz[1] = z0[1];
    return 0;    
}

int rK4(double* z0, double* z, double h, int (*system)(double*, double*)){
    veccpy(z, z0);
    //k1
    double k1[SYSDIM]; system(z,k1);
    
    for(int i = 1; i < SYSDIM; i++){
        k1[i] *= h;
    }
    //now we have k1
    
    //get correct evaluation vector for k2
    z[0] = z0[0] + 0.5*h;
    for(int i = 1; i < SYSDIM; i++){
        z[i] = z0[i] + 0.5*k1[i];
    }
   
    //k2
    double k2[SYSDIM]; system(z, k2);
    for(int i = 1; i < SYSDIM; i++){
        k2[i] *= h;
    }
    //now we have k2
    
    //get correct evaluation vector for k3
    z[0] = z0[0] + 0.5*h;
    for(int i = 1; i < SYSDIM; i++){
        z[i] = z0[i] + 0.5*k2[i];
    }
    //k3
    double k3[SYSDIM];  system(z, k3);
    for(int i = 1; i < SYSDIM; i++){
        k3[i] *= h;
    }
    //now we have k3
    
    //get correct evaluation vector for k4
    z[0] = z0[0] + h;
    for(int i = 1; i < SYSDIM; i++){
        z[i] = z0[i] + k3[i];
    }
    
    //k4
    double k4[SYSDIM]; system(z, k4);
    for(int i = 1; i < SYSDIM; i++){
        k4[i] *= h;
    }
    //now we have k4
    
    z[0] = z0[0] + h;
    for(int i = 1; i < SYSDIM; i++){
        z[i] = z0[i] + (1/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    }
    
    return 0;

}


