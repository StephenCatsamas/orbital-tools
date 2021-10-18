#include <math.h>
#include "util.h"
#include "integrator.h"



int rK4(double* z0, double* z, double h, int (*system)(double*, double*)){
    veccpy(z, z0);
    //k1
    double k1[SYSDIM]; 
    system(z,k1);
    
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
    double k2[SYSDIM]; 
    system(z, k2);
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
    double k3[SYSDIM];  
    system(z, k3);
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
    double k4[SYSDIM]; 
    system(z, k4);
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


int stepper(std::vector<std::array<double,SYSDIM>>& path, double t_stop, double* error, int (*system)(double*, double*)){
    
    static double hdyn = 0.1;
    
    std::array<double,SYSDIM> zc;
    std::array<double,SYSDIM> z0;

    double loc_sf[SYSDIM];
    double loc_h1[SYSDIM];
    double loc_h2[SYSDIM];
    double delta[SYSDIM];
    double ztmp[SYSDIM];

    veccpy(z0, path.back());

    double error_norm = 1E-9;
    
    while(t_stop+z0[0] - path.back()[0] > 0){
        bool got_accuracy = false;
        while(not got_accuracy){
            veccpy(ztmp,path.back());
            rK4(ztmp, loc_sf, 2*hdyn, system);
            rK4(ztmp, loc_h1, hdyn, system);
            rK4(loc_h1, loc_h2, hdyn, system);
            
            for(int i = 1; i < SYSDIM; i++){
                delta[i] = loc_sf[i] - loc_h2[i];
            }
            
            double delta_norm = 0;
            for(int i = 1; i< SYSDIM; i++){
                double v = fabs(delta[i]);
                if(v > delta_norm){
                    delta_norm = v;
                }
            }
            
            if(delta_norm <= error_norm){
                if(delta_norm != 0){
                    hdyn = 0.9 * hdyn * pow(error_norm/delta_norm, 0.2);//pow very costly
                }
                got_accuracy = true;
            }else{
                hdyn = 0.9 * hdyn * pow(error_norm/delta_norm, 0.25);//pow very costly
                got_accuracy = false;
            }  
        }
        if(t_stop+z0[0] - loc_h2[0] < 0){
            hdyn = (t_stop+z0[0] - path.back()[0]);
            veccpy(ztmp,path.back());
            rK4(ztmp, loc_h2, hdyn, system);
        }
        
        veccpy(zc, loc_h2);//converts to std::array
        path.push_back(zc);
    }    
    return 0;
}