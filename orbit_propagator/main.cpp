#include <vector>
#include <array>
#include <stdio.h>
#include <math.h>
#include "util.h"
#include "integrator.h"
#include "main.h"


int main(int argc, char **argv){
	printf("Orbits!\n");
    
    double* error = NULL;

    std::vector<std::array<double,SYSDIM>> path;
    std::array<double, SYSDIM> z0 = {0,2E6,0,0,0,1.8E3,0};
    path.push_back(z0);
        
    stepper(path, landed, error, gravsys);

    write_meta();
    write_path(path);
    printf("Orbited!\n");

	return 0;
}

int write_meta(void){
    FILE* fp = fopen("tmp/met.dat", "w");

    fprintf(fp, "BODY, mass, radius, landing_altitude, rotational_period\n");
    fprintf(fp, "%s, %f, %f, %f, %f\n", 
                (NAME(BODY)), 
                BODY.mass, 
                BODY.radius, 
                BODY.landing_altitude, 
                BODY.rotational_period);

    fclose(fp);

    return 0;
}

int write_path(std::vector<std::array<double,SYSDIM>>& path){
#ifndef _DEBUG
    FILE* fp = fopen("tmp/vis.dat", "w");
    for(long long unsigned int i = 0; i < path.size(); i++){
    for(int j = 0; j < SYSDIM; j++){
        fprintf(fp, "%f,", path[i][j]);
    }
        fprintf(fp, "\n");
    }
#else
    FILE* fp = fopen("NUL", "w");
    fprintf(fp, "%f,", path[0][0]);
#endif
    fclose(fp);
    return 0;
}

int gravsys(double* r0, double* dr0){
    double m = 1;
    double M = BODY.mass;
    //object position at zero;
    
    //unpack
    const double t = r0[0];
    const double_v3 r = {r0[1], r0[2], r0[3]};
    const double_v3 v = {r0[4], r0[5], r0[6]};
    
    const double R = r.mag();
    const double V = v.mag();
    
    
    //dunamds
    double thrust = 1.815;
    
    double_v3 Ft;
    if(V > 1){
        Ft = -thrust*v/V;
    }else{
        Ft = thrust*r/R;
    }
    
    double gf = (-G*M*m/(R*R*R));
    double_v3 Fg = gf*r;
    
    double_v3 F = Ft + Fg;
    
    double_v3 dr = v;
    double_v3 dv = F/m;
    
    //pack
    vecpack(dr0, 0, dr, dv);
    return 0;
}

int expsys(double* z0, double* dz){
    dz[1] = -z0[2];
    dz[2] = z0[1];
    return 0;    
}


