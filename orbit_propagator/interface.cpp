#include <math.h>
#include <stdio.h>
#include <vector>
#include "interface.h"

//we use the apoasis here as it is more likely to be visible when suborbital
std::array<double, SYSDIM> orbit_to_position(double apoapsis, double periapsis, double sea_altitude, double inclination, bool accending){
    apoapsis += BODY.radius;
    periapsis += BODY.radius;
    
    double mu = G*BODY.mass;//graviational parameter
    double E = -mu/(apoapsis + periapsis);//specific orbital energy
    double Lsq = 2*mu*(apoapsis*periapsis)/(apoapsis + periapsis);//specific angular momenum squared
    
    double r_height = sea_altitude+BODY.radius;
    double vsq = 2*(E + mu/(r_height));
    double v_mag = sqrt(vsq);
    
    double sinsq = Lsq/(r_height*r_height*vsq);//sqaured sin of angle between radius and velocity
    double cossq = 1-sinsq;
    double sin = sqrt(sinsq);
    double cos = sqrt(cossq);
    //TODO implement inclination
    
    if(accending){cos = cos; PRINTFLT(cos);}
    else{cos = -1*cos;PRINTFLT(cos);}
    
    
    double t = 0;
    double_v3 r = {r_height,0,0};
    double_v3 v = {v_mag*cos,v_mag*sin,0};
    double mass = 0;
    std::array<double, SYSDIM> loc = {t, r.x,r.y,r.z,v.x,v.y,v.z,mass};
    return loc;
}


int write_meta(void){
    FILE* fp = fopen("tmp/met.dat", "w");

    fprintf(fp, "BODY, mass, radius, landing_altitude, rotational_period, ω_x, ω_y, ω_z\n");
    fprintf(fp, "%s, %e, %e, %e, %e, %e, %e, %e\n", 
                (NAME(BODY)), 
                BODY.mass, 
                BODY.radius, 
                BODY.landing_altitude, 
                BODY.rotational_period,
                BODY.w.x,
                BODY.w.y,
                BODY.w.z);

    fclose(fp);

    return 0;
}

int write_path(std::vector<std::array<double,SYSDIM>>& path){
#ifndef _DEBUG
    FILE* fp = fopen("tmp/vis.dat", "w");
    for(long long unsigned int i = 0; i < path.size(); i++){
    for(int j = 0; j < SYSDIM; j++){
        fprintf(fp, "%e,", path[i][j]);
    }
        fprintf(fp, "\n");
    }
#else
    FILE* fp = fopen("NUL", "w");
    fprintf(fp, "%e,", path[0][0]);
#endif
    fclose(fp);
    return 0;
}