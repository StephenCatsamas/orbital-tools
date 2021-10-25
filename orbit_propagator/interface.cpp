#include "interface.h"

//we use the apoasis here as it is more likely to be visible when suborbital
std::array<double, SYSDIM> orbit_to_position(double apoapsis, double sea_altitude, double inclination, bool accending){
    
    
    
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