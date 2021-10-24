#include <vector>
#include <array>
#include <stdio.h>
#include <math.h>
#include "util.h"
#include "integrator.h"
#include "main.h"

double thrust = 0.01;

int main(int argc, char **argv){
	printf("Orbits!\n");
    BODY.init2();
    

    std::vector<std::array<double,SYSDIM>> path;
    std::array<double, SYSDIM> z0 = {0,2E6,0,0,0,1.4E3,-1.4E3};
    path.push_back(z0);
    
    // set_t_stop(100);
    // stepper(path, time_stop, NULL, gravsys);
    bool sol = solve_BVP(path);
     
    printf("Has solution: %s\n", sol?"true":"false");
    printf("At thrust: %f\n", thrust);

    write_meta();
    write_path(path);
    printf("Orbited!\n");

	return 0;
}

int solve_BVP(std::vector<std::array<double,SYSDIM>>& path){
    double* error = NULL;
    int max_iter = 50;
    
    thrust = TWR_to_thrust(2, BODY.radius+BODY.landing_altitude);
    double lThrust = TWR_to_thrust(1, BODY.radius+BODY.landing_altitude);
    double hThrust = TWR_to_thrust(1E3, BODY.radius+BODY.landing_altitude);
    
    bool sol = false;
    for(int i = 0; i< max_iter; i++){   
        PRINTFLT(thrust);
        
        path.resize(1);
        stepper(path, landed, error, gravsys);
        if(has_sol(path)){
            sol = true;
            break;
        }else{
            double height_error = get_error(path);
            bisect(height_error, &thrust, &hThrust, &lThrust);
            
        }   
    }
    return sol;
}

int gravsys(double* r0, double* dr0){
    double m = 1;
    double M = BODY.mass;
    //object position at zero;
    
    //unpack
    const double_v3 r = vec_unpack_r(r0);
    const double_v3 vs = cross(BODY.w, r);//surface velocity at landing height
    
    const double t = vec_unpack_t(r0);
    const double_v3 v = (vec_unpack_v(r0) - vs);//relative surface velocity
    
    const double R = r.mag();
    const double V = v.mag();

    //dunamds
 
    double_v3 Ft;
    if(V > 1){
        Ft = -thrust*v/V;
    }else{
        Ft = thrust*r/R;
    }

    
    
    double gf = (-G*M*m/(R*R*R));
    double_v3 Fg = gf*r;
    
    double_v3 F = Ft + Fg;
    
    double_v3 dr = v + vs;//proper motion
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
