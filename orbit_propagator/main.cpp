#include <vector>
#include <array>
#include <stdio.h>
#include <math.h>
#include "util.h"
#include "integrator.h"
#include "interface.h"
#include "main.h"

double burn_wait = 1;
double thrust = 1;

int main(int argc, char **argv){
	printf("Orbits!\n");
    
    load_inputs(rocket, BODY);
   
    BODY.init2();
    populate_tolerances(); 

    std::vector<std::array<double,SYSDIM>> path;
    std::vector<std::array<double,AUXDIM>> aux_path;
    std::array<double, SYSDIM> z0 = orbit_to_position(30E3, -5E3, 29E3, 0, false);
    z0[7] = rocket.mass_wet;
    std::array<double, AUXDIM> axz0 = {0,0,0};
    path.push_back(z0);
    aux_path.push_back(axz0);
    // set_t_stop(10000);
    // stepper(path, aux_path, time_stop, NULL, gravsys);
    bool sol = solve_BVP_time(path, aux_path);
     
    printf("Has solution: %s\n", sol?"true":"false");
    printf("At Burn Wait Time: %f\n", burn_wait);

    //calulate path statistics
    get_statistics(path,aux_path);
    
    write_meta();
    write_path(path,aux_path);
    printf("Orbited!\n");

	return 0;
}

int solve_BVP_time(std::vector<std::array<double,SYSDIM>>& path, std::vector<std::array<double,AUXDIM>>& aux_path){
    double* error = NULL;
    int max_iter = 50;
    
    double lwait = 0;
    double hwait = 1E4;
    burn_wait = hwait/2;
    
    bool sol = false;
    for(int i = 0; i< max_iter; i++){   
        PRINTFLT(burn_wait);
        
        path.resize(1);
        aux_path.resize(1);
        stepper(path, aux_path, landed, error, gravsys);
        if(has_sol(path)){
            sol = true;
            break;
        }else{
            double height_error = get_error(path);
            bisect(-height_error, &burn_wait, &hwait, &lwait);
            
        }   
    }
    return sol;
}

int solve_BVP_thrust(std::vector<std::array<double,SYSDIM>>& path, std::vector<std::array<double,AUXDIM>>& aux_path){
    double* error = NULL;
    int max_iter = 50;
    
    thrust = TWR_to_thrust(2, BODY.radius+BODY.landing_altitude);
    double lThrust = 0;
    double hThrust = rocket.thrust;
    
    bool sol = false;
    for(int i = 0; i< max_iter; i++){   
        PRINTFLT(thrust);
        
        path.resize(1);
        aux_path.resize(1);
        stepper(path, aux_path, landed, error, gravsys);
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


int gravsys(std::array<double,SYSDIM>& r0, std::array<double,SYSDIM>& dr0, std::array<double,AUXDIM>& aux_r0){
    double M = BODY.mass;
    double m_dry = rocket.mass_dry;
    const double g = 9.81;//m/s^2
    double exaust_velocity = rocket.ISP*g;
    //object position at zero;
    //unpack
    const double_v3 r = vec_unpack_r(r0);
    const double_v3 vs = cross(BODY.w, r);//surface velocity at landing height
    
    const double t = vec_unpack_t(r0);
    const double_v3 v = (vec_unpack_v(r0) - vs);//relative surface velocity
    double m = vec_unpack_mass(r0);    
    
    const double R = r.mag();
    const double V = v.mag();



    //dunamds
 
    double_v3 Ft;
    if((m < m_dry) or (t < burn_wait)){
        Ft = {0,0,0};
    }else if(V > THRUST_CUT_VEL){
        Ft = -rocket.thrust*v/V;
    }else{
        // Ft = rocket.thrust*r/R;
        Ft = {0,0,0};
    }

    
    
    double gf = (-G*M*m/(R*R*R));
    double_v3 Fg = gf*r;
    
    double_v3 F = Ft + Fg;
    
    double_v3 dr = v + vs;//proper motion
    double_v3 dv = F/m;
    
   
    double dt = 1;//this does nothing as we only need the deriatives of the non-time variables
    double dm = -Ft.mag()/exaust_velocity;
    
    //pack
    vecpack(dr0, dt, dr, dv, dm);
    vecpack(aux_r0, Ft);//log thrust
    return 0;
}

int expsys(double* z0, double* dz){
    dz[1] = -z0[2];
    dz[2] = z0[1];
    return 0;    
}



