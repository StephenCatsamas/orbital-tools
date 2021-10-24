#include <stdio.h>
#include <math.h>
#include <limits>
#include "util.h"
#include "integrator.h"

#define R_LIMIT 10 //maximum height before termination (body radii)
#define T_LIMIT 1E6 //maximum time for simulation (s)


const double e_vr[] = {-3, 3};
const double e_r[] = {50+BODY.radius+BODY.landing_altitude, 150+BODY.radius+BODY.landing_altitude};
double get_error(const std::vector<std::array<double,SYSDIM>>& path){
    double deltH;
    double r_tgt = 0.5*(e_r[0] + e_r[1]); //target height
    double r_min = std::numeric_limits<double>::max();
    for(int i = 0; i < path.size(); i++){
        const double_v3 r = vec_unpack_r(path[i]);
        
        double rr = r.mag();
        if(rr <= r_min){
            r_min = rr;
            deltH = r.mag() - r_tgt;
        }
    }       
    return deltH;
}

int bisect(double v, double* x, double* uB, double* lB){
    if(v>0){
        *uB = *x;
    }else if(v < 0){
        *lB = *x;
    }
    *x = (*uB + *lB) / 2;
    return 0;
}

double t_stop;
std::array<double,SYSDIM> zi;
bool time_stop(const std::array<double,SYSDIM>& z){
    return (t_stop - z[0] > 0);
}

int set_t_stop(const double t){
    t_stop = t;
    return 0;
}

double ccw(double* A, double* B, double* C){
    return (((C[1]-A[1]) * (B[0]-A[0])) > ((B[1]-A[1]) * (C[0]-A[0])));
}
    

bool intersect(double* l1a, double* l1b, double* l2a, double* l2b){
    return ccw(l1a,l2a,l2b) != ccw(l1b,l2a,l2b) and ccw(l1a,l1b,l2a) != ccw(l1a,l1b,l2b);
}

bool has_sol(const std::vector<std::array<double,SYSDIM>>& path){
    const int boxlist[] = {0,1,1,0};
    if(path.size() == 0){return false;}
    for(long long unsigned int i = path.size() - 1; i > 0; i--){
        const double_v3 r = vec_unpack_r(path[i]);
        const double_v3 vs = cross(BODY.w,r);//surface velocity 
        const double_v3 v = (vec_unpack_v(path[i]) - vs);//relative surface velocity
        
        if(ordered(e_r[0], r.r(r) , e_r[1]) and ordered(e_vr[0], v.r(r), e_vr[1])){
            return true;
        }else if (i != 0){
            const double_v3 r2 = vec_unpack_r(path[i-1]);
            const double_v3 v2 = (vec_unpack_v(path[i-1]) - vs);
            double rr1 = r.r(r);
            double rr2 = r2.r(r2);
            double vr1 = v.r(r);
            double vr2 = v2.r(r2);
            
            for(int j = 0; j < 4; j++){
                int irr = boxlist[mod((j),4)];
                int ivr = boxlist[mod((j-1),4)];
                int jrr = boxlist[mod((j+1),4)];
                int jvr = boxlist[mod((j),4)];
                double l1a[] = {rr1, vr1};
                double l1b[] = {rr2, vr2};
                double l2a[] = {e_r[irr], e_vr[ivr]};
                double l2b[] = {e_r[jrr], e_vr[jvr]};
                if(intersect(l1a,l1b, l2a, l2b)){
                    return true;
                }
            }
            
        }
        
        
        
    }
    return 0;
}


bool landed(const std::array<double,SYSDIM>& z){
    double t = z[0];
    double_v3 r = {z[1],z[2],z[3]};
    double_v3 v = {z[4], z[5], z[6]};
    
    //if our radius is less than the landing height
    if(r.mag() < BODY.radius+BODY.landing_altitude){
        // printf("landing termination\n");
        return 0;
    }
    //if we are between the initial height and surface and have a postive velocity
    double_v3 r0 = {zi[1],zi[2],zi[3]};
    if((r.mag() > BODY.radius+BODY.landing_altitude) 
    and (r.mag() < r0.mag()) 
    and (v.r(r) > 0)){
        // printf("re-orbit termination\n");
        return 0;
    }
    if(r.mag() > R_LIMIT*BODY.radius){
        // printf("height termination\n");
        return 0;
    }
    if(t - zi[0] > T_LIMIT){
        printf("time termination\n");
        return 0;
    }
    
    return 1;
}

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


int stepper(std::vector<std::array<double,SYSDIM>>& path, bool (*check_stop)(const std::array<double,SYSDIM>&), double* error, int (*system)(double*, double*)){

    static double hdyn = 0.1;
    
    std::array<double,SYSDIM> zc;
    std::array<double,SYSDIM> z0;

    double loc_sf[SYSDIM];
    double loc_h1[SYSDIM];
    double loc_h2[SYSDIM];
    double delta[SYSDIM];
    double ztmp[SYSDIM];

    veccpy(z0, path.back());
    veccpy(zi, path.back());//save initial state in global variable
    double error_norm = 1E-3;
    
    while(check_stop(path.back())){
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
            
            hdyn = new_step_size(hdyn, delta_norm, error_norm);
            if(delta_norm <= error_norm){got_accuracy = true;}
            
        }
        // if(t_stop+z0[0] - loc_h2[0] < 0){
            // hdyn = (t_stop+z0[0] - path.back()[0]);
            // veccpy(ztmp,path.back());
            // rK4(ztmp, loc_h2, hdyn, system);
        // }
        
        veccpy(zc, loc_h2);//converts to std::array
#ifdef _DEBUG
        if(path.size() > 5E4){
        path.resize(1);
        }
#endif
        path.push_back(zc);
    }


    
    return 0;
}

double new_step_size(double hdyn, double d_norm, double e_norm){
    if(d_norm <= e_norm){
        if(d_norm != 0){
        hdyn *= 0.9 * Q_fqrt(e_norm/d_norm);//pow very costly
        }
    }else{
        hdyn *= 0.9 * Q_qqrt(e_norm/d_norm);
    }  
    return hdyn;
}

//greater than unity only
double Q_fqrt(double n){
    if(n - 1 < 5){
        return 1 + (n-1)/5;
    }else{
        return pow(n, 0.2);
    }
    
    
}

double Q_qqrt(double n){
    float number = (float)(n);
    number = Q_rsqrt(Q_rsqrt(number));
    n = (double)(number);
    return n;
}

//quick inverse square root
float Q_rsqrt(float number){
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y  = number;
	i  = * ( long * ) &y;                       // evil floating point bit level hacking
	i  = 0x5f3759df - ( i >> 1 );                
	y  = * ( float * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration

	return y;
}