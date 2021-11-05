#include <math.h>
#include <stdio.h>
#include <vector>
#include <cstring>
#include "interface.h"

//we use the apoasis here as it is more likely to be visible when suborbital
std::array<double, SYSDIM> orbit_to_position(double apoapsis, double periapsis, double sea_altitude, double inclination, bool ascending){
    //make apoapsis periapsis well ordered;
    if(not ordered(periapsis,apoapsis)){
        double t = apoapsis;
        apoapsis = periapsis;
        periapsis = t;
        printf("Warning: Apoapsis,periapsis backward at %s:%d\n", __FILE__, __LINE__);
    }
    if(sea_altitude < periapsis){
        printf("ERROR: altitude less than periapsis at %s:%d\n\n", __FILE__, __LINE__);
        abort();
    }
    
    apoapsis += BODY.radius;
    periapsis += BODY.radius;
    
    double mu = G*BODY.mass;//graviational parameter
    double E = get_specific_energy(apoapsis,periapsis,mu);
    double Lsq = get_sqr_specific_anglular_momentum(apoapsis, periapsis, mu);
    double esqm1 = get_eccentrity_sqrm1(E, Lsq, mu);//square eccentricity minus 1
    double e = sqrt(esqm1 + 1);
    double s_maj_a = get_s_maj_a(apoapsis, periapsis);
    
    double r_height = sea_altitude+BODY.radius;
    double vsq = 2*get_specific_ke(E,mu,r_height);
    double v_mag = sqrt(vsq);
    
    double sinsq = Lsq/(r_height*r_height*vsq);//squared sin of angle between radius and velocity
    double cossq = 1-sinsq;
    double sin = sqrt(sinsq);
    double cos = sqrt(cossq);
    
    if(ascending){cos = cos;}
    else{cos = -1*cos;}
    
    double t = 0;
    double_v3 r = {r_height,0,0};
    double_v3 v = {v_mag*cos,v_mag*sin,0};

    //angle to SLR
    double sin_slr = get_sin_to_slr(Lsq, mu, esqm1, r_height);
    double slr_angle = asin(sin_slr);
    
    if(ascending){slr_angle = slr_angle;}
    else{slr_angle = -M_PI - slr_angle;}
    

    
    
    
    double_v3 slr_axis = {0,0,1};
    double_v3 inc_axis = {1,0,0};
    
    rotate(r,slr_axis,slr_angle);
    rotate(v,slr_axis,slr_angle);

    
    rotate(r,inc_axis,inclination);
    rotate(v,inc_axis,inclination);

    
    //rotate about body axis of rotation
    
    double mass = 0;
    std::array<double, SYSDIM> loc = {t, r.x,r.y,r.z,v.x,v.y,v.z,mass};
    return loc;
}

double get_s_maj_a(double apoapsis, double periapsis){
    return (apoapsis + periapsis)*0.5;
}

double get_sin_to_slr(double Lsq, double mu, double esqm1, double r){
    double e = sqrt(esqm1 + 1);
    return - (1/e)*((Lsq)/(r*mu) - 1);
}

double get_specific_ke(double E, double mu, double r){
    return (E + mu/(r));
}

double get_specific_energy(double apoapsis, double periapsis, double mu){
    return -mu/(apoapsis + periapsis);
}

double get_sqr_specific_anglular_momentum(double apoapsis, double periapsis, double mu){
    return 2*mu*(apoapsis*periapsis)/(apoapsis + periapsis);
}

double get_eccentrity_sqrm1(double E, double Lsq, double mu){
    return 2*E*Lsq/(mu*mu);
}

int write_meta(void){
    FILE* fp = fopen("tmp/meta.dat", "w");

    fprintf(fp, "BODY, mass, radius, landing_altitude, rotational_period, ω_x, ω_y, ω_z\n");
    fprintf(fp, "%s, %e, %e, %e, %e, %e, %e, %e\n", 
                NAME(BODY), 
                BODY.mass, 
                BODY.radius, 
                BODY.landing_altitude, 
                BODY.rotational_period,
                BODY.w.x,
                BODY.w.y,
                BODY.w.z);

    fprintf(fp, "CRAFT, mass_wet, mass_dry, thrust, ISP\n");
    fprintf(fp, "%s, %e, %e, %e, %e\n", 
                NAME(rocket), 
                rocket.mass_wet, 
                rocket.mass_dry, 
                rocket.thrust, 
                rocket.ISP);

    fclose(fp);  
    return 0;
}

int load_inputs(craft& rocket, body& planet, orbit_param& orbit){
    FILE* fp = fopen("input.cfg", "r");
    
#define BUF_S 2048
    char buff[BUF_S];
    char field_buff[BUF_S];
    
    float mass;
    float radius;
    float landing_altitude;
    float rotational_period;
    
    //line 1
    fgets(buff, BUF_S, fp);
    
    //line 2
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &mass);
    
    //line 3
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &radius);
    
    //line 4
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &landing_altitude);
    
    //line 5
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &rotational_period);
    
    planet.mass = (double)mass;
    planet.radius = (double)radius;
    planet.landing_altitude = (double)landing_altitude;
    planet.rotational_period = (double)rotational_period;
    
    //line 6
    fgets(buff, BUF_S, fp);
    
    float mass_wet;
    float mass_dry;
    float thrust;
    float ISP;
    
    //line 7
    fgets(buff, BUF_S, fp);
    
    //line 8
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &mass_wet);
    
    //line 9
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &mass_dry);
    
    //line 10
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &thrust);
    
    //line 11
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &ISP);
   
    rocket.mass_wet = (double)mass_wet;
    rocket.mass_dry = (double)mass_dry;
    rocket.thrust = (double)thrust;
    rocket.ISP = (double)ISP;
    
    //line 12
    fgets(buff, BUF_S, fp);
    
    float apoapsis;
    float periapsis;
    float sea_altitude;
    float inclination;
    char ascending[256];
    
    //line 13
    fgets(buff, BUF_S, fp);
    
    //line 14
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &apoapsis);
    
    //line 15
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &periapsis);
    
    //line 16
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &sea_altitude);
    
    //line 17
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %e\n", field_buff, &inclination);
    
    //line 18
    fgets(buff, BUF_S, fp);
    sscanf(buff, "%s %s\n", field_buff, ascending);
   
    orbit.apoapsis = (double)apoapsis;
    orbit.periapsis = (double)periapsis;
    orbit.sea_altitude = (double)sea_altitude;
    orbit.inclination = (double)inclination;
    orbit.ascending = strcmp(ascending, "true") ? true : false;
    
    print_inputs(rocket, planet, orbit);

    
    fclose(fp);
    return 0;
}

int print_inputs(craft& rocket, body& planet, orbit_param& orbit){
    
    printf("==========Planet=============\n");
    PRINTFLT(planet.mass);
    PRINTFLT(planet.radius);
    PRINTFLT(planet.landing_altitude);
    PRINTFLT(planet.rotational_period);    
   
    printf("==========Rocket=============\n");
    PRINTFLT(rocket.mass_wet); 
    PRINTFLT(rocket.mass_dry);
    PRINTFLT(rocket.thrust); 
    PRINTFLT(rocket.ISP);
    
    printf("==========Orbit=============\n");
    PRINTFLT(orbit.apoapsis); 
    PRINTFLT(orbit.periapsis);
    PRINTFLT(orbit.sea_altitude); 
    PRINTFLT(orbit.inclination);
    PRINTBOOL(orbit.ascending);
    printf("\n");
    
    return 0;
}

int write_path(std::vector<std::array<double,SYSDIM>>& path, std::vector<std::array<double,AUXDIM>>& aux_path){
#ifndef _DEBUG
    FILE* fp = fopen("tmp/path.dat", "w");
    if(path.size() != aux_path.size()){
        printf("ERROR: path and aux of different lengths %s:%d\n\n", __FILE__, __LINE__);
        abort();
    }
    
    for(size_t i = 0; i < path.size(); i++){
    for(int j = 0; j < SYSDIM; j++){
        fprintf(fp, "%e,", path[i][j]);
    }
    for(int j = 0; j < AUXDIM; j++){
        fprintf(fp, "%e,", aux_path[i][j]);
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