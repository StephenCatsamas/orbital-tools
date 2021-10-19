#ifndef UTIL_H
#define UTIL_H

#include <array>

#define PRINTFLT(var) printf("%s: %f at %d\n", #var, var, __LINE__)
#define PRINTINT(var) printf("%s: %d at %d\n ", #var, var, __LINE__)
#define PRINTDV3(var) printf("%s: (%f,%f,%f) at %d\n ", #var, var.x, var.y, var.z, __LINE__)

#define SYSDIM 7

#define G 6.67E-11



struct quat {
    
    double r;
    double i;
    double j;
    double k;
    
    quat& operator+ (const quat&);
    quat& operator* (const quat&);
    
};

struct double_v3{
    double x;
    double y;
    double z;
    
    double_v3& operator/ (const double);
    double_v3& operator= (const double_v3&);
    bool operator== (const double_v3&);
    
    double mag(void) const;
    
};

inline double_v3 operator+ (double_v3 u, const double_v3& v){ 
    u.x += v.x;    
    u.y += v.y;    
    u.z += v.z;
    return u;
}

inline double_v3 operator- (double_v3 u, const double_v3& v){ 
    u.x -= v.x;    
    u.y -= v.y;    
    u.z -= v.z;
    return u;
}

inline double_v3 operator* (double_v3 v, const double a){ 
    v.x *= a;    
    v.y *= a;    
    v.z *= a;    
    return v;
}

inline double_v3 operator* (const double a, double_v3 v){ 
    return v*a;
}
struct body{
    double mass;//kg
    double radius; //m
    double landing_altitude; //m
    double rotational_period;  //s
};

extern body earth;
extern body moon;


int veccpy(double* a, double* b);
int veccpy(std::array<double,SYSDIM>& a, double* b);
int veccpy(double* a, std::array<double,SYSDIM>& b);
int veccpy(std::array<double,SYSDIM>& a, std::array<double,SYSDIM>& b);
int vecpack(double* a, const double t, const double_v3& b, const double_v3& c);

#endif