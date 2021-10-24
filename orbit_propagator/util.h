#ifndef UTIL_H
#define UTIL_H

#include <array>

#define PRINTFLT(var) printf("%s: %f at %d\n", #var, var, __LINE__)
#define PRINTINT(var) printf("%s: %d at %d\n ", #var, var, __LINE__)
#define PRINTDV3(var) printf("%s: (%f,%f,%f) at %d\n ", #var, var.x, var.y, var.z, __LINE__)

#define STRINGIFY(var) #var
#define NAME(var) STRINGIFY(var)

#define SYSDIM 7
#define BODY moon

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
    double operator* (const double_v3&) const;//dot product
    double_v3& operator= (const double_v3&);
    bool operator== (const double_v3&);
    
    double mag(void) const;
    double_v3 unit(void) const;
    double r(const double_v3& r) const;
    
};

//cross product
double_v3 cross (const double_v3& u, const double_v3& v);

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
    double_v3 w;
    
    void init2(void);//quite dodge
};

extern body earth;
extern body moon;


int veccpy(double* a, double* b);
int veccpy(std::array<double,SYSDIM>& a, double* b);
int veccpy(double* a, std::array<double,SYSDIM>& b);
int veccpy(std::array<double,SYSDIM>& a, const std::array<double,SYSDIM>& b);

int vecpack(double* a, const double t, const double_v3& b, const double_v3& c);

template<typename T>
double vec_unpack_t(const T a){
    return a[0];
}

template<typename T>
double_v3 vec_unpack_r(const T a){
    double_v3 v = {a[1], a[2], a[3]};
    return v;
}

template<typename T>
double_v3 vec_unpack_v(const T a){
    double_v3 v = {a[4], a[5], a[6]};
    return v;
}

template<typename T>
bool ordered(T a){
    return true;
}
template<typename T, typename... Args>
bool ordered(T a, T b, Args... args) {
  return (a <= b) and ordered(b, args...);
}

int mod(const int a, const int b);

double TWR_to_thrust(double TWR, double R);

#endif