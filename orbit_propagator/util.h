#ifndef UTIL_H
#define UTIL_H

#include <array>

#define PRINTFLT(var) printf("%s: %f at %d\n", #var, var, __LINE__)
#define PRINTINT(var) printf("%s: %d at %d\n ", #var, var, __LINE__)
#define PRINTDV3(var) printf("%s: (%f,%f,%f) at %d\n ", #var, var.x, var.y, var.z, __LINE__)

#define STRINGIFY(var) #var
#define NAME(var) STRINGIFY(var)

#define SYSDIM 8 //t,x,y,z,dx,dy,dz,m,
#define AUXDIM 3 //Ftx,Fty,Ftz
#define BODY moon

#define G 6.67E-11


int rotate(double_v3& r,const double_v3& u, const double angle);

struct quat {
    
    double r;
    double i;
    double j;
    double k;
    
    quat& operator+ (const quat&);
    quat& operator* (const quat&);
    quat inv (void) const;
    quat conj (void) const;
    double mag (void) const;
    
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

struct craft{
    double mass_wet;//kg
    double mass_dry; //kg
    double thrust; //kN
    double ISP;  //s

};

extern craft rocket;


template<typename T>
int veccpy(T& a, T& b){
    for(size_t i = 0; i < a.size(); i++){
        a[i] = b[i];
    }
    return 0;
}

int vecpack(std::array<double,SYSDIM>& a, const double t, const double_v3& b, const double_v3& c, const double m);
int vecpack(std::array<double,AUXDIM>& a, const double_v3& b);

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
double vec_unpack_mass(const T a){
    return a[7];
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