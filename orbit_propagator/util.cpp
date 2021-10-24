#include <math.h>
#include "util.h"



int veccpy(double* a, double* b){
    for(int i = 0; i < SYSDIM; i++){
        a[i] = b[i];
    }
    return 0;
}

int veccpy(std::array<double,SYSDIM>& a, double* b){
    for(int i = 0; i < SYSDIM; i++){
        a[i] = b[i];
    }
    return 0;
}

int veccpy(double* a, std::array<double,SYSDIM>& b){
    for(int i = 0; i < SYSDIM; i++){
        a[i] = b[i];
    }
    return 0;
}

int veccpy(std::array<double,SYSDIM>& a, const std::array<double,SYSDIM>& b){
    for(int i = 0; i < SYSDIM; i++){
        a[i] = b[i];
    }
    return 0;
}


quat& quat::operator+ (const quat& q){
    r += q.r;    
    i += q.i;    
    j += q.j;    
    k += q.k;
    return *this;
}
    
quat& quat::operator* (const quat& q){
    double rt = r*q.r - i*q.i - j*q.j - k*q.k;    
    double it = r*q.i + i*q.r - k*q.j + j*q.k;    
    double jt = r*q.j + j*q.r + k*q.i - i*q.k;    
    double kt = r*q.k + k*q.r - j*q.i + i*q.j;
    r = rt;
    i = it;
    j = jt;
    k = kt;
    return *this;
}

//dot product
double double_v3::operator* (const double_v3& r) const{ 
    double inner = x*r.x + y*r.y + z*r.z;
    return inner;
}

double_v3 cross (const double_v3& u, const double_v3& v){ 
    double_v3 w;
    w.x = u.y*v.z - u.z*v.y;    
    w.y = u.z*v.x - u.x*v.z;    
    w.z = u.x*v.y - u.y*v.x;
    return w;
}

double_v3& double_v3::operator/ (const double a){ 
    x /= a;    
    y /= a;    
    z /= a;
    return *this;
}

double_v3& double_v3::operator= (const double_v3& v){ 
    x = v.x;    
    y = v.y;    
    z = v.z;
    return *this;
}

bool double_v3::operator== (const double_v3& v){ 
    if((x == v.x) and (y == v.y) and (z == v.z) ){
        return true;
    }else{
        return false;
    }
}
double double_v3::mag(void) const{
    return sqrt(x*x + y*y + z*z);
}

double_v3 double_v3::unit(void) const{
    double_v3 unit = *this;
    unit = unit/unit.mag();
    return unit;
}

double double_v3::r(const double_v3& r) const{
    double ip = (*this)*r;
    double mag = r.mag();
    return ip/mag;
}

int vecpack(double* a, const double t, const double_v3& b, const double_v3& c){
    a[0] = t;
    a[1] = b.x;
    a[2] = b.y;
    a[3] = b.z;
    a[4] = c.x;
    a[5] = c.y;
    a[6] = c.z;
    return 0;
}


int mod(const int a, const int b){
    if (a >= 0){
        return a%b;
    }else{
        return (a%b + b)%b;
    }
}

double TWR_to_thrust(double TWR, double R){
    double m = 1;
    double thrust = TWR * (G*BODY.mass*m/(R*R));
    return thrust;
}

void body::init2(){
    w.x = 0;
    w.y = 0;
    w.z = 2*M_PI/rotational_period;   
}

struct body earth = {5.972E24, 6.371E6, 0, 86164.1};
struct body moon = {7.342E22, 1.7374E6, 0, 2.3606E6};

