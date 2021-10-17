#include "util.h"

int veccpy(double* a, double* b){
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

double_v3& double_v3::operator+ (const double_v3& v){ 
    x += v.x;    
    y += v.y;    
    z += v.z;
    return *this;
}

double_v3& double_v3::operator- (const double_v3& v){ 
    x -= v.x;    
    y -= v.y;    
    z -= v.z;
    return *this;
}

bool double_v3::operator== (const double_v3& v){ 
    if((x == v.x) and (y == v.y) and (z == v.z) ){
        return true;
    }else{
        return false;
    }
}





