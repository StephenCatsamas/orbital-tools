#include <stdio.h>

#define PRINTFLT(var) printf("%s: %f\n", #var, var)

int divides(int i, int j){
    return i%j;
}


int primes(void){
  for(int i = 1; i < 1E5; i++){
        for(int j = 2; j < i/2 + 1; j++){
            if(divides(i,j) == 0){
                break;
            }
            if(j == i/2){
                printf("Prime! %d\n", i);
            }
        }
    }  
    return 0;
}

struct quat {
    double r;
    double i;
    double j;
    double k;
    
    quat& operator+ (const quat& q){
    r += q.r;    
    i += q.i;    
    j += q.j;    
    k += q.k;
    return *this;
    }
    
    quat& operator* (const quat& q){
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
    
};

int main(int argc, char **argv){
	printf("Quaternions!\n");
    
    //primes();
    struct quat q = {1,2,3,4};
    struct quat p = {-4,-3,-2,-1};
    
    q = q*p;
    
    PRINTFLT(q.r);
    PRINTFLT(q.i);
    PRINTFLT(q.j);
    PRINTFLT(q.k);
    
    
	return 0;
}




