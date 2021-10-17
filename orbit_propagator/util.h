

#define PRINTFLT(var) printf("%s: %f at %d\n", #var, var, __LINE__)
#define PRINTINT(var) printf("%s: %d at %d\n ", #var, var, __LINE__)
#define PRINTDV3(var) printf("%s: (%f,%f,%f) at %d\n ", #var, var.x, var.y, var.z, __LINE__)

#define SYSDIM 2

int veccpy(double* a, double* b);

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
    
    double_v3& operator+ (const double_v3&);
    double_v3& operator- (const double_v3&);
    bool operator== (const double_v3&);
    
};
