from math import *
import matplotlib.pyplot as plt
import numpy as np
import solver

print("====LANDING BURN CALULATOR====")
print("rev: 0.1.0")
print("==============================")

#Minmus
# rs = 60E3;#surface height
# G = 6.67E-11;
# M = 2.6457E19;
# T = 40400;

#Mun
rsea = 200E3;#surface height
rs = 4E3 + rsea;#surface height
G = 6.67E-11;
M = 9.76E20;
T = 138984;

#define initial contitions
z0 = np.array([[11E3+rs,0,0,550/(11E3+rs),0]])#x,y,t
#define initial step size


err = 1E-4
err = np.array([[100*err,100*err,err,err]])

vfs = 2*pi/T;

def print_init_cond():
    print("Planet Radius:", round(rsea,0), "m")
    print("Landing Height:", round(rs,0), "m," , round(rs-rsea,0), "m (ASL)")
    print("Planet Mass:", round(M,3), "kg")
    print("Planet Day Length:", round(T,3), "s")
    print("Initial Height:", round(z0[0,0],0), "m,", round(z0[0,0]-rsea,0), "m (ASL),", round(z0[0,0]-rs,0), "m (AGL)")
    print("Initial Azumith:", round(z0[0,1],3), "rad")
    print("Initial Radial Velocity:", round(z0[0,2],2), "m/s")
    print("Initial Azumith Velocity:", round(z0[0,3]*z0[0,0],3), "m/s")
    vr = z0[0,2]
    vf = z0[0,3]
    r = z0[0,0]
    vo = sqrt(vr*vr + r*vf*r*vf)
    vs = sqrt(vr*vr + (r*vf - rs* vfs)*(r*vf - rs* vfs))
    print("Initial Velocity:", round(vo,1), "m/s (Orbital),", round(vs,1), "m/s (Surface)")
    

def get_TWR(thrust, r = rs):
    m = 1
    TWR = thrust / (G*M*m/(r*r))
    return TWR
def get_thrust(TWR, r = rs):
    m = 1
    thrust = TWR * (G*M*m/(r*r))
    return thrust

#Function: landing_system
#landing burn implementatation of the system function for rK4
#Inputs: 
#    t - time 
#  vec - initial velocity 
#Returns:
#    array of dxdt,dydt, the time derivatives of velocity evaluated at (vec,t)
def landing_config(thrust):
    def landing_system(t,vec):   
        #unpack location vector
        r,f,vr,vf=vec
        
        m = 1;
        #trust
        t = thrust#TWRc*G*M*m/(rs*rs);
        #decompose to retrograde burn
        v = sqrt(vr*vr + (r*vf - rs* vfs)*(r*vf - rs* vfs))##accounting for surface rotation rs*vfs

        if(v > 1):
            tr = -vr/v * t
            tf = -(r*vf - rs* vfs)/v * t
        else:
            tr = t
            tf = 0

        fG = - G*M*m/(r*r)

        #system of equations
        drdt = vr
        dfdt = vf
        dvrdt = (fG + tr)/m + r*vf*vf
        dvfdt = (tf - 2*vr*vf)/(m*r)
        
        #return array of the derivatives
        return np.array([drdt,dfdt,dvrdt,dvfdt])
    return landing_system

def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect
def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)


def main():
    print_init_cond()
    solution = False
    
    h  = 1
    
    thrust = 5
    dthrust = 0.0001

    thrusthigh = 1
    thrustlow = 1
    thrusthighset = False
    thrustlowset = False

    print('=======Solving TWR============')
    while(not solution):
        # solution = True
        allloc = z0
        landed = False;
        
        TWR = get_TWR(thrust)
        print("TWR:", round(TWR,5))
        if(TWR <= 1):
            print("error system will not converge")
            exit()
        system = landing_config(thrust)

        while( not landed):
            #solve the test ODE
            stloc = np.array([allloc[-1]])
            loc = solver.stepper(stloc, 10, system, h, h_adj = True, error = err)
            allloc = np.vstack((allloc,loc))

            r,f,vr,vf,t= allloc[-1]
            
            if((vr > 0) and (abs(vr) >= abs(r*vf - rs*vfs))):
                landed = True
        
        end = [rs+0,rs+100,-3,3]

        for i,vec in enumerate(loc):
            r,f,vr,vf,t= vec
            hit = False
            #case that one of the points lies in the box
            if ((end[0] <= r <= end[1]) and (end[2] <= vr <= end[3])):
                hit = True
                break;
                # print('e')
            else:
            #case that the line segment crosses the box
                if(i != 0):
                    r1 = r;
                    vr1 = vr;
                    r2 = loc[i-1][0]
                    vr2 = loc[i-1][2]
                   
                    if(intersect([r1,vr1],[r2,vr2],[end[1],end[2]],[end[1],end[3]])):
                        hit = True
                        break;
                        # print('a')
                    if(intersect([r1,vr1],[r2,vr2],[end[1],end[3]],[end[0],end[3]])):
                        hit = True
                        break;
                        # print('b')
                    if(intersect([r1,vr1],[r2,vr2],[end[0],end[3]],[end[0],end[2]])):
                        hit = True
                        break;
                        # print('c')
                    if(intersect([r1,vr1],[r2,vr2],[end[0],end[2]],[end[1],end[2]])):
                        hit = True
                        break;
                        # print('d')

        solution = hit
        vrmin = None
        for vec in allloc:
            r,f,vr,vf,t = vec
            
            if(abs(vr) >= abs(r*vf - rs*vfs)):
                if(vrmin == None):
                    vrmin = abs(vr)
                if(abs(vr) <= vrmin):
                    vrmin = abs(vr)
                    deltH = r-rs
                    
        if(not(thrusthighset and thrustlowset)):
            if(deltH > 0):
                 thrusthighset = True
                 thrusthigh = thrust
            else:
                 thrustlowset = True
                 thrustlow = thrust
            TWR = get_TWR(thrust)
            if(not thrusthighset):
                TWR = 1 + 2*(TWR-1)
            elif(not thrustlowset):
                TWR = 1 + 0.5*(TWR-1)
            thrust = get_thrust(TWR)
        if(thrusthighset and thrustlowset):
            if(deltH > 0):
                if(thrust > thrusthigh):
                    print('thrust not bounded')
                thrusthigh = thrust
            else:
                if(thrust < thrustlow):
                    print('thrust not bounded')
                thrustlow = thrust
            thrust = (thrusthigh + thrustlow)/2

        TWR = get_TWR(thrust)
        if(TWR < 1.005):
            thrust = get_thrust(1.005);
        
    print('=======Solution===============')    
    print('Min Velocity:', round(vrmin,0), 'm/s')
    print('Inflection Height:', round(deltH,0), 'm')
    
    fig = plt.figure(figsize=(10,8)) 
    ax = fig.add_subplot(2, 2, 1) 
    ax2 = fig.add_subplot(2, 2, 2) 
    ax3 = fig.add_subplot(2, 2, 3) 
    ax4 = fig.add_subplot(2, 2, 4) 
    ax.plot(allloc[:,1],allloc[:,0])
    ax.set_xlabel('Azumith $\phi$ (rad)');
    ax.set_ylabel('Radius $r$ (m)');
    ax.axhline(rs, color='r')
    ax2.plot(allloc[:,0],allloc[:,2])
    ax2.set_xlabel('Radius $r$ (m)');
    ax2.set_ylabel('Radial Velocity $v_r$ (m/s)');
    ax2.axvline(rs, color='r')
    ax3.plot(allloc[:,4],allloc[:,0])
    ax3.set_xlabel('Time $t$ (s)');
    ax3.set_ylabel('Radius $r$ (m)');
    ax3.axhline(rs, color='r')
    
    vsrfs = [sqrt(vr*vr + (r*vf - rs* vfs)*(r*vf - rs* vfs)) for vr,vf in zip(allloc[:,2],allloc[:,3])]
    
    ax4.plot(allloc[:,0],vsrfs)
    ax4.set_xlabel('Radius $r$ (m)');
    ax4.set_ylabel('Surface Velocity $v_s$ (m/s)');
    ax4.axvline(rs, color='r')
    ax4.invert_xaxis()
    plt.show() 

main();
