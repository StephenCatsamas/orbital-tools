from math import *
import matplotlib.pyplot as plt
import numpy as np
import solver
import csv

print("====LANDING BURN CALULATOR====")
print("rev: 0.2.0")
print("==============================")

G = 6.67E-11;
##step size error
err = 1E-3
err = np.array([[100*err,err,err,err]])



class statestruct():
    def __init__(self):
        self.complete = False;
        self.body = None;
        self.rlh = float('nan');
        self.height = float('nan');
        self.azvel = float('nan');
        self.rsea = float('nan');
        self.rs = float('nan');
        self.M = float('nan');
        self.T = float('nan');
        self.z0 = np.array([[float('nan'),float('nan'),float('nan'),float('nan')]]);
        self.vfs = float('nan');

    def orbit(self, i, v):
        if i == 0:
            self.body = v.strip()
        if i == 1:
            self.rlh = float(v)
        if i == 2:
            self.height = float(v)
        if i == 3:
            self.z0[0,1] = float(v)
        if i == 4:
            self.z0[0,2] = float(v)
        if i == 5:
            self.azvel = float(v)

    def planet(self, i, v):
        if i == 0:
            self.rsea = float(v)
        if i == 1:
            self.M = float(v)
        if i == 2:
            self.T = float(v)

state = statestruct();



def load_init():
    with open('orbit.dat') as csvfile:
        orbit_r = csv.reader(csvfile, delimiter=',', quotechar='|')
        for i,row in enumerate(orbit_r):
            state.orbit(i,row[1])
    with open('planets.dat') as csvfile:
        planet_r = csv.reader(csvfile, delimiter=',', quotechar='|')
        for i,row in enumerate(planet_r):
            if row[0] == state.body:
                for i in range(len(row)-1):
                    state.planet(i , row[i+1]);
                break;
                
                
    state.rs = state.rsea + state.rlh;
    state.vfs = 2*pi/state.T;
    state.z0[0,0] = state.height + state.rsea
    state.z0[0,3] = state.azvel / state.z0[0,0]

    

def round_non(v,n):
    if(isnan(v)):
        return "---"
    else:
        return round(v,n)

def print_int_cond():
    print("1. Planet Radius:", round_non(state.rsea,0), "m")
    print("2. Landing Height:", round_non(state.rs,0), "m," , round_non(state.rs-state.rsea,0), "m (ASL)")
    print("3. Planet Mass:", round_non(state.M,3), "kg")
    print("4. Planet Day Length:", round_non(state.T,3), "s")
    print("5. Initial Height:", round_non(state.z0[0,0],0), "m,", round_non(state.z0[0,0]-state.rsea,0), "m (ASL),", round_non(state.z0[0,0]-state.rs,0), "m (AGL)")
    print("5. Initial Azumith:", round_non(state.z0[0,1],3), "rad")
    print("6. Initial Radial Velocity:", round_non(state.z0[0,2],2), "m/s")
    print("7. Initial Azumith Velocity:", round_non(state.z0[0,3]*state.z0[0,0],3), "m/s")
    vr = state.z0[0,2]
    vf = state.z0[0,3]
    r = state.z0[0,0]
    vo = sqrt(vr*vr + r*vf*r*vf)
    vs = sqrt(vr*vr + (r*vf - state.rs* state.vfs)*(r*vf - state.rs* state.vfs))
    print("8. Initial Velocity:", round_non(vo,1), "m/s (Orbital),", round_non(vs,1), "m/s (Surface)")
    

def check_key(chr):
    keys = [str(i) for i in range(1,9)]
    if chr in keys:
        return True
    else:
        return False
        

def get_TWR(thrust, r = state.rs):
    m = 1
    TWR = thrust / (G*state.M*m/(r*r))
    return TWR
def get_thrust(TWR, r = state.rs):
    m = 1
    thrust = TWR * (G*state.M*m/(r*r))
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
        v = sqrt(vr*vr + (r*vf - state.rs* state.vfs)*(r*vf - state.rs* state.vfs))##accounting for surface rotation rs*vfs

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

    load_init()
    print_int_cond()

    exit();
    
    
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
