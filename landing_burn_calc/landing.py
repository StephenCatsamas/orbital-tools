from math import *
import matplotlib.pyplot as plt
import numpy as np
import solver

print("====LANDING BURN CALULATOR====")
print("rev: 0.0.1")
print("==============================")


# Constants all SI
Dt = 1
Hi = 6860
vi = 162.5
alt = -0#deg
TWR = 1.7
g = 0.491

# A note on coordiantes here +z is up and +x is downrange

z = Hi;
x = 0;

vz = sin(radians(alt))*vi 
vx = cos(radians(alt))*vi


#Function: landing_system
#landing burn implementatation of the system function for rK4
#Inputs: 
#    t - time 
#  vec - initial velocity 
#Returns:
#    array of dxdt,dydt, the time derivatives of velocity evaluated at (vec,t)
def landing_config(TWRc):
    def landing_system(t,vec):   
        #unpack location vector
        x,z,vx,vz=vec
        
        #trust
        t = TWRc*g;
        
        #decompose to retrograde burn
        v = sqrt(vx*vx + vz*vz)

        if(v > 1):
            tz = -vz/v * t
            tx = -vx/v * t
        else:
            tz = t
            tx = 0

        #system of equations
        dxdt = vx
        dzdt = vz
        dvxdt = tx
        dvzdt = -g + tz
        
        #return array of the derivatives
        return np.array([dxdt,dzdt,dvxdt,dvzdt])
    return landing_system

def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect
def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

#define initial contitions
z0 = np.array([[0,3000,157,0,0]])#x,y,t
#define initial step size
h  = 0.1

err = 1E-3

err = np.array([[100*err,100*err,err,err]])

solution = False

TWR = 2
dTWR = 0.0001

TWRhigh = 1
TWRlow = 1
TWRhighset = False
TWRlowset = False

while(not solution):
    solution = True
    allloc = z0
    landed = False;
    
    print(TWR)
    if(TWR <= 1):
        print("error system will not converge")
        exit()
    system = landing_config(TWR)

    while( not landed):
        #solve the test ODE
        stloc = np.array([allloc[-1]])
        loc = solver.stepper(stloc, 10, system, h, h_adj = True, error = err)
        allloc = np.vstack((allloc,loc))

        x,z,vx,vz,t= allloc[-1]
        
        if(vz > 0):
            landed = True
    
    end = [0,100,-3,3]

    for i,vec in enumerate(loc):
        x,z,vx,vz,t= vec
        hit = False
        #case that one of the points lies in the box
        if ((end[0] <= z <= end[1]) and (end[2] <= vz <= end[3])):
            hit = True
            break;
            # print('e')
        else:
        #case that the line segment crosses the box
            if(i != 0):
                z1 = z;
                vz1 = vz;
                z2 = loc[i-1][1]
                vz2 = loc[i-1][3]
               
                if(intersect([z1,vz1],[z2,vz2],[end[1],end[2]],[end[1],end[3]])):
                    hit = True
                    break;
                    # print('a')
                if(intersect([z1,vz1],[z2,vz2],[end[1],end[3]],[end[0],end[3]])):
                    hit = True
                    break;
                    # print('b')
                if(intersect([z1,vz1],[z2,vz2],[end[0],end[3]],[end[0],end[2]])):
                    hit = True
                    break;
                    # print('c')
                if(intersect([z1,vz1],[z2,vz2],[end[0],end[2]],[end[1],end[2]])):
                    hit = True
                    break;
                    # print('d')

    

    solution = hit
    
    vzmin = None
    for vec in allloc:
        x,z,vx,vz,t = vec
        
        if(abs(vz) >= abs(vx)):
            if(vzmin == None):
                vzmin = abs(vz)
            if(abs(vz) < vzmin):
                vzmin = abs(vz)
                deltZ = z
                
    if(not(TWRhighset and TWRlowset)):
        if(deltZ > 0):
             TWRhighset = True
             TWRhigh = TWR
        else:
             TWRlowset = True
             TWRlow = TWR
    
        if(not TWRhighset):
            TWR = 1 + 2*(TWR-1)
        elif(not TWRlowset):
            TWR = 1 + 0.5*(TWR-1)
    if(TWRhighset and TWRlowset):
        if(deltZ > 0):
            if(TWR > TWRhigh):
                print('TWR not bounded')
            TWRhigh = TWR
        else:
            if(TWR < TWRlow):
                print('TWR not bounded')
            TWRlow = TWR
        TWR = (TWRhigh + TWRlow)/2

    
    if(TWR < 1.005):
        TWR = 1.005
    
    

fig = plt.figure() 
ax = fig.add_subplot(1, 2, 1) 
ax2 = fig.add_subplot(1, 2, 2) 
ax.plot(allloc[:,0],allloc[:,1])
ax2.plot(allloc[:,1],allloc[:,3])
ax.axhline(0, color='r')
ax2.axvline(0, color='r')
plt.show()




