from math import *
import matplotlib.pyplot as plt
import numpy as np
import solver

print("====LANDING BURN CALULATOR====")
print("rev: 0.0.1")
print("==============================")

rs = 60E3;#surface height
G = 6.67E-11;
M = 2.6457E19;


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
        r,f,vr,vf=vec
        
        m = 1;
        
        #trust
        t = TWRc*G*M*m/(rs*rs);
        t = 0
        #decompose to retrograde burn
        v = sqrt(vr*vr + r*vf*r*vf)

        if(v > 1):
            tr = -vr/v * t
            tf = -r*vf/v * t
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

#define initial contitions
z0 = np.array([[63E3,0,0,235/63E3,0]])#x,y,t
#define initial step size
h  = 1

err = 1E-3

err = np.array([[100*err,100*err,err,err]])

solution = False

##########
system = landing_config(0)
allloc = solver.stepper(z0, 800000, system, h, h_adj = True, error = err)
############
print(allloc)

fig = plt.figure(figsize=(10,4)) 
ax = fig.add_subplot(1, 2, 1) 
ax2 = fig.add_subplot(1, 2, 2) 
ax.plot(allloc[:,1],allloc[:,0])
ax.set_xlabel('Azumith $\phi$ (rad)');
ax.set_ylabel('Radius $r$ (m)');
ax2.plot(allloc[:,0],allloc[:,2])
ax2.set_xlabel('Radius $r$ (m)');
ax2.set_ylabel('Radial Velocity $v_r$ (m/s)');
# ax.axhline(0, color='r')
# ax2.axvline(0, color='r')
plt.show()




# TWR = 2
# dTWR = 0.0001

# TWRhigh = 1
# TWRlow = 1
# TWRhighset = False
# TWRlowset = False

# while(not solution):
    # solution = True
    # allloc = z0
    # landed = False;
    
    # print(TWR)
    # if(TWR <= 1):
        # print("error system will not converge")
        # exit()
    # system = landing_config(TWR)

    # while( not landed):
        # #solve the test ODE
        # stloc = np.array([allloc[-1]])
        # loc = solver.stepper(stloc, 10, system, h, h_adj = True, error = err)
        # allloc = np.vstack((allloc,loc))

        # x,z,vx,vz,t= allloc[-1]
        
        # if(vz > 0):
            # landed = True
    
    # end = [0,100,-3,3]

    # for i,vec in enumerate(loc):
        # x,z,vx,vz,t= vec
        # hit = False
        # #case that one of the points lies in the box
        # if ((end[0] <= z <= end[1]) and (end[2] <= vz <= end[3])):
            # hit = True
            # break;
            # # print('e')
        # else:
        # #case that the line segment crosses the box
            # if(i != 0):
                # z1 = z;
                # vz1 = vz;
                # z2 = loc[i-1][1]
                # vz2 = loc[i-1][3]
               
                # if(intersect([z1,vz1],[z2,vz2],[end[1],end[2]],[end[1],end[3]])):
                    # hit = True
                    # break;
                    # # print('a')
                # if(intersect([z1,vz1],[z2,vz2],[end[1],end[3]],[end[0],end[3]])):
                    # hit = True
                    # break;
                    # # print('b')
                # if(intersect([z1,vz1],[z2,vz2],[end[0],end[3]],[end[0],end[2]])):
                    # hit = True
                    # break;
                    # # print('c')
                # if(intersect([z1,vz1],[z2,vz2],[end[0],end[2]],[end[1],end[2]])):
                    # hit = True
                    # break;
                    # # print('d')

    

    # solution = hit
    
    # vzmin = None
    # for vec in allloc:
        # x,z,vx,vz,t = vec
        
        # if(abs(vz) >= abs(vx)):
            # if(vzmin == None):
                # vzmin = abs(vz)
            # if(abs(vz) < vzmin):
                # vzmin = abs(vz)
                # deltZ = z
                
    # if(not(TWRhighset and TWRlowset)):
        # if(deltZ > 0):
             # TWRhighset = True
             # TWRhigh = TWR
        # else:
             # TWRlowset = True
             # TWRlow = TWR
    
        # if(not TWRhighset):
            # TWR = 1 + 2*(TWR-1)
        # elif(not TWRlowset):
            # TWR = 1 + 0.5*(TWR-1)
    # if(TWRhighset and TWRlowset):
        # if(deltZ > 0):
            # if(TWR > TWRhigh):
                # print('TWR not bounded')
            # TWRhigh = TWR
        # else:
            # if(TWR < TWRlow):
                # print('TWR not bounded')
            # TWRlow = TWR
        # TWR = (TWRhigh + TWRlow)/2

    
    # if(TWR < 1.005):
        # TWR = 1.005
    
    




