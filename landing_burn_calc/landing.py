from math import *
import matplotlib.pyplot as plt
import solver

print("====LANDING BURN CALULATOR====")
print("rev: 0.0.1")
print("==============================")


# Constants all SI
Dt = 1
Hi = 7000
vi = 160
alt = -0#deg
TWR = 1.67
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
def landing_system(t,vec):   
    #unpack location vector
    vx,vy=vec
    
    #trust
    t = TWR*g;
    
    #decompose to retrograde burn
    v = sqrt(vx*vx + vz*vz)

    if(v > 1):
        tz = -vz/v * t
        tx = -vx/v * t
    else:
        tz = t
        tx = 0

    #get flow velocities from the get_flow function
    dvxdt = -tx
    dvzdt = -g + tz
    
    #return array of the derivatives
    return np.array([dvxdt,dvzdt])

#define initial contitions
z0 = np.array([[1,1,0]])#x,y,t
#define initial step size
h  = 0.1

#solve the test ODE
loc = stepper(z0,10, test_system, h)
loc = stepper(z0,10, test_system, h)


fig = plt.figure() 
ax = fig.add_subplot(1, 1, 1) 

ax.plot(xout,zout)
plt.axhline(0, color='r')
plt.show()




