from math import *
import matplotlib.pyplot as plt

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

zout = []
xout = []

while((vz < 1) and (z > 0)):
    
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
    
    vz += -g*Dt + tz
    vx += tx
    
    z += vz*Dt 
    x += vx*Dt
    
    zout.append(z)
    xout.append(x)

print(len(xout))
print(v)
print(z)

fig = plt.figure() 
ax = fig.add_subplot(1, 1, 1) 

ax.plot(xout,zout)
plt.axhline(0, color='r')
plt.show()




#Function: rK4
#Inputs:
# t0 - initial time
# y0 - initial position
# h - step size
#system - function which defines the system to solve, when given a time and position sh
ould return the x,y time derivatives
#Returns:
# array of [x',y',t'], where x',y' are a step of the system at time t0 + h
def rK4(t0, y0, h, system):
    # Runge-Kutta 4th Order Function (primer)

    #coefficients calulated as per RK4 algorithm
    k1 = h*system(t0, y0)
    k2 = h*system(t0 + 0.5*h, y0 + 0.5*k1)
    k3 = h*system(t0 + 0.5*h, y0 + 0.5*k2)
    k4 = h*system(t0 + h, y0 + k3)

    #increment solution by change found
    y = y0 + (1/6)*(k1+2*k2+2*k3+k4);
    t = np.array([t0 + h])

    #return the new state vector of the system
    return np.concatenate((y, t))

#Function: stepper
#Integrates the ODE specificed by system and returns position vectors
#Inputs:
# z0 - initial state of system [[x0,y0,t0]]
# stop - ammount of time to simulate the system for
#system - function which defines the system to solve. Defined as in rk4 function
# h - initial step size
# h_adj - if True then stepper will adjust time step dynamically, else it will keep the
timestep constant at h.
# bound - if larger than zero will print a warming if solution is further from origin t
han value
#Returns:
# array [[x_i,y_i,t_i]] of solution to system over time from t0 to t0+stop.
def stepper(z0, stop, system, h = 0.01, h_adj = True, bound = 0):
    do_dyn_step = h_adj
    ##set first entry of output array to initial conditions
    loc= z0

    t0 = z0[0,2]

    #define array to give maximum acceptible error
    error = np.array([[1E-9,1E-9]])

    #find abs max norm of vector for later comparison (\Delta_0)
    merror = np.max(np.abs(error))

    #set initial dynamic step size to initial step size
    hdyn = h

    #loc[-1][2] is the latest step time, we want to itterate untill we have gotten to t
    his time
    while(stop+t0 - loc[-1][2] > 0):

    try:
    #enable or disable dynamic step size
    if(do_dyn_step):
    #we now need to loop over one step untill we find the appropriate step
    size
    #for our set error, for most steps this should only take one itteration
    #we define this flag which to loop until we have found the right step s
    ize
    retry_flag = True
    while(retry_flag):
    #take one double step in the solution
    loc1=rK4(loc[-1][2], loc[-1][0:2], 2*hdyn, system)
    #take two steps in the solution
    loc_hf=rK4(loc[-1][2], loc[-1][0:2], hdyn, system)
    loc2=rK4(loc_hf[2], loc_hf[0:2], hdyn, system)
    #note here that loc is an array which defines a state vector, loc
    [i][0:2] defines the x,y position
    #and loc[i][2] is the time
    #take the difference of the double step and 2 steps solution positi
    ons as an estimate for the error (\Delta_1)
    delta = loc1[0:2] - loc2[0:2]



