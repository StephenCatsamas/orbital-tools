import numpy as np
import math

#Function: rK4 
#Inputs: 
#    t0 - initial time 
#    y0 - initial position 
#    h  - step size 
#system - function which defines the system to solve, when given a time and position should return the x,y time derivatives
#Returns:
#   array of [x',y',t'], where x',y' are a step of the system at time t0 + h
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
#    z0 - initial state of system [[x0,y0,t0]]
#  stop - ammount of time to simulate the system for 
#system - function which defines the system to solve. Defined as in rk4 function
#     h - initial step size 
# h_adj - if True then stepper will adjust time step dynamically, else it will keep the timestep constant at h.
# bound - if larger than zero will print a warming if solution is further from origin than value  
#Returns:
#    array [[x_i,y_i,t_i]] of solution to system over time from t0 to t0+stop.
def stepper(z0, stop, system, h = 0.01, h_adj = True, bound = 0, error = 0):
    do_dyn_step = h_adj
    ##set first entry of output array to initial conditions
    loc= z0
    
    t0 = z0[0,-1]
 
    #define array to give maximum acceptible error
    # error = np.array([[1E-9,1E-9,1E-9,1E-9]])
    
    #find abs max norm of vector for later comparison (\Delta_0)
    merror = np.max(np.abs(error))
    
    #set initial dynamic step size to initial step size
    hdyn = h
    
    #loc[-1][-1] is the latest step time, we want to itterate untill we have gotten to this time
    while(stop+t0 - loc[-1][-1] > 0):
    
        try:
            #enable or disable dynamic step size
            if(do_dyn_step):
                #we now need to loop over one step untill we find the appropriate step size
                #for our set error, for most steps this should only take one itteration
                #we define this flag which to loop until we have found the right step size
                retry_flag = True

                while(retry_flag):

                    #take one double step in the solution
                    loc1=rK4(loc[-1][-1], loc[-1][0:-1], 2*hdyn, system)

                    #take two steps in the solution
                    loc_hf=rK4(loc[-1][-1], loc[-1][0:-1], hdyn, system)
                    loc2=rK4(loc_hf[-1], loc_hf[0:-1], hdyn, system)

                    #note here that loc is an array which defines a state vector, loc[i][0:2] defines the x,y position
                    #and loc[i][2] is the time

                    #take the difference of the double step and 2 steps solution positions as an estimate for the error (\Delta_1)
                    delta = loc1[0:-1] - loc2[0:-1]

                    #find the abs max norm of the error
                    mdelta = np.max(np.abs(delta))


                    #this if-statement is the dynamic step size adjustment
                    if mdelta <= merror:
                        #if the error of the step is less than the threshold we can go to the next step with a larger step size
                        if (mdelta!=0):
                            hdyn = 0.9 * hdyn * math.pow(merror/mdelta, 0.2)
                        retry_flag = False
                    else:
                        #otherwise we have more error than our threshold, so we should choose a smaller step size and try again
                        hdyn = 0.9 * hdyn * math.pow(merror/mdelta, 0.25)
                        retry_flag = True
                    #after we have exited the loop we have a solutuion with our error threshold,

            else:
                    #if not using dynmic step size then we just need to
                    #take a step in the solution
                    loc2=rK4(loc[-1][-1], loc[-1][0:-1], hdyn, system)
        except ValueError:
            break
            
        #here we add a warning if the solution has a position outside of the given bounds,
        #this is useful when calulating lyapunov exponents
        if (bound > 0):
            if (np.linalg.norm(z) > bound):
                print('Out of bounds warning')
        
        #if the most recent step oversteps the endpoint, we will should not append this point
        #and set the step size to the distance to the endpoint 
        if (stop+t0 - loc2[-1] < 0):
                hdyn = stop+t0 - loc[-1][-1]
                do_dyn_step = False
                continue
                
        #another error that we should check for is if the step size is zero
        #this is generally caused by a floating point math error
        #it is hard to handle this error, so we will just exit and warn the user
        if (loc[-1][-1] == loc2[-1]):
            print('Loop error')
            break
        
        
        #Add last step to the output array
        loc=np.vstack((loc,loc2))
        
        
        
    #after we have simulate the desired for the desired time we will return the output array    
    return loc
    
