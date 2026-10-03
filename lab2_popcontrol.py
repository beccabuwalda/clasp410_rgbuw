'''
This code was used for the lab Population models, chaos, ordinary differential equations
and so much more
follow this code to model the Lokta-Volterra equations of ecology and population...
as well as simple weather model to simulate deterministic chaos
'''
#importing necessary packages
import numpy as np
import matplotlib
matplotlib.use('TkAgg') #setting up back end
import matplotlib.pyplot as plt
plt.ion() #this makes interactive plots activate
import scipy as sp
from scipy.integrate import ode #used later for ODE45

'''
The twist!!
Orcas consider seals to be both competitors and prey. In this situation, we are separating the competitor and predator-prey relationships for easier conceptualization
in reality, a more realistic model would contain a combined predator-prey and competitor functions into one model
'''


'''
#hand-coded competition model: Euler Method
'''
orca_init = 0.3 #predator initial population density
seal_init = 0.6 #prey initial population density
N_init = [orca_init, seal_init]
dt = 1.0 #time in years
a = 1 #constant
b = 2 #constant
c = 1 #constant
d = 3 #constant
#t = our time steps

#here, i am defining my function arctic that is describing the functions of competition dynamics
def arctic(N):
    a = 1 #constant
    b = 2 #constant
    c = 1 #constant
    d = 3 #constant 
    N1,N2 =N[0],N[1] #we have made a vector of our species 1 and species 2 population densities, so it can return as one value as opposed to having to sparate it
    dN1dt = a*N1*(1-N1)-b*N1*N2 #changes of population density of species 1, orca, for competition
    dN2dt = c*N2*(1-N2)-d*N1*N2 #changes of population density of species 2, seal, for competition
    return np.array([dN1dt,dN2dt]) #report the changes of pop density after a time step with a width size defined by dt
dNdt=arctic(N_init) #the initial change in pop densities over first time step in our func artic

#for later for figures
plt.figure(figsize=(8,8))

#Euler method loop of Lokta-Voltera Competition Model
N = np.copy(N_init)
for t in range(1):
    #N_init = np.copy(N_curr) #in our loop, we change our initial condition to be the next step, as we learned form the euler method
    time_history=[0.0] #time bin to graph later each step rather than it being erased
    popden_history=[N]
    dNdt=arctic(N)
    N_next = N + dt * dNdt
    N = np.copy(N_next)
 #plotting complete trajectory for current dt
    plt.plot(time_history, popden_history, label=f'dt ={dt}')
#Plot aesthetics
plt.xlabel('Time (Years)')
plt.ylabel('Population Density (kg/km^3)')
plt.title('Population Density: Competition of Orcas and Seals')
plt.legend(loc='upper right')
plt.ioff()
plt.show()

'''
#hand-coded predatory-prey model: Euler Method
orca_init = 0.3 #predator initial population density
seal_init = 0.6 #prey initial population density
N_init = [orca_init, seal_init]
dt = 0.5 #time in years
a = 1 #constant
b = 2 #constant
c = 1 #constant
d = 3 #constant
#t = our time steps
def arctic(N):
    a = 1 #constant
    b = 2 #constant
    c = 1 #constant
    d = 3 #constant 
    N1,N2 =N[0],N[1]
    dN1dt = a*N1*(1-N1)-b*N1*N2 #changes of population density of 
    dN2dt = -c*N2*(1-N2)+d*N1*N2
    return [dN1dt,dN2dt]
dNdt=arctic(N_init)


#rk45 predator-prey model



#rk45 competition model




#
'''