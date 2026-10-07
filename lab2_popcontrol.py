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
import scipy as sp
from scipy.integrate import ode #used later for ODE45
from scipy.integrate import solve_ivp

'''
The twist!!
Orcas consider seals to be both competitors and prey. In this situation, we are separating the competitor and predator-prey relationships for easier conceptualization
in reality, a more realistic model would contain a combined predator-prey and competitor functions into one model.
'''


'''
#hand-coded competition model: Euler Method
🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭
'''

'''
to validate the code is running the correct way, perform the following tests by changing line's 45 and 46:
1. Species 1 goes extinct (N1=0 and N2=1)
Set orca_init=0
Set seal_init=1
2. Species 2 goes extinct (N1=0, N2=0)
Set orca_init=1
Set seal_init=0
3.Both species go extinct (N1=0 and N2=0)
Set orca_init=0
Set seal_init=0
4. The species coexist in an uneasy truce by setting N1 and N2 equal to constants
Set orca_init = 0.20
Set seal_init = 0.40
'''
#setting intial conditions
orca_init = 0.30 #predator initial population density
seal_init = 0.60  #prey initial population density
N_init = [orca_init, seal_init]
dt = 1 #time step width in years
num_steps = 100


#here, i am defining my function arctic that is describing the functions of competition dynamics
def arctic_comp(N):
    '''
    the function artic_comp does as follows
    defines constants from Lotka-Volterra Model 
    a constant  reproduction rate of species 1
    b constant scales impacts of Species 2 on Species 1
    c constant reproduction rate
    d constant scales impact of Species 1 on Species 2
    defines the species density of species 1 and species 2 as a vector
    calculates the dN1dt and dN2dt, i.e., changes in pop den related to time step width using L-V competition equations
    returns those values as an array
    Define dNdt of the arctic_comp function as our initial den values
    '''
    a = 1 #constant
    b = 2 #constant
    c = 1 #constant
    d = 3 #constant 
    N1,N2 =N[0],N[1] #we have made a vector of our species 1 and species 2 population densities, so it can return as one value as opposed to having to sparate it
    dN1dt = a*N1*(1-N1)-b*N1*N2 #changes of population density of species 1, orca, for competition
    dN2dt = c*N2*(1-N2)-d*N1*N2 #changes of population density of species 2, seal, for competition
    return np.array([dN1dt,dN2dt]) #report the changes of pop density after a time step with a width size defined by dt
dNdt=arctic_comp(N_init) #the initial change in pop densities over first time step in our func artic

#prepare figure
plt.figure(figsize=(8,8))

#storing history of data outside of loop
N=np.copy(N_init)
time_history=[0.0] #time bin to graph later each step rather than it being erased
popden_history=[np.copy(N)]

#Euler method loop of Lokta-Voltera Competition Model
for t in range(num_steps):
    #N_init = np.copy(N_curr) #in our loop, we change our initial condition to be the next step, as we learned form the euler method
    dNdt=arctic_comp(N)
    N_next = N + dt * dNdt

    # Change window view to only look at y values 0 to 1
    #N_next = np.maximum(0, N + dt * dNdt)
    #therefore, a species density at 0 stays at 0 and doesn't go beyond, aka it is extinct

    N = np.copy(N_next)

    #append new data into our history
    time_history.append((t+1)*dt)
    popden_history.append(np.copy(N))

#convert the history into an array for easy slicing
popden_history = np.array(popden_history)

#plotting complete trajectory for current dt
plt.plot(time_history, popden_history[:,0], label="Orca (N1)", color="blue")
plt.plot(time_history,popden_history[:,1], label="Seal (N2)", color="orange", linestyle='dashed')
#Plot aesthetics
plt.xlabel('Time (Years)')
plt.ylabel('Normalized Population Density (1 = Carrying Capacity')
plt.title('Population Density: Competition of Orcas and Seals')
plt.legend(loc='upper right')
plt.grid(True)
plt.show()


'''
#hand-coded predatory-prey model: Euler Method
🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭
'''


orca_init = 0.3 #predator initial population density
seal_init = 0.6 #prey initial population density
N_init = [seal_init, orca_init]
dt = .05 #time in years
num_steps = 1000

def arctic_pred(N):
    '''
    the function arctic_pred does as follows
    defines constants from Lotka-Volterra Model 
    a constant growth rate of prey
    b constant predator causing death of prey
    c constant consumed prey leading to predator offspring
    d constant natural death rate of predator (starvation)
    defines the species density of species 1 and species 2 as a vector
    calculates the dN1dt and dN2dt, i.e., changes in pop den related to time step width, using L-V predator/prey equations
    returns those values as an array
    Define dNdt of the arctic_pred function as our initial den values
    '''
    a = 1 #constant
    b = 2 #constant
    c = 1 #constant
    d = 3 #constant 
    N1,N2 =N[0],N[1] #N1 is seal, N2 is orca
    dN1dt =a*N1-b*N1*N2 #prey
    dN2dt =-c*N2+d*N1*N2 #predator
    return np.array([dN1dt,dN2dt]) #report the changes of pop density after a time step with a width size defined by dt
dNdt=arctic_pred(N_init) #the initial change in pop densities over first time step in our func artic
    
#prepare figure
plt.figure(figsize=(8,8))
    
#storing history of data outside of loop
N=np.copy(N_init)
time_history=[0.0] #time bin to graph later each step rather than it being erased
popden_history=[np.copy(N)]
    
#Euler method loop of Lokta-Voltera Competition Model
for t in range(num_steps):
    #N_init = np.copy(N_curr) #in our loop, we change our initial condition to be the next step, as we learned form the euler method
    dNdt=arctic_pred(N)
    N_next = N + dt * dNdt
    
    # Prevent populations from going below 0. This was an important check of code soundness because we can't have negative population densities.
    N_next = np.maximum(0, N + dt * dNdt)
    #therefore, a species density at 0 stays at 0 and doesn't go beyond, aka it is extinct
    
    N = np.copy(N_next)
    
    #append new data into our history
    time_history.append((t+1)*dt)
    popden_history.append(np.copy(N))
    
#convert the history into an array for easy slicing
popden_history = np.array(popden_history)
    
#plotting complete trajectory for current dt
plt.plot(time_history, popden_history[:,0], label="Seal (N1)", color="orange", linestyle="dashed")
plt.plot(time_history,popden_history[:,1], label="Orca (N2)", color="blue")
#Plot aesthetics
plt.xlabel('Time (Years)')
plt.ylabel('Normalized Population Density (1 = Carrying Capacity)')
plt.title('Population Density: Predator/Prey Relationship of Orcas and Seals')
plt.legend(loc='upper right')
plt.grid(True)
plt.show()


#plotting phase diagram of species 1 and species 2 
plt.plot(popden_history[:,0], popden_history[:,1])
plt.title('Eulers Method Phase Space Seal (N1) and Orca (N2) Population Density')
plt.xlabel('Normalized Seal Population Density (N1)')
plt.ylabel('Normalized Orca Population Density (N2)')
plt.show()


'''
Runge Kutta 4 solver
about rk45: a scipy product that can solve for an initial value problem for a system of ODEs
'''

'''
#rk45 competition model
🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭
'''


orca_init = 0.3 #species 1 initial population density
seal_init = 0.6 #species 2 initial population density
N_init = [orca_init, seal_init]
num_steps = 100
t_span = (0, num_steps) #time step width in years
t_eval = np.linspace(0, 50, 500)

def arctic_comp(t, N): #now defining function of t and pop den, since we are not manually deciding dt
    '''
    time is now added since we arent manually deciding dt anymore
    the function artic_comp does as follows
    defines constants from Lotka-Volterra Model 
    a constant  reproduction rate of species 1
    b constant scales impacts of Species 2 on Species 1
    c constant reproduction rate
    d constant scales impact of Species 1 on Species 2
    defines the species density of species 1 and species 2 as a vector
    calculates the dN1dt and dN2dt, i.e., changes in pop den related to time step width using L-V competition equations
    returns those values as an array
    Define dNdt of the arctic_comp function as our initial den values
    '''
    a = 1 #constant
    b = 2 #constant
    c = 1 #constant
    d = 3 #constant 
    K = 1 #carrying capacity
    N1,N2 =N[0],N[1] #orca species 1, seal species 2
    dN1dt = a*N1*(1-N1)-b*N1*N2 #changes of population density of species 1, orca, for competition
    dN2dt = c*N2*(1-N2)-d*N1*N2 #changes of population density of species 2, seal, for competition
    return ([dN1dt,dN2dt]) #report the changes of pop density after a time step with a width size defined by dt


sol = solve_ivp(
    fun=arctic_comp, #function
    t_span=t_span, #time span
    y0=N_init, #initial pop density vector
    method='RK45', #package using
    t_eval=t_eval #eval values at these time steps
)

plt.figure(figsize=(8, 5))
plt.plot(sol.t, sol.y[0], label='Orcas (N1)', color='blue')
plt.plot(sol.t, sol.y[1], label='Seals (N2)', color='orange',  linestyle="dashed")
plt.title('RK45 Arctic Competition Model Simulation')
plt.xlabel('Time (Years)')
plt.ylabel('Normalized Population Density (Carrying Capacity = 1)')
plt.legend()
plt.grid(True)
plt.show()

#plotting phase diagram of species 1 and species 2 
plt.plot(sol.y[0], sol.y[1])
plt.title('RK45 Method Phase Space Competition between Orcas (N1) and Seals (N2)')
plt.xlabel('Normalized Orca Population Density (N1)')
plt.ylabel('Normalized Seal Population Density (N2)')
plt.show()




'''
#rk45 predator-prey model
🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭
'''


orca_init = 0.3 #predator initial population density
seal_init = 0.3 #prey initial population density
N_init = [seal_init, orca_init]
num_steps = 100
t_span = (0, num_steps) #time steps
t_eval = np.linspace(0, 50, 100)

def arctic_pred(t, N):
    '''
    the function now includes time since we are not manually deciding time step anymore with euler method
    the function arctic_pred does as follows
    defines constants from Lotka-Volterra Model 
    a constant growth rate of prey
    b constant predator causing death of prey
    c constant consumed prey leading to predator offspring
    d constant natural death rate of predator (starvation)
    defines the species density of species 1 and species 2 as a vector
    calculates the dN1dt and dN2dt, i.e., changes in pop den related to time step width, using L-V predator/prey equations
    returns those values as an array
    Define dNdt of the arctic_pred function as our initial den values
    '''
    a = 1 #constant
    b = 3 #constant
    c = 1 #constant
    d = 3 #constant 
    N1,N2 =N[0],N[1] #N1 is seal, N2 is orca
    dN1dt =a*N1-b*N1*N2 #prey, seal
    dN2dt =-c*N2+d*N1*N2 #predator, orca
    return ([dN1dt,dN2dt]) #report the changes of pop density after a time step with a width size defined by dt

#here is where we enable using the scipy.integrate import solve_ivp. first. we tell the solution that we to solve_ivp
#then we give it the name of our function (mine is arctic), the t_span must be a range of time iterations (num_steps)
#y0 are the initial values of the populations
#the method we are using si RK45
#we can use t_eval to report pop densities at specific time steps
sol = solve_ivp(
    fun=arctic_pred,
    t_span=t_span,
    y0=N_init,
    method='RK45',
    t_eval=t_eval
)

plt.figure(figsize=(8, 5))
plt.plot(sol.t, sol.y[0], label='Seals (N1)', color='orange',  linestyle="dashed")
plt.plot(sol.t, sol.y[1], label='Orcas (N2)', color='blue' )
plt.title('RK45 Arctic Predator/Prey Model Simulation')
plt.xlabel('Time (Years)')
plt.ylabel('Normalized Population Density (Carrying Capacity = 1)')
plt.legend()
plt.grid(True)
plt.show()

#plotting phase diagram of species 1 and species 2 
plt.plot(sol.y[0], sol.y[1])
plt.xlabel('Normalized Seal Population Density (N1)')
plt.ylabel('Normalized Orca Population Density (N2)')
plt.title('RK45 Phase Space of Predator/Prey Relationship Between Seal (N1) and Orca (N2)')
plt.show()



'''
Lorenz Equation Using RK45
'''
#def lorenz system of equations
def lorenz(t, state, sigma, r, b):
    '''
    variables t time, state as x y z, sigma, r, and b. all given in lab guide
    x, y, z are given variables, here referred to as state that can be manipulated at the initial
    dxdt equation given
    dydt equation given
    dzdt equation given
    return the vector of the change of these values with respect to dt
    '''
    x, y, z = state
    dxdt = sigma * (y-x)
    dydt = x * (r-z) - y
    dzdt = x *y - b * z
    return [dxdt, dydt, dzdt]

#model parameters
sigma = 10.0
b = 8.0/3.0
t_span = (0, 100)
t_eval = np.linspace(0,100,10000) #using denser time points here for smoother plotting

#initial conditions
ic_base = [1.0, 1.0, 1.0] #given in lab guide
ic_perturbed = [1.0 + 1e-5, 1.0 + 1e-5, 1.0 + 1e-5] #given in lab guide

def run_experiment(r_val): #we have two r-values we want to test for, r=23 and r=25
    # Solve for unperturbed base initial condition
    sol_base = solve_ivp(
        lorenz, t_span, ic_base, args=(sigma, r_val, b), 
        method='RK45', t_eval=t_eval, rtol=1e-8, atol=1e-10
    )
    
    # Solve for perturbed initial condition
    sol_perturbed = solve_ivp(
        lorenz, t_span, ic_perturbed, args=(sigma, r_val, b), 
        method='RK45', t_eval=t_eval, rtol=1e-8, atol=1e-10
    )
    # Print numerical results at t = 100
    base_final = sol_base.y[:, -1] #we wantt to be able to compare final x y z 
    pert_final = sol_perturbed.y[:, -1]
    
    print(f"=== Results for r = {r_val} at t = 100 ===")
    print(f"Base Trajectory Final (x, y, z):      {base_final}")
    print(f"Perturbed Trajectory Final (x, y, z): {pert_final}")
    
    return sol_base, sol_perturbed
# Run experiments for both r values, r = 23 and r = 25
sol_base_23, sol_pert_23 = run_experiment(r_val=23) #running model for both base xyz and perturbed xyz
sol_base_25, sol_pert_25 = run_experiment(r_val=25)

#Plotting
fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True) #creating two panel figure with scenarios of r=23 and r=25 with both perturbed and base x y z plotted

# Plot for r = 23, labeling
axes[0].plot(sol_base_23.t, sol_base_23.y[0], label="Base (1, 1, 1)", color="blue")
axes[0].plot(sol_pert_23.t, sol_pert_23.y[0], label="Perturbed (1+1e-5, ...)", color="orange", linestyle="--")
axes[0].set_title("r = 23 (Stable Regime: Trajectories Converge)")
axes[0].set_ylabel("x(t)")
axes[0].legend()
axes[0].grid(True)

# Plot for r = 25, labeling
axes[1].plot(sol_base_25.t, sol_base_25.y[0], label="Base (1, 1, 1)", color="blue")
axes[1].plot(sol_pert_25.t, sol_pert_25.y[0], label="Perturbed (1+1e-5, ...)", color="red", linestyle="--")
axes[1].set_title("r = 25 (Chaotic Regime: Trajectories Diverge Exponentially)")
axes[1].set_xlabel("Time (t)")
axes[1].set_ylabel("x(t)")
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()

# Try changing initial conditions:
ic_list = [
    [1.0, 1.0, 1.0],# Standard baseline
    [-1.0, -1.0, 1.0],# Negative Values
    [0.01, 0.0, 0.0],#at 0
    [50.0, 50.0, 50.0]#super far away from initial state
]

for ic in ic_list:
    sol = solve_ivp(lorenz, (0, 100), ic, args=(10.0, 25.0, 8.0/3.0), method='RK45') #print these values for us to compare how changingin initial x y z changes things for final outcome at t=100
    print(f"IC: {ic}  --> Final State (t=100): x={sol.y[0][-1]:.2f}, y={sol.y[1][-1]:.2f}, z={sol.y[2][-1]:.2f}")
