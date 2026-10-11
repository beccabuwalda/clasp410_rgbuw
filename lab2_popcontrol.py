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
to validate the Euler method competition code is running the correct way, perform the following tests by changing line's 45 and 46:
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
N_init = [orca_init, seal_init] #creating a vector N_init of intial population densities
dt = 1 #time step width in years
num_steps = 100 #number of time steps taken


#here, i am defining my function arctic that is describing the functions of competition dynamics
def arctic_comp(N):
    '''
    the function artic_comp does as follows
    defines constants from Lotka-Volterra Model 
    ==================
    inputs
    ==================
    a constant  reproduction rate of species 1
    b constant scales impacts of Species 2 on Species 1
    c constant reproduction rate of species 2
    d constant scales impact of Species 1 on Species 2
    ==================
    action
    ==================
    defines the species density of species 1 and species 2 as a vector
    calculates the dN1dt and dN2dt, i.e., changes in pop den related to time step width using L-V competition equations, using ODEs
    ==================
    outputs
    ==================
    returns those values as an array
    Define dNdt of the arctic_comp function as our initial den values
    '''
    a = 1 #constant reporudction rate of species 1
    b = 2 #constant scales impact of species 2 on species 1
    c = 1 #constant reproduction rate of species 2
    d = 3 #constant scales impact of species 1 on species 2
    N1,N2 =N[0],N[1] #we have made a vector of our species 1 and species 2 population densities, so it can return as one value as opposed to having to sparate it
    dN1dt = a*N1*(1-N1)-b*N1*N2 #changes of population density of species 1, orca, for competition
    dN2dt = c*N2*(1-N2)-d*N1*N2 #changes of population density of species 2, seal, for competition
    return np.array([dN1dt,dN2dt]) #report the changes of pop density after a time step with a width size defined by dt
dNdt=arctic_comp(N_init) #the initial change in pop densities over first time step in our func artic

#prepare figure, load
plt.figure(figsize=(8,8))

#storing history of data outside of loop
N=np.copy(N_init) #make a copy of N_init to N to use at beginning of loop without overwriting our real initial values
time_history=[0.0] #time bin to graph later each step rather than it being erased
popden_history=[np.copy(N)] #making a copy of our N to store as historical data to plot but not overwrite any of the other code

#Euler method loop of Lokta-Voltera Competition Model
for t in range(num_steps): #for time in the range of the number of steps we have decided
    dNdt=arctic_comp(N) #run the dNdt equation from arctic_comp(N) function
    N_next = N + dt * dNdt #calculate the next pop density value for the next step using ODEs

    N = np.copy(N_next) #make our new N pop den value we just calculates a copy and rename it N to use as our new initial value for our next step

    #append new data into our history
    time_history.append((t+1)*dt) #append our time written as time step times time width
    popden_history.append(np.copy(N)) #append a copy of our N values, aka, the pop den values we are getting at each stepm (N_next) 

#convert the history into an array for easy slicing
popden_history = np.array(popden_history) #turn the popden_history into an array sinze we are dealing with a vector of N1 and N2 values 

#plotting complete trajectory for current t*dt
plt.plot(time_history, popden_history[:,0], label="Orca (N1)", color="blue") #plotting first time history and then calling specifically vector y[0] to look at orcas (N1)
plt.plot(time_history,popden_history[:,1], label="Seal (N2)", color="orange", linestyle='dashed') #plotting second time history and then calling specifically vector y[1] to look at seals (N2)
#color blind friwendly, two colo blind firendly colos, dashed and solid lines

#Plot aesthetics
plt.xlabel('Time (Years)')
plt.ylabel('Normalized Population Density (1 = Carrying Capacity') #include units
plt.title('Population Density: Competition of Orcas and Seals') #include units
plt.legend(loc='upper right')
plt.grid(True)
plt.show()


'''
#hand-coded predatory-prey model: Euler Method
🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭🦭
'''
'''
Debugging code: test for convergence!
we want to make sure that changing out time step width does not change the plot
steps:
plot a series of dt:
1. dt = 0.05 (original)
2. dt = 0.01
3. dt = 0.001
4. dt = 0.0001
At what point does the graph no longer change, i.e, when does the size of dt no longer matter and change the plot?
this = convergence
'''

#set initial conditions
orca_init = 0.3 #predator initial population density
seal_init = 0.6 #prey initial population density
N_init = [seal_init, orca_init] #vector of initial pop densities
dt = .01 #time in years
num_steps = 2000 #number of time steps taken

def arctic_pred(N):
    '''
    the function arctic_pred does as follows
    defines constants from Lotka-Volterra Model
    ==================
    inputs
    ==================
    a constant growth rate of prey
    b constant predator causing death of prey
    c constant consumed prey leading to predator offspring
    d constant natural death rate of predator (starvation)
    defines the species density of species 1 and species 2 as a vector
    ==================
    actions
    ==================
    calculates the dN1dt and dN2dt, i.e., changes in pop den related to time step width, using L-V predator/prey equations
    =================
    outputs
    =================
    returns those values as an array
    Define dNdt of the arctic_pred function as our initial den values
    '''
    a = 1 #constant growth rate of prey
    b = 2 #constant predator killing prey
    c = 1 #constant predator growth rate
    d = 3 #constant predator death rate (starvation)
    N1,N2 =N[0],N[1] #N1 is seal, N2 is orca
    dN1dt =a*N1-b*N1*N2 #prey
    dN2dt =-c*N2+d*N1*N2 #predator
    return np.array([dN1dt,dN2dt]) #report the changes of pop density after a time step with a width size defined by dt
dNdt=arctic_pred(N_init) #the initial change in pop densities over first time step in our func artic_pred
    
#prepare figure, visual
plt.figure(figsize=(8,8))
    
#storing history of data outside of loop
N=np.copy(N_init) #making Our first N a copy of N_init so N_init does not get overwritten in our Euler method loop
time_history=[0.0] #time bin to graph later each step rather than it being erased
popden_history=[np.copy(N)] #making copy of N population den data we are getting to store for later to plot
    
#Euler method loop of Lokta-Voltera Competition Model
for t in range(num_steps): #for time in our number of time steps predetermined
    dNdt=arctic_pred(N) #perform the dNdt function from artic_pred
    N_next = N + dt * dNdt #and calculate the next population densit ung that dNdt calculated and our initial N value at the start of the step
    
    N = np.copy(N_next) #make out new population density N_next our new initial step N for the next loop of time step t in time width dt
    
    #append new data into our history
    time_history.append((t+1)*dt) #save history of time with t+1*dt since python starts at 0 
    popden_history.append(np.copy(N)) #save copy of all pop densities calculated over loop
    
#convert the history into an array for easy slicing
popden_history = np.array(popden_history) #turn vector of pop densities [N1, N2] into an easier to read array for plotting
    
#plotting complete trajectory for current dt
plt.plot(time_history, popden_history[:,0], label="Seal (N1)", color="orange", linestyle="dashed") #first plot species 1 aka seals calling y[0]
plt.plot(time_history,popden_history[:,1], label="Orca (N2)", color="blue") #then plot species 2 aka orcas by calling y[1]
#color blind friwendly, two colo blind firendly colos, dashed and solid lines

#Plot aesthetics
plt.xlabel('Time (Years)') #include units
plt.ylabel('Normalized Population Density (1 = Carrying Capacity)') #include units
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
#defining initial conditions
orca_init = 0.3 #species 1 initial population density
seal_init = 0.6 #species 2 initial population density
N_init = [orca_init, seal_init]
num_steps = 100
t_span = (0, num_steps) #time step width in years
t_eval = np.linspace(0, 50, 100) #explained in more context in RK45 function below

def arctic_comp(t, N): #now defining function of t and pop den, since we are not manually deciding dt
    '''
    function arctic_comp now using RK45 method
    time is now added since we arent manually deciding dt anymore
    the function artic_comp does as follows:
    ===========
    inputs
    ===========
    defines constants from Lotka-Volterra Model 
    a constant  reproduction rate of species 1
    b constant scales impacts of Species 2 on Species 1
    c constant reproduction rate of species 2
    d constant scales impact of Species 1 on Species 2
    ===========
    actions
    ===========
    defines the species density of species 1 and species 2 as a vector
    calculates the dN1dt and dN2dt, i.e., changes in pop den related to time step width using L-V competition equations
    ==========
    outputs
    ==========
    returns those values as a vector
    '''
    a = 1 #constant reproduction rate of species 1
    b = 2 #constant scales impacts of species 2 on species 1
    c = 1 #constant reproduction rate of species 2
    d = 3 #constant scales impact of species 1 on species 2
    N1,N2 =N[0],N[1] #orca species 1, seal species 2
    dN1dt = a*N1*(1-N1)-b*N1*N2 #changes of population density of species 1, orca, for competition
    dN2dt = c*N2*(1-N2)-d*N1*N2 #changes of population density of species 2, seal, for competition
    return ([dN1dt,dN2dt]) #report the changes of pop density after a time step with a width size defined by dt

sol = solve_ivp( #our new function sol that includes artic_comp conditions
    fun=arctic_comp, #function we are calling, arctic_comp
    t_span=t_span, #time span we are wanting to go over
    y0=N_init, #initial pop density vector
    method='RK45', #package using
    t_eval=t_eval #eval values at a specific time step size rather than being adaptive from RK45, here we are designating specific equally spaced intervals
)

plt.figure(figsize=(8, 5))
plt.plot(sol.t, sol.y[0], label='Orcas (N1)', color='blue') #plot call to sol function y[0], orca
plt.plot(sol.t, sol.y[1], label='Seals (N2)', color='orange',  linestyle="dashed") #plot call to sol function y[1], seal
#used color blind firendly colors, and dashed vs solid lines
plt.title('RK45 Arctic Competition Model Simulation')
plt.xlabel('Time (Years)') #include units
plt.ylabel('Normalized Population Density (Carrying Capacity = 1)') #include units
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
#initial conditions 
orca_init = 0.6 #predator initial population density
seal_init = 0.3 #prey initial population density
N_init = [seal_init, orca_init] #vector of initial pop densities
num_steps = 100 #num of time steps we are taking
t_span = (0, num_steps) #time steps in a range for RK45 because it needs a start and end point
t_eval = np.linspace(0, 50, 100) #explained further below in function sol

def arctic_pred(t, N):
    '''
    ***the function now includes time since we are not manually deciding time step anymore with euler method***
    the function arctic_pred does as follows
    defines constants from Lotka-Volterra Model 
    =====================
    inputs
    ====================
    a constant growth rate of prey
    b constant predator causing death of prey
    c constant consumed prey leading to predator offspring
    d constant natural death rate of predator (starvation)
    ======================
    actions
    ======================
    defines the species density of species 1 and species 2 as a vector
    calculates the dN1dt and dN2dt, i.e., changes in pop den related to time step width, using L-V predator/prey equations
    =====================
    outputs
    =====================
    returns those values as a vector
    '''
    a = 1 #constant growth rate of prey
    b = 2 #constant death of prey
    c = 1 #constant predator growth rate
    d = 3 #constant death of predator (starvation)
    N1,N2 =N[0],N[1] #N1 is seal, N2 is orca
    dN1dt =a*N1-b*N1*N2 #prey, seal
    dN2dt =-c*N2+d*N1*N2 #predator, orca
    return ([dN1dt,dN2dt]) #report the changes of pop density after a time step with a width size defined by dt

#here is where we enable using the scipy.integrate import solve_ivp. first. we tell the solution that we to solve_ivp
#then we give it the name of our function (mine is arctic), the t_span must be a range of time iterations (num_steps)
#y0 are the initial values of the populations
#the method we are using si RK45
#we can use t_eval to report pop densities at specific time steps
sol = solve_ivp( #function "sol" using scipi solv_ivp
    fun=arctic_pred, #refer to artic_pred that is just newly defined above to now include time
    t_span=t_span, #repeat across range of time steps 0 to 100
    y0=N_init, #start at initial vector of N pop den values
    method='RK45', #method RK45 used
    t_eval=t_eval #evaluate at a specific even time step size (as defined above) as opposed to being a changing time step size. This is jsut a preference saying I dont want the time step to change on its own under the hood.
)

#plotting
plt.figure(figsize=(8, 5))
plt.plot(sol.t, sol.y[0], label='Seals (N1)', color='orange',  linestyle="dashed") #plotting our prey first by calling the sol function y[0]
plt.plot(sol.t, sol.y[1], label='Orcas (N2)', color='blue' ) #calling the predators by calling the sol function y[1]
plt.title('RK45 Arctic Predator/Prey Model Simulation')
plt.xlabel('Time (Years)') #include units
plt.ylabel('Normalized Population Density (Carrying Capacity = 1)') #include units
plt.legend()
plt.grid(True)
plt.show()

#plotting phase diagram of species 1 and species 2 
plt.plot(sol.y[0], sol.y[1])
plt.xlabel('Normalized Seal Population Density (N1)') #include units
plt.ylabel('Normalized Orca Population Density (N2)') #include units
plt.title('RK45 Phase Space of Predator/Prey Relationship Between Seal (N1) and Orca (N2)')
plt.show()



'''
Lorenz Equation Using RK45
'''

#def lorenz system of equations
def lorenz(t, state, sigma, r, b):
    '''
    =====
    inputs
    =====
    time t
    constant sigma
    constant r (manipulated in this lab)
    constant b
    x, y, z are given variables, here referred together as state that can be manipulated at the initial in cohesion.
    ======
    actions
    ======
    dxdt equation given
    dydt equation given
    dzdt equation given
    ======
    outputs
    ======
    return the vector of the change of these values from the equations above [dxdt, dydt, dzdt]
    '''
    x, y, z = state
    dxdt = sigma * (y-x)
    dydt = x * (r-z) - y
    dzdt = x * y - b * z
    return [dxdt, dydt, dzdt]

#model parameters
sigma = 10.0 #constant given in lab guide
b = 8.0/3.0 #constant gien in lab guide
t_span = (0, 100) #number of time steps
t_eval = np.linspace(0,100,10000) #using denser time points here for smoother plotting, but can manipulate later to see how this changes things

#initial conditions
ic_base = [1.0, 1.0, 1.0] #given in lab guide as initial values of [x, y, z]
ic_perturbed = [1.0 + 1e-5, 1.0 + 1e-5, 1.0 + 1e-5] #given in lab guide as perturbed values of [x, y, z]

def run_experiment(r_val): #we have two r-values we want to test for, r=23 and r=25
    '''
    ========
    inputs
    ========
    an open based experiment where the r-value can be changes to anything without having to redo the whole function.
    sol_base, where we are using RK45 model of unperturbed x,y,z values, sigma, r_val, and b
    sol_perturbed where we are using RK45 model of perturbed x,y,z values, sigma, r_val, and b
    =========
    actions
    =========
    solves for RK45 initial conditions of: 
    r_val (defined later can be flexibly decided)
    base xyz values previously defined initial condition
    perturbed xyz values previously defined initial condition
    ==========
    outputs
    ==========
    numerical results of base and perturbed at t = 100
    RK45 resolved of base and perturbed 
    '''
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
    base_final = sol_base.y[:, -1] #we want to be able to compare final x y z 
    pert_final = sol_perturbed.y[:, -1]
    
    print(f"Results for r = {r_val} at t = 100")
    print(f"Base Trajectory Final (x, y, z):      {base_final}")
    print(f"Perturbed Trajectory Final (x, y, z): {pert_final}")
    
    return sol_base, sol_perturbed
###DEFINE WHAT R VALUES TO RUN IN FUNCTION ABOVE!!!
# Run experiments for both r values, r = 23 and r = 25
sol_base_23, sol_pert_23 = run_experiment(r_val=23) #running model for both base xyz and perturbed xyz r=23
sol_base_25, sol_pert_25 = run_experiment(r_val=25) #running model for both base xyz and perturbed xyz r=25

#Plotting
fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True) #creating two panel figure with scenarios of r=23 and r=25 with both perturbed and base x y z plotted

# Plot for r = 23, labeling
axes[0].plot(sol_base_23.t, sol_base_23.y[0], label="Base (1, 1, 1)", color="blue") #initial conditions
axes[0].plot(sol_pert_23.t, sol_pert_23.y[0], label="Perturbed (1+1e-5, ...)", color="orange", linestyle="--") #perturbed conditions, color blind friendly
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

# Try changing initial conditions: what happens when we change x, y, z?
ic_list = [
    [1.0, 1.0, 1.0],# Standard baseline
    [-1.0, -1.0, 1.0],# Negative Values
    [0.01, 0.0, 0.0],#at 0
    [50.0, 50.0, 50.0]#super far away from initial state of 1, 1, 1
]
#used AI here to help me print this in a nice lined up table for comparison purposes
for ic in ic_list:
    sol = solve_ivp(lorenz, (0, 100), ic, args=(10.0, 25.0, 8.0/3.0), method='RK45') #print these values for us to compare how changingin initial x y z changes things for final outcome at t=100 and r=25
    print(f"IC: {ic}  --> Final State (t=100): x={sol.y[0][-1]:.2f}, y={sol.y[1][-1]:.2f}, z={sol.y[2][-1]:.2f}")
