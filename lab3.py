'''
This code was used for the Energy Balance Atmosphere Mode lab
Follow this code to learn the basics of steady state of Earth
And learn to love linear algebra
Goals:
1. gain experience in developing functions to solve
systems of linear equations of the form Ax = b
2. employ our function to create and solve sets of equations representing energy balance
in an N-layer atmosphere
3. use the N-layer atmosphere model to make predictions about atm characteristics
on Earth and at other planets
'''

#importing necessary packages
import numpy as np
import matplotlib
matplotlib.use('TkAgg') #setting up back end
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import ode 
from scipy.integrate import solve_ivp


'''
---------------------------------------------------------------------------------------------------------------------------
1. BY HAND 3 LAYER ATM MODEL
First, let's build a 3-layer atmosphere model
By-hand, I solved the Flux for the 3-layers of atm and the ground, inputting values here in [A][x]=[b] form
--------------------------------------------------------------------------------------------------------------------------
'''
'''
#constants that remain unchanged despite changing emissivity and num_layers
alpha= 0.33 #alpha constant (albedo) (changes also when we go off of earth later)
S = 1350 #S constant (Solar Flux) for EARTH (changes when we go off of earth later)
sigma = 5.6703*10**-8 #sigma constant (Boltzman constant)
num_layers = 3 #3 atmosphere layers
e = 1 #defining epsilon as 1 makes it a black body

c = [-1, -2, -2, -2] #defining an array that will then become the center diag of my matrix
k1 = [1, 1, 1] #defining diag c+1 
k2 = [1, 1, 1] #defining diag c-1


#matrix number 1 is defining the center diagonal conditions
A = np.diag(c, k = 0) #defining the matrix A with the diag of array c now defined, so central diag is now not 0 at indexing diag k = 0 (center diag)
print(A)

B = np.diag(k1, k = 1) #defining the matrix B with the off diag from center + 1 from array k1
print(B)

C = np.diag(k2, k = -1) #defining the matrix C with the off diag from center -1 from array k2
print(C)

#defining our combined 3-atm layer matrix by combining 
blackbody_3layer = A + B + C
print(blackbody_3layer)

#defining our b matrix that will be multiplied with our blackbody_3layer matrix to determined fluxes
b=np.zeros((4,1))
print(b)
b[0,0]=-(((1-alpha)/4)*S)
b[1,0]=0
b[2,0]=0
b[3,0]=0

#multiplying A and b matrices
solve_3layerblackbody = np.linalg.solve(blackbody_3layer, b)
print(solve_3layerblackbody)

#divide by sigma and take the fourth root to get temperature values in each layer
temps_3layerblackbody = (solve_3layerblackbody/sigma)**(1/4) #temp is in Kelvin
print(temps_3layerblackbody)
#check to make sure temperatures reflect a 3-layer blackbody model 
'''
'''
---------------------------------------------------------------------------------------------------------------------------------------------------------
2. BUILDINING INFINITE ATM LAYER BLACKBODY
building an infinite layer matrix; BLACKBODY version; testing with 0, 1, 2, 3 layers of atmosphere!!
blackbody results here will be used to codecheck completely open model at next step (where emissivity and num_layers both need to be defined)
In blackbody, I only need to leave num_layers open-ended, no manipulating emissivity
---------------------------------------------------------------------------------------------------------------------------------------------------------
'''
'''
def blackbody_NLayermatrix(num_layers = 1, e = 1, alpha = 0.3, S = 1350, sigma = 5.6703e-8):
    '''
'''
    Computes the temperature profile for an N-layer atmospheric blackbody model
    Parameters:
    - Num_layers (int): numbers of atmospheric layers (N), creates a system of size (N+1)x(N+1)
    - e(float): Emissivity (remains as 1 for blackbody)
    - alpha (float): albedo (remains same)
    - sigma (float): Boltzmann constant (remains same)
    Returns:
    - temp (ndarray): temperatures from surface (0) to layer N (in Kelvin)
    - matrix (ndarray): the coefficient matrix used to solve the system
    '''
'''
    #system dimensions (N+1) x (N+1)
    size = num_layers + 1

    #Main diagonal: -1 at ground surface, -2 for all atmospheric layers
    c = [-1.0] + [-2.0]*num_layers

    #Upper and lower diagonals filled with 1s
    k1 = np.ones(num_layers)
    k2 = np.ones(num_layers)

    #constructing the full matrix with made diagonals
    A = np.diag(c, k=0) + np.diag(k1, k=1) + np.diag(k2, k=-1)

    #constructing vector b
    b=np.zeros((size,1))
    b[0,0] = -(((1 - alpha) / 4) * S)

    #solve for sigma * T^4 fluxes
    solve_fluxes_blackbody_NLayers = np.linalg.solve(A, b)

    #converting fluxes to temperatures
    blackbody_Nlayer_temps = (solve_fluxes_blackbody_NLayers / sigma) ** (1/4)
    
    return blackbody_Nlayer_temps, A matrix

#dictionary to store outputs of each N layer model in
results = {}

for N in range(4): #loops through N layers of 0, 1, 2, 3 atmospheric layers
    temps, matrix = blackbody_NLayermatrix(num_layers=N) #calls the function for the current layer count N
    #we flatten our temperature array from column vector (N+1, 1) into a 1D array (N+1) and store both temp and matrix in dictionary
    results[N] = {
        "temps": temps.flatten(), #allows us to index temps[0] instead of temps[0,0] by flattening to 1D
        "matrix": matrix
    }
    print(f"Blackbody {N}-Layer Model")
    print(f"Surface Temp: {temps[0, 0]:.2f} K") #saving our special scenario at ground surface level
    if N > 0: # if any atmopsheric layers exist
        for i in range(1, N+1): #for in these additional layers in range N+1
           print(f"Layer {i} Temp: {temps[i,0]:.2f} K") #print their temps
    
    #print resulting coefficient matrix A
    print("\nMatrix A:")
    print(matrix)

    #create space for new block
    print("-" * 35 + "\n")

'''

'''
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
3. BUILDING OPEN LAYER ATM AND OPEN EMISSIVITY MODEL, INCLUDING NUCLEAR WINTER OPTION
Building a infinite layer matrix; able to manipulate emissivity and num_layers of atmosphere
CODE CHECK:
1. Make a no atm situation of the infinite layers model with epsilon = 0, num_layers = 0 -> ensure output matches blackbody no atm scenario above
2. Make a 1 Layer, 2 layer, and 3 layer atm layer models of the infinite layers model with epsilon = 0, num_layers = 1, 2, 3 -> ensure output matches blackbody scenarios above for 1, 2, 3 layer atm models
MAKE SURE THESE VALUES CHECK OUT BEFORE MAKING N LAYERS LAYERS LARGER!
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
#DO NOT COMMENT OUT THIS FUNCTION!!!!!!!!!!!!!!!!!!!
def open_NLayermatrix(num_layers = 1, e = .33, alpha = 0.3, S = 1350, sigma = 5.6703e-8, nuclear_winter = False):
    '''
    Computes the temperature profile for an N-layer atmospheric blackbody model
    
    Parameters:
    - Num_layers (int): numbers of atmospheric layers (N)
    - e(float): Emissivity (**CAN NOW CHANGE!!**)
    - alpha (float): albedo (remains same)
    - sigma (float): Boltzman constant (remains same)
    -nuclear winter: used later as a part of the homework questions (right now set to false)
    Returns:
    - temp (ndarray): temperatures from surface (0) to layer N (in Kelvin)
    - matrix (ndarray): the coefficient matrix used to solve the system, absorption terms scaled by attenuation factor (1-e)^(distance-1)
    '''
    #system dimensions (N+1) x (N+1)
    size = num_layers + 1

    #initialize empy square matrix with zeros
    A = np.zeros((size, size))

    #now I populate the matrix elements based on their spatial interaction - patterns on the diagonal
    for i in range(size):
        for j in range(size):
            if i == j:
                # Main Diagonal: Surface emits -1, Atmospheric layers emit -2
                A[i, j] = -1 if i == 0 else -2.0
            elif i == 0 and j > 0:
                #Row 0 (Surface) recieves down radiation from layer j through layers
                A[0, j] = 1.0 if j == 1 else ((1-e)**(j-1) if e < 1.0 else 0.0)
            elif j == i + 1 or j == i - 1:
                #immediate adjacent neighbors (1 layer next to it) of any given layer emit e
                A[i, j] = e
            elif j > i + 1: #but if we are further away we need to consider more fractioning of emissivity
                #upper triangle of values: layer i received down radiation (distance greater than 1 layer)
                A[i, j] = e * ((1.0-e)**(j-i-1)) if e<1.0 else 0.0 #when we are cosnidering blackbody this will become negative and mess things up so we just tell it to go to 0
            elif j < i - 1: #moving further away from a layer of interest than 1 step now in the opposite direction to consider more fractioning of emissivity
                #lower triangle: flux from layer j up to layer i (distance greater than 1 layer)
                A[i, j] = e*((1.0-e)**(i-j-1)) if e<1.0 else 0.0 
    
    #constructing vector b: right hand side shortwave solar input
    b=np.zeros((size,1))
    if nuclear_winter: #nuclear winter scenarios of shortwave input: inverted atm where all heating is being abosrbed at top of atm
        #top atm layer (index N) absorbs all incoming solar flux
        #surface (index 0) receives zero direct solar shortwave radiation
        b[num_layers, 0] = -(((-1.0 - alpha) / 4.0)*S)
    else:
        b[0,0] = -(((1 - alpha) / 4) * S) #normal scenarios, energy going to surface layer

    #solve for sigma * T^4 fluxes: take abs value 1st to resolve any NaN numbers
    solve_fluxes_open_NLayers = np.linalg.solve(A, b)

    #converting fluxes to temperatures
    open_Nlayer_temps = (np.abs(solve_fluxes_open_NLayers / sigma) ** (1/4))    
    return open_Nlayer_temps.flatten(), A

#bin to store outputs of each N layer model in
results = {}

'''
##RUN THIS SECTION TO VERIFY THE OPEN EMISSIVITY MODEL RUNS PROPERLY##
#first, lets code check by looping through blackbody versions of 0, 1, 2, and 3 atm layers and comparing to our blackbody model
for N in range(4): #loops through N layers of 0, 1, 2, 3 atmospheric layers
    temps, matrix = open_NLayermatrix(num_layers=N, e=1) #calls the function for the current layer count N
    results[N] = {
        "temps": temps,
        "matrix": matrix
    }
    print(f"Non-Blackbody {N}-Layer Model | Code Check: e = 1")
    print(f"Surface Temp: {temps[0]:.2f} K") #saving our special scenario at ground surface level

    if N > 0: # if any atmopsheric layers exist
        for i in range(1, N+1): #for in these additional layers in range N+1
           print(f"Layer {i} Temp: {temps[i]:.2f} K") #print their temps
    
    #print resulting coefficient matrix A
    print("\nMatrix A:")
    print(matrix)

    #create space for new block
    print("-" * 35 + "\n")


#now that we have verified that we get the same values as our blackbody version when e=1, lets change e and n layers and see what happens
e_values =[0.2, 0.5, 0.77, 0.95]
for e in e_values:
    for N in range(4): #loops through N layers of 0, 1, 2, 3 atmospheric layers
        temps, matrix = open_NLayermatrix(num_layers=N, e=e) #calls the function for the current layer count N
        results[(N, e)] = {
        "temps": temps,
        "matrix": matrix
    }
    print(f"{N}-Layer Model with Emissivity = {e}")
    print(f"Surface Temp: {temps[0]:.2f} K") #saving our special scenario at ground surface level

    if N > 0: # if any atmopsheric layers exist
        for i in range(1, N+1): #for in these additional layers in range N+1
           print(f"Layer {i} Temp: {temps[i]:.2f} K") #print their temps
    
    #print resulting coefficient matrix A
    print("\nMatrix A:")
    print(matrix)

    #create space for new block
    print("-" * 35 + "\n")
'''


'''
-----------------------------------------------------------------------------
EXPERIMENT 1: SINGLE LAYER ATM - SURFACE TEMP VS EMISSIVITY @ 288 K ON EARTH
-----------------------------------------------------------------------------
'''

#Finding emissivity when the surface temperature of earth = 288 K (N=1)
e_range = np.linspace(0,1,100) #testing even spaced range of values of e
surf_temps_exp1=[] #empty bin to keep temps in for experimenting dif emissivities

for e in e_range:
    temps, _ = open_NLayermatrix(num_layers=1, e=e, S=1350) #since we only need to experiment with temp, we throw in dummy variable _ into our tuple (temps, N)
    surf_temps_exp1.append(temps[0]) #append what we find for temps here 

#finding predicted emissivity for surface Temperature = 288 K when N =1
#when N=1, T surface = (1-alpha)*S/(4*sigma*(1-e/2))
# e = 2*(1-T_bare^4/T_surf^4)
T_bare = (((1 - alpha) * S)/(4 * sigma))**(1/4) #temp with no layers, bare earth
target_temp = 288.0 #target temp of earth with 1 layer of emissivity e
e_pred_288 = 2 * (1.0 -(T_bare/target_temp)**4) #what the emissivity needs to be to get to target temp, rearranging 1-layer atm model eq. to solve for emissivity

print("----Experiment 1 Results----")
print(f"Bare-Earth Surface Temperature (e=0: {T_bare:.2f} K")
print(f"Predicted Emissivity (e) for 288 K Earth Surface Temp (N = 1): {e_pred_288:.4f}\n")

#plotting experiment 1 results
plt.figure(figsize=(7,5))
plt.plot(e_range, surf_temps_exp1, color='navy', lw=2, label='1-Layer Model') #plot iterations of emissivity and temperatures made
plt.axhline(288, color='red', linestyle='--', label='Desired Earth Temp (288 K)') #reference line at desired temperature
plt.axvline(e_pred_288, color='green', linestyle=':', label=f'e ≈ {e_pred_288:.3f}') #reference line where predicted emissivity will fall at 288 K
plt.title("Experiment 1: Surface Temperature vs Emissivity (1-Layer Atmosphere) for Earth")
plt.ylabel("Earth Surface Temperature (K)")
plt.xlabel("Emissivity (e)")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend() #legend, different line types, colors to be color blind friendly
plt.tight_layout()
plt.show()
#this plot shows a line of emissivities and temp for 1-layer model to compare to the vertical-horzontal line intersection

'''
------------------------------------------------------------------------------------
EXPERIMENT 2: VARYING N-LAYERS WITH FIXED EMISSIVITY (e = 0.255)
------------------------------------------------------------------------------------
'''
'''
e_fixed = 0.255 #new emissivity given
max_N_layer = 15 #create a random number of N layers to test through to find our desired number of N layers (I doubt it will take more than 15)
N_array = np.arange(0,max_N_layer+1) #np.arange is a function used to create a 1D array containing evenly spaced vals, memory optimized array ready for vector ops
surf_temps_exp2 = [] #storing

for N in N_array: #for each N in our range of possible N-layer models
    temps, _ = open_NLayermatrix(num_layers=N, e=e_fixed, S=1350) #we want to look at temps only, put in dummy _ for tuple (temps, N)
    surf_temps_exp2.append(temps[0]) #then save these temps for each layer for each model scenario

#find number of layers required to reach ~288 K for surface layer if e is set at 0.255
N_target = np.argmin(np.abs(np.array(surf_temps_exp2)-288.0)) #find the minimum number of layers it takes for calc surface temp from adding layers to equal 288 (over or under shooting possible with abs value)

print("----Experiment 2 Results----")
print(f"Target Emissivity: e = {e_fixed}")
for N, T in zip(N_array, surf_temps_exp2): #when looking at temps produced from each N layer model
    print(f"Layers (N = {N: 2d}): Surface Temp = {T:.2f} K")
print(f"Required layers for ~280 K Earth Surface: N = {N_target} layers and Surface Temp = {surf_temps_exp2[N_target]:.2f} K\n")

#making an altitude plot for N_target layers for Earth
surf_exp2_target, _ = open_NLayermatrix(num_layers=N_target, e=e_fixed, S=1350)
altitudes_exp2 = np.arange(N_target + 1) # Index 0 = ground (0 km) 1..N = atmospheric layers

plt.figure(figsize=(6,6))
plt.plot(surf_exp2_target, altitudes_exp2, 'o-', color='darkgreen', lw=2, ms=7)
plt.title(f"Experiment 2: Earth Altitude Profile: N = {N_target} Layers, e = {e_fixed}")
plt.xlabel("Temperature (K)")
plt.ylabel("Atmospheric Level / Layer (0 = Surface)")
plt.yticks(altitudes_exp2, labels=['Surface'] + [f'Layer {i}' for i in range(1, N_target+1)])
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

####for comparison, lets plot N-layers 0-3 with same emissivity###
results = {}
print("----Comparison of N = 0 to 14 Amtospheric Layers (e = 0.255)----")
for N in range(15): #loops through N layers of 0, 1, 2, 3 atmospheric layers
    temps, matrix = open_NLayermatrix(num_layers=N, e=e_fixed) #calls the function for the current layer count N, set e to e_fixed for comparison
    results[N] = {
        "temps": temps,
        "matrix": matrix
    }
    print(f"Non-Blackbody {N}-Layer Model | Comaprison to 4-layer atm with e = 0.255")
    print(f"Surface Temp: {temps[0]:.2f} K") #saving our special scenario at ground surface level

    if N > 0: # if any atmopsheric layers exist
        for i in range(1, N+1): #for in these additional layers in range N+1
           print(f"Layer {i} Temp: {temps[i]:.2f} K") #print their temps
    
    #print resulting coefficient matrix A
    print("\nMatrix A:")
    print(matrix)

    #create space for new block
    print("-" * 35 + "\n")

#altitude plotting
plt.figure(figsize=(6,6))

color = [ #asked AI to give me 15 different colors for 15 dif layers
    '#1f77b4',  # Muted Blue
    '#ff7f0e',  # Safety Orange
    '#2ca02c',  # Cooked Asparagus Green
    '#d62728',  # Brick Red
    '#9467bd',  # Muted Purple
    '#8c564b',  # Chestnut Brown
    '#e377c2',  # Raspberry Pink
    '#7f7f7f',  # Middle Gray
    '#bcbd22',  # Curry Yellow-Green
    '#17becf',  # Ice Blue / Cyan
    '#3366cc',  # Royal Blue
    '#dc3912',  # Crimson Red
    '#ff9900',  # Amber
    '#109618',  # Forest Green
    '#990099',  # Deep Magenta
]
marker = ['o', 's', '^', 'D', 'v']

for N in range(15):
    temps = results[N]["temps"]
    altitudes = np.arange(N+1) # 0 = surface, 1...N = atmsophere
    plt.plot(temps, altitudes, marker=marker[N % 5], color=color[N], lw=2, ms=7, label=f"{N}-Layer Model")

plt.title(f"Earth Altitude Profiles: N = 0-14 Layers (e = {e_fixed})", fontsize = 12)
plt.xlabel("Temperature (K)")
plt.ylabel("Atmsopheric Level / Layer (0 = Surface)")
plt.yticks(np.arange(4), labels=['Surface'] + [f'Layer {i}' for i in range(1,4)])
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper right', frameon=True)
plt.tight_layout()
plt.show()
'''
'''
-------------------------------------------------------------------------
EXPERIMENT 3: VENUS ATMOSPHERIC MODEL, CHANGED S0, EMISSIVITY, T SURFACE
-------------------------------------------------------------------------
'''
'''
S_venus = 2600 #new solar flux given
alpha_venus = 0.77 #looked up the albedo for venus 
target_temp_venus = 700 #given temp for surface of venus
venus_N = 1 #step layer
e_venus = 1 #emissivity on venus
max_N_layer_venus = 100 #max layers on venus we want to test (making big number this time)
N_array_venus = np.arange(0,max_N_layer_venus+1) #np.arange is a function used to create a 1D array containing evenly spaced vals, memory optimized array ready for vector ops
surf_temps_exp3 = []

for N in N_array_venus:
    temps, _ = open_NLayermatrix(num_layers=N, e=e_venus, alpha=alpha_venus, S=S_venus) #find temps for each layer in each scenario model with given parameters with dummy in place of N
    surf_temps_exp3.append(temps[0]) #save those temps

#find number of layers required to reach ~288 K for surface layer if e is set at 0.255
N_target_venus = np.argmin(np.abs(np.array(surf_temps_exp3)-700.0)) #find the minimum number of layers it takes for calc surface temp from adding layers to equal 288 (over or under shooting possible with abs value)

print("----Experiment 3 Results----")
print(f"Target Emissivity: e = {e_venus}")
for N, T in zip(N_array_venus, surf_temps_exp3): #when looking at temps produced from each N layer model
    print(f"Layers (N = {N: 2d}): Surface Temp = {T:.2f} K")
print(f"Required layers for ~700 K Venus Surface: N = {N_target_venus} layers and Surface Temp = {surf_temps_exp3[N_target_venus]:.2f} K\n")

#making an altitude plot for N_target layers for Venus
surf_exp3_target, _ = open_NLayermatrix(num_layers=N_target_venus, e=e_venus, S=S_venus) 
altitudes_exp3 = np.arange(N_target_venus + 1) # Index 0 = ground (0 km) 1..N = atmospheric layers

plt.figure(figsize=(6,6))
plt.plot(surf_exp3_target, altitudes_exp3, 'o-', color='brown', lw=2, ms=7)
plt.title(f"Experiment 3: Venus Altitude Profile: N = {N_target_venus} Layers, e = {e_venus}")
plt.xlabel("Temperature (K)")
plt.ylabel("Atmospheric Level / Layer (0 = Surface)")
custom_labels = [""] *len(altitudes_exp3)
custom_labels[0] = 'Surface'
for i in range(10, len(altitudes_exp3),10):
    custom_labels[i] = f'Layer {i}'
plt.yticks(altitudes_exp3, labels=custom_labels)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
'''

'''
---------------------------------------------------------------------
EXPERIMENT 4: NUCLEAR WINTER SCENARIO, 5 LAYERS, e = 0.5, S0 = 1350
---------------------------------------------------------------------
'''
'''
temps_nw, A_nw = open_NLayermatrix(num_layers=5, e=0.5, alpha=0.3, S=1350, nuclear_winter=True) #changed initial conditions as per lab guide
print(f"Parameters: 5 Atmospheric Layers, e = 0.5, S = 1350 (Top-layer solar absorption)")
print(f"Surface Temperature:{temps_nw[0]:.2f} K")
for i in range(1,6): #for layer print the temp
    print(f"Layer {i} Temperature: {temps_nw[i]:.2f} K")

#Plotting altitude profile of Nuclear Winter 
altitudes_nw = np.arange(6) # Index 0 = ground (0 km) 1..N = atmospheric layers

plt.figure(figsize=(6,6))
plt.plot(temps_nw, altitudes_nw, 's--', color='crimson', lw=2, ms=7, label='Nuclear Winter Profile')
plt.title(f"Experiment 4: Experiment 4: Nuclear Winter Temperature Profile\n(5 Layers, e = 0.5, Top Absorption)")
plt.xlabel("Temperature (K)")
plt.ylabel("Atmospheric Level / Layer (0 = Surface)")
plt.yticks(altitudes_nw, labels=['Surface'] + [f'Layer {i}' for i in range(1, 6)])
plt.grid(True, linestyle='--', alpha=0.6)
plt.axvline(273.15, color='blue', linestyle=':', label='Freezing Point (273.15 K)')
plt.legend()
plt.tight_layout()
plt.show()
'''