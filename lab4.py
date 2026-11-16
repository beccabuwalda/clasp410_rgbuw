'''
Snowball Earth Lab
Try not to cry (update: you did cry)
'''
import numpy as np

'''
Helpful functions for the Snowball Earth lab. Initially Given
'''

def temp_warm(lats_in):
    '''
    Create a temperature profile for modern day "warm" earth.
    Parameters
    ----------
    lats_in : Numpy array
    Array of latitudes in degrees where temperature is required.
    0 corresponds to the south pole, 180 to the north.
    Returns
    -------
    temp : Numpy array
    Temperature in Kelvin.
    '''
    # Set initial temperature curve
    T = np.array([-47, -19, -11, 1, 9, 14, 19, 23, 25, 25,
    23, 19, 14, 9, 1, -11, -19, -47])
    T_warm = T + 273.15 #converting to kelvin so this doesn't mess up my pipeline later with solar radiation
    # Get base grid:
    npoints = T_warm.size
    dlat = 180 / npoints # Latitude spacing.
    lats = np.linspace(dlat/2., 180-dlat/2., npoints) # Lat cell centers.
    # Fit a parabola to the above values
    coeffs = np.polyfit(lats, T_warm, 2)
    # Now, return fitting sampled at "lats".
    temp = coeffs[2] + coeffs[1]*lats_in + coeffs[0] * lats_in**2
    return temp



def insolation(S0, lats):
    '''
    Given a solar constant (`S0`), calculate average annual, longitude-averaged
    insolation values as a function of latitude.
    Insolation is returned at position `lats` in units of W/m^2.
    Parameters
    ----------
    S0 : float
    Solar constant (1370 for typical Earth conditions.)
    lats : Numpy array
    Latitudes to output insolation. Following the grid standards set in
    the diffusion program, polar angle is defined from the south pole.
    In other words, 0 is the south pole, 180 the north.
    Returns
    -------
    insolation : numpy array
    Insolation returned over the input latitudes.
    '''
    # Constants:
    max_tilt = 23.5 # tilt of earth in degrees
    # Create an array to hold insolation:
    insolation = np.zeros(lats.size)
    # Daily rotation of earth reduces solar constant by distributing the sun
    # energy all along a zonal band
    dlong = 0.01 # Use 1/100 of a degree in summing over latitudes
    angle = np.cos(np.pi/180. * np.arange(0, 360, dlong))
    angle[angle < 0] = 0
    total_solar = S0 * angle.sum()
    S0_avg = total_solar / (360/dlong)
    # Accumulate normalized insolation through a year.
    # Start with the spin axis tilt for every day in 1 year:
    tilt = [max_tilt * np.cos(2.0*np.pi*day/365) for day in range(365)]
    # Apply to each latitude zone:
    for i, lat in enumerate(lats):
    # Get solar zenith; do not let it go past 180. Convert to latitude.
        zen = lat - 90. + tilt
        zen[zen > 90] = 90 #cap zenith angle at horizon
        # Use zenith angle to calculate insolation as function of latitude.
        insolation[i] = S0_avg * np.sum(np.cos(np.pi/180. * zen)) / 365.
    # Average over entire year; multiply by S0 amplitude:
    insolation = S0_avg * insolation / 365
    
    return insolation

'''
----------------------------------------------------------------------------------------
PART 1: BUILDING THE SOLVER: CONSIDERS DIFFUSION, SPHERICAL SHAPE, & RADIATIVE FORCING
----------------------------------------------------------------------------------------
'''

'''
1. Begin by building only the basic diffusion solver (Eq. 2)
2. Add in spherical correction (Eq. 3-4)
3. Include radiative forcing  (Eq. 5)
In all scenarios: albedo = 0.3
- plot directions 1-3 together to demonstrate how Latitude vs Temperature evolves when adding these considerations
'''

'''
Part 1a. Basic implicity diffusion solver for heat transfer, 1D (Eq. 2)
'''
def solve_diffusion(npoints = 18, n_years=10000, dt_years=1, lam=100.0, init_temp=None):
    '''
    Parameters:
    npoints: int
        number of latitudinal grid cells/bands (18 in this lab)
    n_years: float
        total integration in years (default is 100000)
    dt_years: float
        Time step size in years (default is 1)
    lam: float
        Thermal diffusivity constant λ in m^2/s (default is 100)
    init_temp : numpy.ndarray or None
        Initial temperature vector in degrees Celcius. Defaults to temp_warm(lats) from given function

    Returns:
    lats: numpy.ndarray
        Latitude cell centers in degree (0 = South Pole, 180 = North Pole)
    T: numpy.ndarray 
        Temperatures once equilibrium profile has been reached across all latitudes, in degrees Celcius
    '''
    R = 6378e3 #earth raidus in meters
    dlat = 180 / npoints #making latitudinal spacing in degrees
    lats = np.linspace(dlat / 2.0, 180 - dlat / 2.0, npoints) #making the center of latitudinal cells
    dy = R * np.radians(dlat) #determining the change in y direction in meters, currently in radians 

    #First i build my matrix K that we solved in class with Neumann boundary conditions
    E = np.diag(np.ones(npoints-1), k=1) #shorthand trick to construct tridiagonal matrix K without writing a for loop
    K = (-2 * np.eye(npoints) + E + E.T)/(dy**2) #change in dy instead of dx implies changes along north-sputh (latiduinal)
    #need to apply special Neumann boundary conditions (zero flux at poles)
    #at i = 0 and i = npoints - 1 (north pole and south pole), the boundary cell reflects adjacent cells, modifying -2 to -1
    ##GHOST NODES AT THE POLE AND HOW THEY WORK!
    '''
    3-point 2nd derivative formula requires a node past poles to measure poles as that "center" value
    Neumann Boundary Condition = zero heat flux, dT/dy = 0 forces heeat to bounce back into adjacent interior cell rather than leaking outside of system
    South Pole (index 0):
    Using a ghost node T{-1} and T{1} path the north and south pole, zero flux requires T{-1} = T{1}
    Substituting T{-1} = T{1} into (T{-1} - 2T0 + T1)/dy^2 gives:
    (T1 - 2*T + T1)/dy^2 ----> (-2*T0 + 2*T1)/dy^2 
    OR
    North Pole (index -1): 
    Neumann boundary condition = zero heat flux, dT/dy = 0 forces heat to bounce back into system rather than leaking outside of system
    Using a ghost node TN past the north pole, zero flux requires TN = T{N-2}
    Substituting TN = T{N-2} into (T{N-2} - @*T{N-1} + TN)/dy^2 gives:
    (T{N-2} - 2*T{N-1} + T{N-2})/dy^2 ---> (-2*T{N-1} + @*T{N-2})/dy^2
    '''
    K[0,0] = -2/(dy**2)
    K[0,1] = 2/(dy**2)
    K[-1,-1] = -2/(dy**2)
    K[-1,-2] = 2/(dy**2)

    #defining implicity operator L = I - λ * Δt * K : this is the time-stepping matrix
    dt = dt_years * 365.25 * 86400 #converting to  years (as defined in above of function) to seconds
    I = np.eye(npoints)
    L = I - lam * dt * K
    L_inv = np.linalg.inv(L)

    #setting intial temperature vector Tj
    if init_temp is None:
        T = temp_warm(lats) #if we dont have an init temp lets just assume warm earth at all lats
    else:
        T = np.copy(init_temp) # making a copy instead of printing over original version

    #stepping through time, forward steps
    #total number of time steps
    total_steps = int(n_years/dt_years)
    for step in range(total_steps):
        T = L_inv @ T #multiply the inverse matrix by the current temperature vector T,
        #and overwrite T with the updated temperature values for the next time step
    #Tnew = L-1 @ Told
    return lats, T

lats, Td = solve_diffusion()

'''
Part 1b. 1D Diffusive Solver with Spherical Area Correction Eq. 3 (Spherical Correction considers the first order difference of initial temperature and the correction on each latitundinal band space cross section Axz)
'''
def solve_spherediffusion(npoints =18, n_years=10000, dt_years=1.0, lam=100.0, init_temp=None):
    ''' 
    Parameters:
    npoints: int
        Number of lat grid cells (18 for this lab)
    n_years: float
        total integration of time in years (10000 is default)
    dt_years: float
        time step size in years
    lam: float
        thermal diffusivity constant
    init_temp: numpy,ndarray
        initial temp vector in degrees celcius, defaults to temp_warm(lats)
    
    Returns:
    lats: numpy.ndarray
        latitude cell centers in degrees
    T: numpy.ndarray
        Equilibrium temperatures in lat bands in degrees Celcius
    '''
    R = 6378e3 #earth raidus in meters
    dlat = 180 / npoints #making latitudinal spacing
    lats = np.linspace(dlat / 2.0, 180 - dlat / 2.0, npoints) #making the center of latitudinal cells
    dy = R * np.radians(dlat) #determining the change in y direction in meters, currently in radians 

    #First i build my Second Derivative matrix K that we solved in class
    E = np.diag(np.ones(npoints-1), k=1) #shorthand trick to construct tridiagonal matrix K without writing a for loop
    K = (-2 * np.eye(npoints) + E + E.T)/(dy**2) #change in dy instead of dx implies changes along north-sputh (latiduinal)
    #need to apply special Neumann boundary conditions (zero flux at poles)
    #at i = 0 and i = npoints - 1 (north pole and south pole), the boundary cell reflects adjacent cells, modifying -2 to -1
    
    ##GHOST NODES AT THE POLE AND HOW THEY WORK!
    '''
    3-point 2nd derivative formula requires a node past poles to measure poles as that "center" value
    Neumann Boundary Condition = zero heat flux, dT/dy = 0 forces heeat to bounce back into adjacent interior cell rather than leaking outside of system
    South Pole (index 0):
    Using a ghost node T{-1} and T{1} path the north and south pole, zero flux requires T{-1} = T{1}
    Substituting T{-1} = T{1} into (T{-1} - 2T0 + T1)/dy^2 gives:
    (T1 - 2*T + T1)/dy^2 ----> (-2*T0 + 2*T1)/dy^2 
    OR
    North Pole (index -1): 
    Neumann boundary condition = zero heat flux, dT/dy = 0 forces heat to bounce back into system rather than leaking outside of system
    Using a ghost node TN past the north pole, zero flux requires TN = T{N-2}
    Substituting TN = T{N-2} into (T{N-2} - @*T{N-1} + TN)/dy^2 gives:
    (T{N-2} - 2*T{N-1} + T{N-2})/dy^2 ---> (-2*T{N-1} + @*T{N-2})/dy^2
    '''
    #south pole
    K[0,0] = -2/(dy**2) #fixing boundary condition at pole by recalculating with "ghost node" scenario, to make our pole the desired interior point of the lat band
    K[0,1] = 2/(dy**2) 
    #north pole
    K[-1,-1] = -2/(dy**2) #fixing boundary condition at pole by recalculating with "ghost node" scenario, to make our pole the desired interior point of the lat band
    K[-1,-2] = 2/(dy**2)

    #new here, I build my first derivative matrix, B, to solve for spherical correction
    B = (E - E.T) / (2.0 * dy) #creates centered different first derivative matrix, inputting interior diag values
    '''
    Pole Gradient Boundaries for first derivative:
    at both poles, zero flux requires dT/dy = 0 forces heat to bounace back into cell rather than leak out of poles and earth system
    For South Pole (i=0): ghost now T{-1} = T1 --> (T{i+1} - T{i-1})/(2*dy) = 0
    For North Pole (i=-1): ghost node TN = T{N-2} --> (TN - T{N-2})/(2*dy) = 0 
    Setting boundary rows to zero enforces zero spatial gradient across the poles
    '''
    B[0, :] = 0
    B[-1, :] = 0

    #defining latitundinal cross sectional ring area- we need to consider actual spherical shape and how it changes area of bands and correct for that
    theta = np.radians(lats)
    Axz = 2.0 * np.pi * R * np.sin(theta) #remember we are now in radians, this is the total area of the each latitudinal band

    #identity matrix and converting units, implicit operator
    dt = dt_years * 365.25 * 86400 #converting years (as defined in above of function) to seconds
    I = np.eye(npoints)
    L = I - lam * dt * K
    L_inv = np.linalg.inv(L)

    #setting intial temperature vector T
    if init_temp is None:
        T = temp_warm(lats) #if we dont have an init temp lets just assume warm earth at all lats
    else:
        T = np.copy(init_temp) # making a copy instead of printing over original version

    #stepping through time, to find temp at each time step Tj+1, with spherical geometry flux term
    #total number of time steps
    total_steps = int(n_years/dt_years)
    for step in range(total_steps):
        #apply spherical correction term at each time step for new temp
        f_spherical = lam *dt*(B@T) * (B@Axz)/Axz
        T = L_inv @ (T + f_spherical)
    
    return lats, T


lats, Ts = solve_spherediffusion()

'''
Part 1c. Diffusive Solver with Spherical Correction and Radiative Forcing: THE FINAL FUNCTION
'''

def snowballearth(npoints=18, n_years=10000, dt_years=1.0, lam=100.0, init_temp=None, e = 1.0, sigma = 5.67e-8, init_albedo = None, a_ice = 0.6, a_water = 0.3, S = 1370, p = 1020, dz = 50, heatcap = 4.2e6, dynamic_albedo=False, gamma=1):
    '''
    Parameters:
    npoints: int
        Number of latitudinal cells (18)
    n_years: float
        total simulation duration in yrs
    dt_years: float
        time step size in years
    lam: float
        thermal diffusivitiy of heat in m^2/s
    init_temp: numpy.ndarray or None
        initial vector of temps in kelvin, defaults to temp_warm(lats)
    e: float
        surface emissivity (defaults to 1)
    sigma: float
        Stefan-Boltzmann constant in (W/m^2/K^4)
        is 5.67e-8
    a_ice: float
        albedo of ice/snow (0.6)
    a_water: float
        albedo of water/ground (0.3)
    S: float
        Solar constant S0 in W/m^2 (1370)
    p: float
        seawater density in kg/m^3 (1020)
    dz: float
        ocean mixed-layer depth in meters (50)
    heatcap: float
        specific heat capacity of seawater in C in J/kg/K (default is 4.2e6)
    dynamic_albedo: bool
        if this is true, it will apply an ice-albedo threshold of 263.15 deg C. Default is false
    gamma: float
        solar forcing multiplier factor y (1.0)
    init_albedo: float
        initial albedo if dynamic_albedo is not None
        
    Returns:
    lats: nump.ndarray
        latitudinal cell cennters in degrees
    T: nump.ndarray
        equilibrium steady state surface temp profile along lat bands in deg Celcius
    '''
    R = 6378e3 #earth radius in meters
    dlat = 180 / npoints #making latitudinal spacing
    lats = np.linspace(dlat / 2.0, 180 - dlat / 2.0, npoints) #making the center of latitudinal cells
    dy = R * np.radians(dlat) #determining the change in y direction in meters, currently in radians 

    #First i build my Second Derivative matrix K that we solved in class with Neumann boundary conditions
    E = np.diag(np.ones(npoints-1), k=1) #shorthand trick to construct tridiagonal matrix K without writing a for loop
    K = (-2 * np.eye(npoints) + E + E.T)/(dy**2) #change in dy instead of dx implies changes along north-sputh (latiduinal)

    '''
    3-point 2nd derivative formula requires a node past poles to measure poles as that "center" value
    Neumann Boundary Condition = zero heat flux, dT/dy = 0 forces heeat to bounce back into adjacent interior cell rather than leaking outside of system
    South Pole (index 0):
    Using a ghost node T{-1} and T{1} path the north and south pole, zero flux requires T{-1} = T{1}
    Substituting T{-1} = T{1} into (T{-1} - 2T0 + T1)/dy^2 gives:
    (T1 - 2*T + T1)/dy^2 ----> (-2*T0 + 2*T1)/dy^2 
    OR
    North Pole (index -1): 
    Neumann boundary condition = zero heat flux, dT/dy = 0 forces heat to bounce back into system rather than leaking outside of system
    Using a ghost node TN past the north pole, zero flux requires TN = T{N-2}
    Substituting TN = T{N-2} into (T{N-2} - @*T{N-1} + TN)/dy^2 gives:
    (T{N-2} - 2*T{N-1} + T{N-2})/dy^2 ---> (-2*T{N-1} + @*T{N-2})/dy^2
    '''
    K[0,0] = -2/(dy**2)
    K[0,1] = 2/(dy**2)
    K[-1,-1] = -2/(dy**2)
    K[-1,-2] = 2/(dy**2)

    #here, I build my first derivative matrix, B, to solve for spherical correction
    B = (E - E.T) / (2.0 * dy) #creates centered different first derivative matrix, inputting interior diag values

    '''
    Pole Gradient Boundaries for first derivative:
    at both poles, zero flux requires dT/dy = 0 forces heat to bounace back into cell rather than leak out of poles and earth system
    For South Pole (i=0): ghost now T{-1} = T1 --> (T{i+1} - T{i-1})/(2*dy) = 0
    For North Pole (i=-1): ghost node TN = T{N-2} --> (TN - T{N-2})/(2*dy) = 0 
    Setting boundary rows to zero enforces zero spatial gradient across the poles
    '''

    B[0, :] = 0  #south pole has zero flux
    B[-1, :] = 0 #north pole has zero flux

    #defining latitundinal cross section areas
    theta = np.radians(lats)
    Axz = 2.0 * np.pi * R * np.sin(theta) #remember we are now in radians

    #defining implicity operator L = I - λ * Δt * K, time stepping implicit matrix
    dt = dt_years * 365.25 * 86400 #converting years (as defined in above of function) to seconds
    I = np.eye(npoints)
    L = I - lam * dt * K
    L_inv = np.linalg.inv(L)

    #calculating base solar insolation profile
    S_y = insolation(S, lats)

    #setting intial temperature vector Tj
    if init_temp is None:
        T = temp_warm(lats) #if we dont have an init temp lets just assume warm earth at all lats
    else:
        T = np.copy(init_temp) # making a copy instead of printing over original version

    #stepping through time, to find temp at each time step Tj+1
    #total number of time steps
    total_steps = int(n_years/dt_years)
    for step in range(total_steps):
        #determine spatial albedo array at each time step
        if dynamic_albedo:
            if step == 0 and init_albedo is not None:
                albedo = np.full_like(T, init_albedo) #make an array of that intial albedo defined
            else:
                albedo=np.where(T<= 263.15, a_ice, a_water) #if temp is less than or = to 263.15 K, make random select array of albedos of ice and water spatially
        else:
            albedo = np.full_like(T, a_water) #makes an array with only albedo of water

        #spherical correction flux term
        f_spherical = lam * dt * (B @ T) * (B @ Axz)/Axz

        #radiative forcing (converting celsius to kelvin for Stefan-Boltzmann)
        rad_forcing = (dt/(p*heatcap*dz))*(gamma*S_y*(1-albedo) - e * sigma * (T**4))

        #diffusion and implicit step update
        T = L_inv @ (T + f_spherical + rad_forcing)
        
    return lats, T

lats, Tr = snowballearth()

initial_condition = temp_warm(lats)

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 5))


# Convert Kelvin arrays back to Celsius (- 273.15) for plotting
plt.plot(lats, initial_condition - 273.15, label='Initial Condition', linewidth=2)
plt.plot(lats, Td - 273.15, label='Basic Diffusion', linewidth=2)
plt.plot(lats, Ts - 273.15,linestyle= '--', label='Diff. + Spherical Correction', linewidth=2)
plt.plot(lats, Tr - 273.15, linestyle= '--', label='Diff. + SphCorr + Radiative', linewidth=2)

plt.xlabel('Latitude (°)')
plt.ylabel('Temperature (°C)')
plt.title('Evolution of Latitudinal Temperature Profile (Part 1)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

#confirm plot looks correct before moving on

'''
-------------------------------------------------------------------------------------------
PART 2: WARM EARTH EQUILIBRIUM SCENARIO
-------------------------------------------------------------------------------------------
'''

'''
Turn model so that is can reproduce the warm-Earth equilibrium. 
albedo = 0.3 still
diffusivity and emissivity contain mcuch uncertainty (explore each independently):
    - allow diffusivity (λ) to range from 0 to 150
    - allow emissivity (e) to range from 0 to 1
determine their impact on the equilibrium solution
    - pick a value for each (diffusivitiy and emissivity) when used in combination that allows for you to reproduce "warm Earth" curve given by temp_warm()
    - report these findings and use them for the rest of the lab
'''

'''
Part 2a: changing diffusivity (0-150): warm earth equilibrium target
'''
lam_bin = {} #bin to store results in
lam_values = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150] #testing lambda values in 10 increments

for lam in lam_values: #in range of my lambdas im testestimg
    lats, T_lam = snowballearth( #what are the output lats, Temps from my function snowball earth?
            npoints=18,
            n_years=10000,
            dt_years=1.0,
            lam=lam,
            e=1,
            dynamic_albedo=False
        )
    lam_bin[lam] = T_lam #bin those outputs here can call them T_lam

#plot case scenarios to see how changing diffusivity changes temperature across latitudes
plt.figure(figsize=(8, 5))
for lam, T_lam in lam_bin.items():
    plt.plot(lats, T_lam - 273.15, linewidth=2, label=f'Diffusivity = {lam}')
plt.plot(lats, temp_warm(lats)-273.15, 'k--', linewidth=2, label='Warm Earth Target')
plt.xlabel('Latitude (°)')
plt.ylabel('Temperature (°C)')
plt.title('Part 2a: Testing Diffusivity Range (0-150)')
plt.legend(bbox_to_anchor=(1.04, 1), loc='upper left', borderaxespad=0)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


'''
Part 2b: changing emissivity (0-1): warm earth equilbrium target
'''
eps_bin={} #store eps values here
e_values = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1] #emmisivity values to test

for eps in e_values: #in those range of values im testing
    lats, T_eps = snowballearth( #what are the output temperatures at lats
            npoints=18,
            n_years=10000,
            dt_years=1.0,
            lam=lam,
            e=eps,
            dynamic_albedo=False
        )
    eps_bin[eps]=T_eps #store them here 

#and plot each output to compare how emissivity changes the temp across lats
plt.figure(figsize=(8, 5))
for eps, T_eps in eps_bin.items():
    plt.plot(lats, T_eps - 273.15, linewidth=2, label=f'Epsilon = {eps: .2f}')
plt.plot(lats, temp_warm(lats)-273.15, 'k--', linewidth=2.5, label='Warm Earth Target')
plt.xlabel('Latitude (°)')
plt.ylabel('Temperature (°C)')
plt.title('Part 2b: Testing Emissivity Range (0.30-1)')
plt.legend(bbox_to_anchor=(1.04, 1), loc='upper left', borderaxespad=0)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


'''
Part 2c: combined diffusivity and emissivity values for warm earth target: finding optimal values
'''
#narrowing search to finer steps around my first wider step estimate done with this same code block (30, 0.70) that has been deleted to accomodated a finer search now
lam_values = np.arange(20, 75, 5) #test values in 5 increments between 20 and 75
e_values = np.arange(0.68, 0.79, 0.01) #test values in 0.01 increments between 0.68 and 0.79

T_target = temp_warm(lats) #the target temps at lats im trying to get with tuning emissivity and diffusivity

best_rmse = np.inf #we are going to use root mean square error (standard statistical metric that measures the average magntiude of prediction errors), initialized to go to infinity, so when we iterate to find the "best" rmse later on, automatically the first iteration will be better than infinity and be the "record holder", and so forth through iterations until the best match of e and lam is found for temp_warm conditions
best_params = None #setting up match area for the best matched parameters (emissivitiy and diffusivity) to go that give temps closest to warm earth
best_T = None #setting up match area for the lat temperatures for the optimized emissivity and diffusivity to get closest to T_target warm earth conditions (given to us in beginning)

#search over (lam, e)
for lam in lam_values: #values i made
    for e in e_values: #values i made
        lats, T_model = snowballearth(
            npoints=18,
            n_years=10000,
            dt_years=1.0,
            lam=lam,
            e=e,
            dynamic_albedo=False
        )
        rmse = np.sqrt(np.mean((T_model-T_target)**2)) #for each iteration of lam and e and the temps they output, lets calc the rmse

        if rmse < best_rmse: #tracking the best result so far with the smallest rmse, does it beat the previous best?
            best_rmse = rmse #if so, it becomes the new best to beat
            best_params = (lam, e) #save the winning combination of lam and e
            best_T = T_model #and save the temp outputs of this lam and e combination

#printing winning output
print(f"Optimal Diffusivity (lam): {best_params[0]}")
print(f"Optimal Emissivity (e): {best_params[1]:.2f}")
print(f'Lowest RMSE: {best_rmse: .2f} K')

#plotting the results
plt.figure(figsize=(8, 5))
#plotting target temps (given)
plt.plot(lats, T_target - 273.15, label='Target: Warm Earth (temp_warm)', linewidth=2.5, color='black')
#plotting my model output from given e, lam, and T_model
plt.plot(lats, best_T - 273.15, linestyle='--', linewidth=2.5, color='crimson', label=f'Model Fit, diffusivity = {best_params[0]}, emissivity = {best_params[1]:.2f}')

plt.xlabel('Latitude (°)')
plt.ylabel('Temperature (°C)')
plt.title('Part 2: Calibrating Diffusivity and Emissivity to Modern Earth')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


'''
------------------------------------------------------------------------------------------
PART 3: HOW DO INITIAL CONDITIONS AFFECT THE EQUILIBRIUM SOLUTION?
------------------------------------------------------------------------------------------
'''
npoints = 18
dlat= 180/ npoints
lat = np.linspace(dlat/2, 180-dlat/2, npoints)
#new accepted calibrated diffusivity and emissivity constants from Part 2: to be used for remainder of the lab
lam_cal = 0.35
e_cal = 0.72

#here, this is where it becomes important that my snowball earth function accepts init_temp and init_albedo

#create initial start temp conditions
hot_start = np.ones(npoints) * (60+273.15) #initial hot temps in KELVIN
cold_start = np.ones(npoints) * (-60+273.15) #initial cold temps in KELVIN
warm_start = temp_warm(lats)
#else wise we continue using temp_warm(lats) as our intial conditions as defined previously in our function

'''
1. Switch from constant to dynamic albedo
    a. begin with hot Earth (60*C at all locations)- what is equilibrium solution?
    b. repeat for cold Earth (-60*C at all locations)- what is your equilibrium solution?
    c. "flash freeze" the Earth by starting with the Warm Earth solution curve used in Parts 1/2, but now albedo is set to 0.6
        - How do these results compare to each other?
        - What does it tell us about the stability of snowball vs warm Earth solutions?
'''
#Part 3a
lats, T_hot = snowballearth(npoints=npoints, lam=lam_cal, e=e_cal, init_temp=hot_start, dynamic_albedo=True)

#Part 3b
lats, T_cold = snowballearth(npoints = npoints, lam= lam_cal, e= e_cal, init_temp=cold_start, dynamic_albedo=True)

#Part 3c
lats, T_flash = snowballearth(npoints=npoints, lam=lam_cal, e=e_cal, init_temp=warm_start, a_water=0.6, dynamic_albedo=True)

#plotting the three different scenarios
plt.figure(figsize=(8,5))
plt.plot(lats, T_hot - 273.15, color='red', linewidth=2.5, label='3a: Hot Earth Start (+60°C)')
plt.plot(lats, T_cold - 273.15, color='blue', linewidth=2.5, label='3b: Cold Earth Start (-60°C)')
plt.plot(lats, T_flash - 273.15, color='cyan', linestyle='--',linewidth=2.5, label='3c: Flash Freeze (Warm T, ALbedo = 0.6)')

#plotting warm earth reference curve and freezing threshold line
plt.plot(lats, temp_warm(lats) - 273.15, color='black', linestyle=':', linewidth=2, label='Warm Earth Reference')
plt.axhline(-10, color='gray', linestyle='--', linewidth=1, label='Freezing Threshold (-10°C)')\

plt.xlabel('Latitude (°)')
plt.ylabel('Equilibrium Temperature (°C)')
plt.title('Part 3: How Initial Conditions Affect Equilibrium')
plt.legend(bbox_to_anchor=(1.04,1), loc='upper left', borderaxespad=0)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


'''
-------------------------------------------------------------------------------------------
4. EXPLORING IMPACTS OF SOLAR FORCING ON SNOWBALL EARTH
-------------------------------------------------------------------------------------------
'''

'''
Create a solar multipilier γ, which is applied to the insolation term in your solver
1. Starting with cold Earth intial condition, run simulation with γ = 0.4
    - use the result as an initial condition for another simulation, increasing γ by 0.05
    2. Repeat this until γ = 1.4
        - at this turn around point: lower γ in 0.05 incremements back to γ = 0.4
        - each time, use the previous result as the intial condition to the next simulation
3. Plot average global temperature versus γ
Answer the questions:
- Does the snowball Earth hypothesis represent an equilibrium solution that is stable?
- What does this plot tell you about the stability of the different equilibria?
- Given the range for γ and considering the historical variation of Solar flux, do you think snowball Earth is a valid hypothesis?
'''
npoints = 18
dlat = 180.0 / npoints
lats = np.linspace(dlat / 2.0, 180.0 - dlat / 2.0, npoints)
#this question asks for global avg temperature, so i need to take the avg temperature from lats
weights = np.sin(np.radians(lats)) #area weights for global avg temp
current_T = np.ones(npoints) * (-60.0 + 273.15) #assuming snowball earth passed on from previous step
lam_cal = 0.35 #just rewriting down here for myself so i don't have to scroll
e_cal = 0.72 #just rewriting down here for myself so i don't have to scroll

#creating gamma arrays increasing and decreasing (warming and cooling)
gamma_up = np.round(np.arange(0.4, 1.41, 0.05), 2)
gamma_down = np.round(np.arange(1.4, 0.39, -0.05), 2)

#making bins to store global avg temp at each gamma step
T_avg_up = []
T_avg_down = []

#starting with cold earth (-60°C everywhere), gamma steps up from 0.4 to 1.4 in equal intervals
for g in gamma_up:
    lats, current_T = snowballearth(
        npoints=npoints,
        lam=lam_cal,
        e=e_cal,
        init_temp=current_T, #feeding previous equilibrium into this run
        dynamic_albedo=True,
        gamma=g
    )
    global_mean_T_up = np.average(current_T - 273.15, weights=weights) #find avg global temp at each gamma value
    T_avg_up.append(global_mean_T_up) #save temp in bin

#turning around at gamma = 1.4 and stepping back down to 0.4, current temp is already at gamma = 1.4 equilibrium from the end of the loop before
for g in gamma_down:
    lats, current_T= snowballearth(
        npoints=npoints,
        lam=lam_cal,
        e=e_cal,
        init_temp=current_T, #start at temp previous loop ended at equilibrium
        dynamic_albedo=True,
        gamma=g
    )
    global_mean_T_down = np.average(current_T - 273.15, weights=weights)
    T_avg_down.append(global_mean_T_down)

#plotting global temperature vs gamma radiation factor
plt.figure(figsize=(8,5))
plt.plot(gamma_up, T_avg_up, 'b-o', linewidth=2.5, label='Warming Path (gamma: 0.4 up to 1.4)')
plt.plot(gamma_down, T_avg_down, 'r--s', linewidth=2.5, label='Cooling Path (gamma: 1.4 down to 0.4)')
plt.xlabel('Solar Forcing Multiplier, Gamma')
plt.ylabel('Global Average Temperature (°C)')
plt.title('Part 4: Solar Forcing Hysteris Loop (Snowball vs Warm Eath)')
plt.legend(bbox_to_anchor=(1.04,1), loc='upper left', borderaxespad=0)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()