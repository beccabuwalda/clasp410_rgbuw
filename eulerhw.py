#importing necessary packages
import numpy as np
import matplotlib
matplotlib.use('TkAgg') #setting up back end
import matplotlib.pyplot as plt
plt.ion() #this makes interactive plots activated 

Ti = 100 #initial temperature of coffee
Tf = 80 #desired final temperature of our coffee
Ts = 18 #temperature of the environment
k = 0.2 #constant
#list of different time step widths to test
dt_values = [0.5, 0.1, 0.05, 0.01, 0.001] #when we change our size of time step, the estimation will become more accurate of the curve we are trying to predict
tol = 1 #our tolerance degree threshold that we will permit our coffee to be "drinkable"
#i = our time steps

plt.figure(figsize=(8,8))

for dt in dt_values:
    T_curr = float(Ti) #using float instead of np.copy cleanly saves exact numbers
    #binning our time and temperature for plotting:
    time_history=[0.0] #time bin
    temp_history=[T_curr] #temp bin
    for i in range(10000): #number of time steps permitted
        if T_curr < Tf+tol:
            break
        #if threshold has not been reached, apply eulers method
        T_new = T_curr - dt*(k*(T_curr-Ts))
        T_curr = (T_new)
        #saving the time and temp history at each step within the loop
        time_history.append((i+1)*dt) #we add a plus one to consider our initial step i of 0 
        temp_history.append(T_curr)
    #plotting complete trajectory for current dt
    plt.plot(time_history, temp_history, label=f'dt ={dt}')
    print(f"dt = {dt:<5} | Target reached in {time_history[-1]:.2f} time units ({i} steps)")
#Plot aesthetics
plt.axhline(y=Tf, color='r', linestyle='--', label=f'Target ({Tf}°C)')
plt.axhspan(ymin=Tf-tol, ymax=Tf+tol, color='r', alpha=0.2, label=f'Tolerance ±{tol}°C')

plt.xlabel('Time')
plt.ylabel('Temperature')
plt.title('Coffee Cooling Simulation: Practicing Eulers Method')
plt.legend(loc='upper right')
plt.ioff()
plt.show()




