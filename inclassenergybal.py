#importing necessary packages
import numpy as np
import matplotlib
matplotlib.use('TkAgg') #setting up back end
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import ode #used later for ODE45
from scipy.integrate import solve_ivp

S= 2600
a= 0.71
sigma= 5.67*10**-8

T = (((1-a)*S)/(2*sigma))**(1/4)

print(T)