import numpy as np
import scipy 
from scipy import interpolate
#Make grid
dx = 1e3 #resolution for x-direction this is 1 km
dy = 1e3 #Resolution for y-direction this is 1 km
Lx = 20e3 #Length of domain in x-direction
Ly = 10e3 #Length of domain in y-direction

#Discretize in x and y direction
x = np.arange(0,Lx+dx,dx) 
y = np.arange(0,Ly+dy, dy)

well_info = (
    [0, 0, 967], #x, y, hydropotential m^2/s^2
    [0, 3e3, 985],
    [0, 6e3, 1000],
    [0, 8e3, 990],
    [20e3, 1e3, 965],
    [20e3, 4e3, 995],
    [20e3, 9e3, 970]
)

wells_left = [np.array([0,3,6,8])*1e3, np.array([967, 985,1000,990])]
wells_right=[np.array([1,4,9])*1e3, np.array([965, 995, 970])]


potential_left = interpolate.interp1d(wells_left[0], wells_left[1], kind='linear',fill_value="extrapolate") (y)
potential_right = interpolate.interp1d(wells_right[0], wells_right[1], kind='linear', fill_value="extrapolate") (y)

#Added by Makayla!! :)

potential_top = interpolate.interp1d(np.array([x[0], x[-1]]), np.array([potential_left[0], potential_right[0]]), 
                                     kind='linear',fill_value="extrapolate") (x)
potential_bottom = interpolate.interp1d(np.array([x[-1], x[-1]]), np.array([potential_left[-1], potential_right[-1]]),
                                         kind='linear',fill_value="extrapolate") (x)
