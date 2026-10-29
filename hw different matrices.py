
import numpy as np

N = 4

dx = 2/(N-1)

x = np.linspace(-1, 1, N).reshape(-1,1)

u = 1.0-x**2


E1 = np.diag(np.ones(N))

E2 = np.diag(np.ones(N), -1)

E3 = np.diag(np.ones(N), 1)

I = np.eye(N)

E = np.diag(np.ones(N-1),1)

K = (-2*np.eye(N) + E + E.T)/dx**2

print(K)
'''
The output K:
[[-4.5   2.25  0.    0.  ]
 [ 2.25 -4.5   2.25  0.  ]
 [ 0.    2.25 -4.5   2.25]
 [ 0.    0.    2.25 -4.5 ]]
'''

K.dot(u)

print(K.dot(u))

print(u)
'''
The output K.dot(u):
[[ 2.]
 [-2.]
 [-2.]
 [ 2.]]
'''


'''
Modifying to apply boundary conditions
'''
#exact second derivate expected u''(x)
f = -2 * np.ones(N)

K = (-2*np.eye(N) + E + E.T)/dx**2

#changing boundaries manually to what they actually are based on manual calculation
#modifying top
K[0, :]=0
K[0,0]=1
f[0]=0
#modifying bottom
K[-1,:]=0
K[-1,-1]=1
f[-1]=0

'''
now I have a matric K that I can multiply my u by to get an output for like the atmosphere scenario
'''
K.dot(u)
print(K.dot(u)) #this now properly applies my conditions of 0
'''
[[ 0.]
 [-2.]
 [-2.]
 [ 0.]]
'''
U_solver = np.linalg.solve(K, f) #this now gives me the values of u ! backwards solving by knowing a forward step (f, second derivastive)
print(U_solver)
'''
[0.         0.88888889 0.88888889 0.        ]
'''

'''
How do you modify the system of linear equations to satisfy the boundary conditions du/dx(-1)=0 and u(1) = 1
'''
'''
Dirichlet Condition = Fixed Value 
u(a) = number
In this case: u(1) = 1
physical meaning: setting a fixed state at the edge, holding end of string at fixed value
Advice: set the diagonal entry to 1 and the right-hand side f to 1. Directly sets uN=1
'''
'''
Neumann Condition = Derivative (Rate of Change)
the slope or rate of change of u(x) at a boundary du/dx(a) in this case du/dx(-1)=0
physical meaning: setting a flow, rate, flux or slope at the edge
Advice: Use a ghost point or one-sided differences to express the slope in terms of grid values
'''