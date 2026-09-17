#🌳🌳🌳PART 1🌳🌳🌳
'''
This file contains tools and scripts for completing Lab 1 for ClaSP 410.

To reproduce the plots in the lab report, do this...
'''
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

#🌳🌳🌳PART 2🌳🌳🌳
#defining intitial variables and conditions
nx, ny = 5, 5 #number of cells in X and Y direction
# i have set my grid to 5x5 to have "ghost borders" that already are burned
prob_spread = 1 #Chance to spread to adjacent cells
prob_bare = 0 #Chance of cell to start at bare patch
prob_ignite = 0 #Chance of cell to start on fire
Bare = 1
Forested = 2
OnFire = 3
nstep = 6


#🌳🌳🌳PART 3🌳🌳🌳
#create an initial grid, set all values to "2
# type in array to integers only
#note that (nx, ny) is the same as (i, j)
forest = np.zeros((nstep, nx, ny),
dtype=int) + 2 #here, I have defined a new name "forest" as an array of my 5 by 5 array and conditionally set all cells to 2
#colons call all values in that index
print ("initial forest conditions is: ", forest)

'''
another more streamlined way to make the ghost nodes (our exterior parameter of values that will then be ignored but prevent us from having to write edge conditions)
is to create a padding around the forest with the following:
forest_ghostnodes = np.pad(forest, pad_width+1,mode='constant', constant_values=-1)
having -1 as our ghost values sets them apart
'''


#🌳🌳🌳PART 4🌳🌳🌳
'''
something new I learned here: Python uses zero-based indexing (0,0) is the start of an array
the function numpy.random can roll multiple probablities simultaneously
A Better way to loop through random probabilities is by defining a new function that 
can be turned into True/False
'''

# 3D array of whole forest over ranges nstep, ny, nx
#create an array of randomly generated numbers of range [0,1): here we are creating bare patches
for i in range(nx):
    for j in range(ny):
        #roll our dice to see if we get a bare spot
        if np.random.rand() < prob_bare:
            forest[0,j,i] = 1 #0 is bare
        elif np.random.rand() < prob_ignite:
            forest[0,j,i] = 3 #if not a bare patch then ignite it, leaving a prob for a cell to be forested
forest[0,:,0]=1
forest[0,0,:]=1
forest[0,:,nx-1]=1
forest[0,ny-1,:]=1
# all cells are set to two which is the number to represent forested. our ghost borders have been set to one as alrerady burned so they dont disrupt the model

#set the center grid cell to "burning" this is [time, j, i]. This will be deleted when I set my probabilities to not 1's and 0's
forest[0, 2, 2] = 3

print ("updated initial forest conditions:", forest)

'''
#🌳🌳🌳PART 5🌳🌳🌳
#🔥TIME TO ALLOW THE FIRE TO SPREAD🔥
#we want to check the orthogonal neighbors of each cell of the array (up, down, left, right)
#an example is in array syntax, say forest[j,i] to see if it has status 3 (i.e., actively burning). If so, we look at all neighbor cells to spread the fire
#Loop in "x" direction:
#range will run in a loop, we need to tell it to start at spot 1 so the function doesnt loop back around, then giving you the wrong answer
#Loop in "y" direction
#we need to tell python that if the x statement isnt true, to pass and continue to the next step
'''
fig, ax = plt.subplots(1,1)
import matplotlib.pyplot as plt #pyplot is the main plotting tool
from matplotlib.colors import ListedColormap #importing an object for creating own special color map
#  https://matplotlib.org/stable/gallery/color/named_colors.html
curr_forest_cmap = ListedColormap(['tan', 'darkgreen', 'firebrick'])
curr_forest = np.copy(forest[0, :, :])
pred_forest = np.copy(curr_forest)
for k in range(0, nstep):
    for i in range(1,nx-1):
        for j in range(1, ny-1):
            print(j,i,pred_forest[j,i])
            #Check Right
            if curr_forest[j, i]==3:
                if (curr_forest[j, i+1] == 2) and (np.random.rand() < prob_spread):
                    pred_forest[j, i+1] = 3
                print("changing!")
            #Check Left
                if (curr_forest[j, i-1] == 2) and (np.random.rand() < prob_spread):
                    pred_forest[j, i-1] = 3
                print("changing!")
            #Check Up
                if (curr_forest[j+1, i] == 2) and (np.random.rand() < prob_spread):
                    pred_forest[j+1,i] = 3
                print("changing!")
            #Check down
                if (curr_forest[j-1, i] == 2) and (np.random.rand() < prob_spread):
                    pred_forest[j-1, i] = 3
                print("changing!")
                #This is indexing to the original burning cell/s and is now making it bare 
                pred_forest[j,i] = 1
            #print(j,i,pred_forest[j,i])
    curr_forest = np.copy(pred_forest)
    #print(curr_forest)
    print(pred_forest)

    #Part 6: Creating Figures!!!🎨🎨🎨🌳🌳🌳

    #counting (math) the total number of cells in each condition along each iteration k
    n_bare = (curr_forest == 1).sum() #count our tan squares that are burned and barren 
    n_forested = (curr_forest == 2).sum() #count our forest squares that are unburned
    n_burning = (curr_forest == 3).sum() #count our actove fire squares
    total_cells = curr_forest.size #acknowledging all cells

    #we want to convert our square counts into readable percentages
    pct_bare = (n_bare / total_cells) * 100
    pct_forested = (n_forested / total_cells) * 100
    pct_burning = (n_burning / total_cells) * 100

    #Plotting
    ax.clear()  #clear the axis data, NOT whole figure

    title_text = (
        f'Iteration = {k+1:03d}\n'
        f'Bare: {pct_bare:.1f}%  |   ' #the :.f % tells it we only want one decimal place
        f'Forested: {pct_forested:.1f}%  |   '
        f'Burning: {pct_burning:.1f}%'
    )
    ax.set_title(title_text, fontsize=10, loc='left')

    #use pcolor to display the current state
    ax.pcolor(curr_forest, cmap=curr_forest_cmap, vmin=1, vmax=3)
    
    plt.draw()
    plt.pause(1.0)
    
    #we need to tell our system to quit once no more cells are on fire
    nBurn = (curr_forest == 3).sum()
    if nBurn == 0:
        print(f"Burn completed in {k+1} steps")
    #cleanly escape the loop

    k += 1 #advance time step continously
    #keep final window open when the script finished
plt.ioff()
plt.show()
