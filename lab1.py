#🌳🌳🌳PART 1🌳🌳🌳
'''
This file contains tools and scripts for completing Lab 1 for ClaSP 410.

To reproduce the plots in the lab report, do this...
'''
#importing necessary packages
import numpy as np
import matplotlib
matplotlib.use('TkAgg') #setting up back end
import matplotlib.pyplot as plt
plt.ion() #this makes interactive plots activated 

#🌳🌳🌳PART 2🌳🌳🌳
#defining intitial variables and conditions
nx, ny = 70, 70 #number of cells in X and Y direction
# i have set my grid to have "ghost borders" that already are burned/bare
prob_spread = .3 #Chance to spread to adjacent cells
prob_bare = .02 #Chance of cell to start at bare patch
prob_ignite = .3 #Chance of cell to start on fire
Bare = 1 #identificstion for cells that are bare
Forested = 2 #identification for cells forested
OnFire = 3 #identification for cells on fire
nstep = 5 #range of timesteps "k"


#🌳🌳🌳PART 3🌳🌳🌳
#create an initial grid, set all values to "2" aka forested
# type in array to integers only
#note that (nx, ny) is the rsange of all i, j with indexing starting at 0
forest = np.zeros((nstep, nx, ny),
dtype=int) + 2 #here, I have defined a new name "forest" as an array of my 5 by 5 array and conditionally set all cells to 2
#colons call all values in that index
print ("initial forest conditions is: ", forest)

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

#here i am making my "ghost nodes" that will be set to bare from the beginning. 
forest[0,:,0]=1
forest[0,0,:]=1
forest[0,:,nx-1]=1
forest[0,ny-1,:]=1

#set the center grid cell to "burning" this is [time, j, i]. This will be deleted when I set my probabilities to not 1's and 0's
#forest[0, 2, 2] = 3

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
#for making figs later- lets define these
fig, (ax_map, ax_graph) = plt.subplots(1, 2, figsize=(12,5))
import matplotlib.pyplot as plt #pyplot is the main plotting tool
from matplotlib.colors import ListedColormap #importing an object for creating own special color map
#  https://matplotlib.org/stable/gallery/color/named_colors.html
curr_forest_cmap = ListedColormap(['tan', 'darkgreen', 'firebrick'])
#to make a graph we need to store or old data from our loop somewhere so it isnt just tossed out
history_time = []
history_bare=[]
history_forested=[]
history_burning=[]

#for the looping 
curr_forest = np.copy(forest[0, :, :])
pred_forest = np.copy(curr_forest)

for k in range(0, nstep):
    for i in range(1,nx-1): #here the +1 considers the padding of our ghost nodes predetermined as bare
        for j in range(1, ny-1):
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
    curr_forest = np.copy(pred_forest)
    print(curr_forest)
    print(pred_forest)

    #Part 6: Creating Figures!!!🎨🎨🎨🌳🌳🌳
    #still in our time loop very important
    #first, we need to cut out the ghost nodes to not ruin our probabilities and statistics
    active_view=curr_forest[1:-1,1:-1] #this is cutting the array ghost nodes off for our visualization and statistics portion (they are bare)

    #counting (math) the total number of cells in each condition along each iteration k
    n_bare = (active_view == 1).sum() #count our tan squares that are burned and barren 
    n_forested = (active_view == 2).sum() #count our forest squares that are unburned
    n_burning = (active_view == 3).sum() #count our actove fire squares
    total_cells = active_view.size #acknowledging all cells (minus ghost nodes now cut off)

    #we want to convert our square counts into readable percentages
    pct_bare = (n_bare / total_cells) * 100
    pct_forested = (n_forested / total_cells) * 100
    pct_burning = (n_burning / total_cells) * 100

    #make sure each iteration of our loop is captured in % form to then graph
    history_time.append(k+1)
    history_bare.append(pct_bare)
    history_forested.append(pct_forested)
    history_burning.append(pct_burning)

    #Plotting

    # 1. MAP
    ax_map.clear()

    map_title_text = (
        f'Iteration = {k+1:03d}\n'
        f'Bare: {pct_bare:.1f}%  |   ' #the :.f % tells it we only want one decimal place
        f'Forested: {pct_forested:.1f}%  |   '
        f'Burning: {pct_burning:.1f}%'
    )
    ax_map.set_title(map_title_text, fontsize=10, loc='left')

    #use pcolor to display the current state
    ax_map.pcolor(active_view, cmap=curr_forest_cmap, vmin=1, vmax=3)

    ax_graph.plot(history_time, history_bare, color='tan', label='Bare', linewidth=2)
    ax_graph.plot(history_time, history_forested, color='darkgreen', label='Forested', linewidth=2)
    ax_graph.plot(history_time, history_burning, color='firebrick', label='Burning', linewidth=2)

    ax_graph.set_title("Ecosystem Coverage Over Time", fontsize=10)
    ax_graph.set_xlabel("Time Step (Iteration)")
    ax_graph.set_ylabel("Grid Coverage (%)")

    ax_graph.set_ylim(0,101) #since we are reading in percentages
    ax_graph.set_xlim(1, nstep) #will show graph up to however many time steps I have set
    if k == 0: #if you don't tell it to only put the lgend once it will add a new one each time step 
        ax_graph.legend(loc='upper right')

    plt.draw()       # Tells the global pyplot manager to process modifications at each time step
    plt.pause(0.5)   # Pauses execution briefly in between steps for visualization
plt.ioff()
plt.show()