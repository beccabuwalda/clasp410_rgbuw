#🧟🧟🧟PART 1🧟🧟🧟
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

#🧟🧟🧟PART 2🧟🧟🧟
#defining intitial variables and conditions
nx, ny = 70, 70 #number of people in X and Y direction
# i have set my grid to have "ghost borders" that already are zombified. Once a cell is zombified
prob_infected = .35 #chance of person to become infected
prob_transmission = 0.40 #chance of person being healthy
Healthy = 1 #identification for cells forested
Infected = 2 #identification for cells on fire
Immune = 3
nstep = 5 #range of timesteps "k" minutes


#🧟🧟🧟PART 3🧟🧟🧟
#create an initial grid, set all values to "2" aka forested
# type in array to integers only
#note that (nx, ny) is the rsange of all i, j with indexing starting at 0
infec = np.zeros((nstep, nx, ny),
dtype=int) + 1 #here, I have defined a new name "infec" as an arrayand conditionally set all cells to 1 (aka all my people are not sick and zombies yet)
#colons call all values in that index
print ("initial health of crowd", infec)

#🧟🧟🧟PART 4🧟🧟🧟
'''
something new I learned here: Python uses zero-based indexing (0,0) is the start of an array
the function numpy.random can roll multiple probablities simultaneously
A Better way to loop through random probabilities is by defining a new function that 
can be turned into True/False
'''

# 3D array of whole forest over ranges nstep, ny, nx
#create an array of randomly generated numbers of range [0,1): here we are creating probabilities of a person being a zombie or infected
for i in range(nx):
    for j in range(ny):
        #roll our dice to see if we get a zombie
        if np.random.rand() < prob_infected:
            infec[0,j,i] = Infected #2 is infected

#here i am making my "ghost nodes" that will be set to immune to infection from the beginning. 
infec[0,:,0]= Immune
infec[0,0,:]= Immune
infec[0,:,nx-1]= Immune
infec[0,ny-1,:]= Immune

#set the center grid cell to "infected" this is [time, j, i]. This will be deleted when I set my probabilities to not 1's and 0's
#infec[0, 2, 2] = 3

print ("initial crowd health:", infec)

'''
#🌳🌳🌳PART 5🌳🌳🌳
#🔥TIME TO ALLOW THE FIRE TO SPREAD🔥
#we want to check the orthogonal neighbors of each cell of the array (up, down, left, right)
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
curr_infec_cmap = ListedColormap(['goldenrod','firebrick'])
#to make a graph we need to store or old data from our loop somewhere so it isnt just tossed out
history_time = []
history_healthy=[]
history_infected=[]

#for the looping 
curr_infec = np.copy(infec[0, :, :])
pred_infec = np.copy(curr_infec)

for k in range(0, nstep):
    for i in range(1,nx-1): #here the +1 considers the padding of our ghost nodes predetermined as bare
        for j in range(1, ny-1):
            #Check Right
            if curr_infec[j, i]==2:
                if (curr_infec[j, i+1] == 1) and (np.random.rand() < prob_infected):
                    pred_infec[j, i+1] = 2
                print("changing!")
            #Check Left
                if (curr_infec[j, i-1] == 1) and (np.random.rand() < prob_infected):
                    pred_infec[j, i-1] = 2
                print("changing!")
            #Check Up
                if (curr_infec[j+1, i] == 1) and (np.random.rand() < prob_infected):
                    pred_infec[j+1,i] = 2
                print("changing!")
            #Check down
                if (curr_infec[j-1, i] == 1) and (np.random.rand() < prob_infected):
                    pred_infec[j-1, i] = 2
                print("changing!")
                #This is indexing to the original infected cell/s and leaving them infected
                pred_infec[j,i] = 2
    curr_infec = np.copy(pred_infec)
    print(curr_infec)
    print(pred_infec)

    #Part 6: Creating Figures!!!🎨🎨🎨🌳🌳🌳
    #still in our time loop very important
    #first, we need to cut out the ghost nodes to not ruin our probabilities and statistics
    active_view=curr_infec[1:-1,1:-1] #this is cutting the array ghost nodes off for our visualization and statistics portion (they are bare)

    #counting (math) the total number of cells in each condition along each iteration k
    n_healthy = (active_view == 1).sum() #count our healthy people
    n_infected = (active_view == 2).sum() #count our infected people
    total_cells = active_view.size #acknowledging all cells (minus ghost nodes now cut off)

    #we want to convert our square counts into readable percentages
    pct_healthy = (n_healthy / total_cells) * 100
    pct_infected = (n_infected / total_cells) * 100

    #make sure each iteration of our loop is captured in % form to then graph
    history_time.append(k+1)
    history_healthy.append(pct_healthy)
    history_infected.append(pct_infected)

    #Plotting

    # 1. MAP
    ax_map.clear()

    map_title_text = (
        f'Iteration = {k+1:03d}\n'
        f'Healthy: {pct_healthy:.1f}%  |   ' #the :.f % tells it we only want one decimal place
        f'Infected: {pct_infected:.1f}%'
        )
    ax_map.set_title(map_title_text, fontsize=10, loc='left')

    #use pcolor to display the current state
    ax_map.pcolor(active_view, cmap=curr_infec_cmap, vmin=1, vmax=3)

    ax_graph.plot(history_time, history_healthy, color='goldenrod', label='Healthy', linewidth=2)
    ax_graph.plot(history_time, history_infected, color='firebrick', label='Infected', linewidth=2)

    ax_graph.set_title("Spread of Infectious Disease Over Time", fontsize=10)
    ax_graph.set_xlabel("Time Step (Iteration)")
    ax_graph.set_ylabel("Infection Spread (%)")

    ax_graph.set_ylim(0,101) #since we are reading in percentages
    ax_graph.set_xlim(1, nstep) #will show graph up to however many time steps I have set
    if k == 0: #if you don't tell it to only put the lgend once it will add a new one each time step 
        ax_graph.legend(loc='upper right')

    plt.draw()       # Tells the global pyplot manager to process modifications at each time step
    plt.pause(0.5)   # Pauses execution briefly in between steps for visualization
plt.ioff()
plt.show()