# 🧟🧟🧟 PART 1 🧟🧟🧟

'''
This file contains tools and scripts for completing Lab 1 for ClaSP 410.

This model simulates the spread of an infectious disease through a population.
Individuals can be Healthy, Infected, Immune, or Deceased.

To reproduce the plots in the lab report, run this script.
'''
'''
#Diagnostics Run
# Importing necessary packages
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

# 🧟🧟🧟 PART 2 🧟🧟🧟
# Defining initial variables and conditions
nx, ny = 7, 7
# Probability that a healthy person becomes infected when exposed
prob_infected = 1
# Probability that an infected person survives the disease
prob_survive = 0.0
# Probability that a healthy person is vaccinated/immune at the beginning
prob_immune = 0.0
# Probability that a person is initially healthy
prob_healthy = 1 - prob_immune
# Probability that an infected person dies
prob_fatal = 1 - prob_survive
# Identification numbers for each population state
Deceased = 0
Immune = 1
Healthy = 2
Infected = 3
# Number of time steps
nstep = 2


# 🧟🧟🧟 PART 3 🧟🧟🧟

# Create an initial grid.
# Everyone begins as healthy (2).

infec = np.zeros((nstep, nx, ny), dtype=int) + Healthy

print("Initial health of crowd:")
print(infec[0])
#here i am just checking that my array looks correct before continuing

#cet center cell to infected
infec[0,3,3] = Infected
# 🧟🧟🧟 PART 4 🧟🧟🧟

Create the initial population.

Each person has a probability of being vaccinated/immune.
Everyone else begins healthy.

We then introduce infected individuals randomly into
the population.

###DONT use the probs part for the basic debugging simulation

for i in range(nx):
    for j in range(ny):
        # First determine whether the person is immune/vaccinated
        if np.random.rand() < prob_immune:
            infec[0, j, i] = Immune
        else:
            # Otherwise, determine whether they begin infected
            if np.random.rand() < prob_infected:
                infec[0, j, i] = Infected
            else:
                infec[0, j, i] = Healthy

# Create immune "ghost nodes" around the outside of the grid.
# These cells prevent the disease from spreading outside the edges
infec[0, :, 0] = Immune
infec[0, 0, :] = Immune
infec[0, :, nx-1] = Immune
infec[0, ny-1, :] = Immune

print("Initial crowd health:")
print(infec[0])
#double checking to make sure edges are immune so the edges don't randomly spread illness

# 🧟🧟🧟 PART 5 🧟🧟🧟

# TIME TO ALLOW THE DISEASE TO SPREAD


At each time step, we check the four neighbors
If a healthy person is next to an infected person,
there is a probability that the healthy person becomes infected.
Infected individuals also have a probability of surviving
or dying from the disease. Surviving individuals become immune to the illness.
Immune individuals cannot become infected.
Deceased individuals cannot become infected.


# For making figures later
fig, (ax_map, ax_graph) = plt.subplots(1, 2, figsize=(12, 5))
from matplotlib.colors import ListedColormap
# Color map:
# deceased = black
# immune = gold
# healthy = light green
# infected = red
curr_infec_cmap = ListedColormap([
    'black',
    'goldenrod',
    'lightgreen',
    'firebrick'
])


# Store old data from our loop
history_time = []
history_healthy = []
history_infected = []
history_immune = []
history_deceased = []

# Start with the initial population
curr_infec = np.copy(infec[0, :, :])
pred_infec = np.copy(curr_infec)

# Begin time loop
for k in range(0, nstep):
    # Start the next time step with the current population
    pred_infec = np.copy(curr_infec)
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            if curr_infec[j, i] == Healthy:
                #check right neighbor
                if curr_infec[j, i+1] == Infected:
                    if np.random.rand() < prob_infected:
                        pred_infec[j, i] = Infected
                # Check left neighbor
                elif curr_infec[j, i-1] == Infected:
                    if np.random.rand() < prob_infected:
                        pred_infec[j, i] = Infected
                # Check upper neighbor
                elif curr_infec[j+1, i] == Infected:
                    if np.random.rand() < prob_infected:
                        pred_infec[j, i] = Infected
                # Check lower neighbor
                elif curr_infec[j-1, i] == Infected:
                    if np.random.rand() < prob_infected:
                        pred_infec[j, i] = Infected
                #if the cell started the time step as infected, then the next time step they should be deceased or healthy again
            #two scenarios/probabilities of what happens to infected people
            elif curr_infec[j, i] == Infected:
                if np.random.rand() > prob_survive:
                    pred_infec[j, i] = Deceased
    curr_infec = np.copy(pred_infec)

    # Cut out ghost nodes before calculating statistics
    active_view = curr_infec[1:-1, 1:-1]

    # Count each population state
    n_healthy = (active_view == Healthy).sum()
    n_infected = (active_view == Infected).sum()
    n_immune = (active_view == Immune).sum()
    n_deceased = (active_view == Deceased).sum()

    total_cells = active_view.size


    # Convert counts to percentages
    pct_healthy = (n_healthy / total_cells) * 100
    pct_infected = (n_infected / total_cells) * 100
    pct_immune = (n_immune / total_cells) * 100
    pct_deceased = (n_deceased / total_cells) * 100


    # Store history
    history_time.append(k + 1)
    history_healthy.append(pct_healthy)
    history_infected.append(pct_infected)
    history_immune.append(pct_immune)
    history_deceased.append(pct_deceased)


    # Print diagnostics every 10 time steps
    if (k + 1) % 10 == 0:
        print(
            f"Step {k+1}: "
            f"Healthy={pct_healthy:.1f}%, "
            f"Infected={pct_infected:.1f}%, "
            f"Immune={pct_immune:.1f}%, "
            f"Deceased={pct_deceased:.1f}%"
        )


    ax_map.clear()

    map_title_text = (
        f'Simulation of Infectious Disease Spread over Time'
        f'Iteration = {k+1:03d}\n'
        f'Healthy: {pct_healthy:.1f}% | '
        f'Infected: {pct_infected:.1f}% | '
        f'Immune: {pct_immune:.1f}% | '
        f'Deceased: {pct_deceased:.1f}%'
    )

    ax_map.set_title(map_title_text,fontsize=10,loc='left')
    ax_map.pcolor(active_view, cmap=curr_infec_cmap, vmin=0, vmax=3)
    ax_map.set_xlabel("X Position")
    ax_map.set_ylabel("Y Position")

    # GRAPH

    ax_graph.clear()
    ax_graph.plot(history_time,history_healthy,color='lightgreen',label='Healthy',linewidth=2)
    ax_graph.plot(history_time,history_infected,color='firebrick',label='Infected',linewidth=2)
    ax_graph.plot(history_time,history_immune,color='goldenrod',label='Immune',linewidth=2)
    ax_graph.plot(history_time,history_deceased,color='black',label='Deceased',linewidth=2)
    ax_graph.set_title("Spread of Infectious Disease Over Time",fontsize=10)
    ax_graph.set_xlabel("Time Step (Iteration)")
    ax_graph.set_ylabel("Population (%)")
    ax_graph.set_ylim(0, 101)
    ax_graph.set_xlim(1, nstep)
    ax_graph.legend(loc='upper right')
    plt.tight_layout()

    plt.draw()
    plt.pause(0.05)


plt.ioff()
plt.show()
'''





'''
#THE SIMULATION: THE ACTUAL FINAL PRODUCE!!!!
'''

#here, i am changing my initial variables to explore how initial immunity, spreading probability, and survival rates may change the outcome of our infectious disease simulation
# Importing necessary packages
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

# 🧟🧟🧟 PART 2 🧟🧟🧟
# Defining initial variables and conditions
nx, ny = 100, 100
# Probability that a healthy person becomes infected when exposed
prob_infected = 0.45
# Probability that an infected person survives the disease and becomes immun
prob_survive = 0.7
# Probability that a healthy person is vaccinated/immune at the beginning
prob_immune = 0.1
# Probability that a person is initially healthy
prob_healthy = 1 - prob_immune
# Probability that an infected person dies
prob_fatal = 1 - prob_survive
# Identification numbers for each population state
Deceased = 0
Immune = 1
Healthy = 2
Infected = 3
# Number of time steps
nstep = 25


# 🧟🧟🧟 PART 3 🧟🧟🧟

# Create an initial grid.
# Everyone begins as healthy (2).

infec = np.zeros((nstep, nx, ny), dtype=int) + 2

print("Initial health of crowd:")
print(infec[0])
#here i am just checking that my array looks correct before continuing

# 🧟🧟🧟 PART 4 🧟🧟🧟


#Create the initial population.

#Each person has a probability of being vaccinated/immune.
#Everyone else begins healthy.

#We then introduce infected individuals randomly into
#the population.


for i in range(nx):
    for j in range(ny):
        # First determine whether the person is immune/vaccinated
        if np.random.rand() < prob_immune:
            infec[0, j, i] = Immune
        else:
            # Otherwise, determine whether they begin infected
            if np.random.rand() < prob_infected:
                infec[0, j, i] = Infected
            else:
                infec[0, j, i] = Healthy
# Create immune "ghost nodes" around the outside of the grid.
# These cells prevent the disease from spreading outside the edges
infec[0, :, 0] = Immune
infec[0, 0, :] = Immune
infec[0, :, nx-1] = Immune
infec[0, ny-1, :] = Immune

print("Initial crowd health:")
print(infec[0])
#double checking to make sure edges are immune so the edges don't randomly spread illness

# 🧟🧟🧟 PART 5 🧟🧟🧟

# TIME TO ALLOW THE DISEASE TO SPREAD
'''
At each time step, we check the four neighbors
If a healthy person is next to an infected person,
there is a probability that the healthy person becomes infected.
Infected individuals also have a probability of surviving
or dying from the disease.
Immune individuals cannot become infected.
Deceased individuals cannot become infected.
'''
# For making figures later
fig, (ax_map, ax_graph) = plt.subplots(1, 2, figsize=(12, 5))
from matplotlib.colors import ListedColormap
# Color map:
# deceased = black
# immune = gold
# healthy = light green
# infected = red
curr_infec_cmap = ListedColormap([
    'black',
    'goldenrod',
    'lightgreen',
    'firebrick'
])


# Store old data from our loop
history_time = []
history_healthy = []
history_infected = []
history_immune = []
history_deceased = []

# Start with the initial population
curr_infec = np.copy(infec[0, :, :])
pred_infec = np.copy(curr_infec)

# Begin time loop
for k in range(0, nstep): #time steps
    # Start the next time step with the current population
    pred_infec = np.copy(curr_infec)
    for i in range(1, nx-1): #range in our box array
        for j in range(1, ny-1): #range in our box array
            if curr_infec[j, i] == Healthy:
                if curr_infec[j, i+1] == Infected:
                    if np.random.rand() < prob_infected:
                        pred_infec[j, i] = Infected
                # Check left neighbor
                elif curr_infec[j, i-1] == Infected:
                    if np.random.rand() < prob_infected:
                        pred_infec[j, i] = Infected
                # Check upper neighbor
                elif curr_infec[j+1, i] == Infected:
                    if np.random.rand() < prob_infected:
                        pred_infec[j, i] = Infected
                # Check lower neighbor
                elif curr_infec[j-1, i] == Infected:
                    if np.random.rand() < prob_infected:
                        pred_infec[j, i] = Infected
            #now, there are two possible outcomes for infected cells, and we must set probability statements for both scenarios
            elif curr_infec[j, i] == Infected:
                if np.random.rand() > prob_survive:
                    pred_infec[j, i] = Deceased
                if np.random.rand() < prob_survive:
                    pred_infec[j,i] = Immune
    curr_infec = np.copy(pred_infec)

    # Cut out ghost nodes before calculating statistics
    active_view = curr_infec[1:-1, 1:-1]

    # Count each population state
    n_healthy = (active_view == Healthy).sum()
    n_infected = (active_view == Infected).sum()
    n_immune = (active_view == Immune).sum()
    n_deceased = (active_view == Deceased).sum()

    total_cells = active_view.size


    # Convert counts to percentages
    pct_healthy = (n_healthy / total_cells) * 100
    pct_infected = (n_infected / total_cells) * 100
    pct_immune = (n_immune / total_cells) * 100
    pct_deceased = (n_deceased / total_cells) * 100


    # Store history
    history_time.append(k + 1)
    history_healthy.append(pct_healthy)
    history_infected.append(pct_infected)
    history_immune.append(pct_immune)
    history_deceased.append(pct_deceased)


    # Print diagnostics every 10 time steps
    if (k + 1) % 10 == 0:
        print(
            f"Step {k+1}: "
            f"Healthy={pct_healthy:.1f}%, "
            f"Infected={pct_infected:.1f}%, "
            f"Immune={pct_immune:.1f}%, "
            f"Deceased={pct_deceased:.1f}%"
        )


    ax_map.clear()

    map_title_text = (
        f'INCREASING SURVIVAL PROB: Simulation of Infectious Disease Spread over Time'
        f'Iteration = {k+1:03d}\n'
        f'Healthy: {pct_healthy:.1f}% | '
        f'Infected: {pct_infected:.1f}% | '
        f'Immune: {pct_immune:.1f}% | '
        f'Deceased: {pct_deceased:.1f}%'
    )

    ax_map.set_title(map_title_text,fontsize=12,loc='left')
    ax_map.pcolor(active_view, cmap=curr_infec_cmap, vmin=0, vmax=3)
    ax_map.set_xlabel("X Position")
    ax_map.set_ylabel("Y Position")

    # GRAPH

    ax_graph.clear()
    ax_graph.plot(history_time,history_healthy,color='lightgreen',label='Healthy',linewidth=2)
    ax_graph.plot(history_time,history_infected,color='firebrick',label='Infected',linewidth=2)
    ax_graph.plot(history_time,history_immune,color='goldenrod',label='Immune',linewidth=2)
    ax_graph.plot(history_time,history_deceased,color='black',label='Deceased',linewidth=2)
    ax_graph.set_title("Spread of Infectious Disease Over Time",fontsize=10)
    ax_graph.set_xlabel("Time Step (Iteration)")
    ax_graph.set_ylabel("Population (%)")
    ax_graph.set_ylim(0, 101)
    ax_graph.set_xlim(1, nstep)
    ax_graph.legend(loc='upper right')
    plt.tight_layout()

    plt.draw()
    plt.pause(0.05)


plt.ioff()
plt.show()
