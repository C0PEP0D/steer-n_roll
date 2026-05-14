# Code to Evolve and plot the trajectories in 3D
# It calculates evolution for both strategies: Braitenberg and Triangulation
#

from Time_Int_Traj import * # Import evolution function
from BasicData import *  # Nagents, R0, V
from Strategies import * # Import strategies, Braitenberg or Triangulation
from Strategies_Evolution import * # Import strategies Evolution, Braitenberg or Triangulation
from Flow_functions import * # Import flow functions: stokeslet, grads, strains ... 

from pyquaternion import Quaternion
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Realize Evolution
Strategies = [Evolution_Braitenberg]  # Evolution_Braitenberg ,Evolution_Triangulation
StrategyName = "B" # "B", "T"

for strat_i, Strategy in enumerate(Strategies):
	for grad_i, Grad in enumerate(Grads):
		for twist in twists:
			Trajectories_x, Trajectories_y, Trajectories_z, Successes = Evolution(Grad, Strategy, twist)

			fig = plt.figure(figsize=(6, 6))
			ax = fig.add_subplot(111, projection='3d')

			

			for a in range(Nagents):
				ax.plot(Trajectories_x[a], Trajectories_y[a], Trajectories_z[a], color='cadetblue', alpha=.7, linewidth= 3)	
				ax.scatter(Trajectories_x[a][0], Trajectories_y[a][0], Trajectories_z[a][0], color='red', s= 50)	
			ax.scatter(0,0,0, s=100, marker='o', color = 'limegreen')#, label='Target')
			ax.set_xlim([-2, R0])
			ax.set_ylim([-2, 2])
			ax.set_zlim([0., R0])
			xticks = np.linspace(-2, R0, 7)
			yticks = np.linspace(-2, 2, 4)
			zticks = np.linspace(-2, R0, 7)

			empty_labels_x = ["" for i in range(7)]
			empty_labels_y = ["" for i in range(4)]
			empty_labels_z = ["" for i in range(7)]

			ax.set_xticks(xticks, empty_labels_x)
			ax.set_yticks(yticks, empty_labels_y)
			ax.set_zticks(zticks, empty_labels_z)
		

			ax.set_box_aspect(aspect = (1, .2, 1))
			ax.view_init(15, 305)
			# ax.set_title(StrategyName[strat_i])
			
			# plt.savefig("PlotPaper/Traj"+str(twist)+"_"+StrategyName+".png")
			plt.show()

