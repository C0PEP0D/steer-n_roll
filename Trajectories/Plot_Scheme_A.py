# Plot scheme A, plot the scheme A in Fig.3 
#
#

import numpy as np
import matplotlib.pyplot as plt
from BasicData import *

fig = plt.figure()
ax = fig.add_subplot(111)

# Set agent positions
Thetas = np.linspace(0.2, np.pi/2-0.2, Nagents)
phi = np.zeros(Nagents)
AgentPosition = R0*np.array([np.sin(Thetas)*np.cos(phi), np.sin(Thetas)*np.sin(phi), np.cos(Thetas)]).T.reshape((Nagents,3))

# Plots
ax.scatter(0.0, 0.0, s=500, marker='o', color = 'limegreen', zorder = 3)
ax.scatter(AgentPosition[:,0], AgentPosition[:,2], s=50 ,marker='o', c='red')


ax.set_xlabel("")
ax.set_ylabel("")
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim([0, R0 + 5])
ax.set_ylim([0, R0 + 5])
ax.set_aspect('equal',adjustable='box')
plt.savefig("PlotPaper/Scheme_A_Fig3.png", bbox_inches = 'tight')

plt.show()
