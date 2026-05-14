# BasicData - Import basic data needed :
# - R0 : Euclidean distance from target
# - NTrials : number of agents initialized

import numpy as np
from .Flow_functions import *  # Import flow functions: stokeslet, grads, strains ...

Nagents = 4  # NTrials : number of agents initialized
R0 = 30  # Euclidean distance from target
V = 0.7  # Linear speed of the agents

# Determine domain of motion of the agent
RMIN = 1  # Minimal distance with respect to prey
TMAX = 50 * np.pi * (R0)  # Maximum time allowed for the simulation
dt = 0.1  # Time steps for time integration
e3 = np.array([0, 0, 1])  # Orientation of symmetry axis of source flow
AntennaeLength = 0.001  # Distance between sensors
BodyWidth = 0.1  #

# Define information flow
Grads = [GradientStresslet]  # ShearStresslet, GradientStresslet, ShearStokeslet,
GradName = [
    "GradientStresslet"
]  # 'ShearStresslet', 'GradientStresslet', 'ShearStokeslet',

# Twist setting
twists = [0.0, 0.3]  # Imposed twist
