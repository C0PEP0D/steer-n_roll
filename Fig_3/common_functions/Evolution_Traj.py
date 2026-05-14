# Evolution Trajectories
# Integrate numerically the agent displacements
#
#

import numpy as np

from .BasicData import *  # Import basic data needed in simulations
from .Flow_functions import *  # Import all flow functions: gradients, strains ...
from .Geometric_functions import *  # Import functions for rotations, projection, linear applications ...


def Evolution(Grad, Strategy, twist):

    # Evolution function
    # input:
    # - Grad : gradient of the flow
    # - Strategy : stragegy used by the agent [ Braitenberg, Triangulation ] ...
    # - twist : Imposed twist speed
    #
    # output:
    # - Trajectories position as arrays x, y, z
    # - Success: Booleand value 0 (target failed) , 1 (target reached)
    #

    RMIN = 1  # Minimal distance from the source allowed in simualtions
    TMAX = 3000  # Max time allowed in simulations
    dt = 0.1  # Timesteps
    e3 = np.array([0, 0, 1])  # Symmetric axis of the flow
    AntennaeLength = 0.1  # Distance between sensors
    # BodyWidth = .1 # DIstance between sensors along agent body

    # Initialize positions along axis x, y, z
    Trajectories_x = [[] for _ in range(NTrials)]
    Trajectories_y = [[] for _ in range(NTrials)]
    Trajectories_z = [[] for _ in range(NTrials)]

    # Initial position of agents
    Thetas = np.random.uniform(0, np.pi, 6)  # np.array([0.5, 1.5])*np.pi/2
    Phis = np.random.uniform(0, 2 * np.pi, 6)  # np.array([0.2, .4, .6, .8])*2*np.pi
    Successes = np.zeros(NTrials)

    # Numerical integration
    for k in range(NTrials):
        if k < NTrials / 2:
            theta = Thetas[0]
            phi = Phis[k % len(Phis)] + np.pi / 2
        else:
            theta = Thetas[1]
            phi = Phis[k % len(Phis)]

        AgentPosition = R0 * np.array(
            [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
        )
        Begins = np.linalg.norm(np.copy(AgentPosition))
        alpha, beta, gamma = (2 * np.random.rand(3) - 1) * np.pi

        RotationMatrix = Rotation(alpha, np.abs(beta), gamma)
        n = np.dot(RotationMatrix, [1, 0, 0])
        b = np.dot(RotationMatrix, [0, 1, 0])
        t = np.dot(RotationMatrix, [0, 0, 1])
        time = 0

        while np.linalg.norm(AgentPosition) > RMIN and time < TMAX:

            n, b, t, AgentPosition = Strategy(
                n, b, t, AntennaeLength, AgentPosition, Grad, e3, dt, V, twist, 0
            )  # Add BodyWidth for 4 sensors strategy

            Trajectories_x[k].append(AgentPosition[0])
            Trajectories_y[k].append(AgentPosition[1])
            Trajectories_z[k].append(AgentPosition[2])
            time += 1
        Ends = np.linalg.norm(np.copy(AgentPosition))

        if Ends > Begins:
            Successes[k] = 0
        else:
            Successes[k] = 1

    return Trajectories_x, Trajectories_y, Trajectories_z, Successes
