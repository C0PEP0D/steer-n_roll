# EvolutionTraj : Numerical integration in time for angular - linear dynamics

import numpy as np

from .BasicData import *
from .Geometric_functions import Rotation


def Evolution(Grad, Strategy, twist, Positions):
    # Input:
    #   Grad : information detected by the agent, flow gradient or flow strain
    #   Streategy : strategy performed by the agent, Braitenberg or Triangulation
    #   twist : Imposed twist
    #
    # Outcome:
    #   Trajectories_i: Array with positions along axis i in the frame fixed on the target
    #   Success: Array containing 1s, 0s for success/failures of the strategy

    Nagents = Positions.shape[0]

    # Initialize trajectory vectors
    Trajectories_x = [[] for _ in range(Nagents)]
    Trajectories_y = [[] for _ in range(Nagents)]
    Trajectories_z = [[] for _ in range(Nagents)]

    # Initialize success array
    Successes = np.zeros(Nagents)

    # Time integrations
    for a in range(Nagents):

        AgentPosition = Positions[a]
        Begins = np.linalg.norm(np.copy(AgentPosition))

        # Randomly initialize agent's axes {n, b, t} orientation
        RotationMatrix = Rotation(-0.9, -0.6, -0.2)
        n = np.dot(RotationMatrix, [1, 0, 0])
        b = np.dot(RotationMatrix, [0, 1, 0])
        t = np.dot(RotationMatrix, [0, 0, 1])

        time = 0

        while np.linalg.norm(AgentPosition) > 2 * RMIN and time < TMAX:
            time += 1
            n, b, t, AgentPosition = Strategy(
                n, b, t, AntennaeLength, AgentPosition, Grad, e3, dt, V, twist
            )

            Trajectories_x[a].append(AgentPosition[0])
            Trajectories_y[a].append(AgentPosition[1])
            Trajectories_z[a].append(AgentPosition[2])

        Ends = np.linalg.norm(np.copy(AgentPosition))

        if Ends > Begins:
            Successes[a] = 0
        else:
            Successes[a] = 1

    return Trajectories_x, Trajectories_y, Trajectories_z, Successes
