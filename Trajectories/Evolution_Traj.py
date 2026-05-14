import numpy as np
from BasicData import *
from BasicFunctions import *
from Geometric_functions import *

def Evolution(Grad, Strategy, W, twist):
	
	RMIN = 1
	TMAX = 20*np.pi*(R0/V)
	dt = 0.1
	e3 = np.array([0,0,1])
	AntennaeLength = .1
	BodyWidth = .1

	Trajectories_x = [[] for _ in range(NTrials)]
	Trajectories_y = [[] for _ in range(NTrials)]
	Trajectories_z = [[] for _ in range(NTrials)]
	Thetas = np.random.uniform(0,np.pi, 6)#np.array([0.5, 1.5])*np.pi/2
	Phis = np.random.uniform(0,2*np.pi, 6)#np.array([0.2, .4, .6, .8])*2*np.pi
	Successes = np.zeros(NTrials)
	for k in range(NTrials):
		if k < NTrials/2:
			theta = Thetas[0]
			phi = Phis[k%len(Phis)] + np.pi/2
		else:
			theta = Thetas[1]
			phi = Phis[k%len(Phis)]

		AgentPosition = R0*np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
		Begins = np.linalg.norm(np.copy(AgentPosition))
		alpha, beta, gamma = (2 * np.random.rand(3) - 1) * np.pi

		RotationMatrix = Rotation(alpha, np.abs(beta), gamma)
		n = np.dot(RotationMatrix,[1, 0, 0])
		b = np.dot(RotationMatrix,[0, 1, 0])
		t = np.dot(RotationMatrix,[0, 0, 1])
		time = 0
		print("I am here")
		while np.linalg.norm(AgentPosition) > RMIN and time < TMAX:

			n, b, t, AgentPosition = Strategy(n, b, t, AntennaeLength, AgentPosition, Grad, e3, dt, V, W, twist, 0) # Add BodyWidth for 4 sensors strategy
			print(np.linalg.norm(AgentPosition))

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
        	
