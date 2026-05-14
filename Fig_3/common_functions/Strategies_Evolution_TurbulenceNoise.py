# Strategies_Evolution_TurbulenceNoise - Contain strategies implementation: 
# - Evolution_Braitenberg : Apply Braitenberg strategy
# - Evolution_Triangulation : Apply Triangulation strategy
# both with Turbulence noise implemented

import numpy as np
from Strategies import * 
from pyquaternion import Quaternion # Implement 3D rotations with quaternion

# Function to calculate a pure shear flow, parameter gamma
def PureShear(gamma):
	return gamma*np.array([[np.sqrt(6)/6, 0, 0], [0, np.sqrt(6)/6, 0], [0, 0, -np.sqrt(6)/3]])
	
	
def Evolution_Braitenberg(n, b, t, AntennaeLength, AgentPosition, FlowGradient, e3, dt, twist, Intensity):

    # Calculate the re-orientation of axis due to active rotation and twist, Braitenberg strategy
    # Input:
    #    - {n, b, t} Orientation of the axis at instant t
    #    - AntennaeLength : distance between sensors
    #    - AgentPosition : current agent position
    #    - FlowGradient : information available to agent
    #    - e3 : source orientation
    #    - dt : integration timestep
    #    - twist : Imposed twist
    #
    #Output:
    #    - {n, b, t} at the instant t + dt
    #    - AgentPosition at the instant t + dt
    #
    
    # Set sensors position
    xs1 = AgentPosition + AntennaeLength*n
    xs2 = AgentPosition - AntennaeLength*n

    # Generate noise
    TurbulenceNoise_s1 = np.dot(PureShear(np.sqrt(Intensity)), n)
    TurbulenceNoise_s2 = TurbulenceNoise_s1

    # Detected information
    GradientIntensity_s1 = np.linalg.norm(np.dot(FlowGradient(e3, xs1), n)  + TurbulenceNoise_s1)
    GradientIntensity_s2 = np.linalg.norm(np.dot(FlowGradient(e3, xs2), n)  + TurbulenceNoise_s2)
  
    # Rotation around B
    AngularSpeed = Braitenberg_Control(GradientIntensity_s1, GradientIntensity_s2, AntennaeLength, b, dt)
    Norm = np.linalg.norm(AngularSpeed) 
    
    if Norm < 1e-14:
    	new_n = np.copy(n)
    	new_b = np.copy(b)
    	new_t = np.copy(t)
    else:
    	Rotation_b = Quaternion(axis = AngularSpeed, angle = Norm) 
    	new_n = Rotation_b.rotate(n)
    	new_b = Rotation_b.rotate(b)
    	new_t = Rotation_b.rotate(t)

    # Twist around t 
    Rotation_t = Quaternion(axis= new_t, angle = twist * dt)
    new_n2 = Rotation_t.rotate(new_n)
    new_b2 = Rotation_t.rotate(new_b)
    new_t2 = Rotation_t.rotate(new_t)
    
    new_n2 /= np.linalg.norm(new_n2)
    new_b2 /= np.linalg.norm(new_b2)
    new_t2 /= np.linalg.norm(new_t2)
    
    return new_n2, new_b2, new_t2


# Evolution Braitenberg with multiple sensors (4)  
#  
#

def Evolution_Intensity_MS(n, b, t, AntennaeLength, BodyWidth, AgentPosition, FlowGradient, e3, dt, twist, Intensity):

        # Calculate the re-orientation of axis due to active rotation and twist, Braitenberg multiple sensors
    # Input:
    #    - {n, b, t} Orientation of the axis at instant t
    #    - AntennaeLength : distance between sensors
    #    - BodyWidth : distance between for-aft sensors
    #    - AgentPosition : current agent position
    #    - FlowGradient : information available to agent
    #    - e3 : source orientation
    #    - dt : integration timestep
    #    - V : Linear Speed
    #    - twist : Imposed twist
    #
    #Output:
    #    - {n, b, t} at the instant t + dt
    #    - AgentPosition at the instant t + dt
    #
    
    # Set sensors position along n
    xs1_n = AgentPosition + AntennaeLength*n
    xs2_n = AgentPosition - AntennaeLength*n

    # Generate noise
    TurbulenceNoise_s1 = np.dot(PureShear(np.sqrt(Intensity)), n)
    TurbulenceNoise_s2 = TurbulenceNoise_s1

    # Detected information
    GradientIntensity_s1 = np.linalg.norm(np.dot(FlowGradient(e3, xs1), n)  + TurbulenceNoise_s1)
    GradientIntensity_s2 = np.linalg.norm(np.dot(FlowGradient(e3, xs2), n)  + TurbulenceNoise_s2)
  
    # Rotation around b 
    AngularSpeed = Braitenberg_Control(GradientIntensity_s1, GradientIntensity_s2, AntennaeLength, b, dt)
    Norm = np.linalg.norm(AngularSpeed) 
    
    if Norm < 1e-14:
    	new_n = np.copy(n)
    	new_b = np.copy(b)
    	new_t = np.copy(t)
    else:
    	Rotation_b = Quaternion(axis = AngularSpeed, angle = Norm) 
    	new_n = Rotation_b.rotate(n)
    	new_b = Rotation_b.rotate(b)
    	new_t = Rotation_b.rotate(t)
    
    # Set sensors position along b
    xs1_b = AgentPosition + BodyWidth*new_b
    xs2_b = AgentPosition - BodyWidth*new_b
    
    # Generate noise
    TurbulenceNoise_s1 = np.random.normal(size=3, scale=1)*np.sqrt(Intensity)
    TurbulenceNoise_s2 = np.random.normal(size=3, scale=1)*np.sqrt(Intensity)

    # Detected information
    GradientIntensity_s1 = np.linalg.norm(np.dot(FlowGradient(e3, xs1_b), b)  + TurbulenceNoise_s1)
    GradientIntensity_s2 = np.linalg.norm(np.dot(FlowGradient(e3, xs2_b), b)  + TurbulenceNoise_s2)
  
    # Rotation around n
    AngularSpeed = Braitenberg_Control(GradientIntensity_s1_2, GradientIntensity_s2_2, AntennaeLength, -new_n, dt) 
    Norm = np.linalg.norm(AngularSpeed)
    
    if Norm < 1e-14:
    	new_n2 = np.copy(new_n)
    	new_b2 = np.copy(new_b)
    	new_t2 = np.copy(new_t)
    else:
    	Rotation_n = Quaternion(axis = AngularSpeed, angle = Norm) 
    	new_n2 = Rotation_n.rotate(new_n)
    	new_b2 = Rotation_n.rotate(new_b)
    	new_t2 = Rotation_n.rotate(new_t)

    # Twist around t
    Rotation_t = Quaternion(axis= new_t2, angle = twist* dt)
    new_n3 = Rotation_t.rotate(new_n2)
    new_b3 = Rotation_t.rotate(new_b2)
    new_t3 = Rotation_t.rotate(new_t2)
    
    new_n3 /= np.linalg.norm(new_n3)
    new_b3 /= np.linalg.norm(new_b3)
    new_t3 /= np.linalg.norm(new_t3)
    
    return new_n3, new_b3, new_t3

# Evolution Triangulation
#
#
    
def Evolution_Triangulation(n, b, t, AntennaeLength, AgentPosition, FlowGradient, e3, dt, twist, Intensity):
    
    # Calculate the re-orientation of axis due to active rotation and twist, Triangulation strategy
    # Input:
    #    - {n, b, t} Orientation of the axis at instant t
    #    - AntennaeLength : distance between sensors
    #    - AgentPosition : current agent position
    #    - FlowGradient : information available to agent
    #    - e3 : source orientation
    #    - dt : integration timestep
    #    - V : Linear Speed
    #    - twist : Imposed twist
    #
    #Output:
    #    - {n, b, t} at the instant t + dt
    #    - AgentPosition at the instant t + dt
    #
    
    # Set sensors position 
    xs1 = AgentPosition + AntennaeLength*n
    xs2 = AgentPosition - AntennaeLength*n

    # Generate noise
    TurbulenceNoise_s0 = np.dot(PureShear(np.sqrt(Intensity)), n)
    TurbulenceNoise_s1 = TurbulenceNoise_s0
    TurbulenceNoise_s2 = TurbulenceNoise_s0

    # Detected information
    DirectionalGradient_s0 = np.dot(FlowGradient(e3, AgentPosition), n) + TurbulenceNoise_s0
    DirectionalGradient_s1 = np.dot(FlowGradient(e3, xs1), n) + TurbulenceNoise_s1
    DirectionalGradient_s2 = np.dot(FlowGradient(e3, xs2), n) + TurbulenceNoise_s2
    
    # Triangulation inferred direction d
    TargetDirection = Triangulation3DH(n, DirectionalGradient_s0, DirectionalGradient_s2-DirectionalGradient_s1)
   
    AngularSpeed = np.cross(t, TargetDirection)
    AngularSpeed /= np.linalg.norm(AngularSpeed)
    Angle = np.clip(np.dot(t, TargetDirection), a_min = -1, a_max = 1)
    Norm = np.arccos(Angle)*dt

    if Norm < 1e-14:
    	new_n = np.copy(n)
    	new_b = np.copy(b)
    	new_t = np.copy(t)
    else:	
    	Rotation = Quaternion(axis=AngularSpeed, angle=Norm) 
    	new_n = Rotation.rotate(n)
    	new_b = Rotation.rotate(b)
    	new_t = Rotation.rotate(t)


    # Rotation around t 
    Rotation_t = Quaternion(axis= new_t, angle = twist * dt)
    new_n = Rotation_t.rotate(new_n)
    new_b = Rotation_t.rotate(new_b)
    
    new_n /= np.linalg.norm(new_n)
    new_b /= np.linalg.norm(new_b)
    new_t /= np.linalg.norm(new_t)
    
    
    return new_n, new_b, new_t
    
    
