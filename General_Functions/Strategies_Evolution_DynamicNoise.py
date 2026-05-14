# Strategies_Evolution_DynamicNoise - Contain strategies implementation: 
# - Evolution_Braitenberg : Apply Braitenberg strategy
# - Evolution_Triangulation : Apply Triangulation strategy
# both with Dynamic noise implemented

import numpy as np
from Strategies import * # Import "Braitenberg control" and "Triangulation" functions
from pyquaternion import Quaternion # Implement 3D rotations with quaternion

# Function to find the rotation matrix between two vectors in 3D
def Rotv1tov2(vec1, vec2):
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

    return rotation_matrix

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

    # Detected information
    GradientIntensity_s1 = np.linalg.norm(np.dot(FlowGradient(e3, xs1), n))
    GradientIntensity_s2 = np.linalg.norm(np.dot(FlowGradient(e3, xs2), n))
  
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
    
        
   # Dynamic Noise
    if Intensity != 0:
        DynamicNoise = np.random.normal(size=3, scale=1.)*np.sqrt(Intensity)
        new_t_noise = np.copy(new_t + DynamicNoise*np.sqrt(dt))    
        new_t_noise /= np.linalg.norm(new_t_noise)
        RR = Rotv1tov2(new_t, new_t_noise)
        if np.allclose(np.linalg.inv(RR), np.transpose(RR)):
            RotationNoise = Quaternion(matrix = RR)
            new_t = np.copy(new_t_noise)
            new_n = RotationNoise.rotate(new_n)
            new_b = RotationNoise.rotate(new_b)


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

    # Detected information along n
    GradientIntensity_s1 = np.linalg.norm(np.dot(FlowGradient(e3, xs1_n), n))
    GradientIntensity_s2 = np.linalg.norm(np.dot(FlowGradient(e3, xs2_n), n))
  
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
    
    # Detected information along b
    GradientIntensity_s1_2 = np.linalg.norm(np.dot(FlowGradient(e3, xs1_b), new_b))
    GradientIntensity_s2_2 = np.linalg.norm(np.dot(FlowGradient(e3, xs2_b), new_b))
  
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
    
    # Dynamic Noise
    if Intensity != 0:
        DynamicNoise = np.random.normal(size=3, scale=1.)*np.sqrt(Intensity)
        new_t_noise = np.copy(new_t + DynamicNoise*np.sqrt(dt))    
        new_t_noise /= np.linalg.norm(new_t_noise)
        RR = Rotv1tov2(new_t, new_t_noise)
        if np.allclose(np.linalg.inv(RR), np.transpose(RR)):
            RotationNoise = Quaternion(matrix = RR)
            new_t2 = np.copy(new_t_noise)
            new_n2 = RotationNoise.rotate(new_n2)
            new_b2 = RotationNoise.rotate(new_b2)
    

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

    # Detected information
    DirectionalGradient_s1 = np.dot(FlowGradient(e3, xs1), n)
    DirectionalGradient_s2 = np.dot(FlowGradient(e3, xs2), n)
    DirectionalGradient_s0 = np.dot(FlowGradient(e3, AgentPosition), n)
    
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
    
    # Dynamic Noise
    if Intensity != 0:
        DynamicNoise = np.random.normal(size=3, scale=1.)*np.sqrt(Intensity)
        new_t_noise = np.copy(new_t + DynamicNoise*np.sqrt(dt))    
        new_t_noise /= np.linalg.norm(new_t_noise)
        RR = Rotv1tov2(new_t, new_t_noise)
        if np.allclose(np.linalg.inv(RR), np.transpose(RR)):
            RotationNoise = Quaternion(matrix = RR)
            new_t = np.copy(new_t_noise)
            new_n = RotationNoise.rotate(new_n)
            new_b = RotationNoise.rotate(new_b)


    # Rotation around t 
    Rotation_t = Quaternion(axis= new_t, angle = twist * dt)
    new_n = Rotation_t.rotate(new_n)
    new_b = Rotation_t.rotate(new_b)
    
    new_n /= np.linalg.norm(new_n)
    new_b /= np.linalg.norm(new_b)
    new_t /= np.linalg.norm(new_t)
    
    
    return new_n, new_b, new_t
    
    
