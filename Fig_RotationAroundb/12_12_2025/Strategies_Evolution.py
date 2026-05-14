# Strategies_Evolution - Contain strategies implementation: 
# - Evolution_Braitenberg : Apply Braitenberg strategy
# - Evolution_Triangulation : Apply Triangulation strategy

from Strategies import *  # Import "Braitenberg control" and "Triangulation" functions
import numpy as np
from pyquaternion import Quaternion

# Evolution Braitenberg
#
#

def Evolution_Braitenberg(n, b, t, AntennaeLength, AgentPosition, FlowGradient, e3, dt, V, twist):
    # Calculate the re-orientation of axis due to active rotation and twist, Braitenberg strategy
    # Input:
    #    - {n, b, t} Orientation of the axis at instant T
    #    - AntennaeLength : distance between sensors
    #    - AgentPosition : current agent position
    #    - FlowGradient : information available to agent
    #    - e3 : source orientation
    #    - dt : integration timestep
    #    - V : Linear Speed
    #    - twist : Imposed twist
    #
    # Output:
    #    - {n, b, t} at the instant t + dt
    #    - AgentPosition at the instant t + dt
    #

    # Set sensors position
    xs1 = AgentPosition + AntennaeLength * n
    xs2 = AgentPosition - AntennaeLength * n

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
        Rotation_b = Quaternion(axis=AngularSpeed, angle=Norm)
        new_n = Rotation_b.rotate(n)
        new_b = Rotation_b.rotate(b)
        new_t = Rotation_b.rotate(t)

    # Twist around t
    Rotation_t = Quaternion(axis=new_t, angle=twist * dt)
    new_n2 = Rotation_t.rotate(new_n)
    new_b2 = Rotation_t.rotate(new_b)
    new_t2 = Rotation_t.rotate(new_t)

    new_n2 /= np.linalg.norm(new_n2)
    new_b2 /= np.linalg.norm(new_b2)
    new_t2 /= np.linalg.norm(new_t2)

    # Calculate new Agent position
    AgentPosition = AgentPosition + V * dt * new_t2

    return new_n2, new_b2, new_t2, AgentPosition

# Evolution Triangulation
#
#

def Evolution_Triangulation(n, b, t, AntennaeLength, AgentPosition, FlowGradient, e3, dt, V, twist):
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
    # Output:
    #    - {n, b, t} at the instant t + dt
    #    - AgentPosition at the instant t + dt
    #

    # Set sensors position
    xs1 = AgentPosition + AntennaeLength * n
    xs2 = AgentPosition - AntennaeLength * n

    # Detected information
    DirectionalGradient_s1 = np.dot(FlowGradient(e3, xs1), n)
    DirectionalGradient_s2 = np.dot(FlowGradient(e3, xs2), n)

    # Triangulation inferred direction d
    TargetDirection = Triangulation3D(xs1, DirectionalGradient_s1, xs2, DirectionalGradient_s2)
    AngularSpeed = dt * np.cross(t, TargetDirection)
    Norm = np.linalg.norm(AngularSpeed)

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
    Rotation_t = Quaternion(axis=new_t, angle=twist * dt)
    new_n = Rotation_t.rotate(new_n)
    new_b = Rotation_t.rotate(new_b)

    new_n /= np.linalg.norm(new_n)
    new_b /= np.linalg.norm(new_b)
    new_t /= np.linalg.norm(new_t)

    # Calculate new Agent position
    AgentPosition = AgentPosition + dt * new_t * V

    return new_n, new_b, new_t, AgentPosition

# Evolution Braitenberg with multiple sensors (4)
#
#

def Evolution_Intensity_MS(n, b, t, AntennaeLength, BodyWidth, AgentPosition, FlowGradient, e3, dt, V, twist=0):
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
    # Output:
    #    - {n, b, t} at the instant t + dt
    #    - AgentPosition at the instant t + dt
    #

    # Set sensors position along n
    xs1_n = AgentPosition + AntennaeLength * n
    xs2_n = AgentPosition - AntennaeLength * n

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
        Rotation_b = Quaternion(axis=AngularSpeed, angle=Norm)
        new_n = Rotation_b.rotate(n)
        new_b = Rotation_b.rotate(b)
        new_t = Rotation_b.rotate(t)

    # Set sensors position along b
    xs1_b = AgentPosition + BodyWidth * new_b
    xs2_b = AgentPosition - BodyWidth * new_b

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
        Rotation_n = Quaternion(axis=AngularSpeed, angle=Norm)
        new_b2 = Rotation_n.rotate(new_b)
        new_t2 = Rotation_n.rotate(new_t)

    # Twist around t
    Rotation_t = Quaternion(axis=new_t2, angle=twist * dt)
    new_n3 = Rotation_t.rotate(new_n)
    new_b3 = Rotation_t.rotate(new_b2)
    new_t3 = Rotation_t.rotate(new_t2)

    new_n3 /= np.linalg.norm(new_n3)
    new_b3 /= np.linalg.norm(new_b3)
    new_t3 /= np.linalg.norm(new_t3)

    # Calculate new Agent position
    AgentPosition = AgentPosition + dt * V * new_t3

    return new_n3, new_b3, new_t3, AgentPosition


def Evolution_Triangulation_1D(n, b, t, AntennaeLength, AgentPosition, FlowGradient, e3, dt, V, twist):
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
    # Output:
    #    - {n, b, t} at the instant t + dt
    #    - AgentPosition at the instant t + dt
    #

    # Set sensors position
    xs1 = AgentPosition + AntennaeLength * n
    xs2 = AgentPosition - AntennaeLength * n

    # Potentially detected information
    DirectionalGradient_s1 = np.dot(FlowGradient(e3, xs1), n)
    DirectionalGradient_s2 = np.dot(FlowGradient(e3, xs2), n)

    # Detected information, mimiking coplanar motion
    DirectionalGradient_s1 = np.array([DirectionalGradient_s1[0], DirectionalGradient_s1[1], 0])
    DirectionalGradient_s2 = np.array([DirectionalGradient_s2[0], DirectionalGradient_s2[1], 0])
    # Triangulation inferred direction d
    TargetDirection = Triangulation3D(xs1, DirectionalGradient_s1, xs2, DirectionalGradient_s2)
    AngularSpeed = dt * np.cross(t, TargetDirection)
    Norm = np.linalg.norm(AngularSpeed)

    if Norm < 1e-14:
        new_n = np.copy(n)
        new_b = np.copy(b)
        new_t = np.copy(t)
    else:
        Rotation = Quaternion(axis=b, angle=Norm) # the agent can only rotate around b
        new_n = Rotation.rotate(n)
        new_b = Rotation.rotate(b)
        new_t = Rotation.rotate(t)

    # Rotation around t, implemented via the roll
    Rotation_t = Quaternion(axis=new_t, angle=twist * dt)
    new_n = Rotation_t.rotate(new_n)
    new_b = Rotation_t.rotate(new_b)

    new_n /= np.linalg.norm(new_n)
    new_b /= np.linalg.norm(new_b)
    new_t /= np.linalg.norm(new_t)

    # Calculate new Agent position
    AgentPosition = AgentPosition + dt * new_t * V

    return new_n, new_b, new_t, AgentPosition
