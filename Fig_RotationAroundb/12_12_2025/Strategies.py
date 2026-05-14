# Strategies- implementation of angular speed calculation: 
# - Braitenberg_Control : Calculate angular rotation from Braitenberg strategy
# - Triangulation3DH : Calculate inferred direction axis "d", triangulation strategy 

import numpy as np

# Braitenberg_Control strategy
def Braitenberg_Control(Intensity_s1, Intensity_s2, AntennaeLength, RotAxis, dt):
    DeltaIntensity = Intensity_s1 - Intensity_s2
    RatioIntensity  = Intensity_s1 + Intensity_s2
    Speed =  DeltaIntensity/RatioIntensity*dt*RotAxis/AntennaeLength
    return Speed

# Triangulation strategy, non-linear Sign dynamic
def Triangulation3DH(n, d0, diff):
    n_ = np.cross(d0, diff)
    d = np.cross(n, n_)
    sign = np.sign(np.dot(d0, d))
    target = sign * d0

    return target / np.linalg.norm(target)

# Triangulation strategy, non-linear Sign dynamic
def Triangulation3DHNoise(n, d0, diff, Intensity):
    n_ = np.cross(d0, diff)
    d = np.cross(n, n_)
    sign = np.sign(np.dot(d0, d) + np.random.normal(0, Intensity))
    target = sign * d0

    return target / np.linalg.norm(target)


# Triangulation strategy, skew-symmetric nearest point based strategy
# Find out more at https://en.wikipedia.org/wiki/Skew_lines
def Triangulation3D(p1, d1, p2, d2):
    # Line directions d1, d2
    # Starting points p1, p2
    if np.allclose(d1, d2,rtol=1e-8) :
        return d1
    
    n = np.cross(d1,d2)
    n1 = np.cross(d1,n)
    n2 = np.cross(d2,n)
    
    if abs(np.dot(d1,d2)/(np.linalg.norm(d1)*np.linalg.norm(d2))) > 1 - 1e-8:
        return d1
    
    c1 = p1 + ((np.dot((p2-p1),n2)/np.dot(d1,n2)))*d1
    c2 = p2 + ((np.dot((p1-p2),n1)/np.dot(d2,n1)))*d2
    
    IntersectionPoint = 0.5*(c1+c2)
    TargetDirection = -(p1 + p2)/2 + IntersectionPoint
    return TargetDirection/np.linalg.norm(TargetDirection)

# Simplified triangulation
def TriangulationSimple(AgentPosition, gradient):
    signed_gradient = np.sign(np.dot(gradient, - AgentPosition)) * gradient
    return signed_gradient / np.linalg.norm(signed_gradient)



