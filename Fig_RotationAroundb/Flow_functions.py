# Flow_functions: 
# Calculation of flows, flow gradients, flow shear functions
#
#

# In all following functions:
# e: axisymmetric axis of the flow
# x: position where to calculate the flow, reference system centered on the source
import numpy as np

# Stokeslet
#
#

def FlowStokeslet(e, x):
    r = np.linalg.norm(x)
    return (1/r) * (e + np.dot(e,x) * x / r**2)

def StrainStokeslet(e,x):
    r = np.linalg.norm(x)
    return -(1/r**2) * ( -np.dot(x,e) * np.eye(3) / r + 3 * np.dot(x,e) * np.outer(x,x) / r**3)

def GradientStokeslet(e, x):
    r = np.linalg.norm(x) 
    gradflow = np.outer(e,x)/r**3 - np.outer(x,e)/r**3 -np.eye(3)*np.dot(x, e)/r**3 + 3*np.dot(x,e)*np.outer(x, x)/r**5  
    return gradflow  

#Stresslet
#
#

def FlowStresslet(e, x):
    r = np.linalg.norm(x)
    return (1/r**2) * (- x / r + 3 * np.dot(e,x)**2 * x / r**3)

def StrainStresslet(e, x):
    r = np.linalg.norm(x)
    g1 = np.outer(x,x) + np.dot(e,x)**2 * np.eye(3) + np.dot(e,x) * (np.outer(x,e) + np.outer(e,x))
    return -(1/r**3) * (np.eye(3) + 15 * np.dot(e,x)**2*np.outer(x,x)/r**4 - 3 * g1 / r**2)

def GradientStresslet(e, x):
    r = np.linalg.norm(x) 
    gradflow = + np.eye(3)/r**3 - 3*np.outer(x, x)/r**5 - 6*np.dot(x, e)*np.outer(x, e)/r**5 - 3*np.dot(x, e)**2*np.eye(3)/r**5 + 15*np.dot(x, e)**2*np.outer(x, x)/r**7    
    return -gradflow
    

# Quadruplet
#
#

def FlowQuadripole(e, x):
    g1 = 3 * np.dot(e,x) * x + np.dot(e,x)**2 * e
    r = np.linalg.norm(x)
    return (1/r**3) * (e - 3 * g1 / r**2 + 15 * np.dot(e,x)**3 * x / r**4)

def StrainQuadripole(e, x):
    r = np.linalg.norm(x)
    g2 = 3 * np.dot(e,x) * np.eye(3) + 2 * np.dot(e,x) * np.outer(e,e) + 2 * (np.outer(x,e) + np.outer(e,x))
    g3 = np.dot(e,x)**3 * np.eye(3) + 2 * np.dot(e,x)**2 * (np.outer(x,e) + np.outer(e,x)) + 3 * np.dot(e,x) * np.outer(x,x)
    return -(3/r**4) * (g2 / r - 5 * g3 / r**3 + 35 * np.dot(e,x)**3 * np.outer(x,x) / r**5)

def GradientQuadruplet(e, x):
    r = np.linalg.norm(x)     
    gradflow = -3*np.outer(e, x)/r**5 - 6*np.dot(x, e)*np.outer(e, e)/r**5 + 15*np.dot(x, e)**2*np.outer(e, x)/r**7 - 9*np.outer(x, e)/r**5 - 9*np.dot(x, e)*np.eye(3)/r**5 + 45*np.dot(x, e)**2*np.outer(x, e)/r**7 + 45*np.dot(x, e)*np.outer(x, x)/r**7 + 15*np.dot(e, x)**3*np.eye(3)/r**7 - 105*np.dot(x, e)**3*np.outer(x, x)/r**9   
    return gradflow

def TestFlow(e,x):
    return -x/np.linalg.norm(x)**2

