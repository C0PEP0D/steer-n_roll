# Geometric_functions - implementation of rotations: 
#
#

import numpy as np

# Rotation matrix in 3D, Eulerian frame
def Rotation(alpha, beta, gamma):
    c1 = np.cos(alpha)
    s1 = np.sin(alpha)
    c2 = np.cos(beta)
    s2 = np.sin(beta)
    c3 = np.cos(gamma)
    s3 = np.sin(gamma)
    
    return np.array([
        [c1*c3 - c2*s1*s3, -c1*s3 - c2*c3*s1, s1*s2],
        [c3*s1 + c1*c2*s3, c1*c2*c3 - s1*s3, -c1*s2],
        [s2*s3,            c3*s2,             c2]
    ])

