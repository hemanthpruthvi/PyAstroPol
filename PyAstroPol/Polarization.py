import numpy as np
import copy as cp
import random as rd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime as dt
#

# Normalize array of 3d vectors to get unit vectors
def normalize3DVectors(V):
    Vs = V.shape
    # for singe vectors
    if(len(Vs) == 1):
        V = V/np.linalg.norm(V)
        return V
    # for set of vectors
    if (Vs[0] == 3):
        V = V/np.reshape(np.linalg.norm(V, axis=0), newshape=(1,Vs[1]))
    else :
        V = V/np.reshape(np.linalg.norm(V, axis=1), newshape=(Vs[0],1))
    return V

# Obtain Dot product of two arrays of 3d vectors
def dot3DVectors(V1, V2):
    if (V1.shape == V2.shape):
        Vs = V1.shape
        # for single vector
        if (len(Vs) == 1):
            return np.sum(V1*V2)
        # for set of vectors
        if (Vs[0] == 3):
            V1DotV2 = np.reshape(np.sum(V1*V2, axis=0), newshape=(1, Vs[1]))
        else :
            V1DotV2 = np.reshape(np.sum(V1*V2, axis=1), newshape=(Vs[0], 1))
        return V1DotV2
    else:
        print('Error! dim(V1) != dim(V2)')
        return

# Set aspect ratio of 3d axes to 1:1:1
def adjustAspect(Ax, L, x=0.0, y=0.0, z=0.0):
    Ax.set_xlim([x-L/2.0, x+L/2.0])
    Ax.set_ylim([y-L/2.0, y+L/2.0])
    Ax.set_zlim([z-L/2.0, z+L/2.0])
    return

# Polarimeter
def getSystemMuellerMatrix(Source, Surfaces):
    """
    Polarimetry to compute Mueller matrix
    In_1 : [1,  1, 0, 0],  Ex_1 = 1.0 + 0.0j,  Ey_1 = 0.0 + 0.0j
    In_2 : [1, -1, 0, 0],  Ex_1 = 0.0 + 0.0j,  Ey_1 = 1.0 + 0.0j
    In_3 : [1,  0, 1, 0],  Ex_1 = 0.707 + 0.0j,  Ey_1 = 0.707 + 0.0j
    In_4 : [1,  0, 0, 1],  Ex_1 = 0.707 + 0.0j,  Ey_1 = 0.0 + 0.707j
    #
    Out_1 : [M_00+M_10, M_01+M_11, M_02+M_12, M_03+M_13]
    Out_2 : [M_00-M_10, M_01-M_11, M_02-M_12, M_03-M_13]
    Out_3 : [M_00+M_20, M_01+M_21, M_02+M_22, M_03+M_23]
    Out_4 : [M_00+M_30, M_01+M_31, M_02+M_32, M_03+M_33]
    #
    [M_00, M_01, M_02, M_03] = 0.5*(Out_1 + Out_2)
    [M_10, M_11, M_12, M_13] = 0.5*(Out_1 - Out_2)
    [M_20, M_21, M_22, M_23] = Out_3 - 0.5*(Out_1 + Out_2)
    [M_30, M_31, M_32, M_33] = Out_4 - 0.5*(Out_1 + Out_2)
    
    """
    Ex_In = [1.0+0.0j, 0.0+0.0j, np.sqrt(0.5)+0.0j, np.sqrt(0.5)+0.0j]
    Ey_In = [0.0+0.0j, 1.0+0.0j, np.sqrt(0.5)+0.0j, 0.0+np.sqrt(0.5)*1j]
    OutRays = []
    for Ex, Ey in zip(Ex_In, Ey_In):
        Source.createPolarization(Ex, Ey)
        Surfaces[0].propagateRays(Source)
        for i in range(len(Surfaces)-1):
            if (Surfaces[i].Mirror):
                Surfaces[i+1].propagateRays(Surfaces[i].rRays)
            else:
                Surfaces[i+1].propagateRays(Surfaces[i].tRays)
        if(Surfaces[-1].Mirror):
            OutRays.append(Surfaces[-1].rRays)
        else:
            OutRays.append(Surfaces[-1].tRays)        
    Out = []
    M = []
    for O in OutRays:
        ExTemp, EyTemp = np.sum(O.Ex), np.sum(O.Ey)
        Out.append(ExTemp*np.conjugate(ExTemp)+EyTemp*np.conjugate(EyTemp))
        Out.append(ExTemp*np.conjugate(ExTemp)-EyTemp*np.conjugate(EyTemp))
        Out.append(ExTemp*np.conjugate(EyTemp)+EyTemp*np.conjugate(ExTemp))
        Out.append(1j*(EyTemp*np.conjugate(ExTemp)-ExTemp*np.conjugate(EyTemp)))
    Out = np.real(np.array(Out)).reshape((4,4))
    M = np.copy(Out)
    M[0,:] = 0.5*(Out[0,:] + Out[1,:])
    M[1,:] = 0.5*(Out[0,:] - Out[1,:])
    M[2,:] = Out[2,:] - 0.5*(Out[0,:] + Out[1,:])
    M[3,:] = Out[3,:] - 0.5*(Out[0,:] + Out[1,:])
    MNorm = M/M[0,0]
    return MNorm, M[0,0] 