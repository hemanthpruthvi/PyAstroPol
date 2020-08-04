import numpy as np
import copy as cp
import random as rd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime as dt
EPS = 1e-30 
#

# Normalize array of 3d vectors to get unit vectors
def normalize3DVectors(V):
    Vs = V.shape
    # for singe vectors
    if(len(Vs) == 1):
        V = V/np.linalg.norm(V)
        return V
    # for set of vectors
    if (Vs[1] == 3):
        V = V/np.reshape(np.linalg.norm(V, axis=1), newshape=(Vs[0],1))
    else :
        V = V/np.reshape(np.linalg.norm(V, axis=0), newshape=(1,Vs[1]))
    return V

# Obtain Dot product of two arrays of 3d vectors
def dot3DVectors(V1, V2):
    if (V1.shape == V2.shape):
        Vs = V1.shape
        # for single vector
        if (len(Vs) == 1):
            return np.sum(V1*V2)
        # for set of vectors
        if (Vs[1] == 3):
            V1DotV2 = np.reshape(np.sum(V1*V2, axis=1), newshape=(Vs[0], 1))
        else :
            V1DotV2 = np.reshape(np.sum(V1*V2, axis=0), newshape=(1, Vs[1]))
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

# Rotation angles about Y and X axes (in that order) to rotate "Z" axis to V 
def vectorToAngles(V):
    ThetaY = np.arcsin(V[0])
    ThetaX = np.arctan2(-V[1]/np.cos(ThetaY), V[2]/np.cos(ThetaY))
    return np.degrees(ThetaY), np.degrees(ThetaX)

#  X-rotation matrix
def rotateAboutX(ThetaX):
    ThetaX = np.radians(ThetaX)
    R = np.matrix([[1.0, 0.0, 0.0, 0.0],
                    [0.0, np.cos(ThetaX), -np.sin(ThetaX), 0.0], 
                    [0.0, np.sin(ThetaX),  np.cos(ThetaX), 0.0],
                    [0.0, 0.0, 0.0, 1.0]])
    return R

# Y-rotation matrix
def rotateAboutY(ThetaY):
    ThetaY = np.radians(ThetaY)
    R = np.matrix([[ np.cos(ThetaY), 0.0, np.sin(ThetaY), 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [-np.sin(ThetaY), 0.0, np.cos(ThetaY), 0.0],
                    [0.0, 0.0, 0.0, 1.0]])
    return R

#  Z-rotation matrix
def rotateAboutZ(ThetaZ):
    ThetaZ = np.radians(ThetaZ)
    R = np.matrix([[np.cos(ThetaZ), -np.sin(ThetaZ), 0.0, 0.0],
                    [np.sin(ThetaZ),  np.cos(ThetaZ), 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0]])
    return R

# Translation matrix
def translateOrigin(x=0.0, y=0.0, z=0.0):
    T = np.matrix([[1.0, 0.0, 0.0, x],
                    [0.0, 1.0, 0.0, y],
                    [0.0, 0.0, 1.0, z],
                    [0.0, 0.0, 0.0, 1.0]])
    return T

# Apply transformation matrix to set of points
def applyPointTransformation(P, M):
    R, T = M[0:3,0:3], M[0:3,3]
    # for singe point
    if(len(P.shape) == 1):
        Out = R*np.transpose(np.matrix(P)) + T
        return np.array(Out).flatten()
    # for set of points
    if (P.shape[1] == 3):
        Out = R*np.transpose(np.matrix(P)) + T
        return np.array(np.transpose(Out))
    else :
        Out = R*P + T
        return np.array(Out)

# Apply transformation matrix to set of vectors
def applyVectorTransformation(P, M):
    R, T = M[0:3,0:3], M[0:3,3]
    # for singe point
    if(len(P.shape) == 1):
        Out = R*np.transpose(np.matrix(P))
        return np.array(Out).flatten()
    # for set of points
    if (P.shape[1] == 3):
        Out = R*np.transpose(np.matrix(P))
        return np.array(np.transpose(Out))
    else :
        Out = R*P
        return np.array(Out)

# Rotate Mueller matrix
def rotateMueller(M, Theta):
    Theta = np.radians(Theta)
    R = np.matrix([[1.0,  0.0,             0.0,             0.0], 
                   [0.0,  np.cos(2*Theta), np.sin(2*Theta), 0.0], 
                   [0.0, -np.sin(2*Theta), np.cos(2*Theta), 0.0], 
                   [0.0,  0.0,             0.0,             1.0]])
    return R*M