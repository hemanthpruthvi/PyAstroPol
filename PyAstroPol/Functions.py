"""
|  Author : Hemanth Pruthvi
|  File name : Functions.py
|  Package : PyAstroPol
|  Description : Frequently used functions in other codes
"""
import os
import numpy as np
import copy as cp
import random as rd
import matplotlib
# matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime as dt
EPS = 1e-30 

def normalize3DVectors(V):
    """
    |  Normalize array of 3d vectors to get corresponding unit vectors.
    |  Inputs : Array of 3D vectors.
    |  Returns : Array of 3D vectors same size as input.
    """
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

def dot3DVectors(V1, V2):
    """
    |  Dot product of two arrays of 3D vectors, formatted to match code's convention.
    |  Inputs : Array of 3d vectors, Array of 3d vectors.
    |  Returns : Array of 3d vectors same size as input.
    """
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
 
def adjustAspect(Ax, L, x=0.0, y=0.0, z=0.0):
    """
    |  Set aspect ratio of 3d axes to 1:1:1, and specify the space to be displayed.
    |  Inputs : Pyplot Axis, Size of space, X-coordinate of the center, Y-coordinate of the center, Z-coordinate of the center.
    """
    Ax.set_xlim([x-L/2.0, x+L/2.0])
    Ax.set_ylim([y-L/2.0, y+L/2.0])
    Ax.set_zlim([z-L/2.0, z+L/2.0])
    return

def vectorToAngles(V):
    """
    |  Rotation angles about Y and X axes (in that order) to rotate "Z" vector to "V" vector.
    |  Input : Array of vectors.
    |  Returns : Array of Y-rotation angles, Array of X-rotation angles in degrees.
    """
    ThetaY = np.arcsin(V[0])
    ThetaX = np.arctan2(-V[1]/np.cos(ThetaY), V[2]/np.cos(ThetaY))
    return np.degrees(ThetaY), np.degrees(ThetaX)

def getXRotationMatrix(ThetaX):
    """
    |  Calculate 4x4 affine rotation matrix, to rotate about X-axis.
    |  Inputs : X rotation angle in degrees.
    |  Returns : 4x4 transformation matrix.
    """
    ThetaX = np.radians(ThetaX)
    R = np.matrix([[1.0, 0.0, 0.0, 0.0],
                   [0.0, np.cos(ThetaX), -np.sin(ThetaX), 0.0], 
                   [0.0, np.sin(ThetaX),  np.cos(ThetaX), 0.0],
                   [0.0, 0.0, 0.0, 1.0]])
    return R

def getYRotationMatrix(ThetaY):
    """
    |  Calculate 4x4 affine rotation matrix, to rotate about Y-axis.
    |  Inputs : Y rotation angle in degrees.
    |  Returns : 4x4 transformation matrix.
    """
    ThetaY = np.radians(ThetaY)
    R = np.matrix([[ np.cos(ThetaY), 0.0, np.sin(ThetaY), 0.0],
                   [0.0, 1.0, 0.0, 0.0],
                   [-np.sin(ThetaY), 0.0, np.cos(ThetaY), 0.0],
                   [0.0, 0.0, 0.0, 1.0]])
    return R

def getZRotationMatrix(ThetaZ):
    """
    |  Calculate 4x4 affine rotation matrix, to rotate about Z-axis.
    |  Inputs : Z rotation angle in degrees.
    |  Returns : 4x4 transformation matrix.
    """
    ThetaZ = np.radians(ThetaZ)
    R = np.matrix([[np.cos(ThetaZ), -np.sin(ThetaZ), 0.0, 0.0],
                   [np.sin(ThetaZ),  np.cos(ThetaZ), 0.0, 0.0],
                   [0.0, 0.0, 1.0, 0.0],
                   [0.0, 0.0, 0.0, 1.0]])
    return R

def getTranslationMatrix(x=0.0, y=0.0, z=0.0):
    """
    |  Calculate 4x4 affine translation matrix, to move w.r.t. present position.
    |  Inputs : X translation, Y translation and Z translation.
    |  Returns : 4x4 transformation matrix.
    """
    T = np.matrix([[1.0, 0.0, 0.0, x],
                   [0.0, 1.0, 0.0, y],
                   [0.0, 0.0, 1.0, z],
                   [0.0, 0.0, 0.0, 1.0]])
    return T

def applyPointTransformation(P, M):
    """
    |  Apply transformation matrix to set of points.
    |  Inputs : Array of points, and 4x4 transformation matrix.
    |  Returns : Array of points, same size as input.
    """
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

def applyVectorTransformation(V, M):
    """
    |  Apply transformation matrix to set of vectors.
    |  Inputs : Array of vectors, and 4x4 transformation matrix.
    |  Returns : Array of vectors, same size as input.
    """
    R, T = M[0:3,0:3], M[0:3,3]
    # for singe point
    if(len(V.shape) == 1):
        Out = R*np.transpose(np.matrix(V))
        return np.array(Out).flatten()
    # for set of points
    if (V.shape[1] == 3):
        Out = R*np.transpose(np.matrix(V))
        return np.array(np.transpose(Out))
    else :
        Out = R*V
        return np.array(Out)

def MuellerRotationMatrix(Theta):
    """
    |  Compute Mueller matrix for rotation.
    |  Input : Angle in degrees.
    |  Returns : 4x4 Mueller matrix.
    """
    Theta = np.radians(Theta)
    R = np.matrix([[1.0,  0.0,             0.0,             0.0], 
                   [0.0,  np.cos(2*Theta), np.sin(2*Theta), 0.0], 
                   [0.0, -np.sin(2*Theta), np.cos(2*Theta), 0.0], 
                   [0.0,  0.0,             0.0,             1.0]])
    return R

def rotateMuellerMatrix(M, Theta):
    """
    |  Compute Mueller matrix after rotation.
    |  Input : 4x4 Mueller matrix, Rotation angle in degrees.
    |  Returns : 4x4 Mueller matrix.
    """
    Theta = np.radians(Theta)
    R = np.matrix([[1.0,  0.0,             0.0,             0.0], 
                   [0.0,  np.cos(2*Theta), np.sin(2*Theta), 0.0], 
                   [0.0, -np.sin(2*Theta), np.cos(2*Theta), 0.0], 
                   [0.0,  0.0,             0.0,             1.0]])
    return R*M

# Display digits in iPython
def roundOffDisplay(D):
    """
    |  Set the number of decimals to display while printing.
    |  Input : Integer number of decimals.
    """
    floatParam = '{: 0.' + str(D) + 'f}'
    # print(floatParam)
    np.set_printoptions(formatter={'float':floatParam.format})
    return

def formatMaterialFile(Name):
    """
    |  Format the material file '.csv' that is downloaded from refractiveindex.info, to comply with the code. The file should be copied into "Materials" folder.
    |  Input : Material name without extension.
    """
    FileName = '../Materials/' + str(Name) + '.csv'
    f = open(FileName, 'r')
    w, n, k = [], [], []
    f.readlines(1)

    # Replace missing values with -1000
    for i in f:
        temp = i.split(',')
        w.append(float(temp[0]))
        try:
            n.append(float(temp[1]))
        except:
            n.append(-1000)
        try:
            k.append(float(temp[2]))
        except:
            k.append(-1000)   
    f.close()

    # Find indices of missing values
    w, n, k = np.array(w), np.array(n), np.array(k)
    n_ok = np.argwhere(n != -1000).flatten()
    k_ok = np.argwhere(k != -1000).flatten()
    n_nan = np.argwhere(n == -1000).flatten()
    k_nan = np.argwhere(k == -1000).flatten()

    # Interpolate missing values using existing
    w_proper, n_proper = w[n_ok], n[n_ok]
    for i in n_nan:
        n[i] = np.interp(w[i], w_proper, n_proper)
    w_proper, k_proper = w[k_ok], k[k_ok]
    for i in k_nan:
        k[i] = np.interp(w[i], w_proper, k_proper)

    # Get string data ready to be written to file
    Data = ''
    for i in range(len(w)):
        Data += str(w[i])
        Data += ','
        Data += str(n[i])
        Data += ','
        Data += str(k[i])
        Data += '\n'
    # Write the data to the same file
    f = open(FileName, 'w')
    f.write(Data)
    f.close()
    return