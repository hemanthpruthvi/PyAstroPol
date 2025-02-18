"""
|  Author : Hemanth Pruthvi
|  File name : Functions.py
|  Package : PyAstroPol
|  Description : Frequently used functions in other codes
"""
import os
import copy as cp
import random as rd
from datetime import datetime as dt
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
#
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
EPS = 1e-30             # Error
Z0 = 376.7303134126     # Vaccum impedance

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
    
def getAngle(V1, V2):
    """
    |  get angle between the two sets of vectors
    |  Input : Array of vectors V1
    |          Aray of vectors V2
    |  Returns :  Array of angles
    """
    DOT = np.sum(V1*V2, axis=1)
    CROSS = np.linalg.norm(np.cross(V1, V2), axis=1)
    Theta = np.arctan2(CROSS, DOT).reshape((V1.shape[0], 1))
    return Theta
 
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

def getXRotationMatrices(ThetaX):
    """
    |  Calculate 3x3 Local Euler's rotation matrix about X-axis
    |  Input : array of X rotation angles  
    |  Returns : 3x3 matrix with array as elements 
    """
    ThetaX = np.radians(ThetaX)
    R = np.matrix([[1.0, 0.0, 0.0],
                   [0.0, np.cos(ThetaX), -np.sin(ThetaX)], 
                   [0.0, np.sin(ThetaX),  np.cos(ThetaX)]], dtype=object)
    return R

def getYRotationMatrices(ThetaY):
    """
    |  Calculate 3x3 Local Euler's rotation matrix about Y-axis
    |  Input : array of Y rotation angles  
    |  Returns : 3x3 matrix with array as elements 
    """
    ThetaY = np.radians(ThetaY)
    R = np.matrix([[ np.cos(ThetaY), 0.0, np.sin(ThetaY)],
                   [0.0, 1.0, 0.0],
                   [-np.sin(ThetaY), 0.0, np.cos(ThetaY)]], dtype=object)
    return R

def getZRotationMatrices(ThetaZ):
    """
    |  Calculate 3x3 Local Euler's rotation matrix about Z-axis
    |  Input : array of Z rotation angles  
    |  Returns : 3x3 matrix with array as elements 
    """
    ThetaZ = np.radians(ThetaZ)
    R = np.matrix([[np.cos(ThetaZ), -np.sin(ThetaZ), 0.0],
                   [np.sin(ThetaZ),  np.cos(ThetaZ), 0.0],
                   [0.0, 0.0, 1.0]], dtype=object)
    return R

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

def getBerremanFieldMatrixIsotropic(RI, Beta):
    """
    |   Compute 4x4 Berreman field matrix for isotropic material
    |   Input:  Scalar refractive index
    |           Snell's law quantity n*sin(theta)
    |   Output: 4x4 Berreman field matrix
    """
    Alpha = np.sqrt(RI**2-Beta**2)
    gammaP = RI**2/(Alpha*Z0)
    gammaS = -Alpha/Z0
    #
    FieldMatrix = np.matrix(np.zeros([4,4]), dtype=object)
    FieldMatrix[0,0] = FieldMatrix[0,1] = FieldMatrix[2,2] = FieldMatrix[2,3] = 1+0*Beta
    FieldMatrix[0,2] = FieldMatrix[0,3] = FieldMatrix[1,2] = FieldMatrix[1,3] = 0*Beta
    FieldMatrix[2,0] = FieldMatrix[2,1] = FieldMatrix[3,0] = FieldMatrix[3,1] = 0*Beta
    FieldMatrix[1,0] = gammaS
    FieldMatrix[1,1] = -gammaS
    FieldMatrix[3,2] = gammaP
    FieldMatrix[3,3] = -gammaP
    # FieldMatrix = np.matrix([[1,        1,          0,          0],
    #                          [gammaS,   -gammaS,    0,          0],
    #                          [0,        0,          1,          1],
    #                          [0,        0,          gammaP,     -gammaP]])

    return FieldMatrix, Alpha

def getPermittivityTensor(RI, Angles):
    """
    |   Compute relative permittivity tensor for a general trirefringent material
    |   Input:  Refractive index vector i.e., RI along principal axes (1x3)
    |           Euler's angles for local rotation about Z, Y and Z axes in that order (Nx3) 
    |           here, X: direction of propagation, XY: plane of polarization, Z: axis perpendicular to the plane of the wave  
    |   Output: Relative permittivity/refractive index tensor (3x3)
    """
    nx, ny, nz = RI
    Eta, Psi, Xi = Angles
    EtaMatrix = getZRotationMatrices(Eta)
    PsiMatrix = getYRotationMatrices(Psi)
    XiMatrix = getZRotationMatrices(Xi)
    EtaMatrix_ = getZRotationMatrices(-Eta)
    PsiMatrix_ = getYRotationMatrices(-Psi)
    XiMatrix_ = getZRotationMatrices(-Xi)
    #
    Epsilon = np.matrix(np.diag([nx**2, ny**2, nz**2]))
    Epsilon = EtaMatrix*Epsilon*EtaMatrix_
    Epsilon = PsiMatrix*Epsilon*PsiMatrix_
    Epsilon = XiMatrix*Epsilon*XiMatrix_
    return Epsilon

def getBerremanFieldMatrixAnisotropic(Epsilon, Beta):
    """
    |   Compute 4x4 Berreman field matrix for anisotropic material
    |   Input:  Refractive index tensor
    |           Snell's law quantity n*sin(theta)
    |   Output: 4x4 Berreman field matrix
    |           1x4 alpha vector (propagation constants)
    """
    Epsilon_ = np.array(Epsilon)
    Z0=377
    exx, exy, exz = Epsilon_[0]
    eyx, eyy, eyz = Epsilon_[1]
    ezx, ezy, ezz = Epsilon_[2]
    AM = np.array([[-Beta*ezx/ezz,          Z0*(1-Beta**2/ezz),     -Beta*ezy/ezz,                  0.0*Beta],
                            [(exx-ezx**2/ezz)/Z0,    -Beta*ezx/ezz,          (exy-ezx*ezy/ezz)/Z0,           0.0*Beta],          
                            [0*Beta,                      0*Beta,                      0*Beta,                              0*Beta-Z0],
                            [(-exy+ezy*ezx/ezz)/Z0,  Beta*ezy/ezz,           (Beta**2-eyy+ezy**2/ezz)/Z0,    0.0*Beta]])


    FM, A = np.zeros([4,4,len(Beta)]), np.zeros([4,len(Beta)])
    for i in range(len(Beta)):
        A[:,i], FM[:,:,i] = np.linalg.eig(AM[:,:,i])
    #
    Alpha = np.array(A)
    FieldMatrix = np.matrix(np.zeros([4,4]), dtype=object)
    for i in range(4):
        for j in range(4):
            FieldMatrix[i,j] = np.array(FM[i,j,:])
    return FieldMatrix, Alpha

def getBerremanPhaseMatrix(Alpha, Thick):
    """
    |   Compute 4x4 Berreman phase matrix for a given layer 
    |   Input:  4xN (1xN) alpha vector (scalar) for N rays
    |           Thickness of the layer
    |   Output: 4x4 Berreman phase matrix
    """
    Phase = np.exp(-1j*2*np.pi*Thick*Alpha)
    PhaseMatrix = np.matrix(np.eye(4), dtype=object)
    for i in range(4):
        if len(Alpha.shape)==2 : 
            PhaseMatrix[i,i] = Phase[i,:]
            Zeros = np.zeros([Alpha.shape[1],1])
        elif len(Alpha.shape)==1 : 
            PhaseMatrix[i,i] = Phase
            Zeros = np.zeros([len(Alpha),1])
        else : 
            print('Error in computing the phase matrix, check the propagation!')
    for i in range(4):
        for j in range(4):
            PhaseMatrix[i,j] += Zeros
    return PhaseMatrix

def getInverseMatrix(Mat):
    """
    |  Compute inverse of a 4x4 matrix where each of the elements in a 1-d array
    """
    Mat_ = np.zeros([4,4,len(Mat[0,0])])
    for i in range(4):
        for j in range(4):
            Mat_[i,j] = Mat[i,j].flatten()
    IMat_ = 0.0*Mat_
    for i in range(Mat_.shape[2]):
        IMat_[:,:,i] = np.linalg.inv(Mat_[:,:,i])
    IMat = np.matrix(np.eye(4), dtype=object)
    for i in range(4):
        for j in range(4):
            IMat[i,j] = IMat_[i,j]
    return IMat

def getSellmeierRefractiveIndex(Waves, Coeffs):
    """
    |  Generate material files from Sellmeier coefficients for isotropic material
    |  Input: Name of the material
    |         Sellmeier coefficients
    |  
    """
    B1, C1, B2, C2, B3, C3 = Coeffs
    W = 1/Waves**2
    RI = 1 + B1/(1-C1*W) + B2/(1-C2*W) + B3/(1-C3*W)
    RI = np.sqrt(RI+0j)
    return RI

def flattenMatrix(Mat):
    """
    |  Flatten the matrix elements that are supposed to be 1-d arrays
    """
    for i in range(4):
        for j in range(4):
            Mat[i,j] = Mat[i,j].flatten()
    return Mat