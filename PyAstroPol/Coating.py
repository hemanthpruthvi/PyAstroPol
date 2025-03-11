"""
|  Author : Hemanth Pruthvi
|  File name : Coating.py
|  Package : PyAstroPol
|  Description : Coating and related classes
"""

# import numpy as np
# import copy as cp
# import random as rd
# from matplotlib import pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from datetime import datetime as dt
#
from .Functions import *
from .Source import *
from .Surface import *
from .Material import *
#

class Coating():
    """
    |  Effects of multi-layer coating on polarization of transmission and reflections.
    |  For geometric ray tracing purposes, their thickness will be neglected. 
    
    |  Attributes : General
    |  Layers   : String array (N,1)    : Array of names of materials of all the layers
    |  RI       : Complex Float (N,1)   : Array of refractive indiced of all the layers
    |  Thick    : Float (N,1)           : Array of thicknesses of all the layers in microns
    
    |  Attributes : Polarization
    |  rs : FloatComplex(N,1) : Complex coefficient of reflection for s-polarization of all Rays 
    |  rp : FloatComplex(N,1) : Complex coefficient of reflection for p-polarization of all Rays 
    |  ts : FloatComplex(N,1) : Complex coefficient of transmission for s-polarization of all Rays 
    |  tp : FloatComplex(N,1) : Complex coefficient of transmission for p-polarization of all Rays
    """
    def __init__(self, Layers, Thick, Orientation=None):
        if (len(Layers) != len(Thick)):
            print('Error! Check the inputs!')
            return
        self.Layers = np.array(Layers)
        # self.RI = np.copy(self.Layers)
        self.Thick = np.array(Thick)
        self.Wavelength = 0.6328
        self.Orientation = Orientation
        self.loadRefractiveIndex()
        return
    
    def fromFile(self, FileName):
        """
        |  Load the coating from file; the file should contain generic coating data in a certain format
        |  The file should have N rows and 5 columns. Column data should correspond to
        |  Material Name    Layer Thickness     Angle Eta       Angle Psi       Angle Xi
        """
        Data = np.genfromtxt(FileName, delimiter=',', dtype=str)
        self.Layers = Data[:,0]
        self.Wavelength = 0.6328
        self.Stack = np.zeros([Data.shape[0],7])
        self.Stack[:,0] = np.float64(Data[:,1])
        self.Stack[:,4::] = np.float64(Data[:,2::])
        self.loadRefractiveIndex()
    
    def applyToSurface(self, Surf):
        """ 
        |  Apply defined Coating to given surface, and compute the effects on transmission and reflection.
        |  Input : Surface() to which coating it to be applied. 
        """
        self.Wavelength = np.copy(Surf.Wavelength)
        self.loadRefractiveIndex()
        self.iRI, self.sRI = np.copy(Surf.iRI), np.copy(Surf.tRI)
        self.iTheta, self.sTheta = np.copy(Surf.iTheta), np.copy(Surf.tTheta)
        
        # for s-polarization
        iEta, sEta = self.iRI*np.cos(self.iTheta), self.sRI*np.cos(self.sTheta)
        CoatMatrix = np.array([[1,0],[0,1]])
        for RI_j, Thick_j in zip(self.RI, self.Thick):
            # Angle and phase in each coating media
            Theta_j = np.arcsin(self.iRI*np.sin(self.iTheta)/RI_j)
            Delta_j = 2*np.pi*RI_j*Thick_j*np.cos(Theta_j)/self.Wavelength
            Eta_j = RI_j*np.cos(Theta_j)
            # 
            LayerMatrix =  np.array([[np.cos(Delta_j), 1j*np.sin(Delta_j)/Eta_j], 
                                    [1j*Eta_j*np.sin(Delta_j), np.cos(Delta_j)]])
            CoatMatrix = np.transpose(np.matmul(np.transpose(CoatMatrix), np.transpose(LayerMatrix)))
        # Coefficients    
        E = CoatMatrix[0,0] + CoatMatrix[0,1]*sEta
        H = CoatMatrix[1,0] + CoatMatrix[1,1]*sEta
        # self.sTransferMatrix = np.copy(CoatMatrix)
        self.rs, self.ts = (iEta*E - H)/(iEta*E + H), 2*iEta/(iEta*E + H)
        
        # for p-polarization
        iEta, sEta = self.iRI/np.cos(self.iTheta), self.sRI/np.cos(self.sTheta)
        CoatMatrix = np.array([[1,0],[0,1]])
        for RI_j, Thick_j in zip(self.RI, self.Thick):
            # Angle and phase in each coating media
            Theta_j = np.arcsin(self.iRI*np.sin(self.iTheta)/RI_j)
            Delta_j = 2*np.pi*RI_j*Thick_j*np.cos(Theta_j)/self.Wavelength
            Eta_j = RI_j/np.cos(Theta_j)
            #
            LayerMatrix =  np.array([[np.cos(Delta_j), 1j*np.sin(Delta_j)/Eta_j], 
                                    [1j*Eta_j*np.sin(Delta_j), np.cos(Delta_j)]])
            CoatMatrix = np.transpose(np.matmul(np.transpose(CoatMatrix), np.transpose(LayerMatrix)))
        # Coefficients    
        E = CoatMatrix[0,0] + CoatMatrix[0,1]*sEta
        H = CoatMatrix[1,0] + CoatMatrix[1,1]*sEta
        # self.pTransferMatrix = np.copy(CoatMatrix)
        self.rp, self.tp= -(iEta*E - H)/(iEta*E + H), 2*iEta/(iEta*E + H)

        # To consider interference, compute Berreman matrices for further use
        if (Surf.Interference):
            CharMatrix = np.eye(4)
            for l, s in self.Layers, self.Stack:
                if (Material(l).IsIsotropic):
                    FM_ = getBerremanFieldMatrixIsotropic(s[1], Surf.Beta)
                    Alpha_ = np.sqrt(Surf.RI**2-Surf.Beta**2)
                else:
                    Thick, nx, ny, nz, Eta, Psi, Xi = s
                    Eta = 0.0
                    Psi, Xi = vectorToAngles(self.sCosines)
                    Epsilon_ = getPermittivityTensor([nx,ny,nz], [Eta, Psi, Xi])
                    FM_, Alpha_ = getBerremanFieldMatrixAnisotropic(Epsilon_, Surf.Beta)
                PM_ = getBerremanPhaseMatrix(Alpha_, Thick)
                CharMatrix = CharMatrix*FM_*PM_*np.linalg.inv(FM_)
            self.CharMatrix = CharMatrix
        return
    
    def reverseTheCoating(self):
        """
        |  Reverse the coating layers order along with their thicknesses.
        |  This feature might be useful for elements such as coated lens.
        """
        self.Layers = self.Layers[::-1]
        self.RI = self.RI[::-1]
        self.Thick = self.Thick[::-1]
        return
        
    def loadRefractiveIndex(self):
        """
        |  Load refractive index according to given information on layer material and wavelength. 
        """
        self.RI = []
        for i, l in enumerate(self.Layers):
            try:
                RIs = Material(l).getRefractiveIndicesAt(self.Wavelength)
                self.RI.append(RIs[0])
                self.Stack[i,1:4] = RIs
            except:
                self.RI.append(complex(l))
        self.RI = np.array(self.RI)
        return
