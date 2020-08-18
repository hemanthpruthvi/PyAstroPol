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
#

class Coating():
    """
    |  Effects of multi-layer coating on polarization of transmission and reflections.
    |  For geometric ray tracing purposes, their thickness will be neglected. 
    
    |  Attributes : General
    |  Layers   : String array (N,1)    : Array of names of materials of all the layers
    |  RI       : Complex Float (N,1)   : Array of refractive indiced of all the layers
    |  Thick    : Float (N,1)           : Array of thicknesses of all the layers
    
    |  Attributes : Polarization
    |  rs : FloatComplex(N,1) : Complex coefficient of reflection for s-polarization of all Rays 
    |  rp : FloatComplex(N,1) : Complex coefficient of reflection for p-polarization of all Rays 
    |  ts : FloatComplex(N,1) : Complex coefficient of transmission for s-polarization of all Rays 
    |  tp : FloatComplex(N,1) : Complex coefficient of transmission for p-polarization of all Rays
    """
    def __init__(self, RI, Thick):
        if (len(RI) != len(Thick)):
            print('Error! Check the inputs!')
            return
        self.Layers = np.array(RI)
        # self.RI = np.copy(self.Layers)
        self.Thick = np.array(Thick)
        self.Wavelength = 0.6328
        self.iRI = 1.0 + 0.0j
        self.sRI = 1.5 + 0.0j
        self.iTheta = 1.5 + 0.0j
        self.sTheta = 1.5 + 0.0j
        self.sTransferMatrix = np.array([[1,0],[0,1]])
        self.pTransferMatrix = np.array([[1,0],[0,1]])
        self.loadRefractiveIndex()
        return
    
    
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
        for m in self.Layers:
            try:
                self.RI.append(Material(m).getRefractiveIndexAt(self.Wavelength))
            except:
                self.RI.append(np.complex(m))
        self.RI = np.array(self.RI)
        return