"""
Author : Hemanth Pruthvi
File name : Material.py
Package : PyAstroPol
Description : Material class, for the sake of refractive index
"""

import numpy as np
import copy as cp
import random as rd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime as dt
#
from .Functions import *
#

class Material():
    """
    Simple material definition, material file is required
    
    // Attributes : General
    Wavelength : Float(N, 1) : Wavelengths
    RI : ComplexFloat(N, 1) : Refractive index at given wavelengths
    """
    def __init__ (self, MaterialName):
        FileName = '../Materials/' + str(MaterialName) + '.csv'
        try:
            Array = np.loadtxt(FileName, delimiter=',')
        except:
            print('Invalid material name or file! Check both!')
            return
        self.Wavelength = Array[:,0]
        self.RI = Array[:,1] - 1j*Array[:,2]
        return
    
    # Interpolate RI for a given wavelength
    def getRefractiveIndexAt(self, Wave):
        return np.interp(Wave, self.Wavelength, self.RI)
