"""
|  Author : Hemanth Pruthvi
|  File name : Material.py
|  Package : PyAstroPol
|  Description : Material class, for the sake of refractive index
"""
from .Functions import *

class Material():
    """
    |  Simple material definition, data is read from '.csv' material file.
    |  Attributes : General
    |  Wavelength : Float(N, 1) : Wavelengths
    |  RI : ComplexFloat(N, 1) : Refractive index at given wavelengths
    """
    def __init__ (self, MaterialName):
        try:
            FileName = '../Materials/' + str(MaterialName) + '.csv'
        except:
            raise ValueError('Invalid material name or material file!')
            return
        Array = np.loadtxt(FileName, delimiter=',')
        self.Wavelength = Array[:,0]
        self.RI = Array[:,1] - 1j*Array[:,2]
        return
    
    def getRefractiveIndexAt(self, Wave):
        """
        |  Calculate refractive index for a given wavelength by interpolation
        |  Input : Wavelength
        |  Returns : Complex refractive index
        """
        return np.interp(Wave, self.Wavelength, self.RI)
