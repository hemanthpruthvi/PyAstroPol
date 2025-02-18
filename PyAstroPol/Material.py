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
    def __init__ (self, MaterialName, IsIsotropic=True):
        CurDir = os.path.dirname(__file__)
        try:
            FileName = CurDir + '/../Materials/' + str(MaterialName) + '.csv'
        except:
            raise ValueError(r'%s : Invalid material name or material file!'%(FileName))
            return
        Array = np.loadtxt(FileName, delimiter=',')
        self.Label = MaterialName
        self.Wavelength = Array[:,0]
        # Isotropic material
        if (Array.shape[1] == 3):
            self.IsIsotropic = True
            self.IsUniaxial = False
            self.IsBiaxial = False
            self.RI = Array[:,1] - 1j*Array[:,2]
            self.RIo = self.RIe = self.RI1 = self.RI2 = self.RI3 = self.RI
        # Uniaxial material
        elif (Array.shape[1] == 5):
            self.IsIsotropic = False
            self.IsUniaxial = True
            self.IsBiaxial = False
            self.RIo = Array[:,1] - 1j*Array[:,2]            
            self.RIe = Array[:,3] - 1j*Array[:,4]
            self.RI1 = self.RI2 = self.RI = self.RIo
            self.RI3 = self.RIe
        # Biaxial material
        elif (Array.shape[1] == 7):
            self.IsIsotropic = False
            self.IsUniaxial = False
            self.IsBiaxial = True
            self.RI1 = Array[:,1] - 1j*Array[:,2]            
            self.RI2 = Array[:,3] - 1j*Array[:,4]
            self.RI3 = Array[:,5] - 1j*Array[:,6]
            self.RI = self.RI1
        # Isotropic material, file with Sellmeier coefficients
        elif (len(Array.shape) == 1 and Array.shape[0] == 6):
            self.IsIsotropic = True
            self.IsUniaxial = False
            self.IsBiaxial = False
            self.Wavelength = np.arange(0.2, 2.201, 0.001)
            self.RI = getSellmeierRefractiveIndex(self.Wavelength, Array)
            self.RIo = self.RIe = self.RI1 = self.RI2 = self.RI3 = self.RI
        # Uniaxial material, file with Sellmeier coefficients
        elif (Array.shape[1] == 6 and Array.shape[0] == 2):
            self.IsIsotropic = False
            self.IsUniaxial = True
            self.IsBiaxial = False
            self.Wavelength = np.arange(0.2, 2.201, 0.001)
            self.RIo = getSellmeierRefractiveIndex(self.Wavelength, Array[0,:])
            self.RIe = getSellmeierRefractiveIndex(self.Wavelength, Array[1,:])
            self.RI1 = self.RI2 = self.RI = self.RIo
            self.RI3 = self.RIe
        # Biaxial material, file with Sellmeier coefficients
        elif (Array.shape[1] == 6 and Array.shape[0] == 3):
            self.IsIsotropic = False
            self.IsUniaxial = False
            self.IsBiaxial = True
            self.RI1 = getSellmeierRefractiveIndex(self.Wavelength, Array[0,:])            
            self.RI2 = getSellmeierRefractiveIndex(self.Wavelength, Array[1,:])
            self.RI3 = getSellmeierRefractiveIndex(self.Wavelength, Array[2,:])
            self.RI = self.RI1
        else:
            print(r'%s : Invalid material file, please check its contents!'%(FileName))
        return
    
    def getRefractiveIndexAt(self, Wave):
        """
        |  Calculate refractive index for a given wavelength by interpolation
        |  Returns array depending on the type of the material -- scalar for isotropic, 1x2 for birefringent, 1x3 for trirefringent  
        |  Input : Wavelength
        |  Returns : Complex refractive index
        """
        # if (self.IsIsotropic):
        #     return np.interp(Wave, self.Wavelength, self.RI)
        # if (self.IsUniaxial):
        #     return [np.interp(Wave, self.Wavelength, self.RIo), 
        #             np.interp(Wave, self.Wavelength, self.RIe)]
        # if (self.IsBiaxial):
        #     return [np.interp(Wave, self.Wavelength, self.RI1), 
        #             np.interp(Wave, self.Wavelength, self.RI2),
        #             np.interp(Wave, self.Wavelength, self.RI3)]
        return np.interp(Wave, self.Wavelength, self.RI)


    def getRefractiveIndicesAt(self, Wave):
        """
        |  Calculate refractive indices for a given wavelength by interpolation
        |  Always returns 1x3 array of RI along the principal axes
        |  Input : Wavelength
        |  Returns : 1x3 array of complex refractive indices
        """
        if (self.IsIsotropic):
            RI = np.interp(Wave, self.Wavelength, self.RI)
            return [RI, RI, RI]
        if (self.IsUniaxial):
            RIo = np.interp(Wave, self.Wavelength, self.RIo)
            RIe = np.interp(Wave, self.Wavelength, self.RIe)
            return [RIo, RIo, RIe]
        if (self.IsBiaxial):
            RI1 = np.interp(Wave, self.Wavelength, self.RI1)
            RI2 = np.interp(Wave, self.Wavelength, self.RI2)
            RI3 = np.interp(Wave, self.Wavelength, self.RI3)
            return [RI1, RI2, RI3]
