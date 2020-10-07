"""
|  Author : Hemanth Pruthvi
|  File name : System.py
|  Package : PyAstroPol
|  Description : Optical system class
"""

from .Source import *
from .Surface import *
from .Coating import *
from .Elements import *
from .Functions import *
#


class System():
    """
    |  System class it provided to do all the propagation related activites as once.
    |  Once the Source, Elements and Detector are defined, they are to be passed to this calss.
    |  There a provision to provide a seperate Source for the sake of display.
    |  In this way, analysis can be done with large number of rays and fewer rays can be used for display.

    |  Attributes: Objects
    |  Source : Source() : Source for analysis purposes
    |  Components : N x Surface() : Surface/Lens/Component
    |  Detector : Detector() : Detector
    |  dRays : Source() : Source for display puposes
    
    |  Attributes : General
    |  MuellerMatrix : Float(4,4) : Mueller matrix of the whole system
    |  Transmission : Float : Total fractional throughput of the system.
    """
    def __init__(self, Source, Components, Detector, dRays=None):
        self.Source = cp.deepcopy(Source)
        self.Components = cp.deepcopy(Components)
        self.Detector = cp.deepcopy(Detector)
        if (dRays == None) : 
            self.DisplayRays = cp.deepcopy(Source)
        else :
            self.DisplayRays = cp.deepcopy(dRays)
        return
    
    def propagateRays(self):
        """
        |  Propagates the Rays through the System considering a number of effects.
        """
        self.Components[0].propagateRays(self.Source)
        for i in range(len(self.Components)-1):
            if (self.Components[i].Mirror):
                self.Components[i+1].propagateRays(self.Components[i].rRays)
            else:
                self.Components[i+1].propagateRays(self.Components[i].tRays)
        if (self.Components[-1].Mirror):
            self.Detector.propagateRays(self.Components[-1].rRays)
        else:
            self.Detector.propagateRays(self.Components[-1].tRays)
        return
    
    def getSystemMuellerMatrix(self):
        """
        |  Computes the Mueller matrix for the System with following formulation. 

        |  Four minimum Stokes inputs :
        |  In_1 : [1,  1, 0, 0],  Ex_1 = 1.0 + 0.0j,  Ey_1 = 0.0 + 0.0j
        |  In_2 : [1, -1, 0, 0],  Ex_1 = 0.0 + 0.0j,  Ey_1 = 1.0 + 0.0j
        |  In_3 : [1,  0, 1, 0],  Ex_1 = 0.707 + 0.0j,  Ey_1 = 0.707 + 0.0j
        |  In_4 : [1,  0, 0, 1],  Ex_1 = 0.707 + 0.0j,  Ey_1 = 0.0 - 0.707j
        
        |  Four Stokes outputs as linear combination of Mueller matrix elements :
        |  Out_1 : [M_00+M_10, M_01+M_11, M_02+M_12, M_03+M_13]
        |  Out_2 : [M_00-M_10, M_01-M_11, M_02-M_12, M_03-M_13]
        |  Out_3 : [M_00+M_20, M_01+M_21, M_02+M_22, M_03+M_23]
        |  Out_4 : [M_00+M_30, M_01+M_31, M_02+M_32, M_03+M_33]
        
        |  System Mueller matrix elements : 
        |  [M_00, M_01, M_02, M_03] = 0.5*(Out_1 + Out_2)
        |  [M_10, M_11, M_12, M_13] = 0.5*(Out_1 - Out_2)
        |  [M_20, M_21, M_22, M_23] = Out_3 - 0.5*(Out_1 + Out_2)
        |  [M_30, M_31, M_32, M_33] = Out_4 - 0.5*(Out_1 + Out_2)

        |  Returns : 4x4 Mueller matric of the System, Transmitted fractional intensity
        """
        Source, Components = self.Source, self.Components
        Ex_In = [1.0+0.0j, 0.0+0.0j, np.sqrt(0.5)+0.0j, np.sqrt(0.5)+0.0j]
        Ey_In = [0.0+0.0j, 1.0+0.0j, np.sqrt(0.5)+0.0j, 0.0-np.sqrt(0.5)*1j]
        OutRays = []
        for Ex, Ey in zip(Ex_In, Ey_In):
            Source.createPolarization(Ex, Ey)
            Components[0].propagateRays(Source)
            for i in range(len(Components)-1):
                if (Components[i].Mirror):
                    Components[i+1].propagateRays(Components[i].rRays)
                else:
                    Components[i+1].propagateRays(Components[i].tRays)
            if(Components[-1].Mirror):
                OutRays.append(Components[-1].rRays)
            else:
                OutRays.append(Components[-1].tRays)        
        Out = []
        M = []
        for O in OutRays:
            ExTemp, EyTemp = np.sum(O.Ex*O.Mask), np.sum(O.Ey*O.Mask)
            Out.append(ExTemp*np.conjugate(ExTemp)+EyTemp*np.conjugate(EyTemp))
            Out.append(ExTemp*np.conjugate(ExTemp)-EyTemp*np.conjugate(EyTemp))
            Out.append(ExTemp*np.conjugate(EyTemp)+EyTemp*np.conjugate(ExTemp))
            Out.append(1j*(EyTemp*np.conjugate(ExTemp)-ExTemp*np.conjugate(EyTemp)))
        Out = np.real(np.array(Out)).reshape((4,4))
        M = np.copy(Out)
        M[:,0] = 0.5*(Out[0,:] + Out[1,:])
        M[:,1] = 0.5*(Out[0,:] - Out[1,:])
        M[:,2] = Out[2,:] - 0.5*(Out[0,:] + Out[1,:])
        M[:,3] = Out[3,:] - 0.5*(Out[0,:] + Out[1,:])
        self.MuellerMatrix = np.matrix(M/M[0,0])
        self.Transmission = M[0,0]/np.float(Source.NRays)**2
        return self.MuellerMatrix, self.Transmission
        
    def draw(self, Ax, clear=False):
        """
        |  Draw the optical System.
        |  Input : Pyplot axis, option to clear the previous plots.
        """
        Source = cp.copy(self.DisplayRays)
        Components = cp.copy(self.Components)
        Detector = cp.copy(self.Detector) 
        #
        Components[0].propagateRays(Source)
        for i in range(len(Components)-1):
            if (Components[i].Mirror):
                Components[i+1].propagateRays(Components[i].rRays)
            else:
                Components[i+1].propagateRays(Components[i].tRays)
        if (Components[-1].Mirror):
            Detector.propagateRays(Components[-1].rRays)
        else:
            Detector.propagateRays(Components[-1].tRays)
        #
        if (clear) : Ax.clear()
        for Component in Components:
            Component.draw(Ax, color='r', alpha=0.5)
            Component.drawRays(Ax, color='k', alpha=0.7)
        Detector.draw(Ax, color='r', alpha=0.5)
        Detector.drawRays(Ax, color='k', alpha=0.7)
        return
    
    def drawSpotDiagram(self, Ax, **kwargs):
        """
        |  Draw the spot diagram as seen on the detector.
        |  Input : Pyplot axis, kwargs are directly passed to plot function.
        """
        self.Detector.tRays.drawSpotDiagram(Ax, 0, **kwargs)
        return

    def getOutputPositionAngle(self):
        """
        |  Computes the position angle difference between the detector and output rays.
        |  Returns : Position angle in degrees
        """
        if (self.Components[-1].Mirror):
            OutputRays = self.Components[-1].rRays
        else:
            OutputRays = self.Components[-1].tRays
        DOT = np.sum(self.Detector.xAxis*OutputRays.xAxis)
        CROSSTemp = np.cross(self.Detector.xAxis, OutputRays.xAxis)
        CROSS = np.sum(OutputRays.oAxis*CROSSTemp)
        ThetaFinal = np.degrees(np.arctan2(CROSS, DOT))
        return ThetaFinal