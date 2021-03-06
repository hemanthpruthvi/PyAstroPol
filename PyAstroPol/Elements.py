"""
|  Author : Hemanth Pruthvi
|  File name : Elements.py
|  Package : PyAstroPol
|  Description : Detector, Lens and likewise classes of optical elelements
"""

#
from .Source import *
from .Surface import *
from .Functions import *
#
   
class Detector(Surface):
    """
    |  Plane Surface() with no special attributes.
    |  Created for the purpose of easily specifying detector.
    """
    def __init__(self, Dia):
        Surface.__init__(self, Dia)
        self.Mirror = False
        return    

class UncoatedLens():
    """
    |  A simple singlet lens containing two uncoated surfaces

    |  Attributes: General
    |  S1 : Surface() : Front surface
    |  S2 : Surface() : Back surface
    |  Thick : Float : Thickness of the lens (redundant)
    """
    def __init__(self, Dia, Thick, R1=np.inf, R2=np.inf, K1=0.0, K2=0.0, n=1.5+0.0j):
        self.S1 = Surface(Dia, R=R1, K=K1, n1=1.0+0.0j, n2=n)
        self.S2 = Surface(Dia, R=R2, K=K2, n1=n, n2=1.0+0.0j)
        self.S2.translateOrigin(z=Thick)
        self.S1.rRes = self.S2.rRes
        self.S1.thetaRes = self.S2.thetaRes
        self.Mirror = False
        self.Thick = Thick
        return
    
    def propagateRays(self, Rays):
        """
        |  Propagate given Rays through the element.
        |  Specify a Rays() for the input. 
        """
        self.iRays = cp.copy(Rays)
        self.S1.propagateRays(Rays)
        self.S2.propagateRays(self.S1.tRays)
        self.rRays = self.S1.rRays
        self.tRays = self.S2.tRays
        return
    
    def rotateAboutX(self, ThetaX):
        """
        |  Rotate the element about global X-axis relative to present position.
        |  Input : X-rotation angle in degrees
        """
        self.S1.rotateAboutX(ThetaX)
        self.S2.rotateAboutX(ThetaX)
        return
    
    def rotateAboutY(self, ThetaY):
        """
        |  Rotate the element about global Y-axis relative to present position.
        |  Input : Y-rotation angle in degrees
        """
        self.S1.rotateAboutY(ThetaY)
        self.S2.rotateAboutY(ThetaY)
        return
    
    def rotateAboutZ(self, ThetaZ):
        """
        |  Rotate the element about global Z-axis relative to present position.
        |  Input : Z-rotation angle in degrees
        """
        self.S1.rotateAboutZ(ThetaZ)
        self.S2.rotateAboutZ(ThetaZ)
        return
    
    def translateOrigin(self, x=0.0, y=0.0, z=0.0):
        """
        |  Translate the element relative to present position.
        |  Inputs : X translation, Y translation and Z translation.
        """
        self.S1.translateOrigin(x=x, y=y, z=z)
        self.S2.translateOrigin(x=x, y=y, z=z)
        return
    
    def draw(self, Ax, **kwargs):
        """
        |  Draw the element in 3D.
        |  Inputs : Pyplot Axis, kwargs that are directly pass to plot function.
        """
        x1temp = np.reshape(self.S1.X, newshape=(self.S1.thetaRes, self.S1.rRes))[:,-1]
        y1temp = np.reshape(self.S1.Y, newshape=(self.S1.thetaRes, self.S1.rRes))[:,-1]
        z1temp = np.reshape(self.S1.Z, newshape=(self.S1.thetaRes, self.S1.rRes))[:,-1]
        x2temp = np.reshape(self.S2.X, newshape=(self.S2.thetaRes, self.S2.rRes))[:,-1]
        y2temp = np.reshape(self.S2.Y, newshape=(self.S2.thetaRes, self.S2.rRes))[:,-1]
        z2temp = np.reshape(self.S2.Z, newshape=(self.S2.thetaRes, self.S2.rRes))[:,-1]
        x, y, z = [], [], []
        for i in range(self.S1.thetaRes):
            x.append(x1temp[i])
            x.append(x2temp[i])
            y.append(y1temp[i])
            y.append(y2temp[i])
            z.append(z1temp[i])
            z.append(z2temp[i])
        x = np.array(x).reshape((self.S1.thetaRes,2))
        y = np.array(y).reshape((self.S1.thetaRes,2))
        z = np.array(z).reshape((self.S1.thetaRes,2))
        Ax.plot_surface(x, y, z, antialiased=True, **kwargs)
        #
        self.S1.draw(Ax, **kwargs)
        self.S2.draw(Ax, **kwargs)
        return
    
    def drawRays(self, Ax, **kwargs):
        """
        |  Draw incident rays to the Surfaces of the element.
        |  Inputs : Pyplot Axis, kwargs that are directly pass to plot function.
        """
        self.S1.drawRays(Ax, **kwargs)
        self.S2.drawRays(Ax, **kwargs)
        return
    
    
