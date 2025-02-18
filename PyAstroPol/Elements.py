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
    
    
class CompoundPlate():
    def __init__(self, Dia, FileName, Cover='Air', Substrate='Air'):
        self.Aperture = Dia
        # Data = np.genfromtxt(FileName, delimiter=',', dtype=str)
        # self.Layers = Data[:,0]
        # self.Wavelength = 0.6328
        # self.Stack = np.zeros([Data.shape[0],7])
        # self.Stack[:,0] = np.float64(Data[:,1])
        # self.Stack[:,4::] = np.float64(Data[:,2::])
        # self.Cover=Material(Cover)
        # self.Substrate=Material(Substrate)
        # self.loadRefractiveIndex()
        # self.Origin = np.array([0.0,0.0,0.0])
        # self.xAxis = np.array([1.0,0.0,0.0])
        # self.yAxis = np.array([0.0,1.0,0.0])
        # self.oAxis = np.array([0.0,0.0,1.0])
        # self.Mirror = False


    def propagateRays(self, Rays):
        """
        |  Compute reflection and transmission 'matrices' for the compound coating 
        |  Instead of a simple scalar for the coeff
        """
        self.iRays = cp.copy(Rays)
        self.tRays = cp.copy(Rays)
        self.rRays = cp.copy(Rays)
        self.Wavelength = self.iRays.Wavelength
        self.propagatePolarization()
        return
    
    def propagatePolarization(self):
        # # Compute coordinate rotation angles
        # self.sCosines = normalize3DVectors(np.cross(self.iRays.oCosines, self.oAxis))
        # self.sCosines = np.nan_to_num(self.sCosines)
        # DOT = np.sum(self.iRays.xCosines*self.sCosines, axis=1)
        # CROSSTemp = np.cross(self.iRays.xCosines, self.sCosines)
        # CROSS = np.sum(self.iRays.oCosines*CROSSTemp, axis=1)
        # Theta = np.reshape(np.arctan2(CROSS, DOT), newshape=(self.iRays.NRays, 1))
        #
        Theta = 0
        self.computeBerremanCharMatrix()
        A = self.CharMatrix
        self.rsp = np.matrix([[A[0,0], -A[2,0]],
                              [-A[0,2], A[2,2]]])
        self.rsp = self.rsp/np.linalg.det(self.rsp)
        self.tsp = np.matrix([[A[3,2], A[3,0]],
                               [A[1,2], A[1,0]]])
        self.tsp = self.tsp*self.rsp
        #
        # Incident
        Es =  self.iRays.Ex*np.cos(Theta)+self.iRays.Ey*np.sin(Theta)
        Ep = -self.iRays.Ex*np.sin(Theta)+self.iRays.Ey*np.cos(Theta)
        # Reflection
        Es_r, Ep_r = Es*self.rsp[0,0]+Ep*self.rsp[0,1], Es*self.rsp[1,0]+Ep*self.rsp[1,1] 
        self.rRays.Ex =  Es_r*np.cos(-Theta) + Ep_r*np.sin(-Theta)
        self.rRays.Ey = -Es_r*np.sin(-Theta) + Ep_r*np.cos(-Theta)
        self.rRays.xCosines =  self.sCosines*np.cos(-Theta) + self.pCosines_r*np.sin(-Theta)
        self.rRays.yCosines = -self.sCosines*np.sin(-Theta) + self.pCosines_r*np.cos(-Theta)
        self.rRays.xAxis = self.rRays.xCosines[0,:]
        self.rRays.yAxis = self.rRays.yCosines[0,:]
        self.rRays.oAxis = self.rRays.oCosines[0,:]
        # Transmission
        Es_t, Ep_t = Es*self.tsp[0,0]+Ep*self.tsp[0,1], Es*self.tsp[1,0]+Ep*self.tsp[1,1] 
        self.tRays.Ex =  Es_t*np.cos(-Theta) + Ep_t*np.sin(-Theta)
        self.tRays.Ey = -Es_t*np.sin(-Theta) + Ep_t*np.cos(-Theta)
        self.tRays.xCosines =  self.sCosines*np.cos(-Theta) + self.pCosines_t*np.sin(-Theta)
        self.tRays.yCosines = -self.sCosines*np.sin(-Theta) + self.pCosines_t*np.cos(-Theta)
        self.tRays.xAxis = self.tRays.xCosines[0,:]
        self.tRays.yAxis = self.tRays.yCosines[0,:]
        self.tRays.oAxis = self.tRays.oCosines[0,:]
        return

    def computeBerremanCharMatrix(self):
        """
        |  Compute the Berreman's characteristic matrix for a coated compund plate in the air
        |  Input:  Snell's propagation constant
        |  Output: 4x4 Berreman's characteristic matrix
        """
        # self.loadRefractiveIndex()
        # self.computePropagationConstant()
        # FM_Cover = getBerremanFieldMatrixIsotropic(self.Cover.getRefractiveIndexAt(self.Wavelength), self.Beta)
        # FM_Substrate = getBerremanFieldMatrixIsotropic(self.Substrate.getRefractiveIndexAt(self.Wavelength), self.Beta)
        # CharMatrix = np.linalg.inv(FM_Cover)
        # for s in self.Stack:
        #     Thick, nx, ny, nz, Eta, Psi, Xi = s
        #     Epsilon_ = getPermittivityTensor([nx,ny,nz], [Eta, Psi, Xi])
        #     FM_, Alpha_ = getBerremanFieldMatrixAnisotropic(Epsilon_, self.Beta)
        #     PM_ = getBerremanPhaseMatrix(Alpha_, Thick)
        #     CharMatrix = CharMatrix*FM_*PM_*np.linalg.inv(FM_)
        # self.CharMatrix = CharMatrix*FM_Substrate
        for s in self.Surfaces:
            s.propagateRays()



        return

    def loadRefractiveIndex(self):
        """
        |  Load the refractive indices for all the layers for a given wavelength
        """
        for i, l in enumerate(self.Layers):
            RIs = Material(l).getRefractiveIndicesAt(self.Wavelength)
            self.Stack[i,1:4] = RIs

    def computePropagationConstant(self):
        """
        |  Compute Snell's propagation constant for each ray i.e., n*sin(theta)
        |  Input: DC's of the surface normals at the points of incidence
        """
        RI = self.Cover.getRefractiveIndexAt(self.Wavelength)
        self.Beta = RI*np.linalg.norm(np.cross(self.iRays.oAxis, self.oAxis))
        return


