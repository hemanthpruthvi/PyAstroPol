import numpy as np
import copy as cp
import random as rd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime as dt
#
from .Polarization import *
#

class Rays():
    """
    Rays definition with Point-Eiknol formulation.
    
    // Attributes : General
    NRays : Int : Number of rays
    Wavelength : Float : Wavelength in microns
    
    // Attributes : References
    Origin : Float(1,3) : Starting point of cheif-ray
    xAxis : Float(1,3) : DC of supposed x-axis
    yAxis : Float(1,3) : DC of supposed y-axis
    oAxis : Float(1,3) : DC of optical axis i.e., direction of chief ray
    
    // Attributes : Set of points and cosines
    Points : Float(N,3) : Starting points of rays
    xCosines : Float(N,3) : DC of x-polarization direction for all rays
    yCosines : Float(N,3) : DC of y-polarization direction for all rays
    oCosines : Float(N,3) : DC of propagation direction for all rays
    
    // Attributes : Electric field
    Ex : ComplexFloat(N,1) : Electric field for each ray, in x-direction
    Ey : ComplexFloat(N,1) : Electric field for each ray, in y-direction
    
    (x,y,o) are orthonormal basis for a right handed coordinate system
    """
    
    def __init__(self, N):
        #
        self.NRays = np.copy(N)
        #
        self.Origin = np.array([0.0,0.0,0.0])
        self.xAxis = np.array([1.0,0.0,0.0])
        self.yAxis = np.array([0.0,1.0,0.0])
        self.oAxis = np.array([0.0,0.0,1.0])
        #
        self.Points = np.reshape(np.tile(self.Origin, N), newshape=(N,3))
        self.xCosines = np.reshape(np.tile(self.xAxis, N), newshape=(N,3))
        self.yCosines = np.reshape(np.tile(self.yAxis, N), newshape=(N,3))
        self.oCosines = np.reshape(np.tile(self.oAxis, N), newshape=(N,3))
        #
        self.createPolarization(1.0+0.0j, 0.0+0.0j)
        self.Wavelength = 0.6328
        return
    #
    def createPolarization(self, Ex, Ey):
        self.Ex = Ex*np.ones((self.NRays,1))
        self.Ey = Ey*np.ones((self.NRays,1))
        return 
    
    # Graphics
    def drawRays(self, Ax, Length, **kwargs):
        P1 = self.Points
        P2 = P1 + Length*self.oCosines
        for i in range(len(P1)):
            Ax.plot([P1[i,0], P2[i,0]], [P1[i,1], P2[i,1]], [P1[i,2], P2[i,2]], **kwargs)
        return
    #
    def drawXYAxes(self, Ax, LengthR, LengthP, **kwargs): 
        P1 = self.Points + LengthR*self.oCosines - self.xCosines*0
        P2 = self.Points + LengthR*self.oCosines + self.xCosines*LengthP
        for i in range(len(P1)):
            Ax.plot([P1[i,0], P2[i,0]], [P1[i,1], P2[i,1]], [P1[i,2], P2[i,2]], color='y', **kwargs)
        P2 = self.Points + LengthR*self.oCosines + self.yCosines*LengthP
        for i in range(len(P1)):
            Ax.plot([P1[i,0], P2[i,0]], [P1[i,1], P2[i,1]], [P1[i,2], P2[i,2]], color='b', **kwargs)
        return
    # 
    def drawSpotDiagram(self, Ax, Length, **kwargs): 
        k = np.sum((self.Origin + Length*self.oAxis - self.Points)*self.oAxis, axis=1)/np.sum(self.oCosines*self.oAxis, axis=1)
        k = np.abs(k.reshape((len(k),1)))
        x = np.sum((self.Points + k*self.oCosines)*self.xAxis, axis=1)
        y = np.sum((self.Points + k*self.oCosines)*self.yAxis, axis=1)
        Ax.set_aspect('equal')
        Ax.scatter(x, y, **kwargs)
        return

    
class Source(Rays):
    """
    Source definition with Rays() as parent
    
    // Attributes : General
    Type : Enumerated string : "collimated" or "point" sources
    FNum : Float : Focal ratio of the beam, only used for point sources
    Clear : Float : Clear aperture for the distribution of rays, only used in collimated sources
    Random : Bool : Ray distrbution, if not True rays are distributed in ring fashion with quasi-uniform density
    """
    def __init__(self, NRays, Type='collimated', Clear=1.0, FNum=1.0, Random=False):
        Rays.__init__(self, NRays)
        self.Type = np.copy(Type)
        self.Random = np.copy(Random)
        self.Clear = np.copy(Clear)
        self.FNum = np.copy(FNum)

        # Create random distribution of rays
        Rad = 0.5
        if (Random):
            Radii_temp = np.linspace(0, Rad, NRays-1)
            Thetas_temp = np.linspace(0, 2*np.pi, NRays-1)
            Radii = rd.choices(Radii_temp, weights=Radii_temp, k=NRays-1)
            Thetas = rd.choices(Thetas_temp, k=NRays-1)
            
        # Create concentric distribution of rays, like "rings"
        else:
            RayDens = np.sqrt(NRays/np.pi/Rad**2)
            NRings = np.int(Rad*RayDens)
            Rings = np.linspace(0, Rad, NRings+1)
            NThetas = np.zeros(NRings, dtype=np.int)
            Radii, Thetas = np.array([]), np.array([])
            for i in range(NRings-1):
                NThetas[i] = int(2*np.pi*Rings[i+1]*RayDens)
                Radii = np.append(Radii, Rings[i+1]*np.ones(NThetas[i]))
                Thetas = np.append(Thetas, np.linspace(0, 2*np.pi, NThetas[i]+1)[0:-1])
            i += 1
            NThetas[i] = NRays-np.sum(NThetas)-1
            Radii = np.append(Radii, Rings[i+1]*np.ones(NThetas[i]))
            Thetas = np.append(Thetas, np.linspace(0, 2*np.pi, NThetas[i]+1)[0:-1])
        self.Radii = Radii
        self.Thetas = Thetas

        # Point sources
        if (self.Type == 'point'):
            # make vectorial addition of oAxis unit vector and vector for rings
            Radii = np.array(np.transpose([np.array(Radii)*1.0/FNum]))
            Thetas = np.array(np.transpose([Thetas]))
            oCosines = self.xAxis*Radii*np.cos(Thetas) + self.yAxis*Radii*np.sin(Thetas)
            oCosines += self.oAxis
            # 
            self.oCosines[1::,:] = oCosines/np.array(np.transpose([np.linalg.norm(oCosines, axis=1)]))
            self.Points[1::,:] = np.reshape(np.tile(self.Origin, NRays-1), newshape=(NRays-1, 3))
            self.HalfConeAngle = np.degrees(np.arccos(np.sum(self.oCosines*self.oAxis, axis=1).min()))

        # Collimated sources
        else :
            # Repeat the oCosines 
            Radii = np.array(np.transpose([Clear*np.array(Radii)]))
            Thetas = np.array(np.transpose([Thetas]))
            #
            self.Points[1::,:] = self.xAxis*Radii*np.cos(Thetas) + self.yAxis*Radii*np.sin(Thetas)
            self.oCosines[1::,:] = np.reshape(np.tile(self.oAxis, NRays-1), newshape=(NRays-1, 3))
            
        # Now turn for s- and p- polarization axes
        xAxes = np.cross(self.yAxis, self.oCosines)
        self.xCosines = normalize3DVectors(xAxes)
        self.yCosines = np.cross(self.oCosines, self.xCosines)
        return

    # Apply rotation and translation for all points and directions
    def applyTransformation(self, M):
        #
        temp = np.copy(self.Points)
        x, y, z = np.copy(temp[:,0]), np.copy(temp[:,1]), np.copy(temp[:,2])
        temp[:,0] = M[0,0]*x + M[0,1]*y + M[0,2]*z + M[0,3]
        temp[:,1] = M[1,0]*x + M[1,1]*y + M[1,2]*z + M[1,3]
        temp[:,2] = M[2,0]*x + M[2,1]*y + M[2,2]*z + M[2,3]
        self.Points = np.copy(temp)
        #
        temp = np.copy(self.oCosines)
        x, y, z = np.copy(temp[:,0]), np.copy(temp[:,1]), np.copy(temp[:,2])
        temp[:,0] = M[0,0]*x + M[0,1]*y + M[0,2]*z
        temp[:,1] = M[1,0]*x + M[1,1]*y + M[1,2]*z
        temp[:,2] = M[2,0]*x + M[2,1]*y + M[2,2]*z
        self.oCosines = normalize3DVectors(temp)
        #
        temp = np.copy(self.xCosines)
        x, y, z = np.copy(temp[:,0]), np.copy(temp[:,1]), np.copy(temp[:,2])
        temp[:,0] = M[0,0]*x + M[0,1]*y + M[0,2]*z
        temp[:,1] = M[1,0]*x + M[1,1]*y + M[1,2]*z
        temp[:,2] = M[2,0]*x + M[2,1]*y + M[2,2]*z
        self.xCosines = normalize3DVectors(temp)
        #
        temp = np.copy(self.yCosines)
        x, y, z = np.copy(temp[:,0]), np.copy(temp[:,1]), np.copy(temp[:,2])
        temp[:,0] = M[0,0]*x + M[0,1]*y + M[0,2]*z
        temp[:,1] = M[1,0]*x + M[1,1]*y + M[1,2]*z
        temp[:,2] = M[2,0]*x + M[2,1]*y + M[2,2]*z
        self.yCosines = normalize3DVectors(temp)
        #
        temp = np.copy(self.Origin)
        x, y, z = np.copy(temp[0]), np.copy(temp[1]), np.copy(temp[2])
        temp[0] = M[0,0]*x + M[0,1]*y + M[0,2]*z + M[0,3]
        temp[1] = M[1,0]*x + M[1,1]*y + M[1,2]*z + M[1,3]
        temp[2] = M[2,0]*x + M[2,1]*y + M[2,2]*z + M[2,3]
        self.Origin = np.copy(temp)
        #
        temp = np.copy(self.oAxis)
        x, y, z = np.copy(temp[0]), np.copy(temp[1]), np.copy(temp[2])
        temp[0] = M[0,0]*x + M[0,1]*y + M[0,2]*z
        temp[1] = M[1,0]*x + M[1,1]*y + M[1,2]*z
        temp[2] = M[2,0]*x + M[2,1]*y + M[2,2]*z
        self.oAxis = normalize3DVectors(temp)
        #
        temp = np.copy(self.xAxis)
        x, y, z = np.copy(temp[0]), np.copy(temp[1]), np.copy(temp[2])
        temp[0] = M[0,0]*x + M[0,1]*y + M[0,2]*z
        temp[1] = M[1,0]*x + M[1,1]*y + M[1,2]*z
        temp[2] = M[2,0]*x + M[2,1]*y + M[2,2]*z
        self.xAxis = normalize3DVectors(temp)
        #
        temp = np.copy(self.yAxis)
        x, y, z = np.copy(temp[0]), np.copy(temp[1]), np.copy(temp[2])
        temp[0] = M[0,0]*x + M[0,1]*y + M[0,2]*z
        temp[1] = M[1,0]*x + M[1,1]*y + M[1,2]*z
        temp[2] = M[2,0]*x + M[2,1]*y + M[2,2]*z
        self.yAxis = normalize3DVectors(temp)
        return

    # Rotation matrices for the function
    def rotateAboutX(self, ThetaX):
        ThetaX = np.radians(ThetaX)
        R = np.matrix([[1.0, 0.0, 0.0, 0.0],
                       [0.0, np.cos(ThetaX), -np.sin(ThetaX), 0.0], 
                       [0.0, np.sin(ThetaX),  np.cos(ThetaX), 0.0],
                       [0.0, 0.0, 0.0, 1.0]])
        self.applyTransformation(R)
        return
    def rotateAboutY(self, ThetaY):
        ThetaY = np.radians(ThetaY)
        R = np.matrix([[ np.cos(ThetaY), 0.0, np.sin(ThetaY), 0.0],
                       [0.0, 1.0, 0.0, 0.0],
                       [-np.sin(ThetaY), 0.0, np.cos(ThetaY), 0.0],
                       [0.0, 0.0, 0.0, 1.0]])
        self.applyTransformation(R)
        return
    def rotateAboutZ(self, ThetaZ):
        ThetaZ = np.radians(ThetaZ)
        R = np.matrix([[np.cos(ThetaZ), -np.sin(ThetaZ), 0.0, 0.0],
                       [np.sin(ThetaZ),  np.cos(ThetaZ), 0.0, 0.0],
                       [0.0, 0.0, 1.0, 0.0],
                       [0.0, 0.0, 0.0, 1.0]])
        self.applyTransformation(R)
        return

    # Translation matrix for the function
    def translateOrigin(self, x=0.0, y=0.0, z=0.0):
        T = np.matrix([[1.0, 0.0, 0.0, x],
                       [0.0, 1.0, 0.0, y],
                       [0.0, 0.0, 1.0, z],
                       [0.0, 0.0, 0.0, 1.0]])
        self.applyTransformation(T)
        return