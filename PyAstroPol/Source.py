"""
Author : Hemanth Pruthvi
File name : Source.py
Package : PyAstroPol
Description : Rays, Source and likewise classes
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
        self.Mask = np.ones((N,1), dtype=bool)
        self.createPolarization(1.0+0.0j, 0.0+0.0j)
        self.Wavelength = 0.6328
        return
    #
    def createPolarization(self, Ex, Ey):
        self.Ex = Ex*np.ones((self.NRays,1))
        self.Ey = Ey*np.ones((self.NRays,1))
        return 
    
    def computeKurtosis(self):
        k = np.sum((self.Origin + Length*self.oAxis - self.Points)*self.oAxis, axis=1)/np.sum(self.oCosines*self.oAxis, axis=1)
        k = np.abs(k.reshape((len(k),1)))
        x = np.sum((self.Points + k*self.oCosines)*self.xAxis, axis=1)
        y = np.sum((self.Points + k*self.oCosines)*self.yAxis, axis=1)
        self.Kurt = np.sum((x-x[0])**4 + (y-y[0])**4)/np.sum((x-x[0])**2 + (y-y[0])**2)**2
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
    # 
    def getStokesVector(self):
        Ex, Ey = np.sum(self.Ex*self.Mask), np.sum(self.Ey*self.Mask)
        I = np.real(Ex*np.conjugate(Ex) + Ey*np.conjugate(Ey))
        Q = np.real(Ex*np.conjugate(Ex) - Ey*np.conjugate(Ey))
        U = np.real(Ex*np.conjugate(Ey) + Ey*np.conjugate(Ex))
        V = np.real(1j*(Ey*np.conjugate(Ex) - Ex*np.conjugate(Ey)))
        return np.matrix([[np.real(1.0)],[np.real(Q/I)],[np.real(U/I)],[np.real(V/I)]]), I/np.sum(self.Mask)**2
    
class Source(Rays):
    """
    Source definition with Rays() as parent
    
    // Attributes : General
    Type : Enumerated string : "collimated" or "point" sources
    FNum : Float : Focal ratio of the beam, only used for point sources
    Clear : Float : Clear aperture for the distribution of rays, only used in collimated sources
    Random : Bool : Ray distrbution, if not True rays are distributed in ring fashion with quasi-uniform density
    """
    def __init__(self, NRays, Type='Collimated', Clear=1.0, FNum=1.0, Random=False):
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
        if (self.Type == 'Point'):
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
        elif (self.Type == 'Collimated'):
            # Repeat the oCosines 
            Radii = np.array(np.transpose([Clear*np.array(Radii)]))
            Thetas = np.array(np.transpose([Thetas]))
            #
            self.Points[1::,:] = self.xAxis*Radii*np.cos(Thetas) + self.yAxis*Radii*np.sin(Thetas)
            self.oCosines[1::,:] = np.reshape(np.tile(self.oAxis, NRays-1), newshape=(NRays-1, 3))
        
        #Other keywords
        else:
            print('Error! Invalid source type! \n Use either "Point" or "Collimated" (default) !')
            return
        
        # Now turn for s- and p- polarization axes
        xAxes = np.cross(self.yAxis, self.oCosines)
        self.xCosines = normalize3DVectors(xAxes)
        self.yCosines = np.cross(self.oCosines, self.xCosines)
        return

    # Apply rotation and translation for all points and directions
    def applyTransformation(self, M):
        #
        self.Origin = applyPointTransformation(self.Origin, M)
        self.Points = applyPointTransformation(self.Points, M)
        #
        self.oAxis = applyVectorTransformation(self.oAxis, M)
        self.xAxis = applyVectorTransformation(self.xAxis, M)
        self.yAxis = applyVectorTransformation(self.yAxis, M)
        self.oCosines = applyVectorTransformation(self.oCosines, M)
        self.xCosines = applyVectorTransformation(self.xCosines, M)
        self.yCosines = applyVectorTransformation(self.yCosines, M)
        #
        self.oAxis = normalize3DVectors(self.oAxis)
        self.xAxis = normalize3DVectors(self.xAxis)
        self.yAxis = normalize3DVectors(self.yAxis)
        self.oCosines = normalize3DVectors(self.oCosines)
        self.xCosines = normalize3DVectors(self.xCosines)
        self.yCosines = normalize3DVectors(self.yCosines)
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
    #
    def pointToDirection(self, NewNormal):
        OldNormal = np.copy(self.oAxis)
        ThetaY = np.arcsin(OldNormal[0])
        ThetaX = np.arctan2(-OldNormal[1]/np.cos(ThetaY), OldNormal[2]/np.cos(ThetaY))
        self.rotateAboutX(np.degrees(-ThetaX))
        self.rotateAboutY(np.degrees(-ThetaY))
        #
        ThetaY = np.arcsin(NewNormal[0])
        ThetaX = np.arctan2(-NewNormal[1]/np.cos(ThetaY), NewNormal[2]/np.cos(ThetaY))
        self.rotateAboutY(np.degrees(ThetaY))
        self.rotateAboutX(np.degrees(ThetaX))
        return
    #
    def makeOrigin(self, NewOrigin):
        O, NO = self.Origin, NewOrigin
        self.translateOrigin(x=-O[0], y=-O[1], z=-O[2])
        self.translateOrigin(x=NO[0], y=NO[1], z=NO[2])
        return

class AstroSource(Source):
    #
    def __init__ (self, NRays, HA=0.0, Dec=0.0, Lat=30.0, Clear=1.0, Dist=1000.0):
        Source.__init__(self, NRays, Type='Collimated', Clear=Clear, Random=False)
        self.Type = 'Collimated'
        self.Aperture = Clear
        self.HourAngle = HA
        self.Declination = Dec
        self.Latitude = Lat
        self.Distance = Dist
        #
        HA, Dec, Lat = np.radians(HA), np.radians(Dec), np.radians(Lat)
        Position = np.array([np.sin(HA)*np.cos(Dec), 
                             np.cos(HA)*np.cos(Dec)*np.cos(Lat) + np.sin(Dec)*np.sin(Lat), 
                             -np.cos(HA)*np.cos(Dec)*np.sin(Lat) + np.sin(Dec)*np.cos(Lat)])
        ThetaY = np.arcsin(-Position[0])
        ThetaX = np.arctan2(Position[1]/np.cos(ThetaY), -Position[2]/np.cos(ThetaY))
        #
        self.ThetaY, self.ThetaX = np.degrees(ThetaY), np.degrees(ThetaX)
        self.rotateAboutY(self.ThetaY)
        self.rotateAboutX(self.ThetaX)
        #
        self.translateOrigin(x=-Dist*self.oAxis[0], y=-Dist*self.oAxis[1], z=-Dist*self.oAxis[2])
        return
        
        





























