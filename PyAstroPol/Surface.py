import numpy as np
import copy as cp
import random as rd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime as dt
#
from .Functions import *#
#

class Surface():
    """
    General purpose Surface definition
    
    // Attributes : General
    Aperture : Float : Diameter size (only circular is considered) 
    Origin : Float(1,3) : Center of the aperture (circle)
    oAxis : Float(1,3) : DC of Optical axis i.e., normal to surface, at Origin, pointing towards the direction of the incidence
    Curvature : Float : Radius of curvature of the surface
    Conic : Float : Conic constant of the surface
    SurfaceMatrix : Float(4,4) : Matrix representation of coefficients of the surface equation 
    iRI : FloatComplex : Refractive index of the incident medium
    tRI : FloatComplex : Refrative index of the transmitting medium
    iMU : FloatComplex : Magnetic permeability of the incident medium
    tMU : FloatComplex : Magnetic permeability of the transmitting medium
    Mirror : Bool : Is it only reflecting or only transmitting?
    IncidentPoints : Float(N,3) : Points on the surface corresponding to incidence of Rays
    nCosines : Float(N,3) : DC of the normals to the surface at the points of ray incidences 
    tCosines : Float(N,3) : DC of the tangents to the surface at the points of ray incidences
    Mask : Bool(N,1) : Array representing whether given rays are hitting the surface or not
    
    // Attributes : Polarization
    sCosines : Float(N,3) :  DC of the s-polarization direction to the surface for the Rays
    pCosines_i : Float(N,3) : DC of the p-polarization direction to the surface, for the incident rays 
    pCosines_r : Float(N,3) : DC of the p-polarization direction to the surface, for the reflected rays
    pCosines_t : Float(N,3) : DC of the p-polarization direction to the surface, for the transmitted rays
    rs : FloatComplex(N,1) : Complex coefficient of reflection for s-polarization of all Rays 
    rp : FloatComplex(N,1) : Complex coefficient of reflection for p-polarization of all Rays 
    ts : FloatComplex(N,1) : Complex coefficient of transmission for s-polarization of all Rays 
    tp : FloatComplex(N,1) : Complex coefficient of transmission for p-polarization of all Rays
    
    // Attributes : Objects
    iRays : Rays() : Incident rays
    rRays : Rays() : Reflected rays
    tRays : rays() : Transmitted rays
    
    // Attributes : Only for drawing
    X, Y, Z : Float(N,1) : Coordinates of the points on the surface used for rendering the surface
    """
    def __init__(self, Dia, R=np.inf, K=0, Mirror=False, iDia=0.0,
                 OffAxis=False, OffAxDist=0.0, OffAxAz=0.0, 
                 n1=1.0+0.0j, n2=1.5+0.0j):
        self.iRI = np.copy(n1) # Electric permeability/ refractive index of medium-1
        self.tRI = np.copy(n2) # Electric permeability/ refractive index of medium-2
        self.iMU = 1.0 + 0.0j # Magnetic permeability of medium-1
        self.tMU = 1.0 + 0.0j # Magnetic permeability of medium-2
        self.Aperture = Dia
        self.Mirror = Mirror
        self.InnerDia = iDia
        self.Origin = np.array([0.0,0.0,0.0])
        self.xAxis = np.array([1.0,0.0,0.0])
        self.yAxis = np.array([0.0,1.0,0.0])
        self.oAxis = np.array([0.0,0.0,1.0])
        self.Conic, self.Curvature = K, R
        self.SurfaceMatrix = np.matrix([[1.0/R, 0.0, 0.0, 0.0],
                                        [0.0, 1.0/R, 0.0, 0.0],
                                        [0.0, 0.0, (K+1.0)/R, -2.0],
                                        [0.0, 0.0, 0.0, 0.0]])
        
        self.X, self.Y, self.Z = 0.0, 0.0, 0.0
        self.IncidentPoints = np.array([[0.0,0.0,0.0]])
        # Special things to do in case of off-axis surfaces
        if (OffAxis):
            if (self.Conic == 0 or self.Curvature == np.inf):
                print('Error! Cannot make this surface off-axis!')
                self.OffAxis = False
                self.OffAxisDistance = 0.0
                self.OffAxisAzimuth = 0.0
                return
            else:
                self.OffAxis = True
                self.OffAxisDistance = np.copy(OffAxDist)
                self.OffAxisAzimuth = np.copy(OffAxAz)
                X, Y = OffAxDist*np.cos(np.radians(OffAxAz)), OffAxDist*np.sin(np.radians(OffAxAz))
                if (self.Conic == -1):
                    Z = OffAxDist**2/self.Curvature/2.0
                else:
                    Z = (R - (R/np.abs(R))*np.sqrt(R**2-(K+1.0)*OffAxDist**2))/(K+1.0)
                OldOrigin = np.copy(self.Origin)
                NewOrigin = -X*np.array(self.xAxis).flatten() \
                            -Y*np.array(self.yAxis).flatten() \
                            -Z*np.array(self.oAxis).flatten()
                self.translateOrigin(x=NewOrigin[0], y=NewOrigin[1], z=NewOrigin[2])
                self.Origin = np.copy(OldOrigin)
        else:
            self.OffAxis = False
            self.OffAxisDistance = 0.0
            self.OffAxisAzimuth = 0.0
            
        #
        self.renderSurface(5, 13)
        
    # Create points to render surface
    def renderSurface(self, rRes, thetaRes):
        self.rRes, self.thetaRes = np.copy(rRes), np.copy(thetaRes)
        r = np.linspace(self.InnerDia/2.0, self.Aperture/2.0, rRes)
        theta = np.linspace(0, 2*np.pi, thetaRes)
        rMat, thetaMat = np.meshgrid(r,theta)
        x, y = rMat*np.cos(thetaMat), rMat*np.sin(thetaMat)
        self.X, self.Y = np.copy(x.flatten()), np.copy(y.flatten())
        x = x.flatten() + self.OffAxisDistance*np.cos(np.radians(self.OffAxisAzimuth))
        y = y.flatten() + self.OffAxisDistance*np.sin(np.radians(self.OffAxisAzimuth))
        
        a = (self.Conic+1.0)/self.Curvature
        b = -2
        c = (x**2+y**2)/self.Curvature
        if (self.Curvature == np.inf): # Plane
            z = x*0.0
        if (self.Conic == -1 and self.Curvature != np.inf): # Paraboloid
            z = -c/b
        if (self.Conic != -1 and self.Curvature != np.inf): # Other
            z = (-b-np.sqrt(b**2-4*a*c))/2/a
        self.Z = z - z[0]
        return
    # Apply rotation and translation to equation of surface and points on it
    def applyTransformation(self, M):
        #
        self.SurfaceMatrix = np.transpose(np.linalg.inv(M))*self.SurfaceMatrix*np.linalg.inv(M)
        # self.SurfaceMatrix[3,0:3], self.SurfaceMatrix[3,3] = 0.0, 1.0

        #
        x, y, z = np.copy(self.X), np.copy(self.Y), np.copy(self.Z)
        self.X = M[0,0]*x + M[0,1]*y + M[0,2]*z + M[0,3]
        self.Y = M[1,0]*x + M[1,1]*y + M[1,2]*z + M[1,3]
        self.Z = M[2,0]*x + M[2,1]*y + M[2,2]*z + M[2,3]
        #
        self.Origin = applyPointTransformation(self.Origin, M)
        self.IncidentPoints = applyPointTransformation(self.IncidentPoints, M)
        #
        self.oAxis = applyVectorTransformation(self.oAxis, M)
        self.xAxis = applyVectorTransformation(self.xAxis, M)
        self.yAxis = applyVectorTransformation(self.yAxis, M)
        #
        self.oAxis = normalize3DVectors(self.oAxis)
        self.xAxis = normalize3DVectors(self.xAxis)
        self.yAxis = normalize3DVectors(self.yAxis)
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
        
    # Compute points of incidence on the surface
    def computeIncidence(self, Rays):
        C = np.matrix(self.SurfaceMatrix)
        self.iRays = cp.copy(Rays)
        self.rRays = cp.copy(Rays)
        self.tRays = cp.copy(Rays)
        # Propapage incident rays
        I = np.transpose(np.matrix(np.hstack([Rays.oCosines, np.zeros([Rays.NRays,1])])))
        X0 = np.transpose(np.matrix(np.hstack([Rays.Points, np.ones([Rays.NRays,1])])))
        a = np.sum(np.array(I)*np.array(C*I), axis=0)
        b = np.sum(np.array(I)*np.array(C*X0), axis=0) + np.sum(np.array(X0)*np.array(C*I), axis=0) 
        c = np.sum(np.array(X0)*np.array(C*X0), axis=0)
        self.A, self.B, self.C = a, b, c
        self.I, self.X0 = I, X0
        dnan_replace = -c/b
        dplus = (-b+np.sqrt(b**2-4.0*a*c))/2.0/a
        dminus = (-b-np.sqrt(b**2-4.0*a*c))/2.0/a
        dmask = np.isclose(a, EPS)
        dplus = np.nan_to_num(dplus) + np.nan_to_num(dmask*dnan_replace)
        dminus = np.nan_to_num(dminus) + np.nan_to_num(dmask*dnan_replace)
        if (dminus[0] < 0 and dplus[0] < 0):
            print('Error! Verify direction of propagation!')
            self.IncidentPoints += np.nan
            return
            # print(dplus, dminus)
        if (dminus[0] < dplus[0]):
            if (dminus[0] < 0): 
                d = np.copy(dplus)
            else:
                d = np.copy(dminus)
        else:
            if (dplus[0] < 0): 
                d = np.copy(dminus)
            else:
                d = np.copy(dplus)
        self.DPLUS, self.DMINUS = dplus, dminus
        self.DNAN, self.DPMASK = dnan_replace, dmask
        self.D = d
        self.IncidentPoints = np.transpose((np.array(X0) + np.array(I)*d)[0:3,:])
        OriginDistance = np.sqrt(np.sum(np.array(self.IncidentPoints[0,:].flatten()-self.Origin.flatten())**2))
        if (OriginDistance > self.Aperture/2.0) : 
            print('Error! Rays are not hitting the surface!', 'Chief rays incidence : ', OriginDistance)
            #self.IncidentPoints += np.nan
            return
        self.rRays.Points = self.IncidentPoints
        self.rRays.Origin = self.rRays.Points[0,:]
        self.tRays.Points = self.IncidentPoints
        self.tRays.Origin = self.tRays.Points[0,:]
        # Create a mark to check if the rays are hitting the surface aperture
        d = np.linalg.norm(np.cross(self.rRays.Points-self.Origin, self.oAxis), axis=1)
        self.rRays.Mask = ((d<=self.Aperture/2.0) & (d>=self.InnerDia/2.0)).reshape((self.rRays.NRays, 1))
        self.rRays.Mask = self.iRays.Mask & self.rRays.Mask
        self.tRays.Mask = np.copy(self.rRays.Mask)
        self.tRays.oCosines = self.tRays.oCosines*self.tRays.Mask + self.iRays.oCosines*(~self.tRays.Mask)
        return
    # Compute normal to the surface at points of incedence
    def computeNormals(self):
        C = np.matrix(self.SurfaceMatrix)
        # Solution to normal Gradient(Surface) = 0 
        GX, GY, GZ = np.matrix([[1,0,0,0]]), np.matrix([[0,1,0,0]]), np.matrix([[0,0,1,0]])
        SX = GX*C + GX*np.transpose(C)
        SY = GY*C + GY*np.transpose(C)
        SZ = GZ*C + GZ*np.transpose(C)
        nCosines = np.zeros([self.iRays.NRays, 3])
        self.SX, self.SY, self.SZ = SX, SY, SZ
        nCosines[:,0] = SX[0,0]*self.IncidentPoints[:,0] + \
                        SX[0,1]*self.IncidentPoints[:,1] + \
                        SX[0,2]*self.IncidentPoints[:,2] + SX[0,3]
        nCosines[:,1] = SY[0,0]*self.IncidentPoints[:,0] + \
                        SY[0,1]*self.IncidentPoints[:,1] + \
                        SY[0,2]*self.IncidentPoints[:,2] + SY[0,3]
        nCosines[:,2] = SZ[0,0]*self.IncidentPoints[:,0] + \
                        SZ[0,1]*self.IncidentPoints[:,1] + \
                        SZ[0,2]*self.IncidentPoints[:,2] + SZ[0,3]
        self.nCosines = normalize3DVectors(nCosines)
        if (dot3DVectors(self.nCosines[0,:], -self.iRays.oCosines[0,:]) < 0):
            self.nCosines *= -1.0
        # Tangents
        self.iTheta = np.arccos(dot3DVectors(self.nCosines, -self.iRays.oCosines))
        self.iTheta = np.nan_to_num(self.iTheta)
        tCosines = (self.iRays.oCosines + self.nCosines*np.cos(self.iTheta))/np.sin(self.iTheta)
        self.tCosines = np.nan_to_num(tCosines)
        # s-Polarization
        sCosines = normalize3DVectors(np.cross(self.iRays.oCosines, self.nCosines))
        TempMask = np.reshape(np.isnan(np.sum(sCosines, axis=1)), newshape=(self.iRays.NRays,1))
        self.sCosines = np.nan_to_num(sCosines) + self.iRays.xCosines*TempMask
        self.pCosines_i = np.cross(self.iRays.oCosines, self.sCosines)
        return
    # Compute properties of reflected rays
    def propagateReflectedRays(self):
        # Reflection direction
        rCosines = self.nCosines*np.cos(self.iTheta) + self.tCosines*np.sin(self.iTheta)
        self.rRays.oCosines = normalize3DVectors(rCosines)
        # p-Polarization direction
        self.pCosines_r = np.cross(self.rRays.oCosines, self.sCosines)
        return
    # Compute properties of the transmitted rays
    def propagateTransmittedRays(self):
        # Transmission direction
        nCROSSi = np.reshape(np.linalg.norm(np.cross(self.nCosines, self.iRays.oCosines), axis=1), 
                             newshape=(self.iRays.NRays, 1))
        tTheta = np.real(np.arcsin(self.iRI*np.sin(self.iTheta)/self.tRI))
        tCosines = -self.nCosines*np.cos(tTheta) + self.tCosines*np.sin(tTheta)
        self.tRays.oCosines = normalize3DVectors(tCosines)
        self.tTheta = np.copy(tTheta)
        # p-Polarization direction
        self.pCosines_t = np.cross(self.tRays.oCosines, self.sCosines)
        return
    # Compute the state of polarization
    def propagatePolarization(self):
        # Compute coordinate rotation angles
        DOT = np.sum(self.iRays.xCosines*self.sCosines, axis=1)
        CROSSTemp = np.cross(self.iRays.xCosines, self.sCosines)
        CROSS = np.sum(self.iRays.oCosines*CROSSTemp, axis=1)
        Theta = np.reshape(np.arctan2(CROSS, DOT), newshape=(self.iRays.NRays, 1))
        # Coefficients of reflection and transmission
        self.rp = ((self.iRI/self.iMU)*np.cos(self.tTheta)-(self.tRI/self.tMU)*np.cos(self.iTheta)) / \
                    ((self.iRI/self.iMU)*np.cos(self.tTheta)+(self.tRI/self.tMU)*np.cos(self.iTheta))
        self.tp = (2*(self.iRI/self.iMU)*np.cos(self.iTheta)) / \
                    ((self.iRI/self.iMU)*np.cos(self.tTheta)+(self.tRI/self.tMU)*np.cos(self.iTheta))
        self.rs = ((self.iRI/self.iMU)*np.cos(self.iTheta)-(self.tRI/self.tMU)*np.cos(self.tTheta)) / \
                    ((self.iRI/self.iMU)*np.cos(self.iTheta)+(self.tRI/self.tMU)*np.cos(self.tTheta))
        self.ts = (2*(self.iRI/self.iMU)*np.cos(self.iTheta)) / \
                    ((self.iRI/self.iMU)*np.cos(self.iTheta)+(self.tRI/self.tMU)*np.cos(self.tTheta))
        # Incident
        Es =  self.iRays.Ex*np.cos(Theta)+self.iRays.Ey*np.sin(Theta)
        Ep = -self.iRays.Ex*np.sin(Theta)+self.iRays.Ey*np.cos(Theta)
        # Reflection
        Es_r, Ep_r = Es*self.rs, Ep*self.rp
        self.rRays.Ex =  Es_r*np.cos(np.pi+Theta) + Ep_r*np.sin(np.pi+Theta)
        self.rRays.Ey = -Es_r*np.sin(np.pi+Theta) + Ep_r*np.cos(np.pi+Theta)
        self.rRays.xCosines =  self.sCosines*np.cos(np.pi+Theta) + self.pCosines_r*np.sin(np.pi+Theta)
        self.rRays.yCosines = -self.sCosines*np.sin(np.pi+Theta) + self.pCosines_r*np.cos(np.pi+Theta)
        self.rRays.xAxis = self.rRays.xCosines[0,:]
        self.rRays.yAxis = self.rRays.yCosines[0,:]
        self.rRays.oAxis = self.rRays.oCosines[0,:]
        # Transmission
        Es_t, Ep_t = Es*self.ts, Ep*self.tp
        self.tRays.Ex =  Es_t*np.cos(-Theta) + Ep_t*np.sin(-Theta)
        self.tRays.Ey = -Es_t*np.sin(-Theta) + Ep_t*np.cos(-Theta)
        self.tRays.xCosines =  self.sCosines*np.cos(-Theta) + self.pCosines_t*np.sin(-Theta)
        self.tRays.yCosines = -self.sCosines*np.sin(-Theta) + self.pCosines_t*np.cos(-Theta)
        self.tRays.xAxis = self.tRays.xCosines[0,:]
        self.tRays.yAxis = self.tRays.yCosines[0,:]
        self.tRays.oAxis = self.tRays.oCosines[0,:]
        return
    # Coating
    def applyCoating(self, Coat):
        Coat.applyToSurface(self)
        self.rs, self.ts = np.copy(Coat.rs), np.copy(Coat.ts)
        self.rp, self.tp = np.copy(Coat.rp), np.copy(Coat.tp)
        return


    def propagateRays(self, Rays):
        self.computeIncidence(Rays)
        self.computeNormals()
        self.propagateReflectedRays()
        self.propagateTransmittedRays()
        self.propagatePolarization()
        return
    
    # Graphics
    def draw(self, Ax, **kwargs):
        x = np.reshape(self.X, newshape=(self.thetaRes,self.rRes))
        y = np.reshape(self.Y, newshape=(self.thetaRes,self.rRes))
        z = np.reshape(self.Z, newshape=(self.thetaRes,self.rRes))
        Ax.plot_surface(x, y, z, antialiased=True, **kwargs)
        return
    def drawRays(self, Ax, **kwargs):
        P1 = self.rRays.Points
        P2 = self.iRays.Points
        for i in range(len(P1)):
            if(self.iRays.Mask[i]):
                Ax.plot([P1[i,0], P2[i,0]], [P1[i,1], P2[i,1]], [P1[i,2], P2[i,2]], **kwargs)
        return
    def drawSurfaceNormals(self, Ax, Length, **kwargs):
        P0 = self.rRays.Points
        P1 = P0 - Length*0
        P2 = P0 + Length*self.nCosines
        for i in range(len(P1)):
            Ax.plot([P1[i,0], P2[i,0]], [P1[i,1], P2[i,1]], [P1[i,2], P2[i,2]], **kwargs)
        return
    def drawPolarizationDirection(self, Ax, LengthR, **kwargs): 
        P1 = self.rRays.Points 
        P2 = self.rRays.Points + LengthR*self.sCosines
        for i in range(len(P1)):
            Ax.plot([P1[i,0], P2[i,0]], [P1[i,1], P2[i,1]], [P1[i,2], P2[i,2]], 'b', **kwargs)
        P2 = self.rRays.Points + LengthR*self.pCosines_i
        for i in range(len(P1)):
            Ax.plot([P1[i,0], P2[i,0]], [P1[i,1], P2[i,1]], [P1[i,2], P2[i,2]], 'y', **kwargs)
        P2 = self.rRays.Points + LengthR*self.pCosines_r
        for i in range(len(P1)):
            Ax.plot([P1[i,0], P2[i,0]], [P1[i,1], P2[i,1]], [P1[i,2], P2[i,2]], 'g', **kwargs)
        P2 = self.rRays.Points + LengthR*self.pCosines_t
        for i in range(len(P1)):
            Ax.plot([P1[i,0], P2[i,0]], [P1[i,1], P2[i,1]], [P1[i,2], P2[i,2]], 'w', **kwargs)
        return