import numpy as np
import copy as cp
import random as rd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime as dt
#
from .Source import *
from .Surface import *
#

class Coating():
    def __init__(self, RI, Thick):
        if (len(RI) != len(Thick)):
            print('Error! Both inputs must have the same length!')
            return
        self.RI = np.copy(RI)
        self.Thick = np.copy(Thick)
    #
    def applyToSurface(self, Surf):                                
        self.Wavelength = np.copy(Surf.iRays.Wavelength)
        self.iRI, self.sRI = np.copy(Surf.iRI), np.copy(Surf.tRI)
        self.iTheta, self.sTheta = np.copy(Surf.iTheta), np.copy(Surf.tTheta)
        
        # for s-polarization
        iEta, sEta = self.iRI*np.cos(self.iTheta), self.sRI*np.cos(self.sTheta)
        m_00, m_01, m_10, m_11 = 1, 0, 0, 1
        for RI_j, Thick_j in zip(self.RI, self.Thick):
            Theta_j = np.arcsin(self.iRI/RI_j)*np.sin(self.iTheta)
            Eta_j = RI_j*np.cos(Theta_j)
            Delta_j = RI_j*Thick_j*np.cos(Theta_j)*2*np.pi/self.Wavelength
            temp00 = m_00*np.cos(Delta_j) + m_01*1j*Eta_j*np.sin(Delta_j)
            temp01 = m_00*1j*np.sin(Delta_j)/Eta_j + m_01*np.cos(Delta_j)
            temp10 = m_10*np.cos(Delta_j) + m_11*1j*Eta_j*np.sin(Delta_j)
            temp11 = m_10*1j*np.sin(Delta_j)/Eta_j + m_11*np.cos(Delta_j)
            m_00, m_01 = np.copy(temp00), np.copy(temp01)
            m_10, m_11 = np.copy(temp10), np.copy(temp11)
        E, H = m_00 + m_01*sEta, m_10 + m_11*sEta
        self.rs, self.ts = (iEta*E - H)/(iEta*E + H), 2*iEta/(iEta*E + H)
        
        # for p-polarization
        iEta, sEta = self.iRI/np.cos(self.iTheta), self.sRI/np.cos(self.sTheta)
        m_00, m_01, m_10, m_11 = 1, 0, 0, 1
        for RI_j, Thick_j in zip(self.RI, self.Thick):
            Theta_j = np.arcsin(self.iRI/RI_j)*np.sin(self.iTheta)
            Eta_j = RI_j/np.cos(Theta_j)
            Delta_j = RI_j*Thick_j*np.cos(Theta_j)*2*np.pi/self.Wavelength
            temp00 = m_00*np.cos(Delta_j) + m_01*1j*Eta_j*np.sin(Delta_j)
            temp01 = m_00*1j*np.sin(Delta_j)/Eta_j + m_01*np.cos(Delta_j)
            temp10 = m_10*np.cos(Delta_j) + m_11*1j*Eta_j*np.sin(Delta_j)
            temp11 = m_10*1j*np.sin(Delta_j)/Eta_j + m_11*np.cos(Delta_j)
            m_00, m_01 = np.copy(temp00), np.copy(temp01)
            m_10, m_11 = np.copy(temp10), np.copy(temp11)
        E, H = m_00 + m_01*sEta, m_10 + m_11*sEta
        self.rp, self.tp= (iEta*E - H)/(iEta*E + H), 2*iEta/(iEta*E + H)
                                            
        
        