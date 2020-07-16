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

class CoatedMirror(Surface):
    def __ini__(self, Dia, **kwargs):
        Surface.__init__(self, Dia, **kwargs)
        Surface.Mirror = True
        return
    
    def applyCoating(self, Coat):
        Coat.applyToSurface(self)
        self.rs, self.ts = np.copy(Coat.rs), np.copy(Coat.ts)
        self.rp, self.tp = np.copy(Coat.rp), np.copy(Coat.tp)
        return
    
    
class UncoatedLens():
    def __init__(self, Dia, Thick, R1=np.inf, R2=np.inf, K1=0.0, K2=0.0, n=1.5+0.0j):
        self.S1 = Surface(Dia, R=R1, K=K1, n1=1.0+0.0j, n2=n)
        self.S2 = Surface(Dia, R=R2, K=K2, n1=n, n2=1.0+0.0j)
        S2.translateOrigin(z=Thick)