# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 11:29:41 2016

@author: leem
"""

import numpy as np
from header import *

class Element:
    def __init__(self):
        self.xl = np.zeros(1)
        self.xr = np.zeros(1)
        self.BFq = np.zeros((3,nqmax,mxdeg1))
        self.eMat= np.zeros((3,mxdeg1,mxdeg1))

class Mesh:
    def __init__(self,Pa):
        # allocate objects
        self.x = np.zeros(MAXN)        
        self.EL = [ Element() for i in range(MAXN) ]
        
        self.n = Pa.n + 1      # Number of Nodes !!!
       
        h = (Pa.xb-Pa.xa) / Pa.n
        for i in range(0,self.n):
            self.x[i]=Pa.xa + h*i
       
        
#class LaMat:
#    def __init__(self):
#        self.lda = np.zeros(1, dtype=int)
#        self.n = np.zeros(1, dtype=int)
#        self.cpa = np.zeros(1, dtype=int)
#        self.cpc = np.zeros(1, dtype=int)

class FES:
    def __init__(self,Pa):
        self.dim = Pa.n
        self.M = Mesh(Pa)
#        self.G = LaMat()
#        self.S = LaMat()
#    def SysMat(self):
    
#    def asm_Local_Mat(self):
        
    
        
