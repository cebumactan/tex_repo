# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 11:29:41 2016

@author: leem
"""

import numpy as np

from header import *
from myglobals import *

class Element:
    def __init__(self):
        self.xl = np.zeros(1)
        self.xr = np.zeros(1)
        self.BFq = np.zeros((3,nqmax,mxdeg1))
        self.eMat= np.zeros((3,mxdeg1,mxdeg1))

class Mesh:
    def __init__(self):
        # allocate objects
        self.x = np.zeros(MAXN)    
        print "Mesh : x[MAXN] has allocated."
        
        self.EL = [ Element() for i in range(MAXN) ]
        print "Mesh : EL[MAXN] has allocated."
        
        self.N = P.n + 1      # Number of Nodes !!!

        h = (P.xb-P.xa) / P.n
        for i in range(self.N):
            self.x[i]=P.xa + h*i
    
class FES:
    def __init__(self):
        self.dim = P.n     # finite element space dimension
        self.M = Mesh()   # define mesh with Pa.n intervals
        self.asm_GSH_Mat()  # assemble sparse matrix of G, S, H
        
        
    def LagranBFq(self,deg,nderiv,xq):
        print "not implemented"

    def asm_Local_mat(self):
        print "not implemented"
#        n = self.M.n - 1 
#        nderiv = min(P.deg+1,3)
#        e=Element()
#        for k in range(n):
#            e.xl = self.M.x[k]
#            e.xr = self.M.x[k+1]
#            dx = self.M.x[k+1] - self.M.x[k]
#            jleft = k + P.deg
#            
#            for iq in range(Q.nq):
#                xq = dx*Q.x[iq] + self.M.x[k]
#                dBLagran(P.deg,nderiv,xq,jleft)
#                for nd in range(nderiv):
#                    for j in range(P.deg+1):
#                        e.BFq[nd,i,j] = dBF[nd,j]      
        
#        self.G = LaMat()
#        self.S = LaMat()
    def asm_GSH_Mat(self):
        print "FES asmGSH_Mat : not implemented"


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
        
