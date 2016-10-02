# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 18:27:19 2016

@author: leem
"""
import numpy as np
import scipy.sparse as spa

from header import *
from myglobals import *

# constitute system

class Mysystem:
    def __init__(self,V):
        self.u = np.zeros(MAXN)
        print "u has allocated"
        self.f = np.zeros(MAXN)
        print "f has allocated"
        self.buildSysMat()
        self.quad_fBF(V.M)
        
    def eval_f(self,x):
        return 1+x
        
    def quad_fBF(self,M):
#        global Q
        n = M.N
        e = fes.Element()
        
        for k in range(n):
            e = M.EL[k]
            self.f[k] = 0
            dx = M.x[k+1]-M.x[k]
            for iq in range(Q.nq):
                xq = M.x[k] + dx * Q.x[iq]
                dxw = dx * Q.w[iq]
                for j in range(P.deg+1):
                    self.f[k] += dxw * self.eval_f(xq) * e.BFq[0][iq][j]
        print "quad_fBF : not implemented"
    
    def buildSysMat(self):
        print "buildSysMat : not implemented"

    def solve(self):
        print "solve : not implemente"
        
#        mysum = 0
#        for n in range(0,self.nq):
#            polyvalue = 0
#            for k in range(0,np.size(kn)):
#                polyvalue+= kn[k]*( ( (b-a)*self.x[n] + a )**k ) 
#            mysum+=self.w[n]*polyvalue
#        mysum = mysum*(b-a)
#        return mysum