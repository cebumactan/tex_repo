# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 11:29:41 2016

@author: leem
"""

import numpy as np
import matplotlib.pyplot as plt

from header import *
from myglobals import *

class Element(object):
    def __init__(self,x):
        self.xk = x
        self.BFq = np.zeros((3,nqmax,mxdeg1))
        self.eMat= np.zeros((3,mxdeg1,mxdeg1))

    def LagranBFq(self):
        myxq = Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0]
        for iq in range(Q.nq):
            for iN in range(P.deg+1):
                
                # Lagrange shape function of iN type at xq
                val = np.ones(myxq.size)
                for i in range(P.deg+1):
                    if i==iN:
                        continue
                    val *= ((myxq-self.xk[i])/(self.xk[iN]-self.xk[i]))
                self.BFq[0,0:Q.nq,iN] = val
                
                # first derivative of Lagrange shape function of iN type at xq
                val0 = np.zeros(myxq.size)
                for i in range(P.deg+1):
                    if i==iN:
                        continue
                    val1 = np.ones(myxq.size)/(self.xk[iN]-self.xk[i])
                    for j in range(P.deg+1):
                        if (j==i) or (j==iN):
                            continue
                        val1*= ((myxq-self.xk[j])/(self.xk[iN]-self.xk[j]))
                    val0 +=val1
                self.BFq[1,0:Q.nq,iN]=val0  
                
                # second derivative of Lagrange shape function of iN type at xq
                val0 = np.zeros(myxq.size)
                for i in range(P.deg+1):
                    if i==iN:
                        continue
                    for j in range(P.deg+1):
                        if (j==i) or (j==iN):
                            continue
                        val1 = np.ones(myxq.size)/(self.xk[iN]-self.xk[i])/(self.xk[iN]-self.xk[j])
                        for k in range(P.deg+1):
                            if (k==i) or (k==j) or(k==iN):
                                continue
                            val1*= ((myxq-self.xk[k])/(self.xk[iN]-self.xk[k]))
                        val0 +=val1
                self.BFq[2,0:Q.nq,iN]=val0
        
    def eval_LagranB(self,iN):
        # Lagrange shape function of iN type at xq
        x = np.linspace(self.xk[0],self.xk[-1],30)
        val = np.ones(np.size(x))
        for i in range(P.deg+1):
            if i==iN:
                continue
            val *= (x-self.xk[i])/(self.xk[iN]-self.xk[i])
#        plt.clf()
        plt.plot(x,val,Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0],self.BFq[0,0:Q.nq,iN],'*')
#        plt.show()
        return val
    
    def eval_dLagranB(self,iN):
        # Lagrange shape function of iN type at self.xq
        x = np.linspace(self.xk[0],self.xk[-1],30)
        val0 = np.zeros(np.size(x))
        for i in range(P.deg+1):
            if i==iN:
                continue
            val1 = np.ones(np.size(x))/(self.xk[iN]-self.xk[i])
            for j in range(P.deg+1):
                if (j==i) or (j==iN):
                    continue
                val1*= ((x-self.xk[j])/(self.xk[iN]-self.xk[j]))
    
            val0 +=val1
#        plt.clf()
        plt.plot(x,val0,Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0],self.BFq[1,0:Q.nq,iN],'*')
#        plt.show()
        return val0
    
    def eval_ddLagranB(self,iN):
        # Lagrange shape function of iN type at xq
        x = np.linspace(self.xk[0],self.xk[-1],30)
        val0 = np.zeros(np.size(x))
        for i in range(P.deg+1):
            if i==iN:
                continue
            for j in range(P.deg+1):
                if (j==i) or (j==iN):
                    continue
                val1 = np.ones(np.size(x))/(self.xk[iN]-self.xk[i])/(self.xk[iN]-self.xk[j])
                for k in range(P.deg+1):
                    if (k==i) or (k==j) or(k==iN):
                        continue
                    val1*= ((x-self.xk[k])/(self.xk[iN]-self.xk[k]))
                val0 +=val1
#        plt.clf()
        plt.plot(x,val0,Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0],self.BFq[2,0:Q.nq,iN],'*')
#        plt.show()
        return val0

    def asm_eMat(self):
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


class Mesh(object):
    def __init__(self):
        self.N = P.n + 1      # Number of Nodes !!!
        self.x = np.zeros(MAXN)    
        print "Mesh : x[MAXN] has allocated."
       


        h = (P.xb-P.xa) / P.n
        for i in range(self.N):
            self.x[i]=P.xa + h*i
    
class FES:
    def __init__(self):
        self.dim = P.n + 1  # finite element space dimensions = # of nodes
        self.M = Mesh()     # define mesh with Pa.n intervals

        self.EL = [ Element(self.M.x[iEL*P.deg:(iEL+1)*P.deg+1]) for iEL in range(MAXM) ]
        print "Mesh : EL[MAXM] has allocated."
        self.m = P.m         # Number of Elements
        
        for i in range(self.m):
            self.EL[i].LagranBFq()

    def asm_GSH_Mat(self):
        print "FES asmGSH_Mat : not implemented"

# To verify this module, run this script.
if __name__ == "__main__":
    myV = FES()
    for i in range(P.m):
        myV.EL[i].eval_dLagranB(0)
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
        
