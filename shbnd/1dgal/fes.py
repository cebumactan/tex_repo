# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 11:29:41 2016

@author: leem
"""
import scipy.linalg as splin
import numpy as np
import matplotlib.pyplot as plt

from header import *
#from myglobals import *

import systemdef, quad1d

P = systemdef.Param()
Q = quad1d.Quad1d(P.deg)

class Element(object):
    def __init__(self,ix,x):
        self.ixk = ix
        self.xk = x
        self.BFq = np.zeros((3,mxdeg1,nqmax))
        self.eMatC = np.zeros((mxdeg1,mxdeg1))
        self.eMatB = np.zeros((mxdeg1,mxdeg1))
        self.eMatA = np.zeros((mxdeg1,mxdeg1))
        self.eMatG = np.zeros((mxdeg1,mxdeg1))
#        self.eMat= np.zeros((3,3,mxdeg1,mxdeg1))
        
        self.fq = np.zeros(nqmax)
        self.cq = np.zeros(nqmax)
        self.bq = np.zeros(nqmax)
        self.aq = np.zeros(nqmax)

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
                self.BFq[0,iN,0:Q.nq] = val
                
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
                self.BFq[1,iN,0:Q.nq]=val0  
                
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
                self.BFq[2,iN,0:Q.nq]=val0
        
    def fFq(self):
        myxq = Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0]
        for iq in range(Q.nq):
            self.fq[0:Q.nq] = systemdef.eval_f(myxq)
    def cFq(self):
        myxq = Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0]
        for iq in range(Q.nq):
            self.cq[0:Q.nq] = systemdef.eval_c(myxq)

    def bFq(self):
        myxq = Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0]
        for iq in range(Q.nq):
            self.bq[0:Q.nq] = systemdef.eval_b(myxq)
    def aFq(self):
        myxq = Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0]
        for iq in range(Q.nq):
            self.aq[0:Q.nq] = systemdef.eval_a(myxq)
            
    def asm_eMatCBAG(self):
        for iN in range(P.deg+1):
            for jN in range(P.deg+1):
                for iq in range(Q.nq):
                    self.eMatC[iN,jN]+= self.cq[iq]*self.BFq[0,iN,iq]*self.BFq[0,jN,iq]*Q.w[iq]*(self.xk[-1]-self.xk[0])
                    self.eMatB[iN,jN]+= self.bq[iq]*self.BFq[1,iN,iq]*self.BFq[0,jN,iq]*Q.w[iq]*(self.xk[-1]-self.xk[0])
                    self.eMatA[iN,jN]+= self.aq[iq]*self.BFq[1,iN,iq]*self.BFq[1,jN,iq]*Q.w[iq]*(self.xk[-1]-self.xk[0])
                    self.eMatG[iN,jN]+= self.BFq[0,iN,iq]*self.BFq[0,jN,iq]*Q.w[iq]*(self.xk[-1]-self.xk[0])
                
    def asm_eMat(self):
        for iN in range(P.deg+1):
            for jN in range(P.deg+1):
                for ider in range(3):
                    for jder in range(3):
                        for iq in range(Q.nq):
                            self.eMat[ider,jder,iN,jN]+= self.BFq[ider,iN,iq]*self.BFq[jder,jN,iq]*Q.w[iq]*(self.xk[-1]-self.xk[0])
                            
    def eval_LagranB(self,iN):
        # Lagrange shape function of iN type at xq
        x = np.linspace(self.xk[0],self.xk[-1],30)
        val = np.ones(np.size(x))
        for i in range(P.deg+1):
            if i==iN:
                continue
            val *= (x-self.xk[i])/(self.xk[iN]-self.xk[i])
#        plt.clf()
        plt.plot(x,val,Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0],self.BFq[0,iN,0:Q.nq],'*')
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
        plt.plot(x,val0,Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0],self.BFq[1,iN,0:Q.nq],'*')
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
        plt.plot(x,val0,Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0],self.BFq[2,iN,0:Q.nq],'*')
        return val0

class Mesh(object):
    def __init__(self):
        self.N = P.n + 1      # Number of Nodes !!!
        self.x = np.zeros(MAXN)    
        print "MESH ARRAY x[MAXN] is allocated."

        # uniform mesh
        h = (P.xb-P.xa) / P.n
        for i in range(self.N):
            self.x[i]=P.xa + h*i
    
class FES:
    def __init__(self):
        self.dim = P.n + 1  # finite element space dimensions = # of nodes
        self.M = Mesh()     # define mesh with Pa.n intervals

        self.EL = [ Element(np.arange(iEL*P.deg,(iEL+1)*P.deg+1),self.M.x[iEL*P.deg:(iEL+1)*P.deg+1]) for iEL in range(MAXM) ]
        print "ELEMENTS EL[MAXM] are allocated."

        self.noffdiag = P.deg               # number of off-diagonal = P.deg
        self.iu = self.il = self.noffdiag   # upper offdiag = lower offdiag
        
        # system matrices in band stroage and rhs vector
        self.Abn = np.zeros(shape=(self.iu+self.il+1,self.dim))
        self.Bbn = np.zeros(shape=(self.iu+self.il+1,self.dim))
        self.Cbn = np.zeros(shape=(self.iu+self.il+1,self.dim))
        self.Gbn = np.zeros(shape=(self.iu+self.il+1,self.dim))
        self.SYSbn=np.zeros(shape=(self.iu+self.il+1,self.dim))
        self.fn  = np.zeros(self.dim)
        print "BAND MATRICES, RHS VECTOR are allocated."
        
        # this full matrices will be obsolete in future       
        self.CC=np.zeros(shape=(self.dim,self.dim))
        self.BB=np.zeros(shape=(self.dim,self.dim))
        self.AA=np.zeros(shape=(self.dim,self.dim))
        self.GG=np.zeros(shape=(self.dim,self.dim))
        self.SS=np.zeros(shape=(self.dim,self.dim))
        self.ff=np.zeros(self.dim)
        self.mySbn=np.zeros(shape=(self.iu+self.il+1,self.dim))
        
        for i in range(P.m):
            self.EL[i].LagranBFq()
            self.EL[i].fFq()
            self.EL[i].cFq()
            self.EL[i].bFq()
            self.EL[i].aFq()
#            self.EL[i].asm_eMat()
            self.EL[i].asm_eMatCBAG()  
        print "ELEMENTS INITIALIZED: LOCAL BILINEARS STORED."
    
    def asm_rhs(self):
        for i in range(P.m):
            dx = self.EL[i].xk[-1]-self.EL[i].xk[0]
            for iN in range(P.deg+1):
                k = self.EL[i].ixk[0] + iN      # vector index k in 0..dim-1
                for iq in range(Q.nq):
                    self.fn[k]+=self.EL[i].fq[iq]*self.EL[i].BFq[0,iN,iq]*Q.w[iq]*dx
        print "SOURCE VECTOR (DIM) STORED"
                    
    def asm_CBAGMat_band(self):
        for i in range(P.m):
            dx = self.EL[i].xk[-1]-self.EL[i].xk[0]
            for kN in range(P.deg+1):
                K = self.EL[i].ixk[0] + kN 
                Z = self.iu-K
                for jN in range( max(0,kN-self.iu),min(P.deg+1,kN+1+self.il)):
                    J = self.EL[i].ixk[0] + jN
                    
                    self.Cbn[Z+J,K]+=self.EL[i].eMatC[jN,kN]
                    self.CC[J,K]+=self.EL[i].eMatC[jN,kN] 
                    
                    self.Bbn[Z+J,K]+=self.EL[i].eMatB[jN,kN]
                    self.BB[J,K]+=self.EL[i].eMatB[jN,kN] 
                    
                    self.Abn[Z+J,K]+=self.EL[i].eMatA[jN,kN]
                    self.AA[J,K]+=self.EL[i].eMatA[jN,kN] 
                    
                    self.Gbn[Z+J,K]+=self.EL[i].eMatG[jN,kN]
                    self.GG[J,K]+=self.EL[i].eMatG[jN,kN] 
        print "MATRICES (DIMxDIM) in BAND STORAGE STORED"
                    
    def set_BdryCond(self):
        self.ff[:] = self.fn[:]        
 
        self.SS[0,0] = 1.0
        self.ff[0] = P.ua
        self.SS[0,1:1+self.iu] = 0.0
        self.ff[1:1+self.il] -= P.ua*self.SS[1:1+self.il,0]
        self.SS[1:1+self.il,0] = 0.0

        n = self.dim-1

        self.SS[n,n] = 1.0
        self.ff[n] = P.ub
        self.SS[n,n-self.il:n] = 0.0
        self.ff[n-self.iu:n] -= P.ub*self.SS[n-self.iu:n,n]
        self.SS[n-self.iu:n,n] = 0.0

#        self.SS[-1,-1] = 1.0
#        self.fn[-1] = P.ub
#        self.SS[-1,-1-self.il:-1] = 0.0
#        self.fn[-1-self.iu:-1] -= P.ub*self.SS[-1-self.iu:-1,-1]
#        self.SS[-1-self.iu:-1,-1] = 0.0

        for K in range(self.dim):
            Z = self.iu-K
            for J in range( max(0,K-self.iu),min(self.dim,K+1+self.il)):
                self.mySbn[Z+J,K] =self.SS[J,K]
                
        self.SYSbn[self.iu,0] = 1.0
        self.fn[0] = P.ua
        for K in range(1,1+self.iu):
            self.SYSbn[self.iu-K , K] = 0.0        
#        self.SYSbn[(self.iu-1):0:-1 , 1:1+self.iu] = 0.0        
        for J in range(1,1+self.il):
            self.fn[J] -= P.ua*self.SYSbn[self.iu+J,0]
            self.SYSbn[self.iu+J,0]=0.0
#        self.ff[1:1+self.il] -= P.ua*self.SYSbn[(self.iu+1):(self.iu+1+self.il),0]
#        self.SYSbn[(self.iu+1):(self.iu+1+self.il),0]

        self.SYSbn[self.iu,n]=1.0
        self.fn[n] = P.ub
        for K in range(n-self.il,n):
            self.SYSbn[self.iu+n-K , K] = 0.0 
#        self.SYSbn[(self.iu+self.il):self.iu:-1,n] = 0.0
        for J in range(n-self.iu,n):
            self.fn[J] -= P.ub*self.SYSbn[self.iu+J-n,n]
            self.SYSbn[self.iu+J-n,n] = 0.0        
#        self.ff[n-self.iu:n] -= P.ub*self.SYSbn[0:self.iu,n]
#        self.SYSbn[0:self.iu,n] = 0.0
        print "DIRICHILET BDRY COND REFLECTED"




# To verify this module, run this script.
if __name__ == "__main__":
    
    # try construct FES object
    myV = FES()
    
    # check quadrature and evaluation of Lagran basis upto second-derivative
    for i in range(P.m):
        myV.EL[i].eval_LagranB(1)
    
    # check values of source term
    myV.asm_rhs()
#    print np.array_str(myV.fn, precision=2)
    
    # check reaction term (C), advection term(B), diffusion term(A), grammian(G)
    # First, check values in the full matrix form
    # Second, check if band form coincides
    myV.asm_CBAGMat_band()
    myV.SS = myV.AA + myV.BB + myV.CC
    myV.SYSbn = myV.Abn + myV.Bbn + myV.Cbn

#    print np.array_str(myV.CC, precision=2)
#    print np.array_str(myV.Cbn, precision=2,suppress_small=True)
#    print np.array_str(myV.BB, precision=2)
#    print np.array_str(myV.Bbn, precision=2,suppress_small=True)
#    print np.array_str(myV.AA, precision=2)
#    print np.array_str(myV.Abn, precision=2,suppress_small=True)    
#    print np.array_str(myV.GG, precision=2)
#    print np.array_str(myV.Gbn, precision=2,suppress_small=True) 
#    print np.array_str(myV.SS, precision=2)
#    print np.array_str(myV.SYSbn, precision=2,suppress_small=True)

    # check treatment for dirichilet bdry condition
    # First, check values in the full matrix form
    # second, check if band form coincides,
    myV.set_BdryCond()
#    print np.array_str(myV.ff, precision=2)
#    print np.array_str(myV.SS, precision=1)
#    print np.array_str(myV.mySbn, precision=1,suppress_small=True)  

#    print np.array_str(myV.fn, precision=2)
#    print np.array_str(myV.SYSbn, precision=1,suppress_small=True)  
#
#    print np.array_str(myV.SYSbn-myV.mySbn, precision=1,suppress_small=True)   
    
    # 1D poisson    
    un = splin.solve_banded((myV.il,myV.iu),myV.SYSbn,myV.fn)
    print un
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
        
