# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 11:29:41 2016

@author: leem
"""

import numpy as np
import matplotlib.pyplot as plt

from header import *
from myglobals import *

import systemdef

class Element(object):
    def __init__(self,ix,x):
        self.ixk = ix
        self.xk = x
#        self.dx = x[-1]-x[0]
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
        dx = (self.xk[-1]-self.xk[0])
        for iN in range(P.deg+1):
            for jN in range(P.deg+1):
                for iq in range(Q.nq):
                    self.eMatC[iN,jN]+= self.cq[iq]*self.BFq[0,iN,iq]*self.BFq[0,jN,iq]*Q.w[iq]*dx
                    self.eMatB[iN,jN]+= self.bq[iq]*self.BFq[1,iN,iq]*self.BFq[0,jN,iq]*Q.w[iq]*dx
                    self.eMatA[iN,jN]+= self.aq[iq]*self.BFq[1,iN,iq]*self.BFq[1,jN,iq]*Q.w[iq]*dx
                    self.eMatG[iN,jN]+= self.BFq[0,iN,iq]*self.BFq[0,jN,iq]*Q.w[iq]*dx
                
    def asm_eMat(self):
        dx = (self.xk[-1]-self.xk[0])
        for iN in range(P.deg+1):
            for jN in range(P.deg+1):
                for ider in range(3):
                    for jder in range(3):
                        for iq in range(Q.nq):
                            self.eMat[ider,jder,iN,jN]+= self.BFq[ider,iN,iq]*self.BFq[jder,jN,iq]*Q.w[iq]*dx
                            
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
#        plt.clf()
        plt.plot(x,val0,Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0],self.BFq[1,iN,0:Q.nq],'*')
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
        plt.plot(x,val0,Q.x * (self.xk[-1]-self.xk[0]) + self.xk[0],self.BFq[2,iN,0:Q.nq],'*')
#        plt.show()
        return val0

#    def asm_eMat(self):
#        print "not implemented"
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

        # uniform mesh
        h = (P.xb-P.xa) / P.n
        for i in range(self.N):
            self.x[i]=P.xa + h*i
    
class FES:
    def __init__(self):
        self.dim = P.n + 1  # finite element space dimensions = # of nodes
        self.M = Mesh()     # define mesh with Pa.n intervals

        self.EL = [ Element(np.arange(iEL*P.deg,(iEL+1)*P.deg+1),self.M.x[iEL*P.deg:(iEL+1)*P.deg+1]) for iEL in range(MAXM) ]
        print "Mesh : EL[MAXM] has allocated."
        
        for i in range(P.m):
            self.EL[i].LagranBFq()
            self.EL[i].fFq()
            self.EL[i].cFq()
            self.EL[i].bFq()
            self.EL[i].aFq()
#            self.EL[i].asm_eMat()
            self.EL[i].asm_eMatCBAG()            
            
        self.noffdiag = P.deg
        self.u = self.l = self.noffdiag
        # system matrices in band stroage and rhs vector
        self.Abn = np.zeros(shape=(self.u+self.l+1,self.dim))
        self.Bbn = np.zeros(shape=(self.u+self.l+1,self.dim))
        self.Cbn = np.zeros(shape=(self.u+self.l+1,self.dim))
        self.Gbn = np.zeros(shape=(self.u+self.l+1,self.dim))
        self.fn  = np.zeros(self.dim)
        
    
    def asm_rhs(self):
        for i in range(P.m):
            dx = self.EL[i].xk[-1]-self.EL[i].xk[0]
            for iN in range(P.deg+1):
                k = self.EL[i].ixk[0] + iN      # vector index k in 0..dim-1
                for iq in range(Q.nq):
                    self.fn[k]+=self.EL[i].fq[iq]*self.EL[i].BFq[0,iN,iq]*Q.w[iq]*dx
    def asm_CMat_band(self):
        CC=np.zeros(shape=(self.dim,self.dim))
        for i in range(P.m):
            dx = self.EL[i].xk[-1]-self.EL[i].xk[0]
            for kN in range(P.deg+1):
                K = self.EL[i].ixk[0] + kN 
                Z = self.u-K
                for jN in range( max(0,kN-self.u),min(P.deg+1,kN+1+self.l)):
                    J = self.EL[i].ixk[0] + jN
                    self.Cbn[Z+J,K]+=self.EL[i].eMat[0,0,jN,kN]
                    CC[J,K]+=self.EL[i].eMat[0,0,jN,kN]
        print np.array_str(CC, precision=2)
        print np.array_str(self.Cbn, precision=2,suppress_small=True)  
  
    def asm_BMat_band(self):
        BB=np.zeros(shape=(self.dim,self.dim))
        for i in range(P.m):
            dx = self.EL[i].xk[-1]-self.EL[i].xk[0]
            for kN in range(P.deg+1):
                K = self.EL[i].ixk[0] + kN 
                Z = self.u-K
                for jN in range( max(0,kN-self.u),min(P.deg+1,kN+1+self.l)):
                    J = self.EL[i].ixk[0] + jN
                    self.Bbn[Z+J-1,K]+=self.EL[i].eMat[1,0,jN,kN]
                    BB[J,K]+=self.EL[i].eMat[1,0,jN,kN]
        print np.array_str(BB, precision=2)
        print np.array_str(self.Bbn, precision=2,suppress_small=True)  
        
    def asm_AMat_band(self):
        AA=np.zeros(shape=(self.dim,self.dim))
        for i in range(P.m):
            dx = self.EL[i].xk[-1]-self.EL[i].xk[0]
            for kN in range(P.deg+1):
                K = self.EL[i].ixk[0] + kN 
                Z = self.u-K
                for jN in range( max(0,kN-self.u),min(P.deg+1,kN+1+self.l)):
                    J = self.EL[i].ixk[0] + jN
                    self.Abn[Z+J,K]+=self.EL[i].eMat[1,1,jN,kN]
                    AA[J,K]+=self.EL[i].eMat[1,1,jN,kN]
        print np.array_str(AA, precision=2)
        print np.array_str(self.Abn, precision=2,suppress_small=True)        


# To verify this module, run this script.
if __name__ == "__main__":
#    clear_all()
    myV = FES()
    for i in range(P.m):
        myV.EL[i].eval_LagranB(1)
        
    myV.asm_rhs()
#    myV.asm_AMat_band()
#    myV.asm_BMat_band()
#    myV.asm_CMat_band()
#    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
        
