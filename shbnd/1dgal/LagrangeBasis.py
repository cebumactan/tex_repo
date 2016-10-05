# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 19:07:42 2016

@author: leem
"""

#clear_all()

#import param, quad1d
import numpy as np
import matplotlib.pyplot as plt
#import time
from header import *
from myglobals import *


xk = np.arange(P.deg+1)
BFq = np.zeros((3,nqmax,mxdeg1))
eMat= np.zeros((3,mxdeg1,mxdeg1))
myx = Q.x * P.deg

def LagrangeBq(xk):
    for iq in range(Q.nq):
        for iN in range(P.deg+1):
            
            # Lagrange shape function of iN type at xq
            val = np.ones(size(myx))
            for i in range(P.deg+1):
                if i==iN:
                    continue
                val *= ((myx-xk[i])/(xk[iN]-xk[i]))
            BFq[0,0:Q.nq,iN] = val
            
            # first derivative of Lagrange shape function of iN type at xq
            val0 = np.zeros(size(myx))
            for i in range(P.deg+1):
                if i==iN:
                    continue
                val1 = np.ones(size(myx))/(xk[iN]-xk[i])
                for j in range(P.deg+1):
                    if (j==i) or (j==iN):
                        continue
                    val1*= ((myx-xk[j])/(xk[iN]-xk[j]))
                val0 +=val1
            BFq[1,0:Q.nq,iN]=val0  
            
            # second derivative of Lagrange shape function of iN type at xq
            val0 = np.zeros(size(myx))
            for i in range(P.deg+1):
                if i==iN:
                    continue
                for j in range(P.deg+1):
                    if (j==i) or (j==iN):
                        continue
                    val1 = np.ones(size(myx))/(xk[iN]-xk[i])/(xk[iN]-xk[j])
                    for k in range(P.deg+1):
                        if (k==i) or (k==j) or(k==iN):
                            continue
                        val1*= ((myx-xk[k])/(xk[iN]-xk[k]))
                    val0 +=val1
            BFq[2,0:Q.nq,iN]=val0

def eval_LagrangeB(x,xk,iN):
    # Lagrange shape function of iN type at xq
    val = np.ones(size(x))
    for i in range(P.deg+1):
        if i==iN:
            continue
        val *= (x-xk[i])/(xk[iN]-xk[i])
    return val

def eval_dLagrangeB(x,xk,iN):
    # Lagrange shape function of iN type at xq
    val0 = np.zeros(size(x))
    for i in range(P.deg+1):
        if i==iN:
            continue
        val1 = np.ones(size(x))/(xk[iN]-xk[i])
        for j in range(P.deg+1):
            if (j==i) or (j==iN):
                continue
            val1*= ((x-xk[j])/(xk[iN]-xk[j]))

        val0 +=val1
    return val0

def eval_ddLagrangeB(x,xk,iN):
    # Lagrange shape function of iN type at xq
    val0 = np.zeros(size(x))
    for i in range(P.deg+1):
        if i==iN:
            continue
        for j in range(P.deg+1):
            if (j==i) or (j==iN):
                continue
            val1 = np.ones(size(x))/(xk[iN]-xk[i])/(xk[iN]-xk[j])
            for k in range(P.deg+1):
                if (k==i) or (k==j) or(k==iN):
                    continue
                val1*= ((x-xk[k])/(xk[iN]-xk[k]))
            val0 +=val1
    return val0



LagrangeBq(xk)
x= np.arange(0,P.deg+0.01,0.01)


plt.plot(x,np.zeros(size(x)))
#plt.plot(x, eval_LagrangeB(x,xk,0))
#plt.plot(x, 0.25*x*(x-1)*(x-3)*(x-4),'o')
plt.plot(x, eval_dLagrangeB(x,xk,1),'o-')
#plt.plot(x, x**3-6*x**2+9.5*x-3,'.')
#plt.plot(x, eval_ddLagrangeB(x,xk,1),'*-')
#plt.plot(x, 3*x**2-12*x+9.5,'+')
#for iN in range(P.deg+1):
#plt.plot(x,np.zeros(size(x)))
#plt.plot(myx,BFq[0,0:Q.nq,4],'*')
#plt.plot(x, eval_LagrangeB(x,xk,4))

        
#for iN in range(P.deg+1):
#plt.plot(x,np.zeros(size(x)))
#plt.plot(myx,BFq[1,0:Q.nq,4],'*')
#plt.plot(x, eval_dLagrangeB(x,xk,4))


#for iN in range(P.deg+1):
#plt.plot(x,np.zeros(size(x)))
#plt.plot(myx,BFq[2,0:Q.nq,1],'*')
#plt.plot(x, eval_ddLagrangeB(x,xk,1))

    