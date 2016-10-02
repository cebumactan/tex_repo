# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 16:50:23 2016

@author: leem
"""

import param, quad1d

# In this module, define global valriables

#dBF = zeros((3,mxdeg1))
#BFq = zeros((3,nqmax,mxdeg1))
#eMat = zeros((3,mxdeg1,mxdeg1))
#XN = zeros(MAXN)
#Y = zeros(MAXN)

# allocate / initialize objects
P = param.Param()
Q = quad1d.Quad1d(P.deg)

# allocate memory for mesh, elements of size MAXN, and initialize with the parameter
#V = fes.FES(P) 