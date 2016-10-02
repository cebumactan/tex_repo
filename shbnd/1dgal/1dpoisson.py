# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 10:30:51 2016

@author: leem
"""
clear_all()

#################################
from numpy import *

from header import *

##################################



##################################
# define global variables
dBF = zeros((3,mxdeg1))
BFq = zeros((3,nqmax,mxdeg1))
eMat = zeros((3,mxdeg1,mxdeg1))
XN = zeros(MAXN)
Y = zeros(MAXN)

# allocate / initialize objects
P = param.Param()
Q = quad1d.Quad1d(P.deg)

# allocate memory for mesh, elements of size MAXN, and initialize with the parameter
V = FES(P) 
###################################



###################################
# start program

# Initialize
#P.t = 0.0

#V.asm_Local_Mat()

