# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 10:30:51 2016

@author: leem
"""
clear_all()

from numpy import *

from header import *
import quad1d, param



# define global variables
dBF = zeros((3,mxdeg1))
BFq = zeros((3,nqmax,mxdeg1))
eMat = zeros((3,mxdeg1,mxdeg1))
P = param.Param()
XN = zeros(MAXN)
Y = zeros(MAXN)
Q = quad1d.Quad1d(P.deg)
#Q = quad1d.Quad1d(P.deg)




#def THRE(a):
#    return ( abs )


# define global variables
#dbF = zeros(3,mxdeg1)

