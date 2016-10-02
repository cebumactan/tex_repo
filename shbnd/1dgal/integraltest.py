# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 18:48:22 2016

@author: leem
"""

# my mudules
import numpy as np
import clear_all
import quad1d


clear_all.clear_all()   # clear all variables in workspace memory

q = quad1d.Quad1d(11)   # create a quadrature instance with idegree
    
# 23th polynomial 1+x^23
kn = np.zeros(23+1)
kn[0]=1
kn[7]=1
kn[9]=-1
kn[23]=1
#kn = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
a = q.myintegral(kn,-2,3)
print a