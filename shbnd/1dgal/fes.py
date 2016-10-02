# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 11:29:41 2016

@author: leem
"""

from numpy import *
#from header import *

class Element:
    def __init__(self):
        self.xl = zeros(1)
        self.xr = zeros(1)
        self.BFq = zeros((3,nqmax,mxdeg1))
        self.eMat= zeros((3,mxdeg1,mxdeg1))

class Mesh:
    def __init__(self):
        self.n = zeros(1)
        self.x = zeros(1)
        self.EL = Element() 
        
class LaMat:
    def __init__(self):
        self.lda = zeros(1, dtype=int)
        self.n = zeros(1, dtype=int)
        self.cpa = zeros(1, dtype=int)
        self.cpc = zeros(1, dtype=int)

class Fes:
    def __init__(self):
        self.dim = zeros(1, dtype=int)
        self.M = Mesh()
        self.LaMat G = Mat()
        self.LaMat S = Mat()
        
