# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 10:52:56 2016

@author: leem
"""

# DO NOT MODIFY THIS FILE AND DEFINITIONS

# Used to be #define
MAXN    = 20050
DEG     = 1
nqmax   = DEG+1
mxdeg   = DEG
mxdeg1  = DEG+1
mxdeg2  = 3*DEG+2
DEPS    = 1.11022302462515654E-016
SEPS    = 5.96046448E-08
EPS     = 1.0E-12

# Used to be macro
def THRE(a):
    return ( 0.0 if (abs(a)<EPS) else a )
def THRE2(a):
    return ( 0.0 if (abs(a)<EPS*EPS) else a )   
def DTHRE2(a):
    return ( 0.0 if (abs(a)<DEPS*DEPS) else a )
def ETHRE(a):
    return ( EPS if (abs(a)<EPS) else a )   


from fes import *
import quad1d, param