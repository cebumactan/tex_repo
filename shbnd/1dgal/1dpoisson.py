# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 10:30:51 2016

@author: leem
"""
clear_all()


import numpy as np
import matplotlib.pyplot as pl

from header import *
from myglobals import *

# create finite element space
Vn = fes.FES()                 

# assemble system operator, rhs
SYSn = mysystem.Mysystem(Vn)  



#SYS.buildSysMat()
#SYS.buildfquad()
#SYSn.solve()
#pl.plot(Vn.M.x[0:Vn.M.N], SYSn.f[0:Vn.M.N],'o')
#pl.plot(Vn.M.x[0:Vn.M.n], SYSn.u[0:Vn.M.n])


