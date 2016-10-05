# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 10:24:00 2016

@author: leem
"""

import numpy as np

class Param:
    def __init__(self):
        self.n = 10                 # number of intervals
        self.deg = 4                # basis polynomial degree.
        
        self.m = self.n/self.deg    # number of master elements
        self.n = self.deg * self.m  # make n = m * deg (for simplicity)
        print "dimensions n is adjusted to " + str(self.n)
        
        self.xa = -1.0              # domain = [xa, xb]
        self.xb = 1.0
        self.tend = 2.0             # final time
        
def eval_f(x):
    return np.ones(x.size)
        
def eval_c(x):
    return np.ones(x.size)

def eval_b(x):
    return np.ones(x.size)

def eval_a(x):
    return np.ones(x.size)
