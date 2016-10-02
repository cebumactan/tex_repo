# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 19:13:58 2016

@author: leem
"""

#===================================================
# my modules
#===================================================

def clear_all():
    """Clears all the variables from the workspace of the spyder application."""
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue

        del globals()[var]

#===================================================
# end of my modules
#===================================================