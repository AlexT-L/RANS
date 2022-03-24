"""
Description

Expands Field objects from coarser meshes to finer meshes.

Libraries/Modules

bin.Field \n
numpy
"""
from bin.Field import Field
import numpy as np
import scipy as sp
from scipy.interpolate import interpn

def bilinear4way(coarse, fine):
    """Bilinear interpolation function for expanding Fields. Currently does not work.
    
    Args:
    
    coarse:
        Field on a coarser mesh.
    fine:
        Field on a finer mesh.
    """
    # nx, ny, dim = fine.shape # Doesn't work, fine.shape only returns two items.
    # #nx, ny = fine.shape # This works
    nxf, nyf = Field.size(fine)

    fine[:] = 0

    fine[1:nxf:2, 1:nyf:2] = coarse
    
    for i in range(2,nxf,2):
        for j in range(1,nyf,2):
            fine[i  ,j] = .25*fine[i-1,j]  +.75*fine[i+1,j]
            fine[i-1,j] = .75*fine[i-1,j]  +.25*fine[i+1,j]
    
    
    for j in range(2,nyf,2):
        for i in range(1,nxf,2):
            fine[i,  j] = .25*fine[i,j-1]  +.75*fine[i,j+1]
            fine[i,j-1] = .75*fine[i,j-1]  +.25*fine[i,j+1]
