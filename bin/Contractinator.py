"""
Description

Contracts Field objects from finer meshes to coarser meshes.

Libraries/Modules

bin.Field \n
numpy
"""

from bin.Field import Field, sum
import numpy as np

def simple(fine, coarse):
    """Performs a simple contraction where every other value is deleted.
    
    Args:
    
    fine:
        The Field object on the finer grid
    coarse:
        The Field object on the coarser grid
    """
    if type(fine) is not Field or type(coarse) is not Field:
        raise TypeError('Fine or coarse field is not a field')
    
    x_fine, y_fine, dim = 0, 0, 0
    if(len(fine.shape()) == 3):
        x_fine, y_fine, dim = fine.shape()
    if(len(fine.shape()) == 2):
        x_fine, y_fine = fine.shape()
    if(len(fine.shape()) == 1):
        x_fine = fine.shape()


    # get slice indices
    xSlice = range(0,x_fine,2)
    ySlice = range(0,y_fine,2)

    # copy over
    coarse[:,:] = fine[0:x_fine:2, 0:y_fine:2]


def sum4way(fine, coarse):
    """Contracts the Field by summing 4 values into 1. Only setup for 2D Fields (?) so far.
    
    Args:
    
    fine:
        The Field object on the finer grid
    coarse:
        The Field object on the coarser grid
    """
    if type(fine) is not Field or type(coarse) is not Field:
        raise TypeError('Fine or coarse field is not a field')
    
    # Check that fine grid is divisible by 2 in both dims
    nxf, nyf = fine.size()
    # if (x_fine % 2) == 1 or (y_fine % 2) == 1:
    #     raise ValueError('Fine field dimensions do not allow for 4 way sum')

    # Check that dimensions of coarse grid are half of fine grid
    nxc, nyc = coarse.size()
    # if (x_fine / 2) != x_coarse or (y_fine / 2) != y_coarse:
    #     raise ValueError('Coarse grid size different from expected reduction from fine grid')
    
    # get slice indices
    x = np.arange(0,nxf+1,2)
    y = np.arange(0,nyf+1,2)

    for i in range(nxc):
        for j in range(nyc):
            [i1, i2] = x[i:i+2]
            [j1, j2] = y[j:j+2]
            coarse[i,j] = sum(fine[i1:i2, j1:j2])



def conservative4way(fine, coarse, weights=None):
    """Contracts the Field by averaging 4 values into 1 with weighting terms. Only set up for 2D Fields (?) so far.
    
    Args:
    
    fine:
        The Field object on the finer grid
    coarse:
        The Field object on the coarser grid
    """
    if type(fine) is not Field or type(coarse) is not Field:
        raise TypeError('Fine or coarse field is not a field')
    
    if weights is None:
        sum4way(fine, coarse)
        coarse[:] = coarse/4
        return

    # Combines fine grid by summing over 4 blocks 
    # if not isinstance(fine, Field) or not isinstance(coarse, Field):
    #     raise TypeError('Fine or coarse field is not a field')
    
    # Check that fine grid is divisible by 2 in both dims
    nxf, nyf = fine.size()
    # if (x_fine % 2) == 1 or (y_fine % 2) == 1:
    #     raise ValueError('Fine field dimensions do not allow for 4 way sum')

    # Check that dimensions of coarse grid are half of fine grid
    nxc, nyc = coarse.size()
    # if (x_fine / 2) != x_coarse or (y_fine / 2) != y_coarse:
    #     raise ValueError('Coarse grid size different from expected reduction from fine grid')
    
    # get slice indices
    x = np.arange(0,nxf+1,2)
    y = np.arange(0,nyf+1,2)

    for i in range(nxc):
        for j in range(nyc):
            num = sum(fine[x[i]:x[i+1], y[j]:y[j+1]] * weights[x[i]:x[i+1], y[j]:y[j+1]])
            den = sum(weights[x[i]:x[i+1], y[j]:y[j+1]])
            coarse[i,j] = num / den