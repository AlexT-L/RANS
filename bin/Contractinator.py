import numpy as np
import Field

def sum4way(fine, coarse, weights=np.ones(4)):
    # Combines fine grid by summing over 4 blocks 
    if type(fine) != Field or type(coarse) != Field:
        raise TypeError('Fine or coarse field is not a field')
    
    if len(weights) != 4:
        raise ValueError('Weights not of length 4, is: '+len(weights))
    
    # Check that fine grid is divisible by 2 in both dims
    x_dim = np.shape(fine)[0]
    y_dim = np.shape(fine)[1]
    if (x_dim % 2) == 1 or (y_dim % 2) == 1:
        raise ValueError('Fine field dimensions do not allow for 4 way sum')

    # Check that dimensions of coarse grid are half of fine grid
    x_coarse = np.shape(fine)[0]
    y_coarse = np.shape(fine)[1]
    if (x_dim / 2) != x_coarse or (y_dim / 2) != y_coarse:
        raise ValueError('Coarse grid size different from expected reduction from fine grid')
    
    newField = np.zeros(np.shape(coarse), order='F')
    for x in range(x_coarse):
        for y in range(y_coarse):
            flat = fine[2*x:2*x+1,2*y:2*y+1].flatten() # Pick 4 cells around new center
            newField[x, y] = np.sum(flat * weights) / sum(weights) # Multiply by weights and sum
            
    coarse.set_val(newField)


def conservative4way(fine, coarse, weights):
    pass