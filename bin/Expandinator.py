import Field
import numpy as np
from scipy.interpolate import griddata

def billenear4way(coarse, fine):
    # Combines fine grid by summing over 4 blocks 
    if type(fine) != Field or type(coarse) != Field:
        raise TypeError('Fine or coarse field is not a field')
    
    x_coarse = np.shape(coarse)[0]
    x_step = 1.0/x_coarse
    y_coarse = np.shape(coarse)[1]
    y_step = 1.0/y_coarse
    grid_x, grid_y = np.mgrid[0:x_step:1j, 1:y_step:1j]
    
    
    
    newField = np.zeros(np.shape(fine))
    for x in range(x_coarse):
        for y in range(y_coarse):
            # Set corners to coarse grid
            pass
    