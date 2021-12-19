from bin.Field import Field
import numpy as np
from scipy.interpolate import griddata

# def billenear4way(coarse, fine):
#     # Combines fine grid by summing over 4 blocks 
    
#     x_coarse = np.shape(coarse)[0]
#     y_coarse = np.shape(coarse)[1]
    
#     x_fine = np.shape(fine)[0]
#     y_fine = np.shape(fine)[1]
#     grid_x, grid_y = np.mgrid[0:1:1j*x_fine, 0:1:1j*y_fine]
    
#     x_c_line = np.linspace(0, 1, num=x_coarse)
#     y_c_line = np.linspace(0, 1, num=y_coarse)
#     points = np.array([[x_c_line[0], y_c_line[0]]])
#     values = np.array([])
#     for x in range(x_coarse):
#         for y in range(y_coarse):
#             # Get point coords and value at that point
#             if not(x == 0 and y == 0):
#                 points = np.append(points, [[x_c_line[x], y_c_line[y]]], axis=0)
#             values = np.append(values, coarse[x,y])
    
#     temp = griddata(points, values, (grid_x, grid_y), method='linear')
#     for i in range(x_fine):
#         for j in range(y_fine):
#             fine[i,j] = temp[i,j]

def bilinear4way(coarse, fine):
    [nx, ny, dim] = fine.shape()


    for n in range(dim):
        for i in range(2,nx,2):
            for j in range(1,ny,2):
                    fine[i  ,j,n] = .25*fine[i-1,j,n]  +.75*fine[i+1,j,n]
                    fine[i-1,j,n] = .75*fine[i-1,j,n]  +.25*fine[i+1,j,n]
        
        for i in range(1,nx,2):
            for j in range(2,ny,2):
                    fine[i  ,j,n] = .25*fine[i,j-1,n]  +.75*fine[i,j+1,n]
                    fine[i-1,j,n] = .75*fine[i,j-1,n]  +.25*fine[i,j+1,n]
    
    # Gigi did this, not sure if it's  needed, ignore for now
    # if (mode.gt.0) then
    #      do i=1,ie
    #         dw(i,jb,n)  = 0.
    #      end do
    #   end if