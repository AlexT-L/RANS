from bin.Field import Field
import numpy as np
#from scipy.interpolate import griddata

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
    # nx, ny, dim = fine.shape() # Doesn't work, fine.shape() only returns two items.
    # #nx, ny = fine.shape() # This works
    nxf, nyf = fine.size()
    nxc, nyc = coarse.size()
    
    for i in range(2,nxf,2):
        for j in range(1,nyf,2):
                fine[i  ,j] = .25*coarse[i-1,j]  +.75*coarse[i+1,j]
                fine[i-1,j] = .75*coarse[i-1,j]  +.25*coarse[i+1,j]
    
    for i in range(1,nxf,2):
        for j in range(2,nyf,2):
                fine[i  ,j] = .25*coarse[i,j-1,n]  +.75*coarse[i,j+1]
                fine[i-1,j] = .75*coarse[i,j-1,n]  +.25*coarse[i,j+1]
