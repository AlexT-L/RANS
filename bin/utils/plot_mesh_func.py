"""This module plots the c-mesh in physical space

    Libraries/Modules:
        nump\n
        matplotlib.pyplot\n
"""
import numpy as np 
import matplotlib.pyplot as plt

def plot_mesh(self):
    """ Plots the c-mesh in physical space after conformal mapping
        
        Uses the vertices from x to plot the mesh.
    """
    #x-y vertices
    ver = self.fields['x']
    x=ver[:,:,0]
    y=ver[:,:,1]
    #transpose of x-y vertices
    x_t=x.T().get_vals()
    y_t=y.T().get_vals()
    #plot c-mesh
    # plt.plot(x.get_vals(),y.get_vals())
    # plt.plot(x_t,y_t,linewidth="0.5")
    # plt.xlabel("x [chord lengths]",fontsize=15)
    # plt.ylabel("y [chord lengths]",fontsize=15)
    # plt.axis([-0.75,1.50,-0.8,0.8])
    # plt.show()

    return