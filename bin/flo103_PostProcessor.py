"""This would provide postprocessing of results. 

    Libraries/Modules:
        Would use: numpy\n
        Would use: pandas\n
        """
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from bin.NavierStokes import NavierStokes


class flo103_PostProcessor:
    """Not implemented as this time. 
        What it would do: 
        Plot pressure profile and entropy potential lines.
        Plot outputs of wall stress and aerodynamics of airfoil.

    Attributes:
        None currently used.

    Notes: 
    Currently just passed when called. 
    """

    def __init__(self, input):
        """
        Is not used in current implementation. 
        """
        pass
    
    def print_state(self, workspace, state, pressure=None):
        xc = workspace.get_field('xc')
        x = xc[:,:,0]
        y = xc[:,:,1]

        names = ["density", "x-mom", "y-mom", "energy"]
        plt.figure()
        for i in range(4):
            plt.contourf(x,y,state[:,:,i])
            plt.title(names[i])
            plt.colorbar()
            plt.axis([-0.75,1.50,-0.8,0.8])
            plt.show()
        
        # Pressure contours
        if pressure is not None:
            plt.figure()
            plt.contourf(x,y,pressure)
            plt.title("Pressure contours")
            plt.colorbar()
            plt.axis([-0.75,1.50,-0.8,0.8])
            plt.show()
        
    
    def print_grid(self, model, workspace):
        xcpad = workspace.get_field('xc', model.className)
        x = xcpad[:,:,0]
        y = xcpad[:,:,1]
        plt.plot(x[1:-1,1:-1], y[1:-1,1:-1])
        plt.title("xc")
        plt.axis([-0.75,1.50,-0.8,0.8])
        plt.show()
        
        volpad = workspace.get_field('vol', model.className)
        volpad = volpad[2:-2,2:-2]
        vol = workspace.get_field('vol')
        diff = volpad-vol
        x = x[2:-2,2:-2]
        y = y[2:-2,2:-2]
        xc = workspace.get_field('xc')
        x = xc[:,:,0]
        y = xc[:,:,1]
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        # ax.contour3D(x, y, vol, 50, cmap='binary')
        ax.plot_surface(x, y, diff, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('vol')
        ax.view_init(60, 35)
        plt.show()
        print(np.min(vol))
        assert np.max(np.abs(vol-volpad))<1e-15
        assert np.max(np.abs(xcpad[2:-2,2:-2,:]-xc))<1e-15
