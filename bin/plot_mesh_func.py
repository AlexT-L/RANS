import numpy as np 
import matplotlib.pyplot as plt

def plot_mesh(self,ver):
    #x-y vertices
    x=ver[:,:,0]
    y=ver[:,:,1]
    #transpose of x-y vertices
    x_t=x.T 
    y_t=y.T 
    #plot c-mesh
    plt.plot(x,y)
    plt.plot(x_t,y_t,linewidth="0.5")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis([-0.75,1.50,-0.8,0.8])
    plt.show()

    return