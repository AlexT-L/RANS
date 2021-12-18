import numpy as np
from numpy.core.numeric import Infinity
from Field import Field, max
from Input import Input
from flo103_PostProcessor import flo103_PostProcessor
from flo103_ConvergenceChecker import flo103_ConvergenceChecker
from ImplicitEuler import ImplicitEuler
from NS_Airfoil import NS_Airfoil
from AirfoilMap import AirfoilMap
from CellCenterWS import CellCenterWS
from NavierStokes import NavierStokes
from MultiGrid import MultiGrid
import Contractinator as con

if __name__ == '__main__':
    # # Testing fields
    # field = Field((6,2,2))
    # field2 = Field((6,2,2))
    # field[0:2,0:2,0] = 5
    # field2[0:6,0:2,0] = 1
    # weights = Field((12,8,4))
    # # print(field.vals)

    # newPointer = field
    # print(field)

    # newPointer += field2


    a = Field((5,5,5))

    print(isinstance(a, Field))
    print(type(a) is Field)


    
    print(a[:])
    