import numpy as np
from numpy.core.numeric import Infinity
from Field import Field, max
from Input import Input
from flo103_PostProcessor import flo103_PostProcessor
from flo103_ConvergenceChecker import flo103_ConvergenceChecker
from ImplicitEuler import ImplicitEuler
from NS_AirfoilBC import NS_AirfoilBC
from AirfoilMap import AirfoilMap
from CellCenterWS import CellCenterWS
from NavierStokes import NavierStokes
from MultiGrid import MultiGrid
import Contractinator as con

if __name__ == '__main__':
    # Testing fields
    field = Field([6,2],2)
    field2 = Field([6,2],2)
    field[0:2,0:2,0] = 5
    field2[0:6,0:2,0] = 1
    weights = Field([12,8],4)
    # print(field.vals)

    newPointer = field
    print(field)

    newPointer += field2
    # newPointer.store_sum(newPointer, field2)
    print(newPointer)
    print(field)
    
    print("\n")
    field2[:,:,:] = 2
    print(field2)

    print("\n\n\n")

    # storage = field[0:4, 0:2, 0]
    # storage.copy_from(field2[0:4, 0:2, 0])
    # # field[0:4, 0:2, 0].copy_from(field2[0:4, 0:2, 0])
    # # field.copy_from(field2)

    # print("\n\n")
    # print(storage)
    # print(field)


    print("\n")
    num1 = np.zeros((4,4,6))
    num2 = np.ones((4,4,6))
    # newnum = num1[0:2,0:4,0:2]
    # newnum[:,:,:] = 5
    # print(num1)
    # print(newnum)

    test1 = Field((2,4),6,num1)
    test2 = Field((2,4),6,num1)
    num1[:,:,:] = 3

    temp = test1[0:2,1:3,0:6]
    temp[:,:,:] = 4

    print(test2)

    mask = test2 == 4

    print(test1*mask + test2*(1-mask))

    test1[0,:,0] = [0, 1, 2, 3]
    print(test1[0, 4::-1, 0])

    # field2[0:4,0,0] = 4
    # print(field2.vals)

    # print(field.shape())

    # print("sum")
    # print(sum(sum(field[0:2,0:2,0])))

    # for z in range(4):
    #     field[0,0:8,z] = 3
    #     weights[0:12,0,z] = 2
    
    # print(field.vals)

    # print(field.shape())

    # con.conservative4way(field, field2, weights)

    # print(field2.vals)

    # print(field.vals)