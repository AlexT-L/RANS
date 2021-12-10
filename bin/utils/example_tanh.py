import tanhds_fort
import numpy as np


ns = 10
sp1 = 0.2
sp2 = 0.5
s = np.ones(ns,order = 'f') # this order = f is required to wrap fortran files 
status = 1
print(tanhds_fort.__doc__) # prints documentation for function
tanhds_fort.tanhds(sp1,sp2,s,status) 
# we dont need to give ns = len(s) because the wrapper automatically 
# realized it was optional