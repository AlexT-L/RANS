import sys
sys.path.append("..")

import metricinator
import numpy as np
import Field as Field



print(metricinator.__doc__)
# note you don't need to pass ib and jb 
# the function signature looks like this
nx = 512; ny = 64
il = nx  +1
jl = ny  +1
ie = nx  +2
je = ny  +2
ib = nx + 3
jb = ny + 3
itl = 52 # no idea what a realistic value is
rand = np.array(np.random.standard_normal([ib,jb,2]), order = 'f')
x = np.ones([ib,jb,2],order = 'f')+ rand
xc = np.ones([ib,jb,2],order = 'f')
vol = np.ones([ib,jb],order = 'f')
metricinator.metric(il,jl,ie,je,itl,x,xc,vol)
print(vol)