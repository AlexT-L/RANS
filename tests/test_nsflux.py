'''
import sys
sys.path.append("../")

import nsflux_fort
import numpy as np
#from eflux_arr import eflux
from Field import Field

# grab grid related parameter
#G = ws.grid
nx = 4
ny = 10
il = nx+1
jl = ny+1
ie = il+1
je = jl+1
itl = 1
itu = 3
ib = il + 2
jb = jl + 2

# flow related vars
w = Field.create([ib,jb],4) # state
w = np.array(w + 15*np.random.standard_normal([ib,jb,4]),order = 'f')
P = Field.create([ib,jb]) # pressure
lv = Field.create([ib,jb]) # laminar viscocity
ev = Field.create([ib,jb]) # eddy viscocity
vw = Field.create([ib,jb],4) # residuals

# mesh related vars
porI = Field.create([ib,jb],2) # mesh vertices
porI = np.array(porI + 15*np.random.standard_normal([ib,jb,2]),order = 'f')
porJ = Field.create([ib,jb],2) # mesh centers
porJ = np.array(porJ + 15*np.random.standard_normal([ib,jb,2]),order = 'f')
xc = Field.create([ib,jb],2) # mesh vertices
xc = np.array(porI + 15*np.random.standard_normal([ib,jb,2]),order = 'f')
x = Field.create([ib,jb],2) # mesh centers
x = np.array(porJ + 15*np.random.standard_normal([ib,jb,2]),order = 'f')

# solver related vars
fw = Field.create([ib,jb],4)
radI = Field.create([ib,jb],2) # stability I
radJ = Field.create([ib,jb],2) # stability J

gamma = 1.4
rm = 1.2
scal = 1.8
re = 50000
chord = 2.6
prn = 1000
prt = 10000
mode = 1
rfil = 0.8
vis0 = 0.5
rho0 = 1
p0 = 1;h0 = 1;c0 = 1;u0 = 1;v0 = 1;ca= 1;sa = 1; xm = 1; ym = 1; kvis = 1; bc = 1

print(vw[:][:][0])
print(nsflux_fort.__doc__)
# residuals returned in Field dw
nsflux_fort.nsflux(il, jl, ie, je, \
      w, P, lv, ev,  \
      x, xc, \
      vw,
      gamma,rm,scal,re,chord,prn,prt, \
      rfil)

print(vw[:][:][0])
'''