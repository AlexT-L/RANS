from model_funcs import eflux_fort
import numpy as np
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
w = Field([ib,jb],4) # state
w.vals = np.array(w.vals + 15*np.random.standard_normal([ib,jb,4]),order = 'f')
P = Field([ib,jb]) # pressure
rlv = Field([ib,jb]) # laminar viscocity
rev = Field([ib,jb]) # eddy viscocity
dw = Field([ib,jb],4) # residuals

# mesh related vars
porI = Field([ib,jb],2) # mesh vertices
porI.vals = np.array(porI.vals + 15*np.random.standard_normal([ib,jb,2]),order = 'f')
porJ = Field([ib,jb],2) # mesh centers
porJ.vals = np.array(porJ.vals + 15*np.random.standard_normal([ib,jb,2]),order = 'f')
xc = Field([ib,jb],2) # mesh vertices
xc.vals = np.array(porI.vals + 15*np.random.standard_normal([ib,jb,2]),order = 'f')
x = Field([ib,jb],2) # mesh centers
x.vals = np.array(porJ.vals + 15*np.random.standard_normal([ib,jb,2]),order = 'f')

# solver related vars
fw = Field([ib,jb],4)
radI = Field([ib,jb],2) # stability I
radJ = Field([ib,jb],2) # stability J

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

print(dw.vals[:][:][0])
print(eflux_fort.__doc__)
# residuals returned in Field dw
eflux_fort.eflux(w.vals,dw.vals,P.vals,x.vals,porI.vals,il,jl)

print(dw.vals[:][:][0])
