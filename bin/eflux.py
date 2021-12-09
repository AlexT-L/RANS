from Field import Field
from Grid import Grid
from Workspace import Workspace
from NavierStokes import NavierStokes
import numpy as np

def eflux(ws):
    # take a workspace ws and calculate convective fluxes

    G = ws.grd # grab grid
    pU = ws.flds['pU'] # x- momentum
    pV = ws.flds['pV'] # y- momentum
    pE = ws.flds['pE'] # total energy
    rho = ws.flds['rho'] # density
    porJ = ws.flds['porJ'] # porosity
    n = 4 # number of quantities being convected
    P = ws.flds['P']   # pressure
    
    # storing residuals
    drho = np.zeros(np.shape(rho[:][:]))
    dpU = np.zeros(np.shape(pU[:][:]))
    dpV = np.zeros(np.shape(pV[:][:]))
    dpE = np.zeros(np.shape(pE[:][:]))

    # i direction
    fs = np.zeros([G.il, G.jl,n])
    for j in range(1,G.jl-1):
        for i in range(0,G.il-1):

            # normal vector 
            dxy = G.X[i,j,0] - G.X[i,j-1,0]
            dyy = G.X[i,j,1] - G.X[i,j-1,1]
            
            # pressure averaging
            Pa = P[i+1][j] + P[i][j]

            # flux operator
            qsp       = (dyy*pU[i+1,j]  -dxy*pV[i+1,j])/rho[i+1,j]
            qsm       = (dyy*pU[i,j]  - dxy*pV[i,j])/rho[i,j]

            # add up on faces
            fs[i,j,0] = qsp*rho[i+1,j]   + qsm*rho[i,j] # density
            fs[i,j,1] = qsp*pU[i+1,j]  + qsm*pU[i,j]  + dyy*Pa # x - momentum
            fs[i,j,2] = qsp*pV[i+1,j]  + qsm*pV[i,j]  - dxy*Pa # y - momentum
            fs[i,j,3] = qsp*(pE[i+1,j] + P[+1,j]) + qsm*(pE[i,j] + P[i,j]) # energy

    # now add everything up
    for j in range(1,G.jl-1):
        for i in range(1,G.il-1):
            drho[i,j] = fs[i,j,0] -fs[i-1,j,0]
            dpU[i,j] = fs[i,j,1]  -fs[i-1,j,1]
            dpV[i,j] = fs[i,j,2]  -fs[i-1,j,2]
            dpE[i,j] = fs[i,j,3]  -fs[i-1,j,3]


    # j direction
    for j in range(0,G.jl-1):
      for i in range(1,G.il-1):

          # normal vector
         dxx        = G.X[i,j,0]  -G.X[i-1,j,0]
         dyx        = G.X[i,j,1]  -G.X[i-1,j,1]
         # pressure average
         Pa        = P[i,j+1]  +P[i,j]
         # convective operator
         qsp       = porJ[i,j]*(dxx*pV[i,j+1]  - dyx*pU[i,j+1])/rho[i,j+1]    
         qsm       = porJ[i,j]*(dxx*pV[i,j]  - dyx*pU[i,j])/rho[i,j]
        # add up on faces
         fs[i,j,0] = qsp*rho[i,j+1]  +qsm*rho[i,j]
         fs[i,j,1] = qsp*pU[i,j+1]  +qsm*pU[i,j]  - dyx*Pa
         fs[i,j,2] = qsp*pV[i,j+1]  +qsm*pV[i,j]  + dxx*Pa
         fs[i,j,3] = qsp*(pE[i,j+1]  +P[i,j+1]) +qsm*(pE[i,j]  + P[i,j])
                   
    # now add everything up
    for j in range(1,G.jl-1):
        for i in range(1,G.il-1):
            drho[i,j] = drho[i,j] + fs[i,j,0] - fs[i,j-1,0]
            dpU[i,j]  = dpU[i,j]  + fs[i,j,1] - fs[i,j-1,1]
            dpV[i,j]  = dpV[i,j]  + fs[i,j,2] - fs[i,j-1,2]
            dpE[i,j]  = dpE[i,j]  + fs[i,j,3] - fs[i,j-1,3]

    ws.flds['drho'] = drho
    ws.flds['dpU'] = pU
    ws.flds['dpV'] = pV
    ws.flds['dpE'] = dpE
    

#nx = 10
#ny = 5
#x_bound = np.array([-5,5])
#y_bound = np.array([-10,10])
#I = Input(['pU', 'pV', 'P', 'rho','pE','porJ']) 

#G = Grid(x_bound,y_bound,nx,ny)
#G.X[:][:][:] = G.X[:][:][:] + np.random.standard_normal(np.shape(G.X[:][:][:]))
#M = NavierStokes(I)
#W = Workspace(M,G)
#print(W.flds['pU'][1][1])
#eflux(W)
#print(W.flds['drho'])

#wf = np.zeros([G.il + 3, G.jl + 3, 4])
#dwf = np.zeros([G.il + 3, G.jl + 3, 4])
#print(np.shape(G.X[:][:][:]))
#Pf = np.zeros([G.il + 3, G.jl + 3])
#efluxfort.eflux(wf,dwf,Pf,G.X[:][:][:])
#print(W.flds['dpU'][1][1])