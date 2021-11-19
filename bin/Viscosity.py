import numpy as np
from Grid import Grid

class Viscosity():        
#  from subroutine viscf.f
    def compute_viscosity(grid: Grid, params):
          

#     ******************************************************************
#     *                                                                *
#     *   computes viscosity coefficients                              *
#     *                                                                *
#     ******************************************************************


# subroutines / modules "used"
    # dims
    # flo_var, mesh_var
    # flo_param, solv_param, mg_param

real, dimension(ie,je)            :: u,v,astr,rev0

# useful constants
pi        = np.pi
# ckr       = (.062/(2.*pi)**4)
# ckr       = .01915
ckr       = .0256
# cwk       = .225
cwk       = 0.
#  scf       = re/(sqrt(gamma)*rm)
scf       = (scal*re/chord)/(sqrt(gamma)*rm)

# compute the molecular viscosity

for j in range (1,je):
      for i in range (1,ie):
         tt       = p(i,j)/w(i,j,1)*t0
         rlv(i,j) = 1.461e-06*tt*sqrt(tt)/((tt+110.3)*rmu0)


# for laminar flows we are done.
# for turbulent flows we are also done on the coarser grids.

if (kvis.le.1.or.mode.ne.0):
    return
# if we are using the baldwin and lomax model call turbbl and return

      aturb     = 1.
      if (ncyc.gt.25):
        aturb = .5
      if (kturb.eq.1): 
          for j in range (1,je):
              for i in range (1,ie):
                  rev0(i,j) = rev(i,j)


        # call turbbl
        #   to add: 
        # call turb2
        for j in range (1,je):
            for i in range (1,ie):
                    rev(i,j) = aturb*rev(i,j)  +(1.  -aturb)*rev0(i,j)

                return

#  else start the rng algebraic model

for j in range (1,je):
      for i in range (1,ie):
         u(i,j)   = w(i,j,2)/w(i,j,1)
         v(i,j)   = w(i,j,3)/w(i,j,1)


      for i in range (itl+1,itu):
         u(i,1)   = -u(i,2)
         v(i,1)   = -v(i,2)
for j in range (1,jl):
      for i in range (1,il):

         dx13      = xc(i,j,1)   - xc(i+1,j+1,1)
         dy13      = xc(i,j,2)   - xc(i+1,j+1,2)
         dx24      = xc(i+1,j,1) - xc(i,j+1,1)
         dy24      = xc(i+1,j,2) - xc(i,j+1,2)
         du13      = u(i,j) - u(i+1,j+1)
         dv13      = v(i,j) - v(i+1,j+1)
         du24      = u(i+1,j) - u(i,j+1)
         dv24      = v(i+1,j) - v(i,j+1)
         ua        = .25*(u(i,j) + u(i+1,j+1) + u(i+1,j) + u(i,j+1))
         va        = .25*(v(i,j) + v(i+1,j+1) + v(i+1,j) + v(i,j+1))
         dsij      = 1./(dx13*dy24 - dx24*dy13)
         dvdx      =  dsij * (dv13*dy24 - dv24*dy13)
         dudy      = -dsij * (du13*dx24 - du24*dx13)
         dudx      =  dsij * (du13*dy24 - du24*dy13)
         dvdy      = -dsij * (dv13*dx24 - dv24*dx13)
         astr(i,j) = (dudy+dvdx)**2. +2.*(dudx**2  +dvdy**2  -((dudx+dvdy)**2)/3.)


    #   to add: 
    # call delt


    #   do 30 j=2,jl
    #   do 20 i=2,il
    # what do the 30 and 20 do?
    # also they did not have a corresponding end do?
for j in range (2,jl):
      for i in range (2,il):

      xbi       = .5*(x(i-1,1,1)  +x(i,1,1))
      ybi       = .5*(x(i-1,1,2)  +x(i,1,2))
      astra     = .25*(astr(i-1,j-1)  +astr(i-1,j)
     .                +astr(i,j-1)    +astr(i,j))
      if (i.ge.itl.and.i.le.itu+1):
          a3        = 1./(.225*abs(ynot(i)))
          ysci      = sqrt((xc(i,j,1)  -xbi)**2  +(xc(i,j,2)  -ybi)**2)
          ysc       = w(i,2,1)/(ysci*w(i,j,1))
          csc       = 1./(ysc+a3)**2
      elif
         csc       = (cwk*ynot(i))**2


#     set some parameters

      rnul      = rlv(i,j)/w(i,j,1)
      rnut0     = rev(i,j)/w(i,j,1)
      rnul3     = rnul**3
      a11       = ckr*(csc*csc*scf*scf)/rnul**2
      a2        = 75.
      a1        = a11*(astra)
      rnut0     = rnul+rnut0

#     solve for the eddy viscosity

      if (dim(rnut0*a1,a2).eq.0.):
          rev(i,j)  = 0.
          #  go to 20? 
      elif:
          rnut      = sqrt(a1)
      

k      = 0
fac    = a2 - 1.


# 11?
den    = 1./(4.*rnut*rnut*rnut + fac)
rnut1  = rnut - (rnut**4+rnut*fac  -rnut0*rnut0*a1)*den

      if (abs((rnut1  -rnut)).le.1.e-3):
          rev(i,j) = w(i,j,1)*dim(rnut1,rnul)
         # go to 20
      elif
         k      = k  +1
         if (k.gt.200):
            write (6,*) ' iteration not converged ',i,j
            write (6,*) ' rnut = ',rnut,' rnut1 =',rnut1
            rev(i,j)  = w(i,j,1)*dim(rnut1,rnul)
            # go to 20

         rnut   = rnut1
         # go to 11

# those would 'go' here, and then just continue?
#    20 continue
#    30 continue

#     adjust the near wake

ii        = ie

for i in range (2,itl+1):
      ii        = ii  -1
      for j in range (2,jl):
         pex       = -(xc(i,2,1)  -xc(itl+1,2,1))/(20.*dsti(itl+1))
         rev(i,j)  = rev(i,j)  +(rev(itl+1,j)  -rev(i,j))*exp(pex)
         pex       = -(xc(ii,2,1)  -xc(itu,2,1))/(20.*dsti(itu))
         rev(ii,j) = rev(ii,j)  +(rev(itu,j)  -rev(ii,j))*exp(pex)

for i in range (2,il):
         ii        = ib  -i
         rev(i,je) = rev(i,jl)
         rev(i,1)  = rev(ii,2)

for i in range (itl+1,itu):
         if (xc(i,2,1).le.xtran):
            for j in range (1,jl):
               rev(i,j)  = 0
         rev(i,1)  = -rev(i,2)

for j in range (1,je):
         rev(1,j)  = rev(2,j)
         rev(ie,j) = rev(il,j)

      return
