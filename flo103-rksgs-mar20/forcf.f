      subroutine forcf
c
c     ******************************************************************
c     *                                                                *
c     *   calculates the force coefficients                            *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use flo_var
      use mesh_var
      use out_var
c
c     ******************************************************************
c
      use flo_param
      use mg_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: i,j,l
c
c     ******************************************************************
c
      real     :: xa,ya,dx,dy,dcl,dcd
      real     :: gmm,gm1,gpm,scf,clvis,cdvis
      real     :: rho1,ru1,rv1,p1,u1,v1,t1
      real     :: rho2,ru2,rv2,p2,u2,v2,t2
      real     :: dxi,dxj,dyi,dyj,dui,duj,dvi,dvj,dsj
      real     :: dux,duy,dvx,dvy
      real     :: rlva,reva,rmu,rlam,rk,term
      real     :: sigx,sigy,tauxy
c
c     ******************************************************************
c
      real, dimension(ie,je)            :: u,v,t
      real, dimension(il)               :: gs2,gs3,gs4
c
c     ******************************************************************
c
      cl        = 0.
      cd        = 0.
      cm        = 0.
      gmm       = gamma*rm**2

      do i=2,il
         cp(i)     = ((p(i,2)  +p(i,1))/p0  -2.)/gmm
      end do

      do i=itl+1,itu
         dx        = scal*(x(i,1,1)  -x(i-1,1,1))/chord
         dy        = scal*(x(i,1,2)  -x(i-1,1,2))/chord
c
c     this seems to be a bug.
c     when generating mesh xm is based on input xp.
c     when reading mesh xm is based on x & y.
c     i modified forcf.f gmesh.f st. xm is based on x & y.
c
c        xa        = (.5*scal*(x(i,1,1)  +x(i-1,1,1))  -xm)/chord
c        ya        = (.5*scal*(x(i,1,2)  +x(i-1,1,2))  -ym)/chord
c
         xa        = scal*(.5*(x(i,1,1)  +x(i-1,1,1))  -xm)/chord
         ya        = scal*(.5*(x(i,1,2)  +x(i-1,1,2))  -ym)/chord
         dcl       = -cp(i)*dx
         dcd       = cp(i)*dy
         cl        = cl  +dcl
         cd        = cd  +dcd
         cm        = cm  +dcd*ya  -dcl*xa
      end do

      if (kvis.eq.0) then
         clvis     = 0.
         cdvis     = 0.
      else
c
c     evaluation of viscous corrections to cl and cd.
c
      clvis   = 0.
      cdvis   = 0.

      gm1       = gamma - 1.
      gpm       = 2./(gamma*rm**2)
      scf       = sqrt(gamma)*rm/(scal*re/chord)

      do j=1,2
      do i=1,ie
         u(i,j)    = w(i,j,2)/w(i,j,1)
         v(i,j)    = w(i,j,3)/w(i,j,1)
c        t(i,j)     = p(i,j)/(gm1*w(i,j,1))
      end do
      end do

      if (mode.ne.0) then
         do i=itl+1,itu
            rev(i,1)  = -rev(i,2)
            u(i,1)    = -u(i,2)
            v(i,1)    = -v(i,2)
         end do
      end if
c
c     evaluation of viscous flux terms in i and j directions
c
      do 10 i=itl+1,itu

      rho1    = .25*(w(i,1,1)  +w(i+1,1,1)  +w(i,2,1)  +w(i+1,2,1))
      ru1     = .25*(w(i,1,2)  +w(i+1,1,2)  +w(i,2,2)  +w(i+1,2,2))
      rv1     = .25*(w(i,1,3)  +w(i+1,1,3)  +w(i,2,3)  +w(i+1,2,3))
      p1      = .25*(p(i,1)    +p(i+1,1)    +p(i,2)    +p(i+1,2))
      u1      = ru1/rho1
      v1      = rv1/rho1
c     t1      = p1/(gm1*rho1)

      rho2    = .25*(w(i-1,1,1)  +w(i,1,1)  +w(i-1,2,1)  +w(i,2,1))
      ru2     = .25*(w(i-1,1,2)  +w(i,1,2)  +w(i-1,2,2)  +w(i,2,2))
      rv2     = .25*(w(i-1,1,3)  +w(i,1,3)  +w(i-1,2,3)  +w(i,2,3))
      p2      = .25*(p(i-1,1)    +p(i,1)    +p(i-1,2)    +p(i,2))
      u2      = ru2/rho2
      v2      = rv2/rho2
c     t2      = p2/(gm1*rho2)

      dxi       = x(i,1,1)  -x(i-1,1,1)
      dyi       = x(i,1,2)  -x(i-1,1,2)
      dxj       = xc(i,2,1)  -xc(i,1,1)
      dyj       = xc(i,2,2)  -xc(i,1,2)
      dsj       = 1./(dxi*dyj  -dyi*dxj)

      dui       = u1  -u2
      dvi       = v1  -v2
c     dti       = t1  -t2
      duj       = u(i,2)  -u(i,1)
      dvj       = v(i,2)  -v(i,1)
c     dtj       = t(i,2)  -t(i,1)

      dux       = (dui*dyj  -duj*dyi)*dsj
      dvx       = (dvi*dyj  -dvj*dyi)*dsj
c     dtx       = (dti*dyj  -dtj*dyi)*dsj

      duy       = (duj*dxi  -dui*dxj)*dsj
      dvy       = (dvj*dxi  -dvi*dxj)*dsj
c     dty       = (dtj*dxi  -dti*dxj)*dsj

      rlva      = .5*(rlv(i,2)  +rlv(i,1))
      reva      = .5*(rev(i,2)  +rev(i,1))
      rmu       = rlva  +reva
      rk        = gamma*(rlva/prn  +reva/prt)
      rlam      = -2.*rmu/3.

      term      = rlam*(dux  +dvy)
      sigx      = term  +2.*rmu*dux
      sigy      = term  +2.*rmu*dvy
      tauxy     = rmu*(duy  +dvx)
c     qx        = -rk*dtx
c     qy        = -rk*dty

c
c     residuals in "nsflux.f" have extra 2, compensated later in "euler.f"
c     so, for viscous stresses in lines below:                   -2. ==> -1.
c     also, gs2 and gs3 -are forces on the body (not fluid), so: -1. ==>  1.
c
      gs2(i)    = (-sigx*dyi  +tauxy*dxi)*scf*gpm
      gs3(i)    = ( sigy*dxi  -tauxy*dyi)*scf*gpm

c     gs2(i)    = (-sigx*dyi  +tauxy*dxi)*scf
c     gs3(i)    = ( sigy*dxi  -tauxy*dyi)*scf

      cdvis     = cdvis  +gs2(i)
      clvis     = clvis  +gs3(i)

      cf(i)     = tauxy*scf*gpm

   10 continue

      end if

      dcl       = cl*ca  -cd*sa
      cd        = cl*sa  +cd*ca
      cl        = dcl
      cdv       = clvis*sa  +cdvis*ca
      clv       = clvis*ca  -cdvis*sa

c     cd        = cd +clvis*sa  +cdvis*ca
c     cl        = cl +clvis*ca  -cdvis*sa

      return

      end
