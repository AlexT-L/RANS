      subroutine vflux
c
c     *****************************************************************
c     *                                                               *
c     *   calculates the viscous flux                                 *
c     *                                                               *
c     *****************************************************************
c
      use dims
c
c     ******************************************************************
c
      use flo_var
      use mesh_var
      use solv_var
c
c     ******************************************************************
c
      use flo_param
      use solv_param
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
      real     :: gm1,scf,sfil
      real     :: ua,va,rmu,rlam,rk,rlva,reva,dq
      real     :: dx,dy,dxi,dxj,dyi,dyj,dui,duj,dsj
      real     :: sigx,sigy,tauxy
c
c     ******************************************************************
c
      real, dimension(il,jl,3,2)        :: q
      real, dimension(ie,je,3)          :: u
c     real, dimension(il,jl,3)          :: fs
c
c     *****************************************************************
c
c     if (ny.lt.3) return
c
      sfil      = 1.  -rfil
      gm1       = gamma - 1.
c     scf       = rfil*sqrt(gamma)*rm/re
      scf       = rfil*sqrt(gamma)*rm/(scal*re/chord)
c
c     calcluate the velocities and temperature
c
      do j=jp-1,jp+1
      do i=ip-1,ip+1
         u(i,j,1)  = w(i,j,2)/w(i,j,1)
         u(i,j,2)  = w(i,j,3)/w(i,j,1)
         u(i,j,3)  = p(i,j)/(gm1*w(i,j,1))
      end do
      end do

      if (mode.ne.0.and.jp.eq.2) then
         do i=ip-1,ip+1
            if (i,gt.itl.and.i.le.itu) then
               rev(i,1)  = -rev(i,2)
               u(i,1,1)  = -u(i,2,1)
               u(i,1,2)  = -u(i,2,2)
            end if
         end do
      end if
c
c     evaluate the viscous terms
c
      do 20 j=jp-1,jp
      do 10 i=ip-1,ip
c
c     evaluate the velocity and temperature gradients
c
      dxi       = xc(i+1,j+1,1)  -xc(i,j+1,1)  +xc(i+1,j,1)  -xc(i,j,1)
      dxj       = xc(i+1,j+1,1)  +xc(i,j+1,1)  -xc(i+1,j,1)  -xc(i,j,1)
      dyi       = xc(i+1,j+1,2)  -xc(i,j+1,2)  +xc(i+1,j,2)  -xc(i,j,2)
      dyj       = xc(i+1,j+1,2)  +xc(i,j+1,2)  -xc(i+1,j,2)  -xc(i,j,2)
      dsj       = 1./(dxi*dyj  -dyi*dxj)

      do l=1,3
         dui       = u(i+1,j+1,l)  -u(i,j+1,l)  +u(i+1,j,l)  -u(i,j,l)
         duj       = u(i+1,j+1,l)  +u(i,j+1,l)  -u(i+1,j,l)  -u(i,j,l)
         q(i,j,l,1)  = (dui*dyj  -duj*dyi)*dsj
         q(i,j,l,2)  = (duj*dxi  -dui*dxj)*dsj
      end do
c
c     evaluate the viscous stress tensor and dissipation
c
      ua        = u(i+1,j+1,1)  +u(i,j+1,1)  +u(i+1,j,1)  +u(i,j,1)
      va        = u(i+1,j+1,2)  +u(i,j+1,2)  +u(i+1,j,2)  +u(i,j,2)
      reva      = rev(i+1,j+1)  +rev(i,j+1)  +rev(i+1,j)  +rev(i,j)
      rlva      = rlv(i+1,j+1)  +rlv(i,j+1)  +rlv(i+1,j)  +rlv(i,j)
      rmu       = .25*(reva  +rlva)
      rlam      = 2.*rmu/3.
      rk        = gamma*(rlva/prn  +reva/prt)
      dq        = rlam*(q(i,j,1,1)  +q(i,j,2,2))

      q(i,j,1,1)  = rmu*(q(i,j,1,1)  +q(i,j,1,1))  -dq
      q(i,j,2,2)  = rmu*(q(i,j,2,2)  +q(i,j,2,2))  -dq
      q(i,j,1,2)  = rmu*(q(i,j,1,2)  +q(i,j,2,1))
      q(i,j,2,1)  = q(i,j,1,2)
      q(i,j,3,1)  = .25*(rk*q(i,j,3,1)  +ua*q(i,j,1,1)  +va*q(i,j,2,1))
      q(i,j,3,2)  = .25*(rk*q(i,j,3,2)  +ua*q(i,j,1,2)  +va*q(i,j,2,2))

   10 continue
   20 continue
c
c     viscous fluxes in the i direction
c
      j         = jp

      do i=ip-1,ip
         dx        = x(i,j,1)  -x(i,j-1,1)
         dy        = x(i,j,2)  -x(i,j-1,2)
         do l=1,3
            fs(i,j,l) = dy*(q(i,j,l,1)  +q(i,j-1,l,1))
     .                 -dx*(q(i,j,l,2)  +q(i,j-1,l,2))
         end do
      end do

      i         = ip
      do l=1,3
         vw(i,j,l) = sfil*vw(i,j,l)  -scf*(fs(i,j,l)  -fs(i-1,j,l))
      end do
c
c     viscous fluxes in the j direction
c
      do j=jp-1,jp
         dx        = x(i,j,1)  -x(i-1,j,1)
         dy        = x(i,j,2)  -x(i-1,j,2)
         do l=1,3
            fs(i,j,l) = dx*(q(i,j,l,2)  +q(i-1,j,l,2))
     .                 -dy*(q(i,j,l,1)  +q(i-1,j,l,1))
         end do
      end do

      j         = jp
      do l=1,3
         vw(i,j,l) = vw(i,j,l)  -scf*(fs(i,j,l)  -fs(i,j-1,l))
      end do

      return

      end
