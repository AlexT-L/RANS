      subroutine nsflux(il,jl,ie,je,
     & w,p,rlv,rev,
     & x,xc,
     & vw,
     & gamma,rm,scal,re,chord,prn,prt,
     & mode, rfil)
c     
c     arg list above is:
c     dims
c     flow_var
c     mesh_var
c     solv_var
c     flo_param
c     solv_param

c     
c     
c     *****************************************************************
c     *                                                               *
c     *   calculates the viscous flux                                 *
c     *                                                               *
c     *****************************************************************
c
c      use dims
c
c     ******************************************************************
c
c      use flo_var
c      use mesh_var
c      use solv_var
c
c     ******************************************************************
c
c      use flo_param
c      use solv_param
c      use mg_param
c
c     ******************************************************************
c
      implicit none

c
c     ******************************************************************
c
c     input variables
c
c     ******************************************************************
c     from dims
      integer, intent(in) :: il,jl,ie,je
      integer, intent(in) :: itl, itu


c     from flo_var
      real(8), intent(inout), dimension(:,:,:) :: w
      real(8), intent(inout), dimension(:,:) :: p
      real(8), intent(inout), dimension(:,:) :: rlv, rev

c     from mesh_var
      real(8), intent(inout), dimension(:,:,:) :: x, xc

c     from solv_var
      real(8), intent(inout), dimension(:,:,:) :: vw

c     from flo_param
      real(8), intent(in) :: gamma,rm,scal,re,chord,prn,prt

c     from solv_param
      real(8), intent(in)     :: rfil

c     from mg_param
      integer, intent(in)   :: mode



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
      real(8)     :: gm1,scf,sfil
      real(8)     :: ua,va,rmu,rlam,rk,rlva,reva,dq
      real(8)     :: dx,dy,dxi,dxj,dyi,dyj,dui,duj,dsj
c      real(8)     :: sigx,sigy,tauxy (ignoring for now)
c
c     ******************************************************************
c
      real(8), dimension(il,jl,3,2)        :: q
      real(8), dimension(ie,je,3)          :: u
      real(8), dimension(il,jl,3) :: fs
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
      do j=1,je
      do i=1,ie
         u(i,j,1)  = w(i,j,2)/w(i,j,1)
         u(i,j,2)  = w(i,j,3)/w(i,j,1)
         u(i,j,3)  = p(i,j)/(gm1*w(i,j,1))
      end do
      end do

c
c     evaluate the viscous terms
c
      do 20 j=1,jl
      do 10 i=1,il
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
      do 30 j=2,jl

      do i=1,il
         dx        = x(i,j,1)  -x(i,j-1,1)
         dy        = x(i,j,2)  -x(i,j-1,2)
         do l=1,3
            fs(i,j,l) = dy*(q(i,j,l,1)  +q(i,j-1,l,1))
     .                 -dx*(q(i,j,l,2)  +q(i,j-1,l,2))
         end do
      end do

      do i=2,il
         do l=1,3
            vw(i,j,l) = sfil*vw(i,j,l)  -scf*(fs(i,j,l)  -fs(i-1,j,l))
         end do
      end do

   30 continue
c
c     viscous fluxes in the j direction
c
      do j=1,jl
      do i=2,il
         dx        = x(i,j,1)  -x(i-1,j,1)
         dy        = x(i,j,2)  -x(i-1,j,2)
         do l=1,3
            fs(i,j,l) = dx*(q(i,j,l,2)  +q(i-1,j,l,2))
     .                 -dy*(q(i,j,l,1)  +q(i-1,j,l,1))
         end do
      end do
      end do

      do j=2,jl
      do i=2,il
         do l=1,3
            vw(i,j,l) = vw(i,j,l)  -scf*(fs(i,j,l)  -fs(i,j-1,l))
         end do
      end do
      end do

c     (not doing this anymore)
c     save the maximum shear stress at the wall 
c
c      do i=2,il
c         dx        = x(i,1,1)  -x(i-1,1,1)
c         dy        = x(i,1,2)  -x(i-1,1,2)
c         sigx      = .5*(q(i,1,1,1)  +q(i-1,1,1,1))
c         sigy      = .5*(q(i,1,2,2)  +q(i-1,1,2,2))
c         tauxy     = .5*(q(i,1,1,2)  +q(i-1,1,1,2))
c         tw(i)     = (dx*dy*(sigx  -sigy)  -(dx**2  -dy**2)*tauxy)/
c     .               (dx**2  +dy**2)
c      end do

      return

      end
