      subroutine dflux(ny,il,jl,ie,je,ib,jb,
     & w,p,
     & pori, porj,
     & fw, radi, radj,
c    no flo_params used
     & rfil,vis2,vis4)
c     
c     arg list above is:
c
c     dims
c
c     flow_var
c     mesh_var
c     solv_var
c
c     flo_param
c     solv_param
c
c
c
c     ******************************************************************
c     *                                                                *
c     *   upwind biasing by artificial dissipation using               *
c     *   blended first and third order fluxes                         *
c     *                                                                *
c     ******************************************************************
c
c     w(i,j,1)  = density
c     w(i,j,2)  = momentum in x direction
c     w(i,j,3)  = momentum in y direction
c     w(i,j,4)  = total energy
c
c     ******************************************************************
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
c
c     ******************************************************************
c
      implicit none
      
c     ******************************************************************
c
c     input variables
c
c     ******************************************************************
c     from dims
      integer, intent(in) :: ny,il,jl,ie,je,ib,jb

c     from flo_var
      real(8), intent(inout), dimension(0:ib,0:jb,4) :: w
      real(8), intent(inout), dimension(0:ib,0:jb) :: p

c     from mesh_var
      real(8), intent(in), dimension(1:il,2:jl)   :: pori
      real(8), intent(in), dimension(2:il,1:jl)   :: porj

c     from solv_var
      real(8), intent(inout), dimension(0:ib,0:jb,4) :: fw
      real(8), intent(inout), dimension(0:ib,0:jb)   :: radi,radj

c     from solv_param
      real(8), intent(in)     :: rfil, vis2, vis4

c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: i,j,n
c
c     ******************************************************************
c
      real(8)     :: fis2,fis4,sfil,plim,tol,rad
c
c     ******************************************************************
c
      real(8), dimension(0:ib,0:jb)   :: dp
      real(8), dimension(0:ie,0:je)   :: d
      real(8), dimension(il,jl)       :: e,gs
      real(8), dimension(il,jl)       :: dis2,dis4
c
c     ******************************************************************
c
      fis2      = 2.*rfil*vis2
      fis4      = rfil*vis4/16.
      sfil      = 1.  -rfil
      plim      = .001
      tol       = .25
c
c     replace the energy by the enthalpy
c
      do j=0,jb
      do i=0,ib
         w(i,j,4)  = w(i,j,4)  +p(i,j)
      end do
      end do
c
c     dissipation in the i direction
c
      do j=2,jl
         do i=1,ie
            dp(i,j)   = abs((p(i+1,j)  -2.*p(i,j)  +p(i-1,j))/
     .                      (p(i+1,j)  +2.*p(i,j)  +p(i-1,j)  +plim))
         end do
         dp(0,j)   = dp(1,j)
         dp(ib,j)  = dp(ie,j)
      end do

      do j=2,jl
      do i=1,il
         rad       = min(radi(i+1,j),radi(i,j))
         dis2(i,j) = fis2*rad*
     .               min(tol,max(dp(i+2,j),dp(i+1,j),dp(i,j),dp(i-1,j)))
         dis4(i,j) = fis4*rad
         dis4(i,j) = dim(dis4(i,j),dis2(i,j))
      end do
      end do

      do n=1,4

      do j=2,jl
      do i=0,ie
         d(i,j)    = w(i+1,j,n)  -w(i,j,n)
      end do
      end do

      do j=2,jl
      do i=1,il
         e(i,j)    = d(i+1,j)  -2.*d(i,j)  +d(i-1,j)
      end do
      end do

      do j=2,jl
      do i=1,il
         gs(i,j)   = pori(i,j)*(dis2(i,j)*d(i,j)  -dis4(i,j)*e(i,j))
      end do
      end do

      do j=2,jl
      do i=2,il
         fw(i,j,n) = sfil*fw(i,j,n)  +gs(i-1,j)  -gs(i,j)
      end do
      end do

      end do

      if (ny.lt.3) go to 11
c
c     dissipation in the j direction
c
      do j=1,je
      do i=2,il
         dp(i,j)   = abs((p(i,j+1)  -2.*p(i,j)  +p(i,j-1))/
     .                   (p(i,j+1)  +2.*p(i,j)  +p(i,j-1)  +plim))
      end do
      end do

      do i=2,il
         rad         = min(radj(i,2),radj(i,1))
         dis2(i,1)   = fis2*rad*
     .                 min(tol,max(dp(i,3),dp(i,2),dp(i,1)))
         dis4(i,1)   = fis4*rad
         dis4(i,1)   = dim(dis4(i,1),dis2(i,1))
         rad         = min(radj(i,jl+1),radj(i,jl))
         dis2(i,jl)  = fis2*rad*
     .                 min(tol,max(dp(i,jl+1),dp(i,jl),dp(i,jl-1)))
         dis4(i,jl)  = fis4*rad
         dis4(i,jl)  = dim(dis4(i,jl),dis2(i,jl))
      end do

      do j=2,ny
      do i=2,il
         rad       = min(radj(i,j+1),radj(i,j))
         dis2(i,j) = fis2*rad*
     .               min(tol,max(dp(i,j+2),dp(i,j+1),dp(i,j),dp(i,j-1)))
         dis4(i,j) = fis4*rad
         dis4(i,j) = dim(dis4(i,j),dis2(i,j))
      end do
      end do

      do n=1,4

      do j=0,je
      do i=2,il
         d(i,j)    = w(i,j+1,n)  -w(i,j,n)
      end do
      end do

      do j=1,jl
      do i=2,il
         e(i,j)    = d(i,j+1)  -2.*d(i,j)  +d(i,j-1)
      end do
      end do

      do j=1,jl
      do i=2,il
         gs(i,j)   = porj(i,j)*(dis2(i,j)*d(i,j)  -dis4(i,j)*e(i,j))
      end do
      end do

      do j=2,jl
      do i=2,il
         fw(i,j,n) = fw(i,j,n)  +gs(i,j-1)  -gs(i,j)
      end do
      end do

      end do
c
c     replace the enthalpy by the energy
c
   11 do j=0,jb
      do i=0,ib
         w(i,j,4)  = w(i,j,4)  -p(i,j)
      end do
      end do

      return

      end
