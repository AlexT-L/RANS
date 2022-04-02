      subroutine step(il,jl,ie,je,ib,jb,itl,
     & w,p,rlv,rev,
     & x,vol,
     & rfl,rfli,rflj,radi,radj,dtl,dtlc,
     & gamma,rm,re,prn,prt,kvis,
     & iprec, cfl, vt, adis)
c
c     ******************************************************************
c     *                                                                *
c     *   calculates the permissible time step                         *
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
c     use dims
c
c     ******************************************************************
c
c     use flo_var
c     use mesh_var
c     use solv_var
c
c     ******************************************************************
c
c     use flo_param
c     use solv_param
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
      integer, intent(in) :: il,jl,ie,je,ib,jb,itl

c     from flo_var
      real(8), intent(in), dimension(0:ib,0:jb,4) :: w
      real(8), intent(in), dimension(0:ib,0:jb) :: p,rlv,rev

c     from mesh_var
      real(8), intent(in), dimension(1:il,1:jl,2) :: x
      real(8), intent(in), dimension(0:ib,0:jb) :: vol

c     from solv_var
      real(8), intent(inout), dimension(0:ib,0:jb)   :: rfl,rfli,rflj
      real(8), intent(inout), dimension(0:ib,0:jb)   :: radi,radj
      real(8), intent(inout), dimension(0:ib,0:jb)   :: dtl,dtlc

c     from flo_param
      real(8), intent(in)     :: gamma,rm,re,prn,prt,kvis

c     from solv_param
      real(8), intent(in)     :: iprec, cfl, vt, adis

c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: i,j,ii
      integer  :: imin,jmin
c
c     ******************************************************************
c
      real(8)     :: slim,rlim,b,dtmin
      real(8)     :: xx,yx,xy,yy,qs,cc,cs
      real(8)     :: v1,v2,rmu,rk,dsi,dsj,vsi,vsj,dtv
      real(8)     :: dpi,dpj,a
c
c     ******************************************************************
c
      real(8), dimension(0:ib,0:jb)   :: s
c
c     ******************************************************************
c
      slim      = .001
      rlim      = .001
      b         = 4.
      dtmin     = 0.
      imin      = 0
      jmin      = 0
c
c     permissible time step
c
      do j=2,jl
      do i=2,il
         cc        = gamma*p(i,j)/max(w(i,j,1),rlim)
         xx        = .5*(x(i,j,1) -x(i-1,j,1) +x(i,j-1,1) -x(i-1,j-1,1))
         yx        = .5*(x(i,j,2) -x(i-1,j,2) +x(i,j-1,2) -x(i-1,j-1,2))
         xy        = .5*(x(i,j,1) -x(i,j-1,1) +x(i-1,j,1) -x(i-1,j-1,1))
         yy        = .5*(x(i,j,2) -x(i,j-1,2) +x(i-1,j,2) -x(i-1,j-1,2))
         qs        = (yy*w(i,j,2) -xy*w(i,j,3))/w(i,j,1)
         cs        = cc*(xy**2  +yy**2)
         radi(i,j) = abs(qs)  +sqrt(cs)
         qs        = (xx*w(i,j,3)  -yx*w(i,j,2))/w(i,j,1)
         cs        = cc*(xx**2  +yx**2)
         radj(i,j) = abs(qs)  +sqrt(cs)
         dtl(i,j)  = 1./(radi(i,j)  +radj(i,j))
         dtlc(i,j) = radi(i,j)  +radj(i,j)
      end do
      end do
c
c     pressure or entropy switch
c
c     if (kvis.eq.0) then
         do j=0,jb
         do i=0,ib
            s(i,j)    = p(i,j)
         end do
         end do
c     else
c        do j=0,jb
c        do i=0,ib
c           s(i,j)    = p(i,j)/w(i,j,1)**gamma
c        end do
c        end do
c     end if
c
c     adaptive time step
c
      do j=2,jl
      do i=2,il
         dpi       = abs((s(i+1,j)  -2.*s(i,j)  +s(i-1,j))/
     .                   (s(i+1,j)  +2.*s(i,j)  +s(i-1,j)  +slim))
         dpj       = abs((s(i,j+1)  -2.*s(i,j)  +s(i,j-1))/
     .                   (s(i,j+1)  +2.*s(i,j)  +s(i,j-1)  +slim))
         rfl(i,j)  = 1./(1.  +min(dim(cfl,1.),b*(dpi  +dpj)))
      end do
      end do

      if (vt.eq.0.) go to 11
c
c     fixed time step
c
      dtmin     = dtl(2,2)
      do j=2,jl
      do i=2,il
         if (dtl(i,j).le.dtmin) then
            dtmin     = dtl(i,j)
            imin      = i
            jmin      = j
         end if
      end do
      end do

      do j=2,jl
      do i=2,il
         rfl(i,j)  = dtmin/dtl(i,j)
      end do
      end do
c
c     option to rescale the dissipative coefficients
c
   11 do j=2,jl
      do i=2,il
         if (iprec.ne.0) then
c           rfli(i,j) = sqrt(radj(i,j)/radi(i,j))
c           rflj(i,j) = sqrt(radi(i,j)/radj(i,j))
            rfli(i,j) = (radj(i,j)/radi(i,j))**.25
            rflj(i,j) = (radi(i,j)/radj(i,j))**.25
c           rfli(i,j) = 1.
c           rflj(i,j) = 1.
         end if
         a         = (radi(i,j)/radj(i,j))**adis
         radi(i,j) = radi(i,j)*(1.  +1./a)
         radj(i,j) = radj(i,j)*(1.  +a)
         if (iprec.eq.0) then
            rfli(i,j) = radi(i,j)/dtlc(i,j)
            rflj(i,j) = radj(i,j)/dtlc(i,j)
         end if
      end do
      end do
c
c
c     reduce the artificial dissipation for viscous flows
c     and adjust the time step estimate
c
      if (kvis.gt.0) then

      v1        = sqrt(gamma)*rm/re
      v2        = 0.
      if (kvis.gt.1) v2 = v1

      do j=2,jl
      do i=2,il
         rk        = gamma *(v1*rlv(i,j)/prn + v2*rev(i,j)/prt)/w(i,j,1)
         rmu       = (v1*rlv(i,j)+v2*rev(i,j))/w(i,j,1)
         xx        = .5*(x(i,j,1) -x(i-1,j,1) +x(i,j-1,1) -x(i-1,j-1,1))
         yx        = .5*(x(i,j,2) -x(i-1,j,2) +x(i,j-1,2) -x(i-1,j-1,2))
         xy        = .5*(x(i,j,1) -x(i,j-1,1) +x(i-1,j,1) -x(i-1,j-1,1))
         yy        = .5*(x(i,j,2) -x(i,j-1,2) +x(i-1,j,2) -x(i-1,j-1,2))
         dsi       = xy**2  +yy**2
         dsj       = xx**2  +yx**2
         vsi       = (rk*dsi + rmu*sqrt(dsi*dsj)/6.)/vol(i,j)
         vsj       = (rk*dsj + rmu*sqrt(dsi*dsj)/6.)/vol(i,j)
         dtv       = dtlc(i,j)+4.*(vsi+vsj)
         dtl(i,j)  = 1./dtv
         radi(i,j) = dim(radi(i,j),vsi)
         radj(i,j) = dim(radj(i,j),vsj)
      end do
      end do

      end if
c
c     set boundary values at i=1 and i=ie
c
      do j=2,jl
         radi(1,j)   = radi(2,j)
         radi(ie,j)  = radi(il,j)
         rfl(1,j)    = rfl(2,j)
         rfl(ie,j)   = rfl(il,j)
         rfli(1,j)   = rfli(2,j)
         rfli(ie,j)  = rfli(il,j)
         rflj(1,j)   = rflj(2,j)
         rflj(ie,j)  = rflj(il,j)
         dtl(1,j)    = dtl(2,j)
         dtl(ie,j)   = dtl(il,j)
      end do
c
c     set boundary values at j=1 and j=je
c
      do i=1,ie
         radj(i,1)   = radj(i,2)
         radj(i,je)  = radj(i,jl)
         rfl(i,1)    = rfl(i,2)
         rfl(i,je)   = rfl(i,jl)
         rfli(i,1)   = rfli(i,2)
         rfli(i,je)  = rfli(i,jl)
         rflj(i,1)   = rflj(i,2)
         rflj(i,je)  = rflj(i,jl)
         dtl(i,1)    = dtl(i,2)
         dtl(i,je)   = dtl(i,jl)
      end do
c
c     set boundary values along the cut
c
      do i=1,itl
         ii        = ib  -i
         radj(ii,1)  = radj(i,2)
         radj(i,1)   = radj(ii,2)
         rfl(ii,1)   = rfl(i,2)
         rfl(i,1)    = rfl(ii,2)
         rfli(ii,1)  = rfli(i,2)
         rfli(i,1)   = rfli(ii,2)
         rflj(ii,1)  = rflj(i,2)
         rflj(i,1)   = rflj(ii,2)
         dtl(ii,1)   = dtl(i,2)
         dtl(i,1)    = dtl(ii,2)
      end do

      return

      end
