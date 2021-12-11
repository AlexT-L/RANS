      subroutine bcwall(ny, il, ie, ib, itl, itu,
     & w, p, rev,
     & x,
     & rm, sa, kvis,
     & isym)
c
c     ******************************************************************
c     *                                                                *
c     *   wall boundary condition                                      *
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
c     dims
      integer, intent(in)      :: ny, il, ie, ib
      integer, intent(in)         :: itl,itu      
c
c     ******************************************************************
c
c      use flo_var
      real(8), intent(inout), dimension(:,:,:) :: w
      real(8), intent(inout), dimension(:,:) :: p

c      use mesh_var
      real, dimension(:,:,:), intent(in)    , allocatable :: x

c
c     ******************************************************************
c
c      use flo_param
      real      :: rm,sa
      integer, intent(in) :: kvis

c      use solv_param
      integer, intent(in)   :: isym
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
      integer  :: i,i1,ile
c
c     ******************************************************************
c
      real     :: a,b,sxn,xxn,yxn,xy,yy
      real(8)     :: ra,rua,rva,px,py,qs,qt,qn,gxy,qxy
c
c     ******************************************************************
c
      real, dimension(ie)      :: sx,xx,yx,xxx,yxx
c
c     ******************************************************************
c
c     set values below the cut in the c mesh
c
      do i=0,itl
         w(i,1,1)  = w(ib-i,2,1)
         w(i,1,2)  = w(ib-i,2,2)
         w(i,1,3)  = w(ib-i,2,3)
         w(i,1,4)  = w(ib-i,2,4)
         p(i,1)    = p(ib-i,2)
         w(ib-i,1,1) = w(i,2,1)
         w(ib-i,1,2) = w(i,2,2)
         w(ib-i,1,3) = w(i,2,3)
         w(ib-i,1,4) = w(i,2,4)
         p(ib-i,1)   = p(i,2)
      end do

      if (kvis.gt.0) then
c
c     set the viscous no-slip boundary condition
c
      i1        = itl  +1
      do i=i1,itu
         w(i,1,1)  =  w(i,2,1)
         w(i,1,2)  = -w(i,2,2)
         w(i,1,3)  = -w(i,2,3)
         w(i,1,4)  =  w(i,2,4)
         p(i,1)    =  p(i,2)
      end do

      else
c
c     set the euler boundary condition
c
      i1        = itu  +1
      do i=itl,i1
         a         = x(i,1,1)  -x(i-1,1,1)
         b         = x(i,1,2)  -x(i-1,1,2)
         sx(i)     = 1./sqrt(a**2  +b**2)
         xx(i)     = a*sx(i)
         yx(i)     = b*sx(i)
      end do
c
c     calculate the wall curvature
c
      i1        = itl  +1
      do i=i1,itu
         xxx(i)    = .5*(xx(i+1)  -xx(i-1))
         yxx(i)    = .5*(yx(i+1)  -yx(i-1))
      end do
c
c     special treatment of a sharp leading edge
c
      if (isym.eq.2) then
         ile       = ie/2
         a         = x(ile,2,1)  -x(ile,1,1)
         b         = x(ile,2,2)  -x(ile,1,2)
         sxn       = 1./sqrt(a**2  +b**2)
         xxn       = a*sxn
         yxn       = b*sxn
         xxx(ile)    = .5*(xxn  -xx(ile-1))
         yxx(ile)    = .5*(yxn  -yx(ile-1))
         xxx(ile+1)  = .5*(xx(ile+2)  +xxn)
         yxx(ile+1)  = .5*(yx(ile+2)  +yxn)
      end if
c
c     extrapolation using normal pressure gradient at surface
c
      do i=i1,itu
         qt        = xx(i)*w(i,2,2)  +yx(i)*w(i,2,3)
         qn        = yx(i)*w(i,2,2)  -xx(i)*w(i,2,3)
         w(i,1,1)  = w(i,2,1)
         w(i,1,2)  = xx(i)*qt  -yx(i)*qn
         w(i,1,3)  = xx(i)*qn  +yx(i)*qt
         w(i,1,4)  = w(i,2,4)
         ra        = w(i,2,1)  +w(i,1,1)
         rua       = w(i,2,2)  +w(i,1,2)
         rva       = w(i,2,3)  +w(i,1,3)
         px        = .5*(p(i+1,2)  -p(i-1,2))
         xy        = .5*(x(i,2,1)  -x(i,1,1)  +x(i-1,2,1)  -x(i-1,1,1))
         yy        = .5*(x(i,2,2)  -x(i,1,2)  +x(i-1,2,2)  -x(i-1,1,2))
         qs        = (yy*rua  -xy*rva)/ra
         gxy       = xx(i)*xy  +yx(i)*yy
         qxy       = .5*qs*(xxx(i)*rva  -yxx(i)*rua)
         py        = (px*gxy  +qxy)*sx(i)
         if (ny.lt.3) py = p(i,3)  -p(i,2)
         p(i,1)    = dim(p(i,2),py)
         w(i,1,4)  = w(i,2,4)  +p(i,2)  -p(i,1)
      end do

      end if
c
c     update eddy viscosity boundary conditions
c   
      if (mode.ne.0) then
         do i=itl+1,itu
            rev(i,1)  = -rev(i,2)
         end do
      end if

      return

      end
