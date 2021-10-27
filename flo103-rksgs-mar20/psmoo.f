      subroutine psmoo
c
c     ******************************************************************
c     *                                                                *
c     *   implicit smoothing                                           *
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
      use dims
c
c     ******************************************************************
c
      use solv_var
c
c     ******************************************************************
c
      use solv_param
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
      integer  :: i,j
      integer  :: ii,jj,j1,m2
      integer  :: l,m,n
c
c     ******************************************************************
c
      real     :: r,t
c
c     ******************************************************************
c
      real, dimension(il,2*jl)     :: d,e,sfl
      real, dimension(il)          :: ap,am
c
c     ******************************************************************
c
c     smoothing in the i direction
c
      do j=2,jl
         r         = max(csmoop,rfl0*(rfli(1,j)  +rfli(2,j)))
         am(j)     = .25*smoopi*dim(r**2,1.)
         r         = max(csmoop,rfl0*(rfli(2,j)  +rfli(3,j)))
         ap(j)     = .25*smoopi*dim(r**2,1.)
         t         = 1./(1.  +ap(j)  +am(j))
         d(2,j)    = t*ap(j)
         dw(2,j,1) = t*dw(2,j,1)
         dw(2,j,2) = t*dw(2,j,2)
         dw(2,j,3) = t*dw(2,j,3)
         dw(2,j,4) = t*dw(2,j,4)
      end do

      do i=3,il
      do j=2,jl
         am(j)     = ap(j)
         r         = max(csmoop,rfl0*(rfli(i,j)  +rfli(i+1,j)))
         ap(j)     = .25*smoopi*dim(r**2,1.)
         t         = 1./(1.  +ap(j)  +am(j)  -am(j)*d(i-1,j))
         d(i,j)    = t*ap(j)
         dw(i,j,1) = t*(dw(i,j,1)  +am(j)*dw(i-1,j,1))
         dw(i,j,2) = t*(dw(i,j,2)  +am(j)*dw(i-1,j,2))
         dw(i,j,3) = t*(dw(i,j,3)  +am(j)*dw(i-1,j,3))
         dw(i,j,4) = t*(dw(i,j,4)  +am(j)*dw(i-1,j,4))
      end do
      end do

      do i=nx,2,-1
      do j=2,jl
         dw(i,j,1) = dw(i,j,1)  +d(i,j)*dw(i+1,j,1)
         dw(i,j,2) = dw(i,j,2)  +d(i,j)*dw(i+1,j,2)
         dw(i,j,3) = dw(i,j,3)  +d(i,j)*dw(i+1,j,3)
         dw(i,j,4) = dw(i,j,4)  +d(i,j)*dw(i+1,j,4)
      end do
      end do

      if (ny.lt.3)  return
c
c     smoothing in the j direction
c
      do i=itl+1,itu
         r         = max(csmoop,rfl0*(rflj(i,2)  +rflj(i,3)))
         ap(i)     = .25*smoopj*dim(r**2,1.)
         t         = 1./(1.  +ap(i))
         d(i,2)    = t*ap(i)
         dw(i,2,1) = t*dw(i,2,1)
         dw(i,2,2) = t*dw(i,2,2)
         dw(i,2,3) = t*dw(i,2,3)
         dw(i,2,4) = t*dw(i,2,4)
      end do

      do j=3,jl
      do i=itl+1,itu
         am(i)     = ap(i)
         r         = max(csmoop,rfl0*(rflj(i,j)  +rflj(i,j+1)))
         ap(i)     = .25*smoopj*dim(r**2,1.)
         t         = 1./(1.  +ap(i)  +am(i)  -am(i)*d(i,j-1))
         d(i,j)    = t*ap(i)
         dw(i,j,1) = t*(dw(i,j,1)  +am(i)*dw(i,j-1,1))
         dw(i,j,2) = t*(dw(i,j,2)  +am(i)*dw(i,j-1,2))
         dw(i,j,3) = t*(dw(i,j,3)  +am(i)*dw(i,j-1,3))
         dw(i,j,4) = t*(dw(i,j,4)  +am(i)*dw(i,j-1,4))
      end do
      end do

      do j=ny,2,-1
      do i=itl+1,itu
         dw(i,j,1) = dw(i,j,1)  +d(i,j)*dw(i,j+1,1)
         dw(i,j,2) = dw(i,j,2)  +d(i,j)*dw(i,j+1,2)
         dw(i,j,3) = dw(i,j,3)  +d(i,j)*dw(i,j+1,3)
         dw(i,j,4) = dw(i,j,4)  +d(i,j)*dw(i,j+1,4)
      end do
      end do
c
c     smoothing across the wake cut
c
      m2        = jl  +jl

      do 20 n=1,4

      do i=2,itl
         ii        = il  -i  +2
      do j=2,jl
         l         = jl  -j  +2
         m         = jl  +j  -1
         e(i,l)    = dw(i,j,n)
         e(i,m)    = dw(ii,j,n)
         sfl(i,l)  = rflj(i,j)
         sfl(i,m)  = rflj(ii,j)
         d(i,l)    = 0.
         d(i,m)    = 0.
      end do
         d(i,1)    = 0.
         d(i,m2)   = 0.
         sfl(i,1)  = rflj(i,je)
         sfl(i,m2) = rflj(ii,je)
      end do

      j1        = jl  +jl  -1
      do i=2,itl
         r         = max(csmoop,rfl0*(sfl(i,1)  +sfl(i,2)))
         am(i)     = .25*smoopj*dim(r**2,1.)
         r         = max(csmoop,rfl0*(sfl(i,2)  +sfl(i,3)))
         ap(i)     = .25*smoopj*dim(r**2,1.)
         t         = 1./(1.  +ap(i) +am(i))
         d(i,2)    = t*ap(i)
         e(i,2)    = t*e(i,2)
      end do

      do j=3,j1
      do i=2,itl
         am(i)     = ap(i)
         r         = max(csmoop,rfl0*(sfl(i,j)  +sfl(i,j+1)))
         ap(i)     = .25*smoopj*dim(r**2,1.)
         t         = 1./(1.  +ap(i)  +am(i)  -am(i)*d(i,j-1))
         d(i,j)    = t*ap(i)
         e(i,j)    = t*(e(i,j)  +am(i)*e(i,j-1))
      end do
      end do

      do j=j1-1,2,-1
      do i=2,itl
         e(i,j)    = e(i,j)  +d(i,j)*e(i,j+1)
      end do
      end do

      do j=2,jl
         l         = jl  -j  +2
         m         = jl  +j  -1
      do i=2,itl
         ii        = il  -i  +2
         dw(i,j,n)   = e(i,l)
         dw(ii,j,n)  = e(i,m)
      end do
      end do

   20 continue

      return

      end
