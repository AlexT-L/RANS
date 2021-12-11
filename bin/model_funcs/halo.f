      subroutine halo(il, jl, ie, je, ib, jb, itl, itu,
     & w, p,
     & x, vol)
c
c     ******************************************************************
c     *                                                                *
c     *   sets values in the halo surrounding the mesh                 *
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
      integer, intent(in)           :: il, jl, ie, je, ib, jb

      integer, intent(in)           :: itl,itu
c
c     ******************************************************************
c
c      use flo_var
      real, dimension(:,:,:), intent(inout)  :: w
      real, dimension(:,:), intent(inout)    :: p

c      use mesh_var
      real, dimension(:,:,:), intent(inout) :: x
      real, dimension(:,:), intent(inout)   :: vol
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
      integer  :: i,j,ii,i1,i2
c
c     ******************************************************************
c
c     extended values at i=0
c
      do j=1,je
         w(0,j,1)  = w(1,j,1)  +w(1,j,1)  -w(2,j,1)
         w(0,j,2)  = w(1,j,2)  +w(1,j,2)  -w(2,j,2)
         w(0,j,3)  = w(1,j,3)  +w(1,j,3)  -w(2,j,3)
         w(0,j,4)  = w(1,j,4)  +w(1,j,4)  -w(2,j,4)
         p(0,j)    = p(1,j)    +p(1,j)    -p(2,j)
         vol(0,j)  = 0.
      end do
c
c     extended values at i=ib
c
      do j=1,je
         w(ib,j,1)  = w(ie,j,1)  +w(ie,j,1)  -w(il,j,1)
         w(ib,j,2)  = w(ie,j,2)  +w(ie,j,2)  -w(il,j,2)
         w(ib,j,3)  = w(ie,j,3)  +w(ie,j,3)  -w(il,j,3)
         w(ib,j,4)  = w(ie,j,4)  +w(ie,j,4)  -w(il,j,4)
         p(ib,j)    = p(ie,j)    +p(ie,j)    -p(il,j)
         vol(ib,j)  = 0.
      end do
c
c     extended values at j=0
c     allowing for the cut in the c-mesh
c
      do i=0,ib
         ii        = ib  -i
         w(i,0,1)  = w(ii,3,1)
         w(i,0,2)  = w(ii,3,2)
         w(i,0,3)  = w(ii,3,3)
         w(i,0,4)  = w(ii,3,4)
         p(i,0)    = p(ii,3)
         vol(i,1)  = vol(ii,2)
         vol(i,0)  = vol(ii,3)
      end do
c
c     extended values at j=0 over the surface of the wing
c
      i1        = itl  +1
      i2        = itu

      do i=i1,i2
         w(i,0,1)  = w(i,1,1)  +w(i,1,1)  -w(i,2,1)
         w(i,0,2)  = w(i,1,2)  +w(i,1,2)  -w(i,2,2)
         w(i,0,3)  = w(i,1,3)  +w(i,1,3)  -w(i,2,3)
         w(i,0,4)  = w(i,1,4)  +w(i,1,4)  -w(i,2,4)
         p(i,0)    = p(i,1)    +p(i,1)    -p(i,2)
         vol(i,1)  = vol(i,2)
         vol(i,0)  = 0.
      end do
c
c     extended values at j=jb
c
      do i=0,ib
         w(i,jb,1) = w(i,je,1)  +w(i,je,1)  -w(i,jl,1)
         w(i,jb,2) = w(i,je,2)  +w(i,je,2)  -w(i,jl,2)
         w(i,jb,3) = w(i,je,3)  +w(i,je,3)  -w(i,jl,3)
         w(i,jb,4) = w(i,je,4)  +w(i,je,4)  -w(i,jl,4)
         p(i,jb)   = p(i,je)    +p(i,je)    -p(i,jl)
         vol(i,jb) = 0.
      end do

      return

      end
