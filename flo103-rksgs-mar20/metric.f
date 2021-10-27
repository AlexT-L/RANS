      subroutine metric
c
c     ******************************************************************
c     *                                                                *
c     *   cell areas and cell centers (dual mesh)                      *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use mesh_var
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
c
c     ******************************************************************
c
c     set the areas of interior and exterior cells to zero
c
      do j=0,jb
      do i=0,ib
         vol(i,j)  = 0.
         xc(i,j,1) = 0.
         xc(i,j,2) = 0.
      end do
      end do
c
c     calculate the cell areas and cell centers
c
      do j=2,jl
      do i=2,il
         vol(i,j)  = .5*((x(i,j,1)  -x(i-1,j-1,1))
     .                   *(x(i-1,j,2)  -x(i,j-1,2))
     .                  -(x(i,j,2)  -x(i-1,j-1,2))
     .                   *(x(i-1,j,1)  -x(i,j-1,1)))
         xc(i,j,1) = .25*(x(i,j,1)    +x(i-1,j,1)
     .                   +x(i,j-1,1)  +x(i-1,j-1,1))
         xc(i,j,2) = .25*(x(i,j,2)    +x(i-1,j,2)
     .                   +x(i,j-1,2)  +x(i-1,j-1,2))
      end do
         vol(1,j)    = vol(2,j)
         vol(ie,j)   = vol(il,j)
         xc(1,j,1)   = x(1,j,1)   +x(1,j-1,1)   -xc(2,j,1)
         xc(1,j,2)   = x(1,j,2)   +x(1,j-1,2)   -xc(2,j,2)
         xc(ie,j,1)  = x(il,j,1)  +x(il,j-1,1)  -xc(il,j,1)
         xc(ie,j,2)  = x(il,j,2)  +x(il,j-1,2)  -xc(il,j,2)
      end do
c
c     set boudary values at j=1 and j=je
c
      do i=2,il
         vol(i,1)    = vol(i,2)
         vol(i,je)   = vol(i,jl)
         xc(i,1,1)   = x(i,1,1)   +x(i-1,1,1)   -xc(i,2,1)
         xc(i,1,2)   = x(i,1,2)   +x(i-1,1,2)   -xc(i,2,2)
         xc(i,je,1)  = x(i,jl,1)  +x(i-1,jl,1)  -xc(i,jl,1)
         xc(i,je,2)  = x(i,jl,2)  +x(i-1,jl,2)  -xc(i,jl,2)
      end do
c
c     set boundary values along the cut in the c mesh
c
      do i=2,itl
         vol(ib-i,1)   = vol(i,2)
         xc(ib-i,1,1)  = xc(i,2,1)
         xc(ib-i,1,2)  = xc(i,2,2)
         vol(i,1)    = vol(ib-i,2)
         xc(i,1,1)   = xc(ib-i,2,1)
         xc(i,1,2)   = xc(ib-i,2,2)
      end do
c
c     set corner values
c
      vol(1,1)    = vol(2,1)
      vol(1,je)   = vol(2,je)
      vol(ie,1)   = vol(il,1)
      vol(ie,je)  = vol(il,je)

      xc(1,1,1)   = xc(2,1,1)    +xc(1,2,1)    -xc(2,2,1)
      xc(1,1,2)   = xc(2,1,2)    +xc(1,2,2)    -xc(2,2,2)
      xc(ie,1,1)  = xc(il,1,1)   +xc(ie,2,1)   -xc(il,2,1)
      xc(ie,1,2)  = xc(il,1,2)   +xc(ie,2,2)   -xc(il,2,2)
      xc(1,je,1)  = xc(2,je,1)   +xc(1,jl,1)   -xc(2,jl,1)
      xc(1,je,2)  = xc(2,je,2)   +xc(1,jl,2)   -xc(2,jl,2)
      xc(ie,je,1) = xc(il,je,1)  +xc(ie,jl,1)  -xc(il,jl,1)
      xc(ie,je,2) = xc(il,je,2)  +xc(ie,jl,2)  -xc(il,jl,2)

      return

      end
