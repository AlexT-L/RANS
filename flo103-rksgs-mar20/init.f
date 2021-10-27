      subroutine init
c
c     ******************************************************************
c     *                                                                *
c     *   initial flow field                                           *
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
      use flo_var
c
c     ******************************************************************
c
      use flo_param
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
c     initialize the flow variables at free stream values
c
      do j=0,jb
      do i=0,ib
         w(i,j,1)  = rho0
         w(i,j,2)  = rho0*u0
         w(i,j,3)  = rho0*v0
         w(i,j,4)  = rho0*h0  -p0
         p(i,j)    = p0
         rlv(i,j)  = 0.
         rev(i,j)  = 0.
      end do
      end do

      return

      end
