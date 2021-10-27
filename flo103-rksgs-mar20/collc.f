      subroutine collc (ww,wwr,rlvc,revc)
c
c     ******************************************************************
c     *                                                                *
c     *   collects the residuals and transfers the solution            *
c     *   to a coarser mesh                                            *
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
      use dimsc
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
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      real, dimension(0:iib,0:jjb,4)       :: ww
      real, dimension(iie,jje,4)           :: wwr
      real, dimension(iie,jje)             :: rlvc,revc
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: i,j,ii,jj,n
c
c     ******************************************************************
c
      real, dimension(0:ib,0:jb)           :: rw
c
c     ******************************************************************
c
c     transfer the residuals to the coarse mesh
c
      do n=1,4

         jj        = 1
      do j=2,jl,2
         jj        = jj  +1
         ii        = 1
      do i=2,il,2
         ii        = ii  +1
         wwr(ii,jj,n)  = dw(i,j,n)     +dw(i+1,j,n)
     .                  +dw(i,j+1,n)   +dw(i+1,j+1,n)
      end do
      end do

      end do
c
c     transfer the flow variables to the coarse mesh
c     using area weighting to preserve conservation
c     of mass,momentum and energy
c
      do n=1,4

      do j=0,jb
      do i=0,ib
         rw(i,j)   = vol(i,j)*w(i,j,n)
      end do
      end do

         jj        = 0
      do j=0,je,2
         jj        = jj  +1
         ii        = 0
      do i=0,ie,2
         ii        = ii  +1
         ww(ii,jj,n)   = (rw(i,j)     +rw(i+1,j)
     .                   +rw(i,j+1)   +rw(i+1,j+1))/
     .                   (vol(i,j)    +vol(i+1,j)
     .                   +vol(i,j+1)  +vol(i+1,j+1))
      end do
      end do

      end do

      if (kvis.gt.0) then
c
c     transfer the molecular and turbulent viscosity to the coarse grid
c
         jj        = 1
      do j=2,jl,2
         jj        = jj  +1
         ii        = 1
      do i=2,il,2
         ii        = ii  +1
         rlvc(ii,jj) = .25*(rlv(i,j)    +rlv(i+1,j)
     .                     +rlv(i,j+1)  +rlv(i+1,j+1))
         revc(ii,jj) = .25*(rev(i,j)    +rev(i+1,j)
     .                     +rev(i,j+1)  +rev(i+1,j+1))
      end do
      end do
c
c     set the boundary values at i=1 and i=ie
c
         jj        = 1
      do j=2,jl,2
         jj        = jj  +1
         rlvc(1,jj)    = .5*(rlv(1,j)   +rlv(1,j+1))
         rlvc(iie,jj)  = .5*(rlv(ie,j)  +rlv(ie,j+1))
         revc(1,jj)    = .5*(rev(1,j)   +rev(1,j+1))
         revc(iie,jj)  = .5*(rev(ie,j)  +rev(ie,j+1))
      end do
c
c     set the boundary values at j=1 and j=je
c
         ii        = 1
      do i=2,il,2
         ii        = ii  +1
         rlvc(ii,1)    = .5*(rlv(i,1)   +rlv(i+1,1))
         rlvc(ii,jje)  = .5*(rlv(i,je)  +rlv(i+1,je))
         revc(ii,1)    = .5*(rev(i,1)   +rev(i+1,1))
         revc(ii,jje)  = .5*(rev(i,je)  +rev(i+1,je))
      end do

      end if

      return

      end
