      subroutine eflux
c
c     ******************************************************************
c     *                                                                *
c     *   euler flux                                                   *
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
      use mesh_var
      use solv_var
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
      integer  :: i,j,n
c
c     ******************************************************************
c
      real     :: xx,yx,xy,yy
      real     :: qsp,qsm,pa
c
c     ******************************************************************
c
c     real, dimension(il,jl,4) :: fs
c
c     ******************************************************************
c
c     flux in the i direction
c
      do j=2,jl
      do i=1,il
         xy        = x(i,j,1)  -x(i,j-1,1)
         yy        = x(i,j,2)  -x(i,j-1,2)
         pa        = p(i+1,j)  +p(i,j)
         qsp       = (yy*w(i+1,j,2)  -xy*w(i+1,j,3))/w(i+1,j,1)
         qsm       = (yy*w(i,j,2)  -xy*w(i,j,3))/w(i,j,1)
         fs(i,j,1) = qsp*w(i+1,j,1)  +qsm*w(i,j,1)
         fs(i,j,2) = qsp*w(i+1,j,2)  +qsm*w(i,j,2)  +yy*pa
         fs(i,j,3) = qsp*w(i+1,j,3)  +qsm*w(i,j,3)  -xy*pa
         fs(i,j,4) = qsp*(w(i+1,j,4)  +p(i+1,j))
     .              +qsm*(w(i,j,4)    +p(i,j))
      end do
      end do
c
c     accumulate the flux in the i direction
c
      do n=1,4
         do j=2,jl
         do i=2,il
            dw(i,j,n) = fs(i,j,n)  -fs(i-1,j,n)
         end do
         end do
      end do
c
c     flux in the j direction
c
      do j=1,jl
      do i=2,il
         xx        = x(i,j,1)  -x(i-1,j,1)
         yx        = x(i,j,2)  -x(i-1,j,2)
         pa        = p(i,j+1)  +p(i,j)
         qsp       = porj(i,j)*(xx*w(i,j+1,3)  -yx*w(i,j+1,2))/
     .                          w(i,j+1,1)
         qsm       = porj(i,j)*(xx*w(i,j,3)  -yx*w(i,j,2))/
     .                          w(i,j,1)
         fs(i,j,1) = qsp*w(i,j+1,1)  +qsm*w(i,j,1)
         fs(i,j,2) = qsp*w(i,j+1,2)  +qsm*w(i,j,2)  -yx*pa
         fs(i,j,3) = qsp*w(i,j+1,3)  +qsm*w(i,j,3)  +xx*pa
         fs(i,j,4) = qsp*(w(i,j+1,4)  +p(i,j+1))
     .              +qsm*(w(i,j,4)    +p(i,j))
      end do
      end do
c
c     accumulate the flux in the j direction
c
      do n=1,4
         do j=2,jl
         do i=2,il
            dw(i,j,n) = dw(i,j,n)  +fs(i,j,n)  -fs(i,j-1,n)
         end do
         end do
      end do

      return

      end
