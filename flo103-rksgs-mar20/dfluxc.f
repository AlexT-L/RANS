      subroutine dfluxc
c
c     ******************************************************************
c     *                                                                *
c     *   artificial dissipation on coarse grids                       *
c     *   first order flux scaled to the spectral radius               *
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
      integer  :: i,j,n
c
c     ******************************************************************
c
      real     :: fis0,sfil
c
c     ******************************************************************
c
c     real, dimension(il,jl,4)     :: fs
      real, dimension(il,jl)       :: dis
c
c     ******************************************************************
c
      fis0      = rfil*abs(vis0)/8.
      sfil      = 1.  -rfil
c
c     dissipation in the i direction
c
      do j=2,jl

         do i=1,il
c           dis(i,j)  = fis0*(radi(i+1,j)  +radi(i,j))
            dis(i,j)  = fis0*min(radi(i+1,j),radi(i,j))
         end do

         do n=1,4
            do i=1,il
               fs(i,j,n) = dis(i,j)*(w(i+1,j,n)  -w(i,j,n))
            end do
         end do

         do i=1,il
            fs(i,j,4) = fs(i,j,4)  +dis(i,j)*(p(i+1,j)  -p(i,j))
         end do

         do n=1,4
            do i=2,il
               fw(i,j,n) = sfil*fw(i,j,n)  -fs(i,j,n) +fs(i-1,j,n)
            end do
         end do

      end do

      if (ny.lt.3) return
c
c     dissipation in the j direction
c
      do j=1,jl
      do i=2,il
c        dis(i,j)  = fis0*(radj(i,j+1)  +radj(i,j))
         dis(i,j)  = fis0*porj(i,j)*min(radj(i,j+1),radj(i,j))
      end do
      end do

      do n=1,4
         do j=1,jl
         do i=2,il
            fs(i,j,n) = dis(i,j)*(w(i,j+1,n)  -w(i,j,n))
         end do
         end do
      end do

      do j=1,jl
      do i=2,il
         fs(i,j,4) = fs(i,j,4)  +dis(i,j)*(p(i,j+1)  -p(i,j))
      end do
      end do

      do n=1,4
         do j=2,jl
         do i=2,il
            fw(i,j,n) = fw(i,j,n)  -fs(i,j,n)  +fs(i,j-1,n)
         end do
         end do
      end do

      return

      end
