      subroutine dfluxc(ny,il,jl,
     & w,p,
     & porj,
     & fw, radi, radj,
c    no flo_params used
     & rfil,vis0)
c
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
      integer, intent(in) :: ny,il,jl

c     from flo_var
      real(8), intent(inout), dimension(:,:,:) :: w
      real(8), intent(inout), dimension(:,:) :: p

c     from mesh_var
      real(8), intent(in), dimension(:,:)   :: porj

c     from solv_var
      real(8), intent(inout), dimension(:,:,:) :: fw
      real(8), intent(inout), dimension(:,:)   :: radi,radj

c     from solv_param
      real(8), intent(in)     :: rfil, vis0
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
      real(8)     :: fis0,sfil
c
c     ******************************************************************
c
      real(8), dimension(il,jl,4)     :: fs
      real(8), dimension(il,jl)       :: dis
c
c     ******************************************************************
c
c      write(*,)
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
