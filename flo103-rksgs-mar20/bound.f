      subroutine bound
c
c     ******************************************************************
c     *                                                                *
c     *   sets flags indicating wall and far field boundaries          *
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
      use geo_param
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
      integer  :: ile,ite,i1,i2
c
c     ******************************************************************
c
      ite       = .5000005*xte*(il  -1)
      ile       = il/2  +1
      itl       = ile  -ite
      itu       = ile  +ite
c
c     set the porosity and far field masks to unity
c
      do j=1,jl
      do i=1,il
         pori(i,j) = 1.
         porj(i,j) = 1.
         fint(i,j) = 1.
      end do
      end do
c
c     flag the wall and far field points at the i boundaries
c
      do j=1,jl
         fint(1,j)   = 0.
         fint(il,j)  = 0.
      end do
c
c     flag the wall and far field points at the j boundaries
c
      do i=1,il
         fint(i,jl)  = 0.
      end do

      i1        = itl  +1
      i2        = itu

      do i=i1,i2
         porj(i,1)   = 0.
      end do

      return

      end
