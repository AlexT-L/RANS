      subroutine wplot
c
c     ******************************************************************
c     *                                                                *
c     *   dumps files for plot3d graphic  program                      *
c     *   it makes a conversion to single precision,                   *
c     *   and strips the convex binary                                 *
c     *   for compatibility with sgi iris                              *
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
      integer  :: ix,iw
      integer  :: i,j,n
c
c     ******************************************************************
c
      real*4    xi(idm,jdm,2),wi(idn,jdn,4),conv,srm,sal,sre,stime
c
c     ******************************************************************
c
      data stime/0.0/,conv/4.0/
c
c     ******************************************************************
c
      ix        = 21
      iw        = 22

      open(unit=ix,form='unformatted')
      rewind ix

      do 10 n=1,2
         do 10 j=1,jl
         do 10 i=1,il
            xi(i,j,n) = sngl(x(i,j,n))/conv
         end do
         end do
      end do

      write (ix) il,jl,
     .           ((xi(i,j,1),i=1,il),j=1,jl),
     .           ((xi(i,j,2),i=1,il),j=1,jl)

      istat     = system('tail +5c < fort.21 > xyz.bin')

      close (ix,status='delete')

      open (unit=iw,form='unformatted')
      rewind iw

      fac       = 1./sqrt(gamma)
      do j=1,je
      do i=1,ie
         wi(i,j,1) = sngl(w(i,j,1))/conv
         wi(i,j,2) = sngl(w(i,j,2)*fac)/conv
         wi(i,j,3) = sngl(w(i,j,3)*fac)/conv
         wi(i,j,4) = sngl(w(i,j,4)*fac*fac)/conv
      end do
      end do
c
c     shift w to the cell corners
c
      do n=1,4

         do j=1,jl
         do i=1,il
            wi(i,j,n) = .5*(wi(i,j,n)+wi(i+1,j,n))
         end do
         end do

         do j=1,jl
         do i=1,il
            wi(i,j,n) = .5*(wi(i,j,n)+wi(i,j+1,n))
         end do
         end do

      end do

      srm       = sngl(rm) /conv
      sal       = sngl(al) /conv
      sre       = sngl(re) /conv
      stime     = sngl(time) /conv

      write (iw) il,jl,
     .          srm,sal,sre,stime,
     .          (((wi(i,j,n),i=1,il,1),j=1,jl,1),n=1,4,1)

      istat     = system('tail +5c < fort.22 > flo.bin')

      close (iw,status='delete')

      return

      end
