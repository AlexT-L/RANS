      subroutine shift(ww,ww1,wwr,xx,xxc,volc,rlvc,revc)
c
c     ******************************************************************
c     *                                                                *
c     *   exchanges data for the next mesh                             *
c     *   with data in common blocks                                   *
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
      use mg_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      real, dimension(0:iib,0:jjb,4)          :: ww
      real, dimension(iie,jje,4)              :: ww1,wwr
      real, dimension(iie,jje)                :: rlvc,revc
      real, dimension(iil,jjl,2)              :: xx
      real, dimension(0:iib,0:jjb,2)          :: xxc
      real, dimension(0:iib,0:jjb)            :: volc
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
      real     :: a
c
c     ******************************************************************
c
      do n=1,2
         do j=1,jjl
         do i=1,iil
            a         = x(i,j,n)
            x(i,j,n)  = xx(i,j,n)
            xx(i,j,n) = a
         end do
         end do
      end do

      if (mode.lt.0) return

      do n=1,2
         do j=1,jjb
         do i=1,iib
            a           = xc(i,j,n)
            xc(i,j,n)   = xxc(i,j,n)
            xxc(i,j,n)  = a
         end do
         end do
      end do

      do j=1,jjb
      do i=1,iib
         a         = vol(i,j)
         vol(i,j)  = volc(i,j)
         volc(i,j) = a
      end do
      end do

      do n=1,4
         do j=1,jjb
         do i=1,iib
            a         = w(i,j,n)
            w(i,j,n)  = ww(i,j,n)
            ww(i,j,n) = a
         end do
         end do
      end do

      do n=1,4
         do j=1,jje
         do i=1,iie
            a           = w1(i,j,n)
            w1(i,j,n)   = ww1(i,j,n)
            ww1(i,j,n)  = a
         end do
         end do
      end do

      do n=1,4
         do j=1,jje
         do i=1,iie
            a           = wr(i,j,n)
            wr(i,j,n)   = wwr(i,j,n)
            wwr(i,j,n)  = a
         end do
         end do
      end do

      if (kvis.gt.0) then
         do j=1,jje
         do i=1,iie
            a         = rlv(i,j)
            rlv(i,j)  = rlvc(i,j)
            rlvc(i,j) = a
         end do
         end do
      end if

      if (kvis.gt.1) then
         do j=1,jje
         do i=1,iie
            a         = rev(i,j)
            rev(i,j)  = revc(i,j)
            revc(i,j) = a
         end do
         end do
      end if

      return

      end
