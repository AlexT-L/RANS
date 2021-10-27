      subroutine addw  (ww,ww1)
c
c     ******************************************************************
c     *                                                                *
c     *   interpolates the corrections to a finer mesh                 *
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
      use solv_var
c
c     ******************************************************************
c
      use solv_param
      use mg_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      real, dimension(0:iib,0:jjb,4)       :: ww
      real, dimension(iie,jje,4)           :: ww1
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
c     loop over the flow variables
c
      do 10 n=1,4

      if (mode.eq.0) then
c
c     transfer to a new fine mesh
c     set dw(i,j,n) equal to the existing solution
c
         jj        = 0
      do j=1,jb,2
         jj        = jj  +1
         ii        = 0
      do i=1,ib,2
         ii        = ii  +1
         dw(i,j,n) = ww(ii,jj,n)
      end do
      end do

      else
c
c     transfer to a finer mesh in the multigrid cycle
c     set dw(i,j,n) equal to the total change in the solution
c     on the coarse mesh during the cycle
c
         jj        = 0
      do j=1,jb,2
         jj        = jj  +1
         ii        = 0
      do i=1,ib,2
         ii        = ii  +1
         dw(i,j,n) = ww(ii,jj,n)  -ww1(ii,jj,n)
      end do
      end do

      end if
c
c     transfer the corrections to the finer mesh
c     using bilinear interpolation
c
      if (mode.gt.0) then
         do j=1,jb,2
            dw(1,j,n)   = 0.
            dw(ib,j,n)  = 0.
         end do
      end if

      do i=2,ie,2
      do j=1,jb,2
         dw(i,j,n)   = .25*dw(i-1,j,n)  +.75*dw(i+1,j,n)
         dw(i-1,j,n) = .75*dw(i-1,j,n)  +.25*dw(i+1,j,n)
      end do
      end do

      if (mode.gt.0) then
         do i=1,ie
            dw(i,jb,n)  = 0.
         end do
      end if

      do j=2,je,2
      do i=1,ie
         dw(i,j,n)   = .25*dw(i,j-1,n)  +.75*dw(i,j+1,n)
         dw(i,j-1,n) = .75*dw(i,j-1,n)  +.25*dw(i,j+1,n)
      end do
      end do

   10 continue

      if (mode.eq.0) then
c
c     transfer to a new fine mesh
c     set the solution equal to dw(i,j,n)
c
      do n=1,4
         do j=1,je
         do i=1,ie
            w(i,j,n)    = dw(i,j,n)
         end do
         end do
      end do

      else
c
c     transfer to a finer mesh in the multigrid cycle
c     add dw(i,j,n) to the solution
c     do not update the far field values
c
c     option to smooth the corrections
c
      if (fadd.gt.0.) then
         rfl0      = 0.
         csmoop    = 2.*fadd
         call psmoo
      end if
c
c     add the corrections to the solution
c
      do n=1,4
         do j=1,je
         do i=1,ie
            w(i,j,n)    = w(i,j,n)  +dw(i,j,n)
         end do
         end do
      end do

      end if

      return

      end
