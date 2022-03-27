      subroutine delt(ny,il,jl,ib,jb,w,ynot,dsti,x,xc)
c
c     *****************************************************************
c     *                                                               *
c     *   calculates the boundary layer thickness                     *
c     *                                                               *
c     *****************************************************************
c
c     use dims
c
c     ******************************************************************
c
c     use flo_var
c     use mesh_var
c     use solv_var
c
c     ******************************************************************
c
c     use flo_param
c     use solv_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
c     input variables
c
c     ******************************************************************
c     from dims
      integer, intent(in) :: ny,il,jl,ib,jb

c     from flo_var
      real(8), intent(in), dimension(0:ib,0:jb,4) :: w
      real(8), intent(inout), dimension(0:ib) :: ynot,dsti

c     from mesh_var
      real(8), intent(in), dimension(1:il,1:jl,2) :: x
      real(8), intent(inout), dimension(0:ib,0:jb,2) :: xc
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
c     ******************************************************************
c
      integer :: i,j,js,jmax,lend,lbig,jse
      real(8) :: xy,yy,si,dsi,cdu,fx,uinf,ra,fc,xbi,ybi,ycorr
      real(8), dimension(0:jb)               :: qs,dn,ut
      real(8), dimension(0:ib)               :: ssmax
      integer, external :: idmax
c
c     ******************************************************************
c
      logical locke
c
c     ******************************************************************
c
c     js        = 2*jl/3  -2
      js        = .75*(ny  -4)

      do 20 i=2,il

      qs(1)     = 0.
      ut(1)     = 0.
      j         = js
      xy        = .5*(x(i,j,1)  -x(i,j-1,1)
     .               +x(i-1,j,1)  -x(i-1,j-1,1))
      yy        = .5*(x(i,j,2)  -x(i,j-1,2)
     .               +x(i-1,j,2)  -x(i-1,j-1,2))
      qs(j)     = (yy*w(i,j,2)  -xy*w(i,j,3))/(w(i,j,1))
      si        = sign(real(1.,8),qs(js))

      do j=2,js
         xy        = .5*(x(i,j,1)  -x(i,j-1,1)
     .                  +x(i-1,j,1)  -x(i-1,j-1,1))
         yy        = .5*(x(i,j,2)  -x(i,j-1,2)
     .                  +x(i-1,j,2)  -x(i-1,j-1,2))
         dsi       = 1./sqrt(xy**2  +yy**2)
         qs(j)     = si*(yy*w(i,j,2)  -xy*w(i,j,3))
         dn(j)     = 1./dsi
         ut(j)     = qs(j)*dsi
      end do

      dsti(i)   = 0.
      ynot(i)   = 0.
      ssmax(i)  = 0.
      cdu       = .98
      jmax      = idmax(js,ut,1)
      fx        = .6*ut(jmax)
      lend      = 2
      lbig      = 2
      locke     = .false.

      do  j=3,js
         if ( ut(j-1) .lt. 0. .and. ut(j) .ge. .0) then
            lbig      = j
         end if
         if (.not. locke) then
            if ( ut(j-1) .ge. cdu*ut(j) .and. ut(j) .gt. fx) then
               locke     = .true.
               lend      = j
            end if
         end if
      end do

      uinf      = 1./ut(lend)
      do j=lbig,lend
         dsti (i)  = dsti(i) + (ut(lend)*dn(j) - qs(j))
         ra        = w(i,lend,1)/w(i,j,1)
         ssmax(i)  = ssmax(i) + ra*ut(j)*uinf*(dn(j)-qs(j)*uinf)
      end do

      if (.true.) then

      dsti(i)   = dsti(i)*uinf
      dsti (i)  = max(dsti(i),1.e-6)
      ra        = w(i,lend,1)/w(i,2,1)
      ssmax(i)  = max(ssmax(i),ra*qs(2)*uinf)
      fc        = .95*ut(lend)

      jse       = lend
      do j=2,jse
         lend      = j  -1
         if (ut(j).gt.fc) go to 11
      end do

   11 xbi       = .5*(x(i,1,1)+x(i-1,1,1))
      ybi       = .5*(x(i,1,2)+x(i-1,1,2))
      ycorr     = sqrt((xc(i,lend,1) - xbi)**2+(xc(i,lend,2)-ybi)**2)
      ynot(i)   = 1.5*(ycorr  +dn(lend)*
     .                (fc  -ut(lend))/(ut(lend+1)  -ut(lend)))

      end if
   20 continue

      return

      end subroutine


      integer function idmax(n,sx,incx)
c
c     ******************************************************************
c     *                                                                *
c     *   find the index of element having max value                   *
c     *                                                                *
c     ******************************************************************
c
      integer :: i,incx,ix,n
c
c     *****************************************************************
c
      real(8) :: smax
      real(8), dimension(0:n-1) :: sx
c
c     *****************************************************************
c
      idmax     = 0
      if (n.lt.1) return

      idmax     = 1
      if (n.eq.1) return

      if (incx.eq.1) go to 11
c
c      code for increment not equal to 1
c
      ix        = 1
      smax      = sx(1)
      ix        = ix + incx
      do i=2,n
         if (sx(ix).gt.smax) then
            idmax     = i
            smax      = sx(ix)
         end if
         ix        = ix + incx
      end do

      return
c
c      code for increment equal to 1
c
   11 smax = sx(1)
      do i=2,n
         if (sx(i).gt.smax)then
            idmax     = i
            smax      = sx(i)
         end if
      end do

      return

      end function
