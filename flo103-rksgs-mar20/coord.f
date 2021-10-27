      subroutine  coord
c
c     ******************************************************************
c     *                                                                *
c     *   coordinate stretching functions                              *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use geo_var
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
      integer  :: i,j,ile
c
c     ******************************************************************
c
      real     :: pi
      real     :: dx,dy
      real     :: a,b,c,d
      real     :: px,bp
      real     :: a2,a3
c
c     ******************************************************************
c
      ile       = il/2  +1

      pi        = 4.*atan(1.)

      xlim      = xte*boundx
      dx        = 2.*boundx/float(il  -1)

      px        = pi/xlim
      bp        = bunch/px

      a2        = 3.*ylim1  -4.*ylim2
      a3        = 2.*ylim1  -3.*ylim2

      do i=1,il

      d         = float(i  -ile)*dx
      if (abs(d).le.xlim) d = d  +bp*sin(px*d)

      if (abs(d).gt.xlim) then
         b         = 1.
         if (d.lt.0.) b = -1.
         a         = 1.  -((d  -b*xlim)/(1.  -xlim))**2
         c         = a**ax
         d         = b*xlim  +(1.  -bunch)*(d  -b*xlim)/c
      end if

      a0(i)     = d

      d         = abs(d/xlim)
      if (d.ge.1.) then
         a1(i)     = ylim2*xlim/abs(a0(i))
      else
         a1(i)     = ylim1  -d*d*(a2  -a3*d)
      end if

c     a1(i)     = 1.

      d         = a0(i)*a1(i)
      write (6,1006) i,a0(i),a1(i),d
 1006 format(1x,'i',i6,2x,'a0,a1,d',3f12.4)

      end do

      dy        = boundy/float(jl  -1)

      do j=1,jl

      d         = float(j  -1)*dy
      a         = 1.  -d*d
      c         = a**ay
      b0(j)     = sy*d/c

      end do

      return

      end
