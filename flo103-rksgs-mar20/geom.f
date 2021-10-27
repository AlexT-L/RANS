      subroutine  geom
c
c     ******************************************************************
c     *                                                                *
c     *   square root  mapping of the profile to a slit                *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
c     ******************************************************************
c
      use in_var
      use geo_var
      use mesh_var
c
c     ******************************************************************
c
      use flo_param
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
      integer  :: i,i1,i2
      integer  :: ind
c
c     ******************************************************************
c
      real     :: pi,angl1,angl2,t1,t2
      real     :: xa,ya,u,v,r,angl
      real     :: x0,y0,x1,ang,a,b,d
c
c     ******************************************************************
c
      real, dimension(nn)      :: xs,ys,d1,d2,d3
c
c     ******************************************************************
c
      pi        = 4.*atan(1.)

      scal      = .50001*xlim**2/(xn(nn)  -xsing)
      angl      = atan(slopt)
      angl1     = atan2((yn(1)  -ysing),(xn(1)  -xsing))
      angl2     = atan2((yn(nn)  -ysing),(xn(nn)  -xsing))
      angl1     = angl  -.5*(angl1  -trail)
      angl2     = angl  -.5*(angl2  +trail)
      t1        = tan(angl1)
      t2        = tan(angl2)

      angl      = pi  +pi
      u         = 1.
      v         = 0.

      do i=1,nn
         xa        = xn(i)  -xsing
         ya        = yn(i)  -ysing
         angl      = angl  +atan2((u*ya  -v*xa),(u*xa  +v*ya))
         r         = scal*sqrt(xa**2  +ya**2)
         u         = xa
         v         = ya
         r         = sqrt(2.*r)
         xs(i)     = r*cos(.5*angl)
         ys(i)     = r*sin(.5*angl)
      end do

      scal      = 1./scal

      call splif (1,nn,xs,ys,d1,d2,d3,1,t1,1,t2,0,0.,ind)
      call intpl (itl,itu,a0,s0,1,nn,xs,ys,d1,d2,d3,0)

      x1        = .125*(xlim)**2
      angl      = pi  +pi
      u         = 1.
      v         = 0.
      i1        = 1
      i2        = itl  -1
      i         = itl

   11 r         = sqrt(a0(i)**2  +s0(i)**2)
      ang       = 2.*atan2(s0(i),a0(i))
      r         = .5*r**2
      x0        = r*cos(ang)
      y0        = r*sin(ang)
      a         = slopt*(x0  -x1)
      b         = 1./(x0  -x1)

      do i=i1,i2
         xa        = .5*a0(i)**2
         d         = b*(xa  -x1)
         ya        = y0  +a*log(d)/d
         angl      = angl  +atan2((u*ya  -v*xa),(u*xa  +v*ya))
         r         = sqrt(xa**2  +ya**2)
         u         = xa
         v         = ya
         r         = sqrt(2.*r)
         s0(i)     = r*sin(.5*angl)
      end do

      if (i2.eq.il) return

      angl      = 0.
      u         = 1.
      v         = 0.
      i1        = itu  +1
      i2        = il
      i         = itu

      go to 11

      end
