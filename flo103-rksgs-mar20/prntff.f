      subroutine prntff(ncyc,time,lprnt)
c
c     ******************************************************************
c     *                                                                *
c     *   prints the flow field                                        *
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
      integer  :: ncyc,lprnt
c
c     ******************************************************************
c
      real     :: time
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: iwrit
      integer  :: i,j,jj
c
c     ******************************************************************
c
      real     :: u,v,qq,h,fm,s
      real     :: c1,hfac
c
c     ******************************************************************
c
      iwrit     = 6

      write (iwrit,10)
   10 format(1x,'flow field'/
     .       3x,' ncyc ','    time    ')

      write (iwrit,12) ncyc,time
   12 format(1x,i6,f12.4)

      do i=2,ie,lprnt

      write (iwrit,22)
   22 format(3x,'  i ','  j ',
     .          '   rho   ','    u    ','    v    ',
     .          '    h    ','    p    ','   mach  ',
     .          '    s    ',' eddy vis')

      do j=2,je

      u         = w(i,j,2)/w(i,j,1)
      v         = w(i,j,3)/w(i,j,1)
      qq        = u*u  +v*v
      h         = (w(i,j,4)  +p(i,j))/w(i,j,1)  -h0
      fm        = sqrt((w(i,j,1)*qq)/(gamma*p(i,j)))
      s         = (p(i,j)/p0)*(rho0/w(i,j,1))**gamma  -1.
      write (iwrit,24) i,j,w(i,j,1),u,v,h,p(i,j),fm,s,rev(i,j)
   24 format(1x,2i4,8f10.4)

      end do
      end do

      jj        = je  -2
      do j=2,je,jj

      write (iwrit,22)

      do i=1,ie

      u         = w(i,j,2)/w(i,j,1)
      v         = w(i,j,3)/w(i,j,1)
      qq        = u*u  +v*v
      h         = (w(i,j,4)  +p(i,j))/w(i,j,1)  -h0
      fm        = sqrt((w(i,j,1)*qq)/(gamma*p(i,j)))
      s         = (p(i,j)/p0)*(rho0/w(i,j,1))**gamma  -1.
      write (iwrit,24) i,j,w(i,j,1),u,v,h,p(i,j),fm,s,rev(i,j)

      end do
      end do

      if (kvis.eq.0) return
c
c     turbulent parameters
c
      c1        = 2./(sqrt(gamma) *rm * re)
      j         = 2

      write (iwrit,42)
   42 format(3x,'  i ','  j ',
     .          '   tauw  ','  dstar  ','  delta  ',
     .          '  theta  ','  h-fac  ','    cf   ')

      do i=itl+1,itu

      cfric(i)  = c1*tw(i)
      if (ssmax(i).gt.0.) hfac = dsti(i)/ssmax(i)
      write(iwrit,44) i,j,tw(i),dsti(i),ynot(i),ssmax(i),hfac,cfric(i)
   44 format(1x,2i4,8e10.3)

      end do

      return

      end
