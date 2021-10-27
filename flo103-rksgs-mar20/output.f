      subroutine output
c
c     ******************************************************************
c     *                                                                *
c     *   outputs aerodynamic and force data                           *
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
      use out_var
c
c     ******************************************************************
c
      use flo_param
      use solv_param
c
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      common/tit/ title
c
c     ******************************************************************
c
      character(80) :: title
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: iwrit,iout
      integer  :: i,j,l,n
      integer  :: i1,i2,j1,j2
c
c     ******************************************************************
c
      real     :: rad,al,gmm
c
c     ******************************************************************
c
      iwrit     = 6

      rad       = 45./atan(1.)
      al        = alpha*rad
      gmm       = .5*gamma*rm**2
c
c     print the flow field
c
      if (iprnt.gt.0)
     .call prntff(ncyc,0.,lprnt)

      write (iwrit,600)
      write (iwrit,12)
   12 format(1x,'Section characteristics'/
     .       3x,' Mach number',' Ang attack ',' Reynolds no',
     .          '     CL     ','     CD     ','     CM     ')
      write (iwrit,610) rm,al,re,cl,cd,cm

      write (iwrit,14)
   14 format(1x,'Section characteristics'/
     .       3x,'     CLv    ','     CDv    ')
      write (iwrit,611) clv,cdv
c
c     print-plot the pressure distribution on the standard output
c
c     if (iprnt.gt.0)
c    .call cplot
c
c     plot the mesh
c
      call grid(title)
         i1        = 320
         i2        = 324
         j1        = 1
         j2        =  38
      call zoomgrid(title,i1,i2,j1,j2)

      do i=1,il
         xp(i)     = scal*x(i,1,1)
         yp(i)     = scal*x(i,1,2)
      end do
c
c     print the skin friction coefficient
c
      do i=itl+1,itu
         write (6,1006) i,xp(i),yp(i),cp(i),cf(i)
      end do

 1006 format(1x,'i',i6,2x,'xp,yp,Cp,Cf',4f12.4)
c
c     plot the pressure distribution
c
      call graph (title,rm,al,re,cl,cd,cm,clv,cdv,xp,yp,cp,rtrms)
      call lgraph

      return

  600 format(1x)
  610 format(1x,2f12.4,f12.2,3f12.4)
  611 format(1x,2f12.4)

      end
