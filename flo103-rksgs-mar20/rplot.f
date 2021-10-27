      subroutine  rplot(nres,res,fsup,count,title,rm,al,nx,ny)
c
c     ******************************************************************
c     *                                                                *
c     *   plots the convergence rate                                   *
c     *                                                                *
c     ******************************************************************
c
      implicit none
c
c     ******************************************************************
c
      common/psp/ npgrid,npcp,npcnv
c
c     ******************************************************************
c
      integer  :: nres,nx,ny
      integer  :: npgrid,npcp,npcnv
c
c     ******************************************************************
c
      real     :: rm,al
      real     :: res,fsup,count
c
c     ******************************************************************
c
      character(80) :: title
c
c     ******************************************************************
c
      dimension   res(*),fsup(*),count(*)
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: i
c
c     ******************************************************************
c
      real     :: rate,rmin,rmax,res1,fsupf
      real     :: yscal,yint,ylo,yhi
      real     :: xscal,xmax,xint
      real     :: xp,yp
c
c     ******************************************************************
c
      character(80) :: b
c
c     ******************************************************************
c
      rate      = (res(nres)/res(1))**(1./count(nres))
c
c     open a new page
c
c     open  (18,file='CNV.PLOT',access='append')
      open  (18,file='CNV.PLOT',position='append')

      npcnv = npcnv +1
      call initpl(2.,1.25,1.,npcnv)
c
c     write the subtitles
c
      write (b,12)  title
   12 format(a80)
      call symbol(0.,0.,.14,b,0.,80)
      write (b,14)  rm,al
   14 format('Mach ',f11.3,4x,'Alpha',f11.3)
      call symbol(0.,-.25,.14,b,0.,36)
      write (b,16)  res(1),res(nres)
   16 format('Resid1 ',e9.3,4x,'Resid2 ',e9.3)
      call symbol(0.,-.5,.14,b,0.,36)
      write (b,18)  count(nres),rate
   18 format('Work ',f11.2,4x,'Rate ',f11.4)
      call symbol(0.,-.75,.14,b,0.,36)
      write (b,20)  nx,ny
   20 format('Grid   ',i4,'x',i4)
      call symbol(0.,-1.,.14,b,0.,16)
c
c     calculate the vertical scale
c     depending on the magnitude of the reduction in the error
c
      rmin      = 0.
      rmax      = 0.
      res1      = res(1)
      fsupf     = fsup(nres)
      if (abs(fsupf).gt.1.e-06)  fsupf = 1./fsupf

      do i=1,nres
         fsup(i)   = fsup(i)*fsupf
         res(i)    = log(res(i)/res1)
         rmax      = max(rmax,res(i))
         rmin      = min(rmin,res(i))
      end do

      yscal     = 1./log(10.)
      yint      = 2.
      if (yscal*rmin.lt.-12.) yint = 4.
      ylo       = -6.*yint
      yhi       = 2.*yint
      yscal     = yscal/yint
c
c     calculate the horizontal scale
c     depending on the number of cycles
c
      xmax      = 300.
      xint      = 50.

      if (count(nres).gt.300.) then
         xmax = 600.
         xint = 100.
      end if

      if (count(nres).gt.600.) then
         xmax = 1200.
         xint = 200.
      end if

      if (count(nres).gt.1200.) then
         xmax = 3000.
         xint = 500.
      end if

      if (count(nres).gt.3000.) then
         xmax = 6000.
         xint = 1000.
      end if

      xscal     = 1./xint
c
c     draw the axes
c
      call xaxis(0.,xmax,xint,-.5,1.,6.,-1,-1,'Work',4,12,0)
      call yaxis(ylo,yhi,yint,-.5,1.,8.,-1,-1,'Log(error)',10,12,0)
      call yaxis(-.2,1.4,.2,5.5,1.,8.,1,1,'Nsup',4,12,0)
c
c     plot the convergence history
c
      xp        = xscal*count(1)  -.5
      yp        = 7.
      call plot(xp,yp,3)

      do i=2,nres
         xp        = xscal*count(i)  -.5
         yp        = min(2.,yscal*res(i))  +7.
         call plot(xp,yp,2)
      end do
c
c     plot the number of supersonic points
c
      xp        = xscal*count(1) - .5
      yp        = 5.*fsup(1) + 2.
      call plot(xp,yp,3)

      do  i=2,nres
         xp        = xscal*count(i) - .5
         yp        = 5.*fsup(i) + 2.
         call plot(xp,yp,2)
      end do
c
c     close the plot file
c
      call endplt
      endfile 18
      close (18)

      return

      end
