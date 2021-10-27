      subroutine  graph(title,rm,al,re,cl,cd,cm,clv,cdv,xp,yp,cp,ares)
c
c     ******************************************************************
c     *                                                                *
c     *   plots the profile and streamwise pressure distribution       *
c     *                                                                *
c     ******************************************************************
c
      use dims
c
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
      integer  :: npgrid,npcp,npcnv
c
c     ******************************************************************
c
      real     :: rm,al,re,cl,cd,cm,clv,cdv,ares
      real     :: xp,yp,cp
c
c     ******************************************************************
c
      character(80) :: title
c
c     ******************************************************************
c
      dimension   xp(*),yp(*),cp(*)
c
c     ******************************************************************
c
c     local variables
c
c     ******************************************************************
c
      integer  :: i,imax
c
c     ******************************************************************
c
      real     :: xpt,ypt
      real     :: xmax,xmin,ymin,cpmax,scale
      real     :: cpc
c
c     ******************************************************************
c
      character(80) :: b
c
c     ******************************************************************
c
c     open a new page
c
c     open  (18,file='CP.PLOT',access='append')
      open  (18,file='CP.PLOT',position='append')

      npcp = npcp +1
      call initpl(2.,1.25,1.,npcp)
c
c     write the subtitles
c
      write (b,12)  title
   12 format(a80)
      call symbol(0.,-.25,.20,b,0.,80)
      write (b,14)  rm,al,re
   14 format('Mach ',f7.3,4x,'Alpha',f7.3,4x,'Re ',e9.3)
      call symbol(0.,-.5,.14,b,0.,44)
      write (b,16)  cl,cd,cm,clv,cdv
   16 format('CL ',f7.4,4x,'CD ',f7.4,4x,'CM ',f7.4,4x,
     .       'CLv ',f7.4,4x,'CDv ',f7.4)
      call symbol(0.,-.75,.14,b,0.,68)
c
c     calculate the scale
c
      xmax      = xp(itl)
      xmin      = xp(itl)
      ymin      = yp(itl)

      do i=itl,itu
        xmax      = max(xp(i),xmax)
        xmin      = min(xp(i),xmin)
        ymin      = min(yp(i),ymin)
      end do

      scale     = 5./(xmax  -xmin)
c
c     draw the profile
c
      xpt       = scale*(xp(itl)  -xmin)
      ypt       = scale*(yp(itl)  -ymin)  +.25
      call plot(xpt,ypt,3)

      do i=itl,itu
         xpt       = scale*(xp(i)  -xmin)
         ypt       = scale*(yp(i)  -ymin)  +.25
         call plot(xpt,ypt,2)
      end do
c
c     find the leading edge stagnation point
c
      cpmax     = 0.
      imax      = itl  +1
      do i=itl+1,itu
         if (cp(i).gt.cpmax) then
            cpmax     = cp(i)
            imax      = i
         end if
      end do
c
c     draw the vertical axis
c
      call yaxis(1.2,-2.,-.4,-.5,1.,8.,-1,-1,'Cp',2,10.,0)
      cpc       = -10.
      if (rm.gt.0.)
     .cpc       = (((5.  +rm**2)/6.)**3.5  -1.)/(.7*rm**2)
      xpt       = -.5
      ypt       = 4.  -2.5*cpc
c     if (cpc.ge.-2.0) call symbol(xpt,ypt,.36,'-',0.,1)

c
c     plot the lower surface pressures
c
      do i=itl+1,imax
         if (cp(i).ge.-2.0) then
            xpt       = scale*(.5*(xp(i)  +xp(i-1))  -xmin)
            ypt       = 4.  -2.5*cp(i)
            call symbol(xpt,ypt,.18,'+',45.,1)
         end if
      end do
c
c     plot the upper surface pressures
c
      do i=imax,itu
         if (cp(i).ge.-2.0) then
            xpt       = scale*(.5*(xp(i)  +xp(i-1))  -xmin)
            ypt       = 4.  -2.5*cp(i)
            call symbol(xpt,ypt,.18,'+',0.,1)
         end if
      end do
c
c     close the plot
c
      call endplt
      endfile 18
      close (18)

      return

      end
