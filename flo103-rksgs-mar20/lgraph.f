      subroutine  lgraph
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
      common/tit/ title
c
c     ******************************************************************
c
      character(80) :: title
      character(80) :: b
c
c     ******************************************************************
c
c     real x2(ie,je),y2(ie,je),qt(ie,je)
      real q1,qep
c
c     ******************************************************************
c
      integer iwatch_cp,npage
      data iwatch_cp,npage/0,0/
      save iwatch_cp,npage
c
c     ******************************************************************
c
      rad       = 45./atan(1.)
      al        = alpha*rad
      i1        = itl  +1
c
c     open a new page
c
      open  (18,file='LCP.PLOT')
      call initpl(2.,1.25,1.,1)
c
c     initialize the left window
c
      call upperport(0)
c
c     write the subtitles
c
      write (b,12)  title
   12 format(a80)
      call symbol(-1.,9.2,.18,b,-90.,80)
      write (b,14)  rm,al
   14 format('Mach ',f7.3,4x,'Alpha',f7.3)
      call symbol(-1.2,9.2,.16,b,-90.,28)
      write (b,16)  cl,cd,cm
   16 format('CL   ',f7.4,4x,'CD   ',f7.4,4x,'CM   ',f7.4)
      call symbol(-1.4,9.2,.16,b,-90.,44)
      write (b,18)  nx,ny,ncyc,rtrms
   18 format('Grid ',i5,'x',i5,4x,'Ncyc  ',i4,4x,'Res  ',e9.3)
      call symbol(-1.6,9.2,.16,b,-90.,48)

      call translate(-0.5,8.4)
      call scaling(0.68,.68)
      call rotate(-90.)
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
      imax      = i1

      do i=i1,itu
         if (cp(i).gt.cpmax) then
            cpmax     = cp(i)
            imax      = i
         end if
      end do
c
c     draw the vertical axis
c
      call yaxis(1.2,-2.,-.4,-.5,1.,8.,-1,-1,'Cp',2,40.,0)

      cpc       = -10.
      if (rm.gt.0.)
     .cpc       = (((5.  +rm**2)/6.)**3.5  -1.)/(.7*rm**2)
      xpt       = -.5
      ypt       = 4.  -2.5*cpc
c     if (cpc.ge.-2.0) call symbol(xpt,ypt,.36,'-',0.,1)

c
c     plot the lower surface pressures
c
      do i=i1,imax
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
c     initialize the right window
c
      call upperport(1)

      call lowerport(0)

      call translate(0.,3.)
      call scaling(.7,.7)
      call rotate(-90.)
c
c     draw the profile
c
      xpt       = scale*(xp(itl)  -xmin)
      ypt       = scale*(yp(itl)  -ymin)  +.5
      call plot(xpt,ypt,3)

      do i=itl,itu
         xpt       = scale*(xp(i)  -xmin)
         ypt       = scale*(yp(i)  -ymin)  +.5
         call plot(xpt,ypt,2)
      end do
c
c     calculate the local mach number
c
      qmin      = 1.e16
      qmax      = -qmin

      do i=1,il
        do j=1,jl
          xq(i,j) = scale*(scal*x(i,j,1)  -xmin)
          yq(i,j) = scale*(scal*x(i,j,2)  -ymin)  +.5
          pa      = .25*(p(i,j) +p(i+1,j) +p(i,j+1) +p(i+1,j+1))
          ra      = .25*(w(i,j,1) +w(i+1,j,1) +w(i,j+1,1) +w(i+1,j+1,1))
          rua     = .25*(w(i,j,2) +w(i+1,j,2) +w(i,j+1,2) +w(i+1,j+1,2))
          rva     = .25*(w(i,j,3) +w(i+1,j,3) +w(i,j+1,3) +w(i+1,j+1,3))
          qm      = sqrt((rua*rua +rva*rva)/(gamma*pa*ra))
          qmin    = min(qmin,qm)
          qmax    = max(qmax,qm)
          qc(i,j) = qm
        end do
      end do
c
c     define intervals for the mach contours
c
      nct       = 40
      hinc      = 2./3./(nct-1)
      qep       = 0.0001*(qmax -qmin)/real(nct-1)
c
c     plot the mach contours
c
      do nc = 1, nct
c        huet    = hinc*(nct-nc)
c        call pshue(huet)
         q1      = qmin + (qmax -qmin)*real(nc-1)/real(nct-1)
c        call pscont(x2,y2,qt,ie,je,1,il,1,jl,q1,qep)
         call pscont(ie,je,1,il,1,jl,q1,qep)
      end do

      call lowerport(1)
c
c     close the plot
c
      call endplt
      endfile 18
      close (18)
c
c     launch or refresh ghostview
c
      npage     = npage  +1
      if (iwatch_cp.eq.0) then
c        call system("gv -watch -seascape LCP.PLOT &")
c        call system("gv -watch -safer -quiet LCP.PLOT &")
c        call system("gv -watch -orientation=seascape LCP.PLOT &")
         call system("gv --orientation=seascape --watch LCP.PLOT &")
         iwatch_cp = 1
         call system("cp -f LCP.PLOT HCP.PLOT")
      else
c        call system("killall -hup ghostview")
c        open(180,file='HCP.PLOT',access='append')
         open(180,file='HCP.PLOT',position='append')
         write(180,'(a,2i4)') '%%page:',npage,npage
         close(180)
         call system("tail +7 LCP.PLOT >> HCP.PLOT")
c        call system("tail --lines=+7 LCP.PLOT >> HCP.PLOT")
      end if

      return

      end
