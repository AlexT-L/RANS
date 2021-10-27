c
c===============================================================================
c
      subroutine initpl(xinch,yinch,scale,npage)
c
c===============================================================================
c
c-----print prologue for postscript file
c
c=======================================
c
      common /mic/ nen
c
c=======================================
c
      nen = 0
      if (npage.eq.1) then
      write(18,1001) 
      write(18,1002)
      write(18,1003)
      write(18,1004)
      write(18,1005)
      end if
      write(18,1011) npage,npage
      write(18,1012) 
      write(18,1013) 
      write(18,1014) 
      write(18,1015) 
      write(18,1016) 
      write(18,1017) 
      write(18,1018) xinch,yinch,scale,scale
      write(18,1020)
      write(18,1021)
      write(18,1022)
      write(18,1023)
      write(18,1024)
      write(18,1025)
 1001 format('%!PS-Adobe-2.0')
 1002 format('% definition operators')
 1003 format('/bdef {bind def} bind def')
 1004 format('% page state control')
 1005 format('/bpage {/pgsv save def} bdef')
 1011 format('%%Page: ',i3,1x,i3)
 1012 format('/inch {72 mul} def')
 1013 format('/fni /Helvetica-Oblique findfont 8 scalefont def')
 1014 format('/fnb /Helvetica-Oblique findfont 12 scalefont def')
 1015 format('fni setfont')
 1016 format('0 setlinewidth')
 1017 format('/shwb { moveto fnb setfont show} def')
 1018 format(2((f5.2),' inch'), ' translate ',
     .       2((f5.2),' inch'), ' scale')
 1020 format('/upperport {newpath',
     .       '            -2     4.26 moveto  6.47  4.26 lineto',
     .       '             6.47  9.73 lineto -2     9.73 lineto',
     .       '            closepath } def')
 1021 format('/lowerport {newpath',
     .       '            -2    -1.24 moveto  6.47 -1.24 lineto',
     .       '             6.47  4.24 lineto -2     4.24 lineto',
     .       '            closepath } def')
 1022 format('/lowerleft {newpath',
     .       '            -2    -1.24 moveto  2.23 -1.24 lineto',
     .       '             2.23  4.24 lineto -2     4.24 lineto',
     .       '            closepath } def')
 1023 format('/lowerright {newpath',
     .       '             2.24 -1.24 moveto  6.47 -1.24 lineto',
     .       '             6.47  4.24 lineto  2.24  4.24 lineto',
     .       '            closepath } def')
 1024 format('/upperleft {newpath',
     .       '            -2     4.26 moveto  2.23  4.26 lineto',
     .       '             2.23  9.73 lineto -2     9.73 lineto',
     .       '            closepath } def')
 1025 format('/upperright {newpath',
     .       '             2.24  4.26 moveto  6.47  4.26 lineto',
     .       '             6.47  9.73 lineto  2.24  9.73 lineto',
     .       '            closepath } def')
      return
      end
c
c===============================================================================
c
      subroutine plot(x,y,ipen)
c
c===============================================================================
c
c-----Postscript equivalent
c
c=======================================
c
      common/mic/ nen
c
c=======================================
c
      if(ipen .eq. 3) then
      nen = nen + 1
      write(18,1) x,y
      return
      end if
      if(ipen .eq. 2)  then
      nen = nen + 1
      write(18,2) x,y
      if(nen.ge.500) then
      write(18,3) 
      write(18,1) x,y
      nen = 0
      end if
      end if
      return
    1 format(2f10.4,'  moveto')
    2 format(2f10.4,'  lineto')
    3 format('stroke')
      end
c
c===============================================================================
c
      subroutine symbol(xo,yo,chasiz,flabel,rot,nchar)
c
c===============================================================================
c
c-----This subroutine writes a text string about the current position
c
c-----List of external parameters
c 
c     flabel  : String to be written
c     chasiz  : Size of characters in string (ie. font size)
c     nchar   : Number of characters in string
c     rot     : Angle of rotation from horizontal axis in degrees
c               about the current position
c
c=======================================
c
c-----Postscript equivalent
c
c=======================================
c
      character(80) flabel
      data sscale/.5/
c
c=======================================
c
      ang      = acos (-1.0)/180. * rot
      xdir     = cos (ang)
      ydir     = sin (ang)
      cent     = - chasiz *nchar/2
      centx    = xdir * cent
      centy    = ydir * cent
c
c     call text (flabel)
c
      icentr = 0
      if ( nchar .eq. 1 .or. ichar(flabel(1:1)) .eq. 92) then
      icentr = 1
      end if
      scrpsz   = chasiz
      rrot     = -rot
      rcentx   = -centx
      rcenty   = -centy
      if (ichar(flabel(1:1)) .eq. 92 ) then
      write(18,2940) scrpsz
      else
      write(18,2930) scrpsz
      end if
      if (icentr .eq. 1) then
      write(18,2931) xo,yo      
      write(18,2932) rot
      write(18,2934) flabel(1:nchar)
      write(18,2933) flabel(1:nchar)
      write(18,2932) rrot
      else   
      write(18,2931) xo,yo      
      write(18,2932) rot
      write(18,2933) flabel(1:nchar)
      write(18,2932) rrot
      end if
 2930 format('/Times-Roman findfont ',f10.5,'  scalefont setfont')
 2940 format('/Symbol findfont ',f10.5,'  scalefont setfont')
 2931 format(2f10.5,' moveto')
 2932 format(f10.5,' rotate')
 2933 format('(', a, ') show')
 2934 format('(', a, ') stringwidth pop -.5 mul dup rmoveto ')
      return
      end
c
c===============================================================================
c
      subroutine xaxis (xmin,xmax,xint,xorig,yorig,xlong,
     .                  itic,ilabel,xtitle,ntitle,size,iend)
c
c===============================================================================
c
c-----This subroutine draws a horizontal axis and returns a scale factor
c-----for plotting variables on this axis.
c-----Absolute size of tics,labels and titles is set by size
c-----(Normalized NDC coord)
c-----Relative size of tics,labels and titles are set by internal parameters
c-----which may be adjusted
c
c-----List of external parameters
c
c     xmin,  xmax,  xint  : min max and intervals of x along axis
c     xorig, yorig, xlong : origin and length of x axis in window coordinates
c     itic = -1  ----->     tics on negative y side of axis
c     itic =  0  ----->     no tics
c     itic =  1  ----->     tics on positive y side of axis
c     itic =  2  ----->     tics on both sides of axis
c     ilabel=-1  ----->     x-coord labels on negative y side of axis
c     ilabel= 0  ----->     no x-coord labels
c     ilabel= 1  ----->     x-coord labels on positive y side of axis
c     ilabel=-+2 ----->     same as above except x-coord labels are integers
c     xtitle              : Title for x axis
c     ntitle              : Number of characters in Title (Omit title if =0)
c     size                : Size of tics,labels and title in NDC coordinates
c     iend = -2  ----->     omit all x-coord labels except for at ends of axis
c     iend = -1  ----->     omit x-coord label at axis origin
c     iend =  0  ----->     write all x-coord labels
c     iend =  1  ----->     omit x-coord label at end of axis
c     iend =  2  ----->     omit x-coord label at both ends of axis
c
c=======================================
c
c-----Postscript equivalent
c
c=======================================
c
      common/xax/ xax1,xax2,xax3,xax4
      character(80) flabel,xtitle
c
c=======================================
c
      xax1     = xmin
      xax2     = xmax
      xax3     = xlong
      xax4     = xorig
c
c-----List of axis internal parameters
c     (dist is distance between labels and axis)
c     (relsiz is relative size of title characters versus label characters)
c
      siz      = size/160.
c     call mapwld(0., 0.,a1,b1)
c     call mapwld(siz,0.,a2,b2)
c     ticsiz   = a2-a1
      ticsiz   = 8./72.
      chasiz   = 1.25*ticsiz
      dist     = 3.5*ticsiz
      relsiz   = 1.5
      if (ilabel .lt.  0) dist = -dist 
      if (dist/ticsiz .lt. 0) dist = dist - abs(dist)/dist *ticsiz
c
c-----Draw the axis
c
      xend     = xorig + xlong
      call moveab (xorig,yorig)
      call lineab (xend, yorig)
c
c-----Prepare parameters for labeling the axis
c
      itik     = itic + 2
      itik     = min(itik,3)
      go to (11,30,12) itik
   11 tic      = -ticsiz
      go to 15
   12 tic      =  ticsiz
   15 xpos     = xorig
c
      xscal    = xlong / (xmax - xmin)
      dx       = xint * xscal
      n        = (xlong/dx) + 1.01
c
c-----Loop to tic and Label the axis
c
      do 20 i=1,n
c
c-----Draw tics
c
      call moveab (xpos, yorig)
      call lineab (xpos, yorig+tic)
      if (itic .eq. 2) call lineab (xpos,yorig-tic)
c
c-----Move to center location of label and write x-value using mlabel
c
      if (i .eq. 1  .and.  iend .eq.-1) go to 19
      if (i .eq. 1  .and.  iend .eq. 2) go to 19
      if (i .eq. n  .and.  iend .gt. 0) go to 19
      if (i .gt. 1  .and.  i .lt. n  .and. iend .eq. -2) go to 19
      if (ilabel .eq. 0) go to 19
c
      call moveab (xpos, yorig+dist)
      if (abs(ilabel) .eq. 2) then
          ixlabl   =  xmin + (i-1)*xint
          write(flabel,21) ixlabl
   21     format(i5)
          nchar    = 5
      else
          xlabl   =  xmin + (i-1)*xint
          write(flabel,22) xlabl
   22     format(f7.2)
          nchar    = 7
      end if
      call mlabel(flabel,chasiz,nchar,0.,1)
c
   19 xpos     = xpos + dx
   20 continue
c
   30 continue
c
c-----Write title of axis
c
      if (ntitle .eq. 0) go to 40
      xmid    = (xorig + xend)/2.0
      csiz    = chasiz * relsiz
      call moveab(xmid, yorig + (1+relsiz)*dist)
      call mlabel(xtitle,csiz,ntitle,0.,1)
   40 return
      end
c
c===============================================================================
c
      subroutine yaxis (ymin,ymax,yint,xorig,yorig,ylong,
     .                  itic,ilabel,ytitle,ntitle,size,iend)
c
c===============================================================================
c
c-----This subroutine draws a horizontal axis and returns a scale factor
c-----for plotting variables on this axis.
c-----Absolute size of tics,labels and titles is set by size
c-----(Normalized NDC coord)
c-----Relative size of tics, labels and titles are set by internal parameters
c-----which may be adjusted
c
c-----List of external parameters
c
c     ymin,  ymax,  yint  : min max and intervals of y along axis
c     xorig, yorig, ylong : origin and length of y axis in window coordinates
c     itic = -1  ----->     tics on negative x side of axis
c     itic =  0  ----->     no tics
c     itic =  1  ----->     tics on positive x side of axis
c     itic =  2  ----->     tics on both sides of axis
c     ilabel=-1  ----->     y-coord labels on negative x side of axis
c     ilabel= 0  ----->     no y-coord labels
c     ilabel= 1  ----->     y-coord labels on positive x side of axis
c     ilabel=-+2 ----->     same as above except x-coord labels are integers
c     ytitle              : Title for y axis
c     ntitle              : Number of characters in Title (omit title if =0)
c     size                : Size of tics,labels and title in NDC coordinates
c     iend = -2  ----->     omit all x-coord labels except for at ends of axis
c     iend = -1  ----->     omit x-coord label at axis origin
c     iend =  0  ----->     write all x-coord labels
c     iend =  1  ----->     omit x-coord label at end of axis
c     iend =  2  ----->     omit x-coord label at both ends of axis
c     kdigit                flag to control format of y axis labels
c=======================================
c
c-----Postscript equivalent
c
c=======================================
c
      common/yax/ yax1,yax2,yax3,yax4
      character(80) flabel,ytitle
c
c=======================================
c
      yax1     = ymin
      yax2     = ymax
      yax3     = ylong
      yax4     = yorig
c
c-----List of axis internal parameters
c     (dist is distance between labels and axis)
c     (relsiz is the relative size of title characters versus label characters)
c
      siz      = size/160.
c     call mapwld(0., 0.,a1,b1)
c     call mapwld(0.,siz,a2,b2)
c     ticsiz   = b2-b1
      ticsiz   = 8./72.
      chasiz   = 1.25*ticsiz
      dist     = 3.5*ticsiz
      relsiz   = 1.5
      if (ilabel .lt. 0) dist = -dist
      if (dist/ticsiz .lt. 0) dist = dist - abs(dist)/dist *ticsiz
c
c-----Draw the axis
c
      yend     = yorig + ylong
      call moveab (xorig,yorig)
      call lineab (xorig,yend )
c
c-----Prepare parameters for labeling the axis
c
      itik     = itic + 2
      itik     = min(itik,3)
      go to (11,36,12) itik
   11 tic      = -ticsiz
      go to 15
   12 tic      =  ticsiz
   15 ypos     = yorig
c
      yscal    = ylong / (ymax - ymin)
      dy       = yint * yscal
      n        = (ylong/dy) + 1.01
c
c-----Loop to tic and Label the axis
c
      do 34 i=1,n
c
c-----Draw tics
c
      call moveab (xorig    , ypos)
      call lineab (xorig+tic, ypos)
      if (itic .eq. 2) call lineab(xorig-tic,ypos)
c
c-----Move to center location of label and write y-value using mlabel
c
      if (i .eq. 1  .and.  iend .eq.-1) go to 32
      if (i .eq. 1  .and.  iend .eq. 2) go to 32
      if (i .eq. n  .and.  iend .gt. 0) go to 32
      if (i .gt. 1  .and.  i .lt. n  .and. iend .eq. -2) go to 32
      if (ilabel .eq. 0) go to 32
c
      call moveab (xorig+dist, ypos)
      if (abs(ilabel) .eq. 2) then
          iylabl   = ymin + (i-1)*yint
          write(flabel,21) iylabl
   21     format(i5)
          nchar    = 5
      else
          ylabel   = ymin + (i-1)*yint
          write(flabel,23) ylabel
c  23     format(e7.1)
c  23     format(f6.1)
c  23     format(f6.2)
   23     format(f7.3)
          nchar    = 7
      end if
      call mlabel(flabel,chasiz,nchar,90.,1)
c
   32 ypos     = ypos + dy
   34 continue
c
   36 continue
c
c-----Write title of axis
c
      if (ntitle .eq. 0) go to 40
      ymid    = (yorig + yend)/2.0
      csiz    = chasiz * relsiz
      call moveab(xorig + (1+relsiz)*dist, ymid)
      call mlabel(ytitle,csiz,ntitle,90.,1)
   40 return
      end
c
c===============================================================================
c
      subroutine mlabel(flabel,chasiz,nchar,rot,icentr)
c
c===============================================================================
c
c-----This subroutine writes a text string about the current position
c
c-----List of external parameters
c 
c     flabel  : String to be written
c     chasiz  : Size of characters in string (ie. font size)
c     nchar   : Number of characters in string
c     rot     : Angle of rotation from horizontal axis in degrees
c               about the current position
c     icentr  : Center string about current position if icentr=1
c               otherwise write string beginning at current position
c
c=======================================
c
c-----Postscript equivalent
c
c=======================================
c
      character(80) flabel
      data sscale/.5/
c
c=======================================
c
      ang      = acos (-1.0)/180. * rot
      xdir     = cos (ang)
      ydir     = sin (ang)
      cent     = - chasiz *nchar/2
      centx    = xdir * cent
      centy    = ydir * cent
c
c     call text (flabel)
c
      scrpsz   = chasiz
      rrot     = -rot
      rcentx   = -centx
      rcenty   = -centy
      write(18,2930) scrpsz
c.jcv if (icentr .eq. 1) then
c.jcv write(18,2931) centx,centy
c.jcv end if
      write(18,2932) rot
      if(icentr .eq. 1 ) write(18,2934) flabel(1:nchar)
      write(18,2933) flabel(1:nchar)
      write(18,2932) rrot
      if (icentr .eq. 1) then
      write(18,2931) rcentx,rcenty
      end if
 2930 format('/Times-Roman findfont ' ,f10.5,'  scalefont setfont' )
 2931 format(2f10.5, ' rmoveto' )
 2932 format(f10.5, ' rotate' )
 2933 format('(',a,') show' )
 2934 format('(',a,') stringwidth pop -.5 mul 0. rmoveto' )
      return
      end
c
c===============================================================================
c
      subroutine moveab(x,y)
c
c===============================================================================
c
      write(18,2921) x,y
 2921 format(2f10.5,' moveto')
      return
      end
c
c===============================================================================
c
      subroutine lineab(x,y)
c
c===============================================================================
c
      common /strc/ ncst
c
c=======================================
c
      write(18,2922) x,y
      if (ncst .lt. 100) return
      ncst = 0
      write(18,2923)
 2922 format(2f10.5,' lineto')
 2923 format('stroke')
      return
      end
c
c===============================================================================
c
      subroutine endplt
c
c===============================================================================
c
      write(18,10)
      write(18,11)
   10 format('stroke')
   11 format('showpage')
      return
      end

      subroutine upperport(n)
      if (n.eq.0) then
        write(18,10)
   10   format('gsave'/'upperport gsave stroke grestore clip')
      else
        write(18,20)
   20   format('stroke'/'grestore')
      end if
      return
      end

      subroutine lowerport(n)
      if (n.eq.0) then
        write(18,10)
   10   format('gsave'/'lowerport gsave stroke grestore clip')
      else
        write(18,20)
   20   format('stroke'/'grestore')
      end if
      return
      end

      subroutine lowerleft(n)
      if (n.eq.0) then
        write(18,10)
   10   format('gsave'/'lowerleft gsave stroke grestore clip')
      else
        write(18,20)
   20   format('stroke'/'grestore')
      end if
      return
      end

      subroutine lowerright(n)
      if (n.eq.0) then
        write(18,10)
   10   format('gsave'/'lowerright gsave stroke grestore clip')
      else
        write(18,20)
   20   format('stroke'/'grestore')
      end if
      return
      end

      subroutine upperleft(n)
      if (n.eq.0) then
        write(18,10)
   10   format('gsave'/'upperleft gsave stroke grestore clip')
      else
        write(18,20)
   20   format('stroke'/'grestore')
      end if
      return
      end

      subroutine upperright(n)
      if (n.eq.0) then
        write(18,10)
   10   format('gsave'/'upperright gsave stroke grestore clip')
      else
        write(18,20)
   20   format('stroke'/'grestore')
      end if
      return
      end

      subroutine rotate(ang)
      write(18,20) ang
   20   format(f12.4,' rotate')
      return
      end

      subroutine translate(x,y)
      write(18,20) x,y
   20   format(2f12.4,' translate')
      return
      end

      subroutine scaling(x,y)
      write(18,20) x,y
   20   format(2f12.4,' scale')
      return
      end

      subroutine stroke
      write(18,20)
   20   format('stroke')
      return
      end

      subroutine setdash(n,dsh)
      real dsh(n)
      write(18,20) (dsh(i),i=1,n)
   20   format('[',f12.4,'] 0 setdash')
      return
      end

      subroutine setsolid
      write(18,20)
   20   format('[] 0 setdash')
      return
      end

      subroutine setrgbcolor(r,g,b)
      write(18,20) r,g,b
   20   format(3f12.4,' setrgbcolor')
      return
      end
