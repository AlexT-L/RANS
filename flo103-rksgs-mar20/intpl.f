      subroutine  intpl(mi,ni,si,fi,m,n,s,f,fp,fpp,fppp,mode)
c
c     ******************************************************************
c     *                                                                *
c     *   interpolation of spline using local taylor series            *
c     *                                                                *
c     ******************************************************************
c
c     interpolation of cubic spline f(i) over s(i) from i = m to n
c     on to si(i) from i = mi to ni.
c     m and mi may be larger or smaller than n and ni,
c     while s(i) and si(i) may be increasing or decreasing,
c     as long as they continue in the same direction.
c     adds correction for piecewise constant fourth deriviative
c     if mode greater than 0.
c
c     ******************************************************************
c
c     programmed by antony jameson july 1970 for the ibm 1130 computer.
c     revised december 2001 for fortran 95 compliance.
c
c     ******************************************************************
c
      dimension   si(*),fi(*),s(*),f(*),fp(*),fpp(*),fppp(*)
c
c     ******************************************************************
c
      k         = sign(1,(n  -m))

      i         = m
      min       = mi
      nin       = ni

      d         = s(n)  -s(m)
      if (d*(si(ni)  -si(mi)).lt.0.) then
         min       = ni
         nin       = mi
      end if

      ki        = 0
      if (min.ne.nin) ki = sign(1,(nin  -min))

      ii        = min  -ki
      c         = 0.
      if (mode.gt.0) c = 1.

   11 ii        = ii  +ki
      ss        = si(ii)

   21 i         = i  +k
      if (i.ne.n.and.d*(s(i)  -ss).le.0.) go to 21

      j         = i
      i         = i  -k
      ss        = ss  -s(i)
      fpppp     = c*(fppp(j)  -fppp(i))/(s(j)  -s(i))
      ff        = fppp(i)  +.25*ss*fpppp
      ff        = fpp(i)  +ss*ff/3.
      ff        = fp(i)  +.5*ss*ff
      fi(ii)    = f(i)  +ss*ff

      if (ii.ne.nin) go to 11

      return

      end
