      subroutine  splif(m,n,s,f,fp,fpp,fppp,km,vm,kn,vn,mode,fqm,ind)
c
c     ******************************************************************
c     *                                                                *
c     *   spline fit - jameson                                         *
c     *                                                                *
c     ******************************************************************
c
c     cubic spline fit of f(i) over s(i) from i = m to n.
c     m may be larger or smaller than n,
c     while s(i) may be increasing or decreasing,
c     as long as it continues in the same direction.
c     fp, fpp and fppp are the first, second and third derivatives of f.
c     vm and vn are end values at i = m and n
c     of the first, second or third derivatives of f,
c     depending on whether km or kn = 1,2 or 3.
c     the integral of f is placed in fppp if mode is greater than 0.
c     ind is set to zero if the data is illegal.
c
c     ******************************************************************
c
c     programmed by antony jameson july 1970 for the ibm 1130 computer.
c     revised december 2001 for fortran 95 compliance.
c
c     ******************************************************************
c
      dimension   s(*),f(*),fp(*),fpp(*),fppp(*)
c
c     ******************************************************************
c
      ind       = 0
      if (m.eq.n) return
      k         = sign(1,(n  -m))

      i         = m
      j         = m  +k
      ds        = s(j)  -s(i)
      d         = ds
      if (ds.eq.0.) return

      df        = (f(j)  -f(i))/ds

      if (km.le.1) then
         u         = .5
         v         = 3.*(df  -vm)/ds
      else if (km.eq.2) then
         u         = 0.
         v         = vm
      else
         u         = -1.
         v         = -ds*vm
      end if

      go to 21

   11 i         = j
      j         = j  +k
      ds        = s(j)  -s(i)
      if (d*ds.le.0.) return

      df        = (f(j)  -f(i))/ds
      b         = 1./(ds  +ds  +u)
      u         = b*ds
      v         = b*(6.*df  -v)

   21 fp(i)     = u
      fpp(i)    = v
      u         = (2.  -u)*ds
      v         = 6.*df  +ds*v

      if (j.ne.n) go to 11

      if (kn.le.1) then
         v         = (6.*vn  -v)/u
      else if (kn.eq.2) then
         v         = vn
      else
         v         = (ds*vn  +fpp(i))/(1.  +fp(i))
      end if

      b         = v
      d         = ds

   31 ds        = s(j)  -s(i)
      u         = fpp(i)  -fp(i)*v
      fppp(i)   = (v  -u)/ds
      fpp(i)    = u
      fp(i)     = (f(j)  -f(i))/ds  -ds*(v  +u  +u)/6.
      v         = u
      j         = i
      i         = i  -k

      if (j.ne.m) go to 31

      i         = n  -k
      fppp(n)   = fppp(i)
      fpp(n)    = b
      fp(n)     = df  +d*(fpp(i)  +b  +b)/6.
      ind       = 1

      if (mode.le.0) return
      fppp(j)   = fqm
      v         = fpp(j)

   41 i         = j
      j         = j  +k
      ds        = s(j)  -s(i)
      u         = fpp(j)
      fppp(j)   = fppp(i)  +.5*ds*(f(i)  +f(j)  -ds*ds*(u  +v)/12.)
      v         = u

      if (j.ne.n) go to 41

      return

      end
